#ifndef lint
static const char	RCSid[] = "$Id: ambcomp.c,v 2.86 2021/02/17 01:29:22 greg Exp $";
#endif
/*
 * Routines to compute "ambient" values using Monte Carlo
 *
 *  Hessian calculations based on "Practical Hessian-Based Error Control
 *	for Irradiance Caching" by Schwarzhaupt, Wann Jensen, & Jarosz
 *	from ACM SIGGRAPH Asia 2012 conference proceedings.
 *
 *  Added book-keeping optimization to avoid calculations that would
 *	cancel due to traversal both directions on edges that are adjacent
 *	to same-valued triangles.  This cuts about half of Hessian math.
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include  "ray.h"
#include  "ambient.h"
#include  "random.h"

#ifndef MINADIV
#define MINADIV		7	/* minimum # divisions in each dimension */
#endif

extern void		SDsquare2disk(double ds[2], double seedx, double seedy);

typedef struct {
	COLOR	v;		/* hemisphere sample value */
	float	d;		/* reciprocal distance */
	FVECT	p;		/* intersection point */
} AMBSAMP;		/* sample value */

typedef struct {
	RAY	*rp;		/* originating ray sample */
	int	ns;		/* number of samples per axis */
	int	sampOK;		/* acquired full sample set? */
	COLOR	acoef;		/* division contribution coefficient */
	double	acol[3];	/* accumulated color */
	FVECT	ux, uy;		/* tangent axis unit vectors */
	AMBSAMP	sa[1];		/* sample array (extends struct) */
}  AMBHEMI;		/* ambient sample hemisphere */

#define AI(h,i,j)	((i)*(h)->ns + (j))
#define ambsam(h,i,j)	(h)->sa[AI(h,i,j)]

typedef struct {
	FVECT	r_i, r_i1, e_i, rcp, rI2_eJ2;
	double	I1, I2;
} FFTRI;		/* vectors and coefficients for Hessian calculation */


static int
ambcollision(				/* proposed direciton collides? */
	AMBHEMI	*hp,
	int	i,
	int	j,
	FVECT	dv
)
{
	double	cos_thresh;
	int	ii, jj;
					/* min. spacing = 1/4th division */
	cos_thresh = (PI/4.)/(double)hp->ns;
	cos_thresh = 1. - .5*cos_thresh*cos_thresh;
					/* check existing neighbors */
	for (ii = i-1; ii <= i+1; ii++) {
		if (ii < 0) continue;
		if (ii >= hp->ns) break;
		for (jj = j-1; jj <= j+1; jj++) {
			AMBSAMP	*ap;
			FVECT	avec;
			double	dprod;
			if (jj < 0) continue;
			if (jj >= hp->ns) break;
			if ((ii==i) & (jj==j)) continue;
			ap = &ambsam(hp,ii,jj);
			if (ap->d <= .5/FHUGE)
				continue;	/* no one home */
			VSUB(avec, ap->p, hp->rp->rop);
			dprod = DOT(avec, dv);
			if (dprod >= cos_thresh*VLEN(avec))
				return(1);	/* collision */
		}
	}
	return(0);			/* nothing to worry about */
}


static int
ambsample(				/* initial ambient division sample */
	AMBHEMI	*hp,
	int	i,
	int	j,
	int	n
)
{
	AMBSAMP	*ap = &ambsam(hp,i,j);
	RAY	ar;
	int	hlist[3], ii;
	double	spt[2], zd;
					/* generate hemispherical sample */
					/* ambient coefficient for weight */
	if (ambacc > FTINY)
		setcolor(ar.rcoef, AVGREFL, AVGREFL, AVGREFL);
	else
		copycolor(ar.rcoef, hp->acoef);
	if (rayorigin(&ar, AMBIENT, hp->rp, ar.rcoef) < 0)
		return(0);
	if (ambacc > FTINY) {
		multcolor(ar.rcoef, hp->acoef);
		scalecolor(ar.rcoef, 1./AVGREFL);
	}
	hlist[0] = hp->rp->rno;
	hlist[1] = j;
	hlist[2] = i;
	multisamp(spt, 2, urand(ilhash(hlist,3)+n));
resample:
	SDsquare2disk(spt, (j+spt[1])/hp->ns, (i+spt[0])/hp->ns);
	zd = sqrt(1. - spt[0]*spt[0] - spt[1]*spt[1]);
	for (ii = 3; ii--; )
		ar.rdir[ii] =	spt[0]*hp->ux[ii] +
				spt[1]*hp->uy[ii] +
				zd*hp->rp->ron[ii];
	checknorm(ar.rdir);
					/* avoid coincident samples */
	if (!n && ambcollision(hp, i, j, ar.rdir)) {
		spt[0] = frandom(); spt[1] = frandom();
		goto resample;		/* reject this sample */
	}
	dimlist[ndims++] = AI(hp,i,j) + 90171;
	rayvalue(&ar);			/* evaluate ray */
	ndims--;
	zd = raydistance(&ar);
	if (zd <= FTINY)
		return(0);		/* should never happen */
	multcolor(ar.rcol, ar.rcoef);	/* apply coefficient */
	if (zd*ap->d < 1.0)		/* new/closer distance? */
		ap->d = 1.0/zd;
	if (!n) {			/* record first vertex & value */
		if (zd > 10.0*thescene.cusize + 1000.)
			zd = 10.0*thescene.cusize + 1000.;
		VSUM(ap->p, ar.rorg, ar.rdir, zd);
		copycolor(ap->v, ar.rcol);
	} else {			/* else update recorded value */
		hp->acol[RED] -= colval(ap->v,RED);
		hp->acol[GRN] -= colval(ap->v,GRN);
		hp->acol[BLU] -= colval(ap->v,BLU);
		zd = 1.0/(double)(n+1);
		scalecolor(ar.rcol, zd);
		zd *= (double)n;
		scalecolor(ap->v, zd);
		addcolor(ap->v, ar.rcol);
	}
	addcolor(hp->acol, ap->v);	/* add to our sum */
	return(1);
}


/* Estimate variance based on ambient division differences */
static float *
getambdiffs(AMBHEMI *hp)
{
	const double	normf = 1./bright(hp->acoef);
	float	*earr = (float *)calloc(hp->ns*hp->ns, sizeof(float));
	float	*ep;
	AMBSAMP	*ap;
	double	b, b1, d2;
	int	i, j;

	if (earr == NULL)		/* out of memory? */
		return(NULL);
					/* sum squared neighbor diffs */
	for (ap = hp->sa, ep = earr, i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j++, ap++, ep++) {
		b = bright(ap[0].v);
		if (i) {		/* from above */
			b1 = bright(ap[-hp->ns].v);
			d2 = b - b1;
			d2 *= d2*normf/(b + b1);
			ep[0] += d2;
			ep[-hp->ns] += d2;
		}
		if (!j) continue;
					/* from behind */
		b1 = bright(ap[-1].v);
		d2 = b - b1;
		d2 *= d2*normf/(b + b1);
		ep[0] += d2;
		ep[-1] += d2;
		if (!i) continue;
					/* diagonal */
		b1 = bright(ap[-hp->ns-1].v);
		d2 = b - b1;
		d2 *= d2*normf/(b + b1);
		ep[0] += d2;
		ep[-hp->ns-1] += d2;
	    }
					/* correct for number of neighbors */
	earr[0] *= 8./3.;
	earr[hp->ns-1] *= 8./3.;
	earr[(hp->ns-1)*hp->ns] *= 8./3.;
	earr[(hp->ns-1)*hp->ns + hp->ns-1] *= 8./3.;
	for (i = 1; i < hp->ns-1; i++) {
		earr[i*hp->ns] *= 8./5.;
		earr[i*hp->ns + hp->ns-1] *= 8./5.;
	}
	for (j = 1; j < hp->ns-1; j++) {
		earr[j] *= 8./5.;
		earr[(hp->ns-1)*hp->ns + j] *= 8./5.;
	}
	return(earr);
}


/* Perform super-sampling on hemisphere (introduces bias) */
static void
ambsupersamp(AMBHEMI *hp, int cnt)
{
	float	*earr = getambdiffs(hp);
	double	e2rem = 0;
	float	*ep;
	int	i, j, n, nss;

	if (earr == NULL)		/* just skip calc. if no memory */
		return;
					/* accumulate estimated variances */
	for (ep = earr + hp->ns*hp->ns; ep > earr; )
		e2rem += *--ep;
	ep = earr;			/* perform super-sampling */
	for (i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j++) {
		if (e2rem <= FTINY)
			goto done;	/* nothing left to do */
		nss = *ep/e2rem*cnt + frandom();
		for (n = 1; n <= nss && ambsample(hp,i,j,n); n++)
			if (!--cnt) goto done;
		e2rem -= *ep++;		/* update remainder */
	}
done:
	free(earr);
}


static AMBHEMI *
samp_hemi(				/* sample indirect hemisphere */
	COLOR	rcol,
	RAY	*r,
	double	wt
)
{
	AMBHEMI	*hp;
	double	d;
	int	n, i, j;
					/* insignificance check */
	if (bright(rcol) <= FTINY)
		return(NULL);
					/* set number of divisions */
	if (ambacc <= FTINY &&
			wt > (d = 0.8*intens(rcol)*r->rweight/(ambdiv*minweight)))
		wt = d;			/* avoid ray termination */
	n = sqrt(ambdiv * wt) + 0.5;
	i = 1 + (MINADIV-1)*(ambacc > FTINY);
	if (n < i)			/* use minimum number of samples? */
		n = i;
					/* allocate sampling array */
	hp = (AMBHEMI *)malloc(sizeof(AMBHEMI) + sizeof(AMBSAMP)*(n*n - 1));
	if (hp == NULL)
		error(SYSTEM, "out of memory in samp_hemi");
	hp->rp = r;
	hp->ns = n;
	hp->acol[RED] = hp->acol[GRN] = hp->acol[BLU] = 0.0;
	memset(hp->sa, 0, sizeof(AMBSAMP)*n*n);
	hp->sampOK = 0;
					/* assign coefficient */
	copycolor(hp->acoef, rcol);
	d = 1.0/(n*n);
	scalecolor(hp->acoef, d);
					/* make tangent plane axes */
	if (!getperpendicular(hp->ux, r->ron, 1))
		error(CONSISTENCY, "bad ray direction in samp_hemi");
	VCROSS(hp->uy, r->ron, hp->ux);
					/* sample divisions */
	for (i = hp->ns; i--; )
	    for (j = hp->ns; j--; )
		hp->sampOK += ambsample(hp, i, j, 0);
	copycolor(rcol, hp->acol);
	if (!hp->sampOK) {		/* utter failure? */
		free(hp);
		return(NULL);
	}
	if (hp->sampOK < hp->ns*hp->ns) {
		hp->sampOK *= -1;	/* soft failure */
		return(hp);
	}
	if (hp->sampOK <= MINADIV*MINADIV)
		return(hp);		/* don't bother super-sampling */
	n = ambssamp*wt + 0.5;
	if (n > 8) {			/* perform super-sampling? */
		ambsupersamp(hp, n);
		copycolor(rcol, hp->acol);
	}
	return(hp);			/* all is well */
}


/* Return brightness of farthest ambient sample */
static double
back_ambval(AMBHEMI *hp, const int n1, const int n2, const int n3)
{
	if (hp->sa[n1].d <= hp->sa[n2].d) {
		if (hp->sa[n1].d <= hp->sa[n3].d)
			return(colval(hp->sa[n1].v,CIEY));
		return(colval(hp->sa[n3].v,CIEY));
	}
	if (hp->sa[n2].d <= hp->sa[n3].d)
		return(colval(hp->sa[n2].v,CIEY));
	return(colval(hp->sa[n3].v,CIEY));
}


/* Compute vectors and coefficients for Hessian/gradient calcs */
static void
comp_fftri(FFTRI *ftp, AMBHEMI *hp, const int n0, const int n1)
{
	double	rdot_cp, dot_e, dot_er, rdot_r, rdot_r1, J2;
	int	ii;

	VSUB(ftp->r_i, hp->sa[n0].p, hp->rp->rop);
	VSUB(ftp->r_i1, hp->sa[n1].p, hp->rp->rop);
	VSUB(ftp->e_i, hp->sa[n1].p, hp->sa[n0].p);
	VCROSS(ftp->rcp, ftp->r_i, ftp->r_i1);
	rdot_cp = 1.0/DOT(ftp->rcp,ftp->rcp);
	dot_e = DOT(ftp->e_i,ftp->e_i);
	dot_er = DOT(ftp->e_i, ftp->r_i);
	rdot_r = 1.0/DOT(ftp->r_i,ftp->r_i);
	rdot_r1 = 1.0/DOT(ftp->r_i1,ftp->r_i1);
	ftp->I1 = acos( DOT(ftp->r_i, ftp->r_i1) * sqrt(rdot_r*rdot_r1) ) *
			sqrt( rdot_cp );
	ftp->I2 = ( DOT(ftp->e_i, ftp->r_i1)*rdot_r1 - dot_er*rdot_r +
			dot_e*ftp->I1 )*0.5*rdot_cp;
	J2 =  ( 0.5*(rdot_r - rdot_r1) - dot_er*ftp->I2 ) / dot_e;
	for (ii = 3; ii--; )
		ftp->rI2_eJ2[ii] = ftp->I2*ftp->r_i[ii] + J2*ftp->e_i[ii];
}


/* Compose 3x3 matrix from two vectors */
static void
compose_matrix(FVECT mat[3], FVECT va, FVECT vb)
{
	mat[0][0] = 2.0*va[0]*vb[0];
	mat[1][1] = 2.0*va[1]*vb[1];
	mat[2][2] = 2.0*va[2]*vb[2];
	mat[0][1] = mat[1][0] = va[0]*vb[1] + va[1]*vb[0];
	mat[0][2] = mat[2][0] = va[0]*vb[2] + va[2]*vb[0];
	mat[1][2] = mat[2][1] = va[1]*vb[2] + va[2]*vb[1];
}


/* Compute partial 3x3 Hessian matrix for edge */
static void
comp_hessian(FVECT hess[3], FFTRI *ftp, FVECT nrm)
{
	FVECT	ncp;
	FVECT	m1[3], m2[3], m3[3], m4[3];
	double	d1, d2, d3, d4;
	double	I3, J3, K3;
	int	i, j;
					/* compute intermediate coefficients */
	d1 = 1.0/DOT(ftp->r_i,ftp->r_i);
	d2 = 1.0/DOT(ftp->r_i1,ftp->r_i1);
	d3 = 1.0/DOT(ftp->e_i,ftp->e_i);
	d4 = DOT(ftp->e_i, ftp->r_i);
	I3 = ( DOT(ftp->e_i, ftp->r_i1)*d2*d2 - d4*d1*d1 + 3.0/d3*ftp->I2 )
			/ ( 4.0*DOT(ftp->rcp,ftp->rcp) );
	J3 = 0.25*d3*(d1*d1 - d2*d2) - d4*d3*I3;
	K3 = d3*(ftp->I2 - I3/d1 - 2.0*d4*J3);
					/* intermediate matrices */
	VCROSS(ncp, nrm, ftp->e_i);
	compose_matrix(m1, ncp, ftp->rI2_eJ2);
	compose_matrix(m2, ftp->r_i, ftp->r_i);
	compose_matrix(m3, ftp->e_i, ftp->e_i);
	compose_matrix(m4, ftp->r_i, ftp->e_i);
	d1 = DOT(nrm, ftp->rcp);
	d2 = -d1*ftp->I2;
	d1 *= 2.0;
	for (i = 3; i--; )		/* final matrix sum */
	    for (j = 3; j--; ) {
		hess[i][j] = m1[i][j] + d1*( I3*m2[i][j] + K3*m3[i][j] +
						2.0*J3*m4[i][j] );
		hess[i][j] += d2*(i==j);
		hess[i][j] *= -1.0/PI;
	    }
}


/* Reverse hessian calculation result for edge in other direction */
static void
rev_hessian(FVECT hess[3])
{
	int	i;

	for (i = 3; i--; ) {
		hess[i][0] = -hess[i][0];
		hess[i][1] = -hess[i][1];
		hess[i][2] = -hess[i][2];
	}
}


/* Add to radiometric Hessian from the given triangle */
static void
add2hessian(FVECT hess[3], FVECT ehess1[3],
		FVECT ehess2[3], FVECT ehess3[3], double v)
{
	int	i, j;

	for (i = 3; i--; )
	    for (j = 3; j--; )
		hess[i][j] += v*( ehess1[i][j] + ehess2[i][j] + ehess3[i][j] );
}


/* Compute partial displacement form factor gradient for edge */
static void
comp_gradient(FVECT grad, FFTRI *ftp, FVECT nrm)
{
	FVECT	ncp;
	double	f1;
	int	i;

	f1 = 2.0*DOT(nrm, ftp->rcp);
	VCROSS(ncp, nrm, ftp->e_i);
	for (i = 3; i--; )
		grad[i] = (0.5/PI)*( ftp->I1*ncp[i] + f1*ftp->rI2_eJ2[i] );
}


/* Reverse gradient calculation result for edge in other direction */
static void
rev_gradient(FVECT grad)
{
	grad[0] = -grad[0];
	grad[1] = -grad[1];
	grad[2] = -grad[2];
}


/* Add to displacement gradient from the given triangle */
static void
add2gradient(FVECT grad, FVECT egrad1, FVECT egrad2, FVECT egrad3, double v)
{
	int	i;

	for (i = 3; i--; )
		grad[i] += v*( egrad1[i] + egrad2[i] + egrad3[i] );
}


/* Compute anisotropic radii and eigenvector directions */
static void
eigenvectors(FVECT uv[2], float ra[2], FVECT hessian[3])
{
	double	hess2[2][2];
	FVECT	a, b;
	double	evalue[2], slope1, xmag1;
	int	i;
					/* project Hessian to sample plane */
	for (i = 3; i--; ) {
		a[i] = DOT(hessian[i], uv[0]);
		b[i] = DOT(hessian[i], uv[1]);
	}
	hess2[0][0] = DOT(uv[0], a);
	hess2[0][1] = DOT(uv[0], b);
	hess2[1][0] = DOT(uv[1], a);
	hess2[1][1] = DOT(uv[1], b);
					/* compute eigenvalue(s) */
	i = quadratic(evalue, 1.0, -hess2[0][0]-hess2[1][1],
			hess2[0][0]*hess2[1][1]-hess2[0][1]*hess2[1][0]);
	if (i == 1)			/* double-root (circle) */
		evalue[1] = evalue[0];
	if (!i || ((evalue[0] = fabs(evalue[0])) <= FTINY*FTINY) |
			((evalue[1] = fabs(evalue[1])) <= FTINY*FTINY) ) {
		ra[0] = ra[1] = maxarad;
		return;
	}
	if (evalue[0] > evalue[1]) {
		ra[0] = sqrt(sqrt(4.0/evalue[0]));
		ra[1] = sqrt(sqrt(4.0/evalue[1]));
		slope1 = evalue[1];
	} else {
		ra[0] = sqrt(sqrt(4.0/evalue[1]));
		ra[1] = sqrt(sqrt(4.0/evalue[0]));
		slope1 = evalue[0];
	}
					/* compute unit eigenvectors */
	if (fabs(hess2[0][1]) <= FTINY)
		return;			/* uv OK as is */
	slope1 = (slope1 - hess2[0][0]) / hess2[0][1];
	xmag1 = sqrt(1.0/(1.0 + slope1*slope1));
	for (i = 3; i--; ) {
		b[i] = xmag1*uv[0][i] + slope1*xmag1*uv[1][i];
		a[i] = slope1*xmag1*uv[0][i] - xmag1*uv[1][i];
	}
	VCOPY(uv[0], a);
	VCOPY(uv[1], b);
}


static void
ambHessian(				/* anisotropic radii & pos. gradient */
	AMBHEMI	*hp,
	FVECT	uv[2],			/* returned */
	float	ra[2],			/* returned (optional) */
	float	pg[2]			/* returned (optional) */
)
{
	static char	memerrmsg[] = "out of memory in ambHessian()";
	FVECT		(*hessrow)[3] = NULL;
	FVECT		*gradrow = NULL;
	FVECT		hessian[3];
	FVECT		gradient;
	FFTRI		fftr;
	int		i, j;
					/* be sure to assign unit vectors */
	VCOPY(uv[0], hp->ux);
	VCOPY(uv[1], hp->uy);
			/* clock-wise vertex traversal from sample POV */
	if (ra != NULL) {		/* initialize Hessian row buffer */
		hessrow = (FVECT (*)[3])malloc(sizeof(FVECT)*3*(hp->ns-1));
		if (hessrow == NULL)
			error(SYSTEM, memerrmsg);
		memset(hessian, 0, sizeof(hessian));
	} else if (pg == NULL)		/* bogus call? */
		return;
	if (pg != NULL) {		/* initialize form factor row buffer */
		gradrow = (FVECT *)malloc(sizeof(FVECT)*(hp->ns-1));
		if (gradrow == NULL)
			error(SYSTEM, memerrmsg);
		memset(gradient, 0, sizeof(gradient));
	}
					/* compute first row of edges */
	for (j = 0; j < hp->ns-1; j++) {
		comp_fftri(&fftr, hp, AI(hp,0,j), AI(hp,0,j+1));
		if (hessrow != NULL)
			comp_hessian(hessrow[j], &fftr, hp->rp->ron);
		if (gradrow != NULL)
			comp_gradient(gradrow[j], &fftr, hp->rp->ron);
	}
					/* sum each row of triangles */
	for (i = 0; i < hp->ns-1; i++) {
	    FVECT	hesscol[3];	/* compute first vertical edge */
	    FVECT	gradcol;
	    comp_fftri(&fftr, hp, AI(hp,i,0), AI(hp,i+1,0));
	    if (hessrow != NULL)
		comp_hessian(hesscol, &fftr, hp->rp->ron);
	    if (gradrow != NULL)
		comp_gradient(gradcol, &fftr, hp->rp->ron);
	    for (j = 0; j < hp->ns-1; j++) {
		FVECT	hessdia[3];	/* compute triangle contributions */
		FVECT	graddia;
		double	backg;
		backg = back_ambval(hp, AI(hp,i,j),
					AI(hp,i,j+1), AI(hp,i+1,j));
					/* diagonal (inner) edge */
		comp_fftri(&fftr, hp, AI(hp,i,j+1), AI(hp,i+1,j));
		if (hessrow != NULL) {
		    comp_hessian(hessdia, &fftr, hp->rp->ron);
		    rev_hessian(hesscol);
		    add2hessian(hessian, hessrow[j], hessdia, hesscol, backg);
		}
		if (gradrow != NULL) {
		    comp_gradient(graddia, &fftr, hp->rp->ron);
		    rev_gradient(gradcol);
		    add2gradient(gradient, gradrow[j], graddia, gradcol, backg);
		}
					/* initialize edge in next row */
		comp_fftri(&fftr, hp, AI(hp,i+1,j+1), AI(hp,i+1,j));
		if (hessrow != NULL)
		    comp_hessian(hessrow[j], &fftr, hp->rp->ron);
		if (gradrow != NULL)
		    comp_gradient(gradrow[j], &fftr, hp->rp->ron);
					/* new column edge & paired triangle */
		backg = back_ambval(hp, AI(hp,i+1,j+1),
					AI(hp,i+1,j), AI(hp,i,j+1));
		comp_fftri(&fftr, hp, AI(hp,i,j+1), AI(hp,i+1,j+1));
		if (hessrow != NULL) {
		    comp_hessian(hesscol, &fftr, hp->rp->ron);
		    rev_hessian(hessdia);
		    add2hessian(hessian, hessrow[j], hessdia, hesscol, backg);
		    if (i < hp->ns-2)
			rev_hessian(hessrow[j]);
		}
		if (gradrow != NULL) {
		    comp_gradient(gradcol, &fftr, hp->rp->ron);
		    rev_gradient(graddia);
		    add2gradient(gradient, gradrow[j], graddia, gradcol, backg);
		    if (i < hp->ns-2)
			rev_gradient(gradrow[j]);
		}
	    }
	}
					/* release row buffers */
	if (hessrow != NULL) free(hessrow);
	if (gradrow != NULL) free(gradrow);
	
	if (ra != NULL)			/* extract eigenvectors & radii */
		eigenvectors(uv, ra, hessian);
	if (pg != NULL) {		/* tangential position gradient */
		pg[0] = DOT(gradient, uv[0]);
		pg[1] = DOT(gradient, uv[1]);
	}
}


/* Compute direction gradient from a hemispherical sampling */
static void
ambdirgrad(AMBHEMI *hp, FVECT uv[2], float dg[2])
{
	AMBSAMP	*ap;
	double	dgsum[2];
	int	n;
	FVECT	vd;
	double	gfact;

	dgsum[0] = dgsum[1] = 0.0;	/* sum values times -tan(theta) */
	for (ap = hp->sa, n = hp->ns*hp->ns; n--; ap++) {
					/* use vector for azimuth + 90deg */
		VSUB(vd, ap->p, hp->rp->rop);
					/* brightness over cosine factor */
		gfact = colval(ap->v,CIEY) / DOT(hp->rp->ron, vd);
					/* sine = proj_radius/vd_length */
		dgsum[0] -= DOT(uv[1], vd) * gfact;
		dgsum[1] += DOT(uv[0], vd) * gfact;
	}
	dg[0] = dgsum[0] / (hp->ns*hp->ns);
	dg[1] = dgsum[1] / (hp->ns*hp->ns);
}


/* Compute potential light leak direction flags for cache value */
static uint32
ambcorral(AMBHEMI *hp, FVECT uv[2], const double r0, const double r1)
{
	const double	max_d = 1.0/(minarad*ambacc + 0.001);
	const double	ang_res = 0.5*PI/hp->ns;
	const double	ang_step = ang_res/((int)(16/PI*ang_res) + 1.01);
	double		avg_d = 0;
	uint32		flgs = 0;
	FVECT		vec;
	double		u, v;
	double		ang, a1;
	int		i, j;
					/* don't bother for a few samples */
	if (hp->ns < 8)
		return(0);
					/* check distances overhead */
	for (i = hp->ns*3/4; i-- > hp->ns>>2; )
	    for (j = hp->ns*3/4; j-- > hp->ns>>2; )
		avg_d += ambsam(hp,i,j).d;
	avg_d *= 4.0/(hp->ns*hp->ns);
	if (avg_d*r0 >= 1.0)		/* ceiling too low for corral? */
		return(0);
	if (avg_d >= max_d)		/* insurance */
		return(0);
					/* else circle around perimeter */
	for (i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j += !i|(i==hp->ns-1) ? 1 : hp->ns-1) {
		AMBSAMP	*ap = &ambsam(hp,i,j);
		if ((ap->d <= FTINY) | (ap->d >= max_d))
			continue;	/* too far or too near */
		VSUB(vec, ap->p, hp->rp->rop);
		u = DOT(vec, uv[0]);
		v = DOT(vec, uv[1]);
		if ((r0*r0*u*u + r1*r1*v*v) * ap->d*ap->d <= u*u + v*v)
			continue;	/* occluder outside ellipse */
		ang = atan2a(v, u);	/* else set direction flags */
		for (a1 = ang-ang_res; a1 <= ang+ang_res; a1 += ang_step)
			flgs |= 1L<<(int)(16/PI*(a1 + 2.*PI*(a1 < 0)));
	    }
	return(flgs);
}


int
doambient(				/* compute ambient component */
	COLOR	rcol,			/* input/output color */
	RAY	*r,
	double	wt,
	FVECT	uv[2],			/* returned (optional) */
	float	ra[2],			/* returned (optional) */
	float	pg[2],			/* returned (optional) */
	float	dg[2],			/* returned (optional) */
	uint32	*crlp			/* returned (optional) */
)
{
	AMBHEMI	*hp = samp_hemi(rcol, r, wt);
	FVECT	my_uv[2];
	double	d, K;
	AMBSAMP	*ap;
	int	i;
					/* clear return values */
	if (uv != NULL)
		memset(uv, 0, sizeof(FVECT)*2);
	if (ra != NULL)
		ra[0] = ra[1] = 0.0;
	if (pg != NULL)
		pg[0] = pg[1] = 0.0;
	if (dg != NULL)
		dg[0] = dg[1] = 0.0;
	if (crlp != NULL)
		*crlp = 0;
	if (hp == NULL)			/* sampling falure? */
		return(0);

	if ((ra == NULL) & (pg == NULL) & (dg == NULL) ||
			(hp->sampOK < 0) | (hp->ns < MINADIV)) {
		free(hp);		/* Hessian not requested/possible */
		return(-1);		/* value-only return value */
	}
	if ((d = bright(rcol)) > FTINY) {	/* normalize Y values */
		d = 0.99*(hp->ns*hp->ns)/d;
		K = 0.01;
	} else {			/* or fall back on geometric Hessian */
		K = 1.0;
		pg = NULL;
		dg = NULL;
		crlp = NULL;
	}
	ap = hp->sa;			/* relative Y channel from here on... */
	for (i = hp->ns*hp->ns; i--; ap++)
		colval(ap->v,CIEY) = bright(ap->v)*d + K;

	if (uv == NULL)			/* make sure we have axis pointers */
		uv = my_uv;
					/* compute radii & pos. gradient */
	ambHessian(hp, uv, ra, pg);

	if (dg != NULL)			/* compute direction gradient */
		ambdirgrad(hp, uv, dg);

	if (ra != NULL) {		/* scale/clamp radii */
		if (pg != NULL) {
			if (ra[0]*(d = fabs(pg[0])) > 1.0)
				ra[0] = 1.0/d;
			if (ra[1]*(d = fabs(pg[1])) > 1.0)
				ra[1] = 1.0/d;
			if (ra[0] > ra[1])
				ra[0] = ra[1];
		}
		if (ra[0] < minarad) {
			ra[0] = minarad;
			if (ra[1] < minarad)
				ra[1] = minarad;
		}
		ra[0] *= d = 1.0/sqrt(wt);
		if ((ra[1] *= d) > 2.0*ra[0])
			ra[1] = 2.0*ra[0];
		if (ra[1] > maxarad) {
			ra[1] = maxarad;
			if (ra[0] > maxarad)
				ra[0] = maxarad;
		}
					/* flag encroached directions */
		if (crlp != NULL)
			*crlp = ambcorral(hp, uv, ra[0]*ambacc, ra[1]*ambacc);
		if (pg != NULL) {	/* cap gradient if necessary */
			d = pg[0]*pg[0]*ra[0]*ra[0] + pg[1]*pg[1]*ra[1]*ra[1];
			if (d > 1.0) {
				d = 1.0/sqrt(d);
				pg[0] *= d;
				pg[1] *= d;
			}
		}
	}
	free(hp);			/* clean up and return */
	return(1);
}

