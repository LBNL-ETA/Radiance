#ifndef lint
static const char	RCSid[] = "$Id: ambcomp.c,v 2.50 2014/05/07 16:02:26 greg Exp $";
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

#ifdef NEWAMB

extern void		SDsquare2disk(double ds[2], double seedx, double seedy);

				/* vertex direction bit positions */
#define	VDB_xy	0
#define VDB_y	01
#define VDB_x	02
#define VDB_Xy	03
#define VDB_xY	04
#define VDB_X	05
#define VDB_Y	06
#define VDB_XY	07
				/* get opposite vertex direction bit */
#define VDB_OPP(f)	(~(f) & 07)
				/* adjacent triangle vertex flags */
static const int  adjacent_trifl[8] = {
			0,			/* forbidden diagonal */
			1<<VDB_x|1<<VDB_y|1<<VDB_Xy,
			1<<VDB_y|1<<VDB_x|1<<VDB_xY,
			1<<VDB_y|1<<VDB_Xy|1<<VDB_X,
			1<<VDB_x|1<<VDB_xY|1<<VDB_Y,
			1<<VDB_Xy|1<<VDB_X|1<<VDB_Y,
			1<<VDB_xY|1<<VDB_Y|1<<VDB_X,
			0,			/* forbidden diagonal */
		};

typedef struct {
	COLOR	v;		/* hemisphere sample value */
	float	d;		/* reciprocal distance (1/rt) */
	FVECT	p;		/* intersection point */
} AMBSAMP;		/* sample value */

typedef struct {
	RAY	*rp;		/* originating ray sample */
	FVECT	ux, uy;		/* tangent axis unit vectors */
	int	ns;		/* number of samples per axis */
	COLOR	acoef;		/* division contribution coefficient */
	AMBSAMP	sa[1];		/* sample array (extends struct) */
}  AMBHEMI;		/* ambient sample hemisphere */

#define ambndx(h,i,j)	((i)*(h)->ns + (j))
#define ambsam(h,i,j)	(h)->sa[ambndx(h,i,j)]

typedef struct {
	FVECT	r_i, r_i1, e_i, rcp, rI2_eJ2;
	double	I1, I2;
	int	valid;
} FFTRI;		/* vectors and coefficients for Hessian calculation */


/* Get index for adjacent vertex */
static int
adjacent_verti(AMBHEMI *hp, int i, int j, int dbit)
{
	int	i0 = i*hp->ns + j;

	switch (dbit) {
	case VDB_y:	return(i0 - hp->ns);
	case VDB_x:	return(i0 - 1);
	case VDB_Xy:	return(i0 - hp->ns + 1);
	case VDB_xY:	return(i0 + hp->ns - 1);
	case VDB_X:	return(i0 + 1);
	case VDB_Y:	return(i0 + hp->ns);
				/* the following should never occur */
	case VDB_xy:	return(i0 - hp->ns - 1);
	case VDB_XY:	return(i0 + hp->ns + 1);
	}
	return(-1);
}


/* Get vertex direction bit for the opposite edge to complete triangle */
static int
vdb_edge(int db1, int db2)
{
	switch (db1) {
	case VDB_x:	return(db2==VDB_y ? VDB_Xy : VDB_Y);
	case VDB_y:	return(db2==VDB_x ? VDB_xY : VDB_X);
	case VDB_X:	return(db2==VDB_Xy ? VDB_y : VDB_xY);
	case VDB_Y:	return(db2==VDB_xY ? VDB_x : VDB_Xy);
	case VDB_xY:	return(db2==VDB_x ? VDB_y : VDB_X);
	case VDB_Xy:	return(db2==VDB_y ? VDB_x : VDB_Y);
	}
	error(INTERNAL, "forbidden diagonal in vdb_edge()");
	return(-1);
}


static AMBHEMI *
inithemi(			/* initialize sampling hemisphere */
	COLOR	ac,
	RAY	*r,
	double	wt
)
{
	AMBHEMI	*hp;
	double	d;
	int	n, i;
					/* set number of divisions */
	if (ambacc <= FTINY &&
			wt > (d = 0.8*intens(ac)*r->rweight/(ambdiv*minweight)))
		wt = d;			/* avoid ray termination */
	n = sqrt(ambdiv * wt) + 0.5;
	i = 1 + 5*(ambacc > FTINY);	/* minimum number of samples */
	if (n < i)
		n = i;
					/* allocate sampling array */
	hp = (AMBHEMI *)malloc(sizeof(AMBHEMI) + sizeof(AMBSAMP)*(n*n - 1));
	if (hp == NULL)
		return(NULL);
	hp->rp = r;
	hp->ns = n;
					/* assign coefficient */
	copycolor(hp->acoef, ac);
	d = 1.0/(n*n);
	scalecolor(hp->acoef, d);
					/* make tangent plane axes */
	hp->uy[0] = 0.5 - frandom();
	hp->uy[1] = 0.5 - frandom();
	hp->uy[2] = 0.5 - frandom();
	for (i = 3; i--; )
		if ((-0.6 < r->ron[i]) & (r->ron[i] < 0.6))
			break;
	if (i < 0)
		error(CONSISTENCY, "bad ray direction in inithemi");
	hp->uy[i] = 1.0;
	VCROSS(hp->ux, hp->uy, r->ron);
	normalize(hp->ux);
	VCROSS(hp->uy, r->ron, hp->ux);
					/* we're ready to sample */
	return(hp);
}


/* Sample ambient division and apply weighting coefficient */
static int
getambsamp(RAY *arp, AMBHEMI *hp, int i, int j, int n)
{
	int	hlist[3], ii;
	double	spt[2], zd;
					/* ambient coefficient for weight */
	if (ambacc > FTINY)
		setcolor(arp->rcoef, AVGREFL, AVGREFL, AVGREFL);
	else
		copycolor(arp->rcoef, hp->acoef);
	if (rayorigin(arp, AMBIENT, hp->rp, arp->rcoef) < 0)
		return(0);
	if (ambacc > FTINY) {
		multcolor(arp->rcoef, hp->acoef);
		scalecolor(arp->rcoef, 1./AVGREFL);
	}
	hlist[0] = hp->rp->rno;
	hlist[1] = j;
	hlist[2] = i;
	multisamp(spt, 2, urand(ilhash(hlist,3)+n));
	if (!n) {			/* avoid border samples for n==0 */
		if ((spt[0] < 0.1) | (spt[0] >= 0.9))
			spt[0] = 0.1 + 0.8*frandom();
		if ((spt[1] < 0.1) | (spt[1] >= 0.9))
			spt[1] = 0.1 + 0.8*frandom();
	}
	SDsquare2disk(spt, (j+spt[1])/hp->ns, (i+spt[0])/hp->ns);
	zd = sqrt(1. - spt[0]*spt[0] - spt[1]*spt[1]);
	for (ii = 3; ii--; )
		arp->rdir[ii] =	spt[0]*hp->ux[ii] +
				spt[1]*hp->uy[ii] +
				zd*hp->rp->ron[ii];
	checknorm(arp->rdir);
	dimlist[ndims++] = ambndx(hp,i,j) + 90171;
	rayvalue(arp);			/* evaluate ray */
	ndims--;			/* apply coefficient */
	multcolor(arp->rcol, arp->rcoef);
	return(1);
}


static AMBSAMP *
ambsample(				/* initial ambient division sample */
	AMBHEMI	*hp,
	int	i,
	int	j
)
{
	AMBSAMP	*ap = &ambsam(hp,i,j);
	RAY	ar;
					/* generate hemispherical sample */
	if (!getambsamp(&ar, hp, i, j, 0) || ar.rt <= FTINY) {
		memset(ap, 0, sizeof(AMBSAMP));
		return(NULL);
	}
	ap->d = 1.0/ar.rt;		/* limit vertex distance */
	if (ar.rt > 10.0*thescene.cusize)
		ar.rt = 10.0*thescene.cusize;
	VSUM(ap->p, ar.rorg, ar.rdir, ar.rt);
	copycolor(ap->v, ar.rcol);
	return(ap);
}


/* Estimate errors based on ambient division differences */
static float *
getambdiffs(AMBHEMI *hp)
{
	float	*earr = (float *)calloc(hp->ns*hp->ns, sizeof(float));
	float	*ep;
	AMBSAMP	*ap;
	double	b, d2;
	int	i, j;

	if (earr == NULL)		/* out of memory? */
		return(NULL);
					/* compute squared neighbor diffs */
	for (ap = hp->sa, ep = earr, i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j++, ap++, ep++) {
		b = bright(ap[0].v);
		if (i) {		/* from above */
			d2 = b - bright(ap[-hp->ns].v);
			d2 *= d2;
			ep[0] += d2;
			ep[-hp->ns] += d2;
		}
		if (j) {		/* from behind */
			d2 = b - bright(ap[-1].v);
			d2 *= d2;
			ep[0] += d2;
			ep[-1] += d2;
		}
	    }
					/* correct for number of neighbors */
	earr[0] *= 2.f;
	earr[hp->ns-1] *= 2.f;
	earr[(hp->ns-1)*hp->ns] *= 2.f;
	earr[(hp->ns-1)*hp->ns + hp->ns-1] *= 2.f;
	for (i = 1; i < hp->ns-1; i++) {
		earr[i*hp->ns] *= 4./3.;
		earr[i*hp->ns + hp->ns-1] *= 4./3.;
	}
	for (j = 1; j < hp->ns-1; j++) {
		earr[j] *= 4./3.;
		earr[(hp->ns-1)*hp->ns + j] *= 4./3.;
	}
	return(earr);
}


/* Perform super-sampling on hemisphere (introduces bias) */
static void
ambsupersamp(double acol[3], AMBHEMI *hp, int cnt)
{
	float	*earr = getambdiffs(hp);
	double	e2sum = 0.0;
	AMBSAMP	*ap;
	RAY	ar;
	double	asum[3];
	float	*ep;
	int	i, j, n;

	if (earr == NULL)		/* just skip calc. if no memory */
		return;
					/* add up estimated variances */
	for (ep = earr + hp->ns*hp->ns; ep-- > earr; )
		e2sum += *ep;
	ep = earr;			/* perform super-sampling */
	for (ap = hp->sa, i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j++, ap++) {
		int	nss = *ep/e2sum*cnt + frandom();
		asum[0] = asum[1] = asum[2] = 0.0;
		for (n = 1; n <= nss; n++) {
			if (!getambsamp(&ar, hp, i, j, n)) {
				nss = n-1;
				break;
			}
			addcolor(asum, ar.rcol);
		}
		if (nss) {		/* update returned ambient value */
			const double	ssf = 1./(nss + 1);
			for (n = 3; n--; )
				acol[n] += ssf*asum[n] +
						(ssf - 1.)*colval(ap->v,n);
		}
		e2sum -= *ep++;		/* update remainders */
		cnt -= nss;
	}
	free(earr);
}


/* Compute vertex flags, indicating farthest in each direction */
static uby8 *
vertex_flags(AMBHEMI *hp)
{
	uby8	*vflags = (uby8 *)calloc(hp->ns*hp->ns, sizeof(uby8));
	uby8	*vf;
	AMBSAMP	*ap;
	int	i, j;

	if (vflags == NULL)
		error(SYSTEM, "out of memory in vertex_flags()");
	vf = vflags;
	ap = hp->sa;		/* compute farthest along first row */
	for (j = 0; j < hp->ns-1; j++, vf++, ap++)
		if (ap[0].d <= ap[1].d)
			vf[0] |= 1<<VDB_X;
		else
			vf[1] |= 1<<VDB_x;
	++vf; ++ap;
				/* flag subsequent rows */
	for (i = 1; i < hp->ns; i++) {
	    for (j = 0; j < hp->ns-1; j++, vf++, ap++) {
		if (ap[0].d <= ap[-hp->ns].d)	/* row before */
			vf[0] |= 1<<VDB_y;
		else
			vf[-hp->ns] |= 1<<VDB_Y;
		if (ap[0].d <= ap[1-hp->ns].d)	/* diagonal we care about */
			vf[0] |= 1<<VDB_Xy;
		else
			vf[1-hp->ns] |= 1<<VDB_xY;
		if (ap[0].d <= ap[1].d)		/* column after */
			vf[0] |= 1<<VDB_X;
		else
			vf[1] |= 1<<VDB_x;
	    }
	    if (ap[0].d <= ap[-hp->ns].d)	/* final column edge */
		vf[0] |= 1<<VDB_y;
	    else
		vf[-hp->ns] |= 1<<VDB_Y;
	    ++vf; ++ap;
	}
	return(vflags);
}


/* Return brightness of farthest ambient sample */
static double
back_ambval(AMBHEMI *hp, int i, int j, int dbit1, int dbit2, const uby8 *vflags)
{
	const int	v0 = ambndx(hp,i,j);
	const int	tflags = (1<<dbit1 | 1<<dbit2);
	int		v1, v2;

	if ((vflags[v0] & tflags) == tflags)	/* is v0 the farthest? */
		return(colval(hp->sa[v0].v,CIEY));
	v1 = adjacent_verti(hp, i, j, dbit1);
	if (vflags[v0] & 1<<dbit2)		/* v1 farthest if v0>v2 */
		return(colval(hp->sa[v1].v,CIEY));
	v2 = adjacent_verti(hp, i, j, dbit2);
	if (vflags[v0] & 1<<dbit1)		/* v2 farthest if v0>v1 */
		return(colval(hp->sa[v2].v,CIEY));
						/* else check if v1>v2 */
	if (vflags[v1] & 1<<vdb_edge(dbit1,dbit2))
		return(colval(hp->sa[v1].v,CIEY));
	return(colval(hp->sa[v2].v,CIEY));
}


/* Compute vectors and coefficients for Hessian/gradient calcs */
static void
comp_fftri(FFTRI *ftp, AMBHEMI *hp, int i, int j, int dbit, const uby8 *vflags)
{
	const int	i0 = ambndx(hp,i,j);
	double		rdot_cp, dot_e, dot_er, rdot_r, rdot_r1, J2;
	int		i1, ii;

	ftp->valid = 0;			/* check if we can skip this edge */
	ii = adjacent_trifl[dbit];
	if ((vflags[i0] & ii) == ii)	/* cancels if vertex used as value */
		return;
	i1 = adjacent_verti(hp, i, j, dbit);
	ii = adjacent_trifl[VDB_OPP(dbit)];
	if ((vflags[i1] & ii) == ii)	/* on either end (for both triangles) */
		return;
					/* else go ahead with calculation */
	VSUB(ftp->r_i, hp->sa[i0].p, hp->rp->rop);
	VSUB(ftp->r_i1, hp->sa[i1].p, hp->rp->rop);
	VSUB(ftp->e_i, hp->sa[i1].p, hp->sa[i0].p);
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
	ftp->valid++;
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

	if (!ftp->valid) {		/* preemptive test */
		memset(hess, 0, sizeof(FVECT)*3);
		return;
	}
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

	if (!ftp->valid) {		/* preemptive test */
		memset(grad, 0, sizeof(FVECT));
		return;
	}
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
static int
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
			((evalue[1] = fabs(evalue[1])) <= FTINY*FTINY) )
		error(INTERNAL, "bad eigenvalue calculation");

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
	uby8		*vflags;
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
					/* get vertex position flags */
	vflags = vertex_flags(hp);
					/* compute first row of edges */
	for (j = 0; j < hp->ns-1; j++) {
		comp_fftri(&fftr, hp, 0, j, VDB_X, vflags);
		if (hessrow != NULL)
			comp_hessian(hessrow[j], &fftr, hp->rp->ron);
		if (gradrow != NULL)
			comp_gradient(gradrow[j], &fftr, hp->rp->ron);
	}
					/* sum each row of triangles */
	for (i = 0; i < hp->ns-1; i++) {
	    FVECT	hesscol[3];	/* compute first vertical edge */
	    FVECT	gradcol;
	    comp_fftri(&fftr, hp, i, 0, VDB_Y, vflags);
	    if (hessrow != NULL)
		comp_hessian(hesscol, &fftr, hp->rp->ron);
	    if (gradrow != NULL)
		comp_gradient(gradcol, &fftr, hp->rp->ron);
	    for (j = 0; j < hp->ns-1; j++) {
		FVECT	hessdia[3];	/* compute triangle contributions */
		FVECT	graddia;
		double	backg;
		backg = back_ambval(hp, i, j, VDB_X, VDB_Y, vflags);
					/* diagonal (inner) edge */
		comp_fftri(&fftr, hp, i, j+1, VDB_xY, vflags);
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
		comp_fftri(&fftr, hp, i+1, j+1, VDB_x, vflags);
		if (hessrow != NULL)
		    comp_hessian(hessrow[j], &fftr, hp->rp->ron);
		if (gradrow != NULL)
		    comp_gradient(gradrow[j], &fftr, hp->rp->ron);
					/* new column edge & paired triangle */
		backg = back_ambval(hp, i+1, j+1, VDB_x, VDB_y, vflags);
		comp_fftri(&fftr, hp, i, j+1, VDB_Y, vflags);
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
	free(vflags);
	
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
	const double	ang_res = 0.5*PI/(hp->ns-1);
	const double	ang_step = ang_res/((int)(16/PI*ang_res) + (1+FTINY));
	uint32		flgs = 0;
	int		i, j;
					/* circle around perimeter */
	for (i = 0; i < hp->ns; i++)
	    for (j = 0; j < hp->ns; j += !i|(i==hp->ns-1) ? 1 : hp->ns-1) {
		AMBSAMP	*ap = &ambsam(hp,i,j);
		FVECT	vec;
		double	u, v;
		double	ang, a1;
		int	abp;
		if ((ap->d <= FTINY) | (ap->d >= max_d))
			continue;	/* too far or too near */
		VSUB(vec, ap->p, hp->rp->rop);
		u = DOT(vec, uv[0]) * ap->d;
		v = DOT(vec, uv[1]) * ap->d;
		if ((r0*r0*u*u + r1*r1*v*v) * ap->d*ap->d <= 1.0)
			continue;	/* occluder outside ellipse */
		ang = atan2a(v, u);	/* else set direction flags */
		for (a1 = ang-.5*ang_res; a1 <= ang+.5*ang_res; a1 += ang_step)
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
	AMBHEMI	*hp = inithemi(rcol, r, wt);
	int	cnt;
	FVECT	my_uv[2];
	double	d, K, acol[3];
	AMBSAMP	*ap;
	int	i, j;
					/* check/initialize */
	if (hp == NULL)
		return(0);
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
					/* sample the hemisphere */
	acol[0] = acol[1] = acol[2] = 0.0;
	cnt = 0;
	for (i = hp->ns; i--; )
		for (j = hp->ns; j--; )
			if ((ap = ambsample(hp, i, j)) != NULL) {
				addcolor(acol, ap->v);
				++cnt;
			}
	if (!cnt) {
		setcolor(rcol, 0.0, 0.0, 0.0);
		free(hp);
		return(0);		/* no valid samples */
	}
	if (cnt < hp->ns*hp->ns) {	/* incomplete sampling? */
		copycolor(rcol, acol);
		free(hp);
		return(-1);		/* return value w/o Hessian */
	}
	cnt = ambssamp*wt + 0.5;	/* perform super-sampling? */
	if (cnt > 8)
		ambsupersamp(acol, hp, cnt);
	copycolor(rcol, acol);		/* final indirect irradiance/PI */
	if ((ra == NULL) & (pg == NULL) & (dg == NULL)) {
		free(hp);
		return(-1);		/* no radius or gradient calc. */
	}
	if ((d = bright(acol)) > FTINY) {	/* normalize Y values */
		d = 0.99*(hp->ns*hp->ns)/d;
		K = 0.01;
	} else {			/* or fall back on geometric Hessian */
		K = 1.0;
		pg = NULL;
		dg = NULL;
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
		ra[0] *= d = 1.0/sqrt(sqrt(wt));
		if ((ra[1] *= d) > 2.0*ra[0])
			ra[1] = 2.0*ra[0];
		if (ra[1] > maxarad) {
			ra[1] = maxarad;
			if (ra[0] > maxarad)
				ra[0] = maxarad;
		}
		if (crlp != NULL)	/* flag encroached directions */
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


#else /* ! NEWAMB */


void
inithemi(			/* initialize sampling hemisphere */
	AMBHEMI  *hp,
	COLOR ac,
	RAY  *r,
	double  wt
)
{
	double	d;
	int  i;
					/* set number of divisions */
	if (ambacc <= FTINY &&
			wt > (d = 0.8*intens(ac)*r->rweight/(ambdiv*minweight)))
		wt = d;			/* avoid ray termination */
	hp->nt = sqrt(ambdiv * wt / PI) + 0.5;
	i = ambacc > FTINY ? 3 : 1;	/* minimum number of samples */
	if (hp->nt < i)
		hp->nt = i;
	hp->np = PI * hp->nt + 0.5;
					/* set number of super-samples */
	hp->ns = ambssamp * wt + 0.5;
					/* assign coefficient */
	copycolor(hp->acoef, ac);
	d = 1.0/(hp->nt*hp->np);
	scalecolor(hp->acoef, d);
					/* make axes */
	VCOPY(hp->uz, r->ron);
	hp->uy[0] = hp->uy[1] = hp->uy[2] = 0.0;
	for (i = 0; i < 3; i++)
		if (hp->uz[i] < 0.6 && hp->uz[i] > -0.6)
			break;
	if (i >= 3)
		error(CONSISTENCY, "bad ray direction in inithemi");
	hp->uy[i] = 1.0;
	fcross(hp->ux, hp->uy, hp->uz);
	normalize(hp->ux);
	fcross(hp->uy, hp->uz, hp->ux);
}


int
divsample(				/* sample a division */
	AMBSAMP  *dp,
	AMBHEMI  *h,
	RAY  *r
)
{
	RAY  ar;
	int  hlist[3];
	double  spt[2];
	double  xd, yd, zd;
	double  b2;
	double  phi;
	int  i;
					/* ambient coefficient for weight */
	if (ambacc > FTINY)
		setcolor(ar.rcoef, AVGREFL, AVGREFL, AVGREFL);
	else
		copycolor(ar.rcoef, h->acoef);
	if (rayorigin(&ar, AMBIENT, r, ar.rcoef) < 0)
		return(-1);
	if (ambacc > FTINY) {
		multcolor(ar.rcoef, h->acoef);
		scalecolor(ar.rcoef, 1./AVGREFL);
	}
	hlist[0] = r->rno;
	hlist[1] = dp->t;
	hlist[2] = dp->p;
	multisamp(spt, 2, urand(ilhash(hlist,3)+dp->n));
	zd = sqrt((dp->t + spt[0])/h->nt);
	phi = 2.0*PI * (dp->p + spt[1])/h->np;
	xd = tcos(phi) * zd;
	yd = tsin(phi) * zd;
	zd = sqrt(1.0 - zd*zd);
	for (i = 0; i < 3; i++)
		ar.rdir[i] =	xd*h->ux[i] +
				yd*h->uy[i] +
				zd*h->uz[i];
	checknorm(ar.rdir);
	dimlist[ndims++] = dp->t*h->np + dp->p + 90171;
	rayvalue(&ar);
	ndims--;
	multcolor(ar.rcol, ar.rcoef);	/* apply coefficient */
	addcolor(dp->v, ar.rcol);
					/* use rt to improve gradient calc */
	if (ar.rt > FTINY && ar.rt < FHUGE)
		dp->r += 1.0/ar.rt;
					/* (re)initialize error */
	if (dp->n++) {
		b2 = bright(dp->v)/dp->n - bright(ar.rcol);
		b2 = b2*b2 + dp->k*((dp->n-1)*(dp->n-1));
		dp->k = b2/(dp->n*dp->n);
	} else
		dp->k = 0.0;
	return(0);
}


static int
ambcmp(					/* decreasing order */
	const void *p1,
	const void *p2
)
{
	const AMBSAMP	*d1 = (const AMBSAMP *)p1;
	const AMBSAMP	*d2 = (const AMBSAMP *)p2;

	if (d1->k < d2->k)
		return(1);
	if (d1->k > d2->k)
		return(-1);
	return(0);
}


static int
ambnorm(				/* standard order */
	const void *p1,
	const void *p2
)
{
	const AMBSAMP	*d1 = (const AMBSAMP *)p1;
	const AMBSAMP	*d2 = (const AMBSAMP *)p2;
	int	c;

	if ( (c = d1->t - d2->t) )
		return(c);
	return(d1->p - d2->p);
}


double
doambient(				/* compute ambient component */
	COLOR  rcol,
	RAY  *r,
	double  wt,
	FVECT  pg,
	FVECT  dg
)
{
	double  b, d=0;
	AMBHEMI  hemi;
	AMBSAMP  *div;
	AMBSAMP  dnew;
	double  acol[3];
	AMBSAMP  *dp;
	double  arad;
	int  divcnt;
	int  i, j;
					/* initialize hemisphere */
	inithemi(&hemi, rcol, r, wt);
	divcnt = hemi.nt * hemi.np;
					/* initialize */
	if (pg != NULL)
		pg[0] = pg[1] = pg[2] = 0.0;
	if (dg != NULL)
		dg[0] = dg[1] = dg[2] = 0.0;
	setcolor(rcol, 0.0, 0.0, 0.0);
	if (divcnt == 0)
		return(0.0);
					/* allocate super-samples */
	if (hemi.ns > 0 || pg != NULL || dg != NULL) {
		div = (AMBSAMP *)malloc(divcnt*sizeof(AMBSAMP));
		if (div == NULL)
			error(SYSTEM, "out of memory in doambient");
	} else
		div = NULL;
					/* sample the divisions */
	arad = 0.0;
	acol[0] = acol[1] = acol[2] = 0.0;
	if ((dp = div) == NULL)
		dp = &dnew;
	divcnt = 0;
	for (i = 0; i < hemi.nt; i++)
		for (j = 0; j < hemi.np; j++) {
			dp->t = i; dp->p = j;
			setcolor(dp->v, 0.0, 0.0, 0.0);
			dp->r = 0.0;
			dp->n = 0;
			if (divsample(dp, &hemi, r) < 0) {
				if (div != NULL)
					dp++;
				continue;
			}
			arad += dp->r;
			divcnt++;
			if (div != NULL)
				dp++;
			else
				addcolor(acol, dp->v);
		}
	if (!divcnt) {
		if (div != NULL)
			free((void *)div);
		return(0.0);		/* no samples taken */
	}
	if (divcnt < hemi.nt*hemi.np) {
		pg = dg = NULL;		/* incomplete sampling */
		hemi.ns = 0;
	} else if (arad > FTINY && divcnt/arad < minarad) {
		hemi.ns = 0;		/* close enough */
	} else if (hemi.ns > 0) {	/* else perform super-sampling? */
		comperrs(div, &hemi);			/* compute errors */
		qsort(div, divcnt, sizeof(AMBSAMP), ambcmp);	/* sort divs */
						/* super-sample */
		for (i = hemi.ns; i > 0; i--) {
			dnew = *div;
			if (divsample(&dnew, &hemi, r) < 0) {
				dp++;
				continue;
			}
			dp = div;		/* reinsert */
			j = divcnt < i ? divcnt : i;
			while (--j > 0 && dnew.k < dp[1].k) {
				*dp = *(dp+1);
				dp++;
			}
			*dp = dnew;
		}
		if (pg != NULL || dg != NULL)	/* restore order */
			qsort(div, divcnt, sizeof(AMBSAMP), ambnorm);
	}
					/* compute returned values */
	if (div != NULL) {
		arad = 0.0;		/* note: divcnt may be < nt*np */
		for (i = hemi.nt*hemi.np, dp = div; i-- > 0; dp++) {
			arad += dp->r;
			if (dp->n > 1) {
				b = 1.0/dp->n;
				scalecolor(dp->v, b);
				dp->r *= b;
				dp->n = 1;
			}
			addcolor(acol, dp->v);
		}
		b = bright(acol);
		if (b > FTINY) {
			b = 1.0/b;	/* compute & normalize gradient(s) */
			if (pg != NULL) {
				posgradient(pg, div, &hemi);
				for (i = 0; i < 3; i++)
					pg[i] *= b;
			}
			if (dg != NULL) {
				dirgradient(dg, div, &hemi);
				for (i = 0; i < 3; i++)
					dg[i] *= b;
			}
		}
		free((void *)div);
	}
	copycolor(rcol, acol);
	if (arad <= FTINY)
		arad = maxarad;
	else
		arad = (divcnt+hemi.ns)/arad;
	if (pg != NULL) {		/* reduce radius if gradient large */
		d = DOT(pg,pg);
		if (d*arad*arad > 1.0)
			arad = 1.0/sqrt(d);
	}
	if (arad < minarad) {
		arad = minarad;
		if (pg != NULL && d*arad*arad > 1.0) {	/* cap gradient */
			d = 1.0/arad/sqrt(d);
			for (i = 0; i < 3; i++)
				pg[i] *= d;
		}
	}
	if ((arad /= sqrt(wt)) > maxarad)
		arad = maxarad;
	return(arad);
}


void
comperrs(			/* compute initial error estimates */
	AMBSAMP  *da,	/* assumes standard ordering */
	AMBHEMI  *hp
)
{
	double  b, b2;
	int  i, j;
	AMBSAMP  *dp;
				/* sum differences from neighbors */
	dp = da;
	for (i = 0; i < hp->nt; i++)
		for (j = 0; j < hp->np; j++) {
#ifdef  DEBUG
			if (dp->t != i || dp->p != j)
				error(CONSISTENCY,
					"division order in comperrs");
#endif
			b = bright(dp[0].v);
			if (i > 0) {		/* from above */
				b2 = bright(dp[-hp->np].v) - b;
				b2 *= b2 * 0.25;
				dp[0].k += b2;
				dp[-hp->np].k += b2;
			}
			if (j > 0) {		/* from behind */
				b2 = bright(dp[-1].v) - b;
				b2 *= b2 * 0.25;
				dp[0].k += b2;
				dp[-1].k += b2;
			} else {		/* around */
				b2 = bright(dp[hp->np-1].v) - b;
				b2 *= b2 * 0.25;
				dp[0].k += b2;
				dp[hp->np-1].k += b2;
			}
			dp++;
		}
				/* divide by number of neighbors */
	dp = da;
	for (j = 0; j < hp->np; j++)		/* top row */
		(dp++)->k *= 1.0/3.0;
	if (hp->nt < 2)
		return;
	for (i = 1; i < hp->nt-1; i++)		/* central region */
		for (j = 0; j < hp->np; j++)
			(dp++)->k *= 0.25;
	for (j = 0; j < hp->np; j++)		/* bottom row */
		(dp++)->k *= 1.0/3.0;
}


void
posgradient(					/* compute position gradient */
	FVECT  gv,
	AMBSAMP  *da,			/* assumes standard ordering */
	AMBHEMI  *hp
)
{
	int  i, j;
	double  nextsine, lastsine, b, d;
	double  mag0, mag1;
	double  phi, cosp, sinp, xd, yd;
	AMBSAMP  *dp;

	xd = yd = 0.0;
	for (j = 0; j < hp->np; j++) {
		dp = da + j;
		mag0 = mag1 = 0.0;
		lastsine = 0.0;
		for (i = 0; i < hp->nt; i++) {
#ifdef  DEBUG
			if (dp->t != i || dp->p != j)
				error(CONSISTENCY,
					"division order in posgradient");
#endif
			b = bright(dp->v);
			if (i > 0) {
				d = dp[-hp->np].r;
				if (dp[0].r > d) d = dp[0].r;
							/* sin(t)*cos(t)^2 */
				d *= lastsine * (1.0 - (double)i/hp->nt);
				mag0 += d*(b - bright(dp[-hp->np].v));
			}
			nextsine = sqrt((double)(i+1)/hp->nt);
			if (j > 0) {
				d = dp[-1].r;
				if (dp[0].r > d) d = dp[0].r;
				mag1 += d * (nextsine - lastsine) *
						(b - bright(dp[-1].v));
			} else {
				d = dp[hp->np-1].r;
				if (dp[0].r > d) d = dp[0].r;
				mag1 += d * (nextsine - lastsine) *
						(b - bright(dp[hp->np-1].v));
			}
			dp += hp->np;
			lastsine = nextsine;
		}
		mag0 *= 2.0*PI / hp->np;
		phi = 2.0*PI * (double)j/hp->np;
		cosp = tcos(phi); sinp = tsin(phi);
		xd += mag0*cosp - mag1*sinp;
		yd += mag0*sinp + mag1*cosp;
	}
	for (i = 0; i < 3; i++)
		gv[i] = (xd*hp->ux[i] + yd*hp->uy[i])*(hp->nt*hp->np)/PI;
}


void
dirgradient(					/* compute direction gradient */
	FVECT  gv,
	AMBSAMP  *da,			/* assumes standard ordering */
	AMBHEMI  *hp
)
{
	int  i, j;
	double  mag;
	double  phi, xd, yd;
	AMBSAMP  *dp;

	xd = yd = 0.0;
	for (j = 0; j < hp->np; j++) {
		dp = da + j;
		mag = 0.0;
		for (i = 0; i < hp->nt; i++) {
#ifdef  DEBUG
			if (dp->t != i || dp->p != j)
				error(CONSISTENCY,
					"division order in dirgradient");
#endif
							/* tan(t) */
			mag += bright(dp->v)/sqrt(hp->nt/(i+.5) - 1.0);
			dp += hp->np;
		}
		phi = 2.0*PI * (j+.5)/hp->np + PI/2.0;
		xd += mag * tcos(phi);
		yd += mag * tsin(phi);
	}
	for (i = 0; i < 3; i++)
		gv[i] = xd*hp->ux[i] + yd*hp->uy[i];
}

#endif	/* ! NEWAMB */
