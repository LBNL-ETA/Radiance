#ifndef lint
static const char	RCSid[] = "$Id: ambcomp.c,v 2.26 2014/04/16 20:32:00 greg Exp $";
#endif
/*
 * Routines to compute "ambient" values using Monte Carlo
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include  "ray.h"
#include  "ambient.h"
#include  "random.h"

#ifdef NEWAMB

extern void		SDsquare2disk(double ds[2], double seedx, double seedy);

typedef struct {
	RAY	*rp;		/* originating ray sample */
	FVECT	ux, uy;		/* tangent axis directions */
	int	ns;		/* number of samples per axis */
	COLOR	acoef;		/* division contribution coefficient */
	struct s_ambsamp {
		COLOR	v;		/* hemisphere sample value */
		float	p[3];		/* intersection point */
	} sa[1];		/* sample array (extends struct) */
}  AMBHEMI;		/* ambient sample hemisphere */

#define ambsamp(h,i,j)	(h)->sa[(i)*(h)->ns + (j)]


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
	i = 1 + 4*(ambacc > FTINY);	/* minimum number of samples */
	if (n < i)
		n = i;
					/* allocate sampling array */
	hp = (AMBHEMI *)malloc(sizeof(AMBHEMI) +
				sizeof(struct s_ambsamp)*(n*n - 1));
	if (hp == NULL)
		return(NULL);
	hp->rp = r;
	hp->ns = n;
					/* assign coefficient */
	copycolor(hp->acoef, ac);
	d = 1.0/(n*n);
	scalecolor(hp->acoef, d);
					/* make tangent axes */
	hp->uy[0] = hp->uy[1] = hp->uy[2] = 0.0;
	for (i = 0; i < 3; i++)
		if (r->rn[i] < 0.6 && r->rn[i] > -0.6)
			break;
	if (i >= 3)
		error(CONSISTENCY, "bad ray direction in inithemi()");
	hp->uy[i] = 1.0;
	VCROSS(hp->ux, hp->uy, r->rn);
	normalize(hp->ux);
	VCROSS(hp->uy, r->rn, hp->ux);
					/* we're ready to sample */
	return(hp);
}


static int
ambsample(				/* sample an ambient direction */
	AMBHEMI	*hp,
	int	i,
	int	j,
)
{
	struct s_ambsamp	*ap = &ambsamp(hp,i,j);
	RAY			ar;
	int			hlist[3];
	double			spt[2], dz;
	int			ii;
					/* ambient coefficient for weight */
	if (ambacc > FTINY)
		setcolor(ar.rcoef, AVGREFL, AVGREFL, AVGREFL);
	else
		copycolor(ar.rcoef, hp->acoef);
	if (rayorigin(&ar, AMBIENT, hp->rp, ar.rcoef) < 0) {
		setcolor(ap->v, 0., 0., 0.);
		ap->r = 0.;
		return(0);		/* no sample taken */
	}
	if (ambacc > FTINY) {
		multcolor(ar.rcoef, hp->acoef);
		scalecolor(ar.rcoef, 1./AVGREFL);
	}
					/* generate hemispherical sample */
	SDsquare2disk(spt, (i+frandom())/hp->ns, (j+frandom())/hp->ns);
	zd = sqrt(1. - spt[0]*spt[0] - spt[1]*spt[1]);
	for (ii = 3; ii--; )
		ar.rdir[ii] =	spt[0]*hp->ux[ii] +
				spt[1]*hp->uy[ii] +
				zd*hp->rp->ron[ii];
	checknorm(ar.rdir);
	dimlist[ndims++] = i*hp->ns + j + 90171;
	rayvalue(&ar);			/* evaluate ray */
	ndims--;
	multcolor(ar.rcol, ar.rcoef);	/* apply coefficient */
	copycolor(ap->v, ar.rcol);
	if (ar.rt > 20.0*maxarad)	/* limit vertex distance */
		ar.rt = 20.0*maxarad;
	VSUM(ap->p, ar.rorg, ar.rdir, ar.rt);
	return(1);
}


static void
ambHessian(				/* anisotropic radii & pos. gradient */
	AMBHEMI	*hp,
	FVECT	uv[2],			/* returned */
	float	ra[2],			/* returned */
	float	pg[2]			/* returned */
)
{
	if (ra != NULL) {		/* compute Hessian-derived radii */
	} else {			/* else copy original tangent axes */
		VCOPY(uv[0], hp->ux);
		VCOPY(uv[1], hp->uy);
	}
	if (pg == NULL)			/* no position gradient requested? */
		return;
}

int
doambient(				/* compute ambient component */
	COLOR	rcol,			/* input/output color */
	RAY	*r,
	double	wt,
	FVECT	uv[2],			/* returned */
	float	ra[2],			/* returned */
	float	pg[2],			/* returned */
	float	dg[2]			/* returned */
)
{
	int			cnt = 0;
	FVECT			my_uv[2];
	AMBHEMI			*hp;
	double			d, acol[3];
	struct s_ambsamp	*ap;
	int			i, j;
					/* initialize */
	if ((hp = inithemi(rcol, r, wt)) == NULL)
		return(0);
	if (uv != NULL)
		memset(uv, 0, sizeof(FVECT)*2);
	if (ra != NULL)
		ra[0] = ra[1] = 0.0;
	if (pg != NULL)
		pg[0] = pg[1] = 0.0;
	if (dg != NULL)
		dg[0] = dg[1] = 0.0;
					/* sample the hemisphere */
	acol[0] = acol[1] = acol[2] = 0.0;
	for (i = hemi.ns; i--; )
		for (j = hemi.ns; j--; )
			if (ambsample(hp, i, j)) {
				ap = &ambsamp(hp,i,j);
				addcolor(acol, ap->v);
				++cnt;
			}
	if (!cnt) {
		setcolor(rcol, 0.0, 0.0, 0.0);
		free(hp);
		return(0);		/* no valid samples */
	}
	d = 1.0 / cnt;			/* final indirect irradiance/PI */
	acol[0] *= d; acol[1] *= d; acol[2] *= d;
	copycolor(rcol, acol);
	if (cnt < hp->ns*hp->ns ||	/* incomplete sampling? */
			(ra == NULL) & (pg == NULL) & (dg == NULL)) {
		free(hp);
		return(-1);		/* no radius or gradient calc. */
	}
	d = 0.01 * bright(rcol);	/* add in 1% before Hessian comp. */
	if (d < FTINY) d = FTINY;
	ap = hp->sa;			/* using Y channel from here on... */
	for (i = hp->ns*hp->ns; i--; ap++)
		colval(ap->v,CIEY) = bright(ap->v) + d;

	if (uv == NULL)			/* make sure we have axis pointers */
		uv = my_uv;
					/* compute radii & pos. gradient */
	ambHessian(hp, uv, ra, pg);
	if (dg != NULL)			/* compute direction gradient */
		ambdirgrad(hp, uv, dg);
	if (ra != NULL) {		/* adjust/clamp radii */
		d = pow(wt, -0.25);
		if ((ra[0] *= d) > maxarad)
			ra[0] = maxarad;
		if ((ra[1] *= d) > 2.0*ra[0])
			ra[1] = 2.0*ra[0];
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
