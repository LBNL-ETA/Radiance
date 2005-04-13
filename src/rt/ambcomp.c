#ifndef lint
static const char	RCSid[] = "$Id: ambcomp.c,v 2.13 2005/04/13 23:00:59 greg Exp $";
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


static int
ambcmp(d1, d2)				/* decreasing order */
AMBSAMP  *d1, *d2;
{
	if (d1->k < d2->k)
		return(1);
	if (d1->k > d2->k)
		return(-1);
	return(0);
}


static int
ambnorm(d1, d2)				/* standard order */
AMBSAMP  *d1, *d2;
{
	register int  c;

	if ( (c = d1->t - d2->t) )
		return(c);
	return(d1->p - d2->p);
}


int
divsample(dp, h, r)			/* sample a division */
register AMBSAMP  *dp;
AMBHEMI  *h;
RAY  *r;
{
	RAY  ar;
	int  hlist[3];
	double  spt[2];
	double  xd, yd, zd;
	double  b2;
	double  phi;
	register int  i;

	if (rayorigin(&ar, r, AMBIENT, AVGREFL) < 0)
		return(-1);
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
	dimlist[ndims++] = dp->t*h->np + dp->p + 90171;
	rayvalue(&ar);
	ndims--;
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


double
doambient(acol, r, wt, pg, dg)		/* compute ambient component */
COLOR  acol;
RAY  *r;
double  wt;
FVECT  pg, dg;
{
	double  b, d;
	AMBHEMI  hemi;
	AMBSAMP  *div;
	AMBSAMP  dnew;
	register AMBSAMP  *dp;
	double  arad;
	int  ndivs, ns;
	register int  i, j;
					/* initialize color */
	setcolor(acol, 0.0, 0.0, 0.0);
					/* initialize hemisphere */
	inithemi(&hemi, r, wt);
	ndivs = hemi.nt * hemi.np;
	if (ndivs == 0)
		return(0.0);
					/* set number of super-samples */
	ns = ambssamp * wt + 0.5;
	if (ns > 0 || pg != NULL || dg != NULL) {
		div = (AMBSAMP *)malloc(ndivs*sizeof(AMBSAMP));
		if (div == NULL)
			error(SYSTEM, "out of memory in doambient");
	} else
		div = NULL;
					/* sample the divisions */
	arad = 0.0;
	if ((dp = div) == NULL)
		dp = &dnew;
	for (i = 0; i < hemi.nt; i++)
		for (j = 0; j < hemi.np; j++) {
			dp->t = i; dp->p = j;
			setcolor(dp->v, 0.0, 0.0, 0.0);
			dp->r = 0.0;
			dp->n = 0;
			if (divsample(dp, &hemi, r) < 0)
				goto oopsy;
			arad += dp->r;
			if (div != NULL)
				dp++;
			else
				addcolor(acol, dp->v);
		}
	if (ns > 0 && arad > FTINY && ndivs/arad < minarad)
		ns = 0;			/* close enough */
	else if (ns > 0) {		/* else perform super-sampling */
		comperrs(div, &hemi);			/* compute errors */
		qsort(div, ndivs, sizeof(AMBSAMP), ambcmp);	/* sort divs */
						/* super-sample */
		for (i = ns; i > 0; i--) {
			dnew = *div;
			if (divsample(&dnew, &hemi, r) < 0)
				goto oopsy;
							/* reinsert */
			dp = div;
			j = ndivs < i ? ndivs : i;
			while (--j > 0 && dnew.k < dp[1].k) {
				*dp = *(dp+1);
				dp++;
			}
			*dp = dnew;
		}
		if (pg != NULL || dg != NULL)	/* restore order */
			qsort(div, ndivs, sizeof(AMBSAMP), ambnorm);
	}
					/* compute returned values */
	if (div != NULL) {
		arad = 0.0;
		for (i = ndivs, dp = div; i-- > 0; dp++) {
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
			b = ndivs/b;
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
		} else {
			if (pg != NULL)
				for (i = 0; i < 3; i++)
					pg[i] = 0.0;
			if (dg != NULL)
				for (i = 0; i < 3; i++)
					dg[i] = 0.0;
		}
		free((void *)div);
	}
	b = 1.0/ndivs;
	scalecolor(acol, b);
	if (arad <= FTINY)
		arad = maxarad;
	else
		arad = (ndivs+ns)/arad;
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
oopsy:
	if (div != NULL)
		free((void *)div);
	return(0.0);
}


void
inithemi(hp, r, wt)		/* initialize sampling hemisphere */
register AMBHEMI  *hp;
RAY  *r;
double  wt;
{
	register int  i;
					/* set number of divisions */
	hp->nt = sqrt(ambdiv * wt / PI) + 0.5;
	i = ambacc > FTINY ? 3 : 1;	/* minimum number of samples */
	if (hp->nt < i)
		hp->nt = i;
	hp->np = PI * hp->nt + 0.5;
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


void
comperrs(da, hp)		/* compute initial error estimates */
AMBSAMP  *da;		/* assumes standard ordering */
register AMBHEMI  *hp;
{
	double  b, b2;
	int  i, j;
	register AMBSAMP  *dp;
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
posgradient(gv, da, hp)				/* compute position gradient */
FVECT  gv;
AMBSAMP  *da;			/* assumes standard ordering */
register AMBHEMI  *hp;
{
	register int  i, j;
	double  nextsine, lastsine, b, d;
	double  mag0, mag1;
	double  phi, cosp, sinp, xd, yd;
	register AMBSAMP  *dp;

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
		gv[i] = (xd*hp->ux[i] + yd*hp->uy[i])/PI;
}


void
dirgradient(gv, da, hp)				/* compute direction gradient */
FVECT  gv;
AMBSAMP  *da;			/* assumes standard ordering */
register AMBHEMI  *hp;
{
	register int  i, j;
	double  mag;
	double  phi, xd, yd;
	register AMBSAMP  *dp;

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
		gv[i] = (xd*hp->ux[i] + yd*hp->uy[i])/(hp->nt*hp->np);
}
