#ifndef lint
static const char	RCSid[] = "$Id: pcond3.c,v 3.12 2003/04/23 00:52:34 greg Exp $";
#endif
/*
 * Routines for computing and applying brightness mapping.
 */

#include "pcond.h"


#define CVRATIO		0.025		/* fraction of samples allowed > env. */

#define LN_10		2.30258509299404568402
#define exp10(x)	exp(LN_10*(x))

double	modhist[HISTRES];		/* modified histogram */
double	mhistot;			/* modified histogram total */
double	cumf[HISTRES+1];		/* cumulative distribution function */


getfixations(fp)		/* load fixation history list */
FILE	*fp;
{
#define	FIXHUNK		128
	RESOLU	fvres;
	int	pos[2];
	register int	px, py, i;
				/* initialize our resolution struct */
	if ((fvres.rt=inpres.rt)&YMAJOR) {
		fvres.xr = fvxr;
		fvres.yr = fvyr;
	} else {
		fvres.xr = fvyr;
		fvres.yr = fvxr;
	}
				/* read each picture position */
	while (fscanf(fp, "%d %d", &pos[0], &pos[1]) == 2) {
				/* convert to closest index in foveal image */
		loc2pix(pos, &fvres,
				(pos[0]+.5)/inpres.xr, (pos[1]+.5)/inpres.yr);
				/* include nine neighborhood samples */
		for (px = pos[0]-1; px <= pos[0]+1; px++) {
			if (px < 0 || px >= fvxr)
				continue;
			for (py = pos[1]-1; py <= pos[1]+1; py++) {
				if (py < 0 || py >= fvyr)
					continue;
				for (i = nfixations; i-- > 0; )
					if (fixlst[i][0] == px &&
							fixlst[i][1] == py)
						break;
				if (i >= 0)
					continue;	/* already there */
				if (nfixations % FIXHUNK == 0) {
					if (nfixations)
						fixlst = (short (*)[2])
							realloc((void *)fixlst,
							(nfixations+FIXHUNK)*
							2*sizeof(short));
					else
						fixlst = (short (*)[2])malloc(
							FIXHUNK*2*sizeof(short)
							);
					if (fixlst == NULL)
						syserror("malloc");
				}
				fixlst[nfixations][0] = px;
				fixlst[nfixations][1] = py;
				nfixations++;
			}
		}
	}
	if (!feof(fp)) {
		fprintf(stderr, "%s: format error reading fixation data\n",
				progname);
		exit(1);
	}
#undef	FIXHUNK
}


gethisto(fp)			/* load precomputed luminance histogram */
FILE	*fp;
{
	double	histo[MAXPREHIST];
	double	histart, histep;
	double	l, b, lastb, w;
	int	n;
	register int	i;
					/* load data */
	for (i = 0; i < MAXPREHIST &&
			fscanf(fp, "%lf %lf", &b, &histo[i]) == 2; i++) {
		if (i > 1 && fabs(b - lastb - histep) > .001) {
			fprintf(stderr,
				"%s: uneven step size in histogram data\n",
					progname);
			exit(1);
		}
		if (i == 1)
			if ((histep = b - (histart = lastb)) <= FTINY) {
				fprintf(stderr,
					"%s: illegal step in histogram data\n",
						progname);
				exit(1);
			}
		lastb = b;
	}
	if (i < 2 || !feof(fp)) {
		fprintf(stderr,
		"%s: format/length error loading histogram (log10L %f at %d)\n",
				progname, b, i);
		exit(1);
	}
	n = i;
	histart *= LN_10;
	histep *= LN_10;
					/* find extrema */
	for (i = 0; i < n && histo[i] <= FTINY; i++)
		;
	bwmin = histart + (i-.001)*histep;
	for (i = n; i-- && histo[i] <= FTINY; )
		;
	bwmax = histart + (i+1.001)*histep;
	if (bwmax > Bl(LMAX))
		bwmax = Bl(LMAX);
	if (bwmin < Bl(LMIN))
		bwmin = Bl(LMIN);
	else				/* duplicate bottom bin */
		bwmin = bwmax - (bwmax-bwmin)*HISTRES/(HISTRES-1);
					/* convert histogram */
	bwavg = 0.; histot = 0.;
	for (i = 0; i < HISTRES; i++)
		bwhist[i] = 0.;
	for (i = 0, b = histart; i < n; i++, b += histep) {
		if (b < bwmin+FTINY) continue;
		if (b >= bwmax-FTINY) break;
		w = histo[i];
		bwavg += w*b;
		bwhist[bwhi(b)] += w;
		histot += w;
	}
	bwavg /= histot;
	if (bwmin > Bl(LMIN)+FTINY) {	/* add false samples at bottom */
		bwhist[1] *= 0.5;
		bwhist[0] += bwhist[1];
	}
}


double
centprob(x, y)			/* center-weighting probability function */
int	x, y;
{
	double	xr, yr, p;
				/* paraboloid, 0 at 90 degrees from center */
	xr = (x - .5*(fvxr-1))/90.;	/* 180 degree fisheye has fv?r == 90 */
	yr = (y - .5*(fvyr-1))/90.;
	p = 1. - xr*xr - yr*yr;
	return(p < 0. ? 0. : p);
}


comphist()			/* create foveal sampling histogram */
{
	double	l, b, w, lwmin, lwmax;
	register int	x, y;
					/* check for precalculated histogram */
	if (what2do&DO_PREHIST)
		return;
	lwmin = 1e10;			/* find extrema */
	lwmax = 0.;
	for (y = 0; y < fvyr; y++)
		for (x = 0; x < fvxr; x++) {
			l = plum(fovscan(y)[x]);
			if (l < lwmin) lwmin = l;
			if (l > lwmax) lwmax = l;
		}
	lwmax *= 1.01;
	if (lwmax > LMAX)
		lwmax = LMAX;
	bwmax = Bl(lwmax);
	if (lwmin < LMIN) {
		lwmin = LMIN;
		bwmin = Bl(LMIN);
	} else {			/* duplicate bottom bin */
		bwmin = bwmax - (bwmax-Bl(lwmin))*HISTRES/(HISTRES-1);
		lwmin = Lb(bwmin);
	}
					/* (re)compute histogram */
	bwavg = 0.;
	histot = 0.;
	for (x = 0; x < HISTRES; x++)
		bwhist[x] = 0.;
					/* global average */
	if (!(what2do&DO_FIXHIST) || fixfrac < 1.-FTINY)
		for (y = 0; y < fvyr; y++)
			for (x = 0; x < fvxr; x++) {
				l = plum(fovscan(y)[x]);
				if (l < lwmin+FTINY) continue;
				if (l >= lwmax-FTINY) continue;
				b = Bl(l);
				w = what2do&DO_CWEIGHT ? centprob(x,y) : 1.;
				bwavg += w*b;
				bwhist[bwhi(b)] += w;
				histot += w;
			}
					/* average fixation points */
	if (what2do&DO_FIXHIST && nfixations > 0) {
		if (histot > FTINY)
			w = fixfrac/(1.-fixfrac)*histot/nfixations;
		else
			w = 1.;
		for (x = 0; x < nfixations; x++) {
			l = plum(fovscan(fixlst[x][1])[fixlst[x][0]]);
			if (l < lwmin+FTINY) continue;
			if (l >= lwmax-FTINY) continue;
			b = Bl(l);
			bwavg += w*b;
			bwhist[bwhi(b)] += w;
			histot += w;
		}
	}
	bwavg /= histot;
	if (lwmin > LMIN+FTINY) {	/* add false samples at bottom */
		bwhist[1] *= 0.5;
		bwhist[0] += bwhist[1];
	}
}


mkcumf()			/* make cumulative distribution function */
{
	register int	i;
	register double	sum;

	mhistot = 0.;		/* compute modified total */
	for (i = 0; i < HISTRES; i++)
		mhistot += modhist[i];

	sum = 0.;		/* compute cumulative function */
	for (i = 0; i < HISTRES; i++) {
		cumf[i] = sum/mhistot;
		sum += modhist[i];
	}
	cumf[HISTRES] = 1.;
}


double
cf(b)				/* return cumulative function at b */
double	b;
{
	double	x;
	register int	i;

	i = x = HISTRES*(b - bwmin)/(bwmax - bwmin);
	x -= (double)i;
	return(cumf[i]*(1.-x) + cumf[i+1]*x);
}


double
BLw(Lw)				/* map world luminance to display brightness */
double	Lw;
{
	double	b;

	if (Lw <= LMIN || (b = Bl(Lw)) <= bwmin+FTINY)
		return(Bldmin);
	if (b >= bwmax-FTINY)
		return(Bldmax);
	return(Bldmin + cf(b)*(Bldmax-Bldmin));
}


double
htcontrs(La)		/* human threshold contrast sensitivity, dL(La) */
double	La;
{
	double	l10La, l10dL;
				/* formula taken from Ferwerda et al. [SG96] */
	if (La < 1.148e-4)
		return(1.38e-3);
	l10La = log10(La);
	if (l10La < -1.44)		/* rod response regime */
		l10dL = pow(.405*l10La + 1.6, 2.18) - 2.86;
	else if (l10La < -.0184)
		l10dL = l10La - .395;
	else if (l10La < 1.9)		/* cone response regime */
		l10dL = pow(.249*l10La + .65, 2.7) - .72;
	else
		l10dL = l10La - 1.255;

	return(exp10(l10dL));
}


double
clampf(Lw)			/* histogram clamping function */
double	Lw;
{
	double	bLw, ratio;

	bLw = BLw(Lw);		/* apply current brightness map */
	ratio = what2do&DO_HSENS ? htcontrs(Lb(bLw))/htcontrs(Lw) : Lb(bLw)/Lw;
	return(ratio/(Lb1(bLw)*(Bldmax-Bldmin)*Bl1(Lw)));
}

double
crfactor(Lw)			/* contrast reduction factor */
double	Lw;
{
	int	i = HISTRES*(Bl(Lw) - bwmin)/(bwmax - bwmin);
	double	bLw, ratio, Tdb;

	if (i <= 0)
		return(1.0);
	if (i >= HISTRES)
		return(1.0);
	bLw = BLw(Lw);
	ratio = what2do&DO_HSENS ? htcontrs(Lb(bLw))/htcontrs(Lw) : Lb(bLw)/Lw;
	Tdb = mhistot * (bwmax - bwmin) / HISTRES;
	return(modhist[i]*Lb1(bLw)*(Bldmax-Bldmin)*Bl1(Lw)/(Tdb*ratio));
}


#if ADJ_VEIL
mkcrfimage()			/* compute contrast reduction factor image */
{
	int	i;
	float	*crfptr;
	COLOR	*fovptr;

	if (crfimg == NULL)
		crfimg = (float *)malloc(fvxr*fvyr*sizeof(float));
	if (crfimg == NULL)
		syserror("malloc");
	crfptr = crfimg;
	fovptr = fovimg;
	for (i = fvxr*fvyr; i--; crfptr++, fovptr++)
		crfptr[0] = crfactor(plum(fovptr[0]));
}
#endif


int
mkbrmap()			/* make dynamic range map */
{
	double	Tdb, b, s;
	double	ceiling, trimmings;
	register int	i;
					/* copy initial histogram */
	bcopy((char *)bwhist, (char *)modhist, sizeof(modhist));
	s = (bwmax - bwmin)/HISTRES;	/* s is delta b */
					/* loop until satisfactory */
	do {
		mkcumf();			/* sync brightness mapping */
		if (mhistot <= histot*CVRATIO)
			return(-1);		/* no compression needed! */
		Tdb = mhistot * s;
		trimmings = 0.;			/* clip to envelope */
		for (i = 0, b = bwmin + .5*s; i < HISTRES; i++, b += s) {
			ceiling = Tdb*clampf(Lb(b));
			if (modhist[i] > ceiling) {
				trimmings += modhist[i] - ceiling;
				modhist[i] = ceiling;
			}
		}
	} while (trimmings > histot*CVRATIO);

#if ADJ_VEIL
	mkcrfimage();			/* contrast reduction image */
#endif

	return(0);			/* we got it */
}


scotscan(scan, xres)		/* apply scotopic color sensitivity loss */
COLOR	*scan;
int	xres;
{
	COLOR	ctmp;
	double	incolor, b, Lw;
	register int	i;

	for (i = 0; i < xres; i++) {
		Lw = plum(scan[i]);
		if (Lw >= TopMesopic)
			incolor = 1.;
		else if (Lw <= BotMesopic)
			incolor = 0.;
		else
			incolor = (Lw - BotMesopic) /
					(TopMesopic - BotMesopic);
		if (incolor < 1.-FTINY) {
			b = (1.-incolor)*slum(scan[i])*inpexp/SWNORM;
			if (lumf == rgblum) b /= WHTEFFICACY;
			setcolor(ctmp, b, b, b);
			if (incolor <= FTINY)
				setcolor(scan[i], 0., 0., 0.);
			else
				scalecolor(scan[i], incolor);
			addcolor(scan[i], ctmp);
		}
	}
}


mapscan(scan, xres)		/* apply tone mapping operator to scanline */
COLOR	*scan;
int	xres;
{
	double	mult, Lw, b;
	register int	x;

	for (x = 0; x < xres; x++) {
		Lw = plum(scan[x]);
		if (Lw < LMIN) {
			setcolor(scan[x], 0., 0., 0.);
			continue;
		}
		b = BLw(Lw);		/* apply brightness mapping */
		mult = (Lb(b) - ldmin)/(ldmax - ldmin)/(Lw*inpexp);
		if (lumf == rgblum) mult *= WHTEFFICACY;
		scalecolor(scan[x], mult);
	}
}


putmapping(fp)			/* put out mapping function */
FILE	*fp;
{
	double	b, s;
	register int	i;
	double	wlum, sf, dlum;

	sf = scalef*inpexp;
	if (lumf == cielum) sf *= WHTEFFICACY;
	s = (bwmax - bwmin)/HISTRES;
	for (i = 0, b = bwmin + .5*s; i < HISTRES; i++, b += s) {
		wlum = Lb(b);
		if (what2do&DO_LINEAR) {
			dlum = sf*wlum;
			if (dlum > ldmax) dlum = ldmax;
			else if (dlum < ldmin) dlum = ldmin;
			fprintf(fp, "%e %e\n", wlum, dlum);
		} else
			fprintf(fp, "%e %e\n", wlum, Lb(BLw(wlum)));
	}
}
