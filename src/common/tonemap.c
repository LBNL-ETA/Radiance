#ifndef lint
static const char	RCSid[] = "$Id: tonemap.c,v 3.9 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Tone mapping functions.
 * See tonemap.h for detailed function descriptions.
 * Added von Kries white-balance calculations 10/01 (GW).
 *
 * Externals declared in tonemap.h
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include	<stdio.h>
#include	<math.h>
#include	"tmprivat.h"
#include	"tmerrmsg.h"

#define	exp10(x)	exp(M_LN10*(x))

struct tmStruct	*tmTop = NULL;		/* current tone mapping stack */

					/* our list of conversion packages */
struct tmPackage	*tmPkg[TM_MAXPKG];
int	tmNumPkgs = 0;			/* number of registered packages */

int	tmLastError;			/* last error incurred by library */
char	*tmLastFunction;		/* error-generating function name */


struct tmStruct *
tmInit(flags, monpri, gamval)		/* initialize new tone mapping */
int	flags;
RGBPRIMP	monpri;
double	gamval;
{
	COLORMAT	cmat;
	register struct tmStruct	*tmnew;
	register int	i;
						/* allocate structure */
	tmnew = (struct tmStruct *)malloc(sizeof(struct tmStruct));
	if (tmnew == NULL)
		return(NULL);

	tmnew->flags = flags & ~TM_F_UNIMPL;
						/* set monitor transform */
	if (monpri == NULL || monpri == stdprims || tmnew->flags & TM_F_BW) {
		tmnew->monpri = stdprims;
		tmnew->clf[RED] = rgb2xyzmat[1][0];
		tmnew->clf[GRN] = rgb2xyzmat[1][1];
		tmnew->clf[BLU] = rgb2xyzmat[1][2];
	} else {
		comprgb2xyzWBmat(cmat, tmnew->monpri=monpri);
		tmnew->clf[RED] = cmat[1][0];
		tmnew->clf[GRN] = cmat[1][1];
		tmnew->clf[BLU] = cmat[1][2];
	}
						/* set gamma value */
	if (gamval < MINGAM)
		tmnew->mongam = DEFGAM;
	else
		tmnew->mongam = gamval;
						/* set color divisors */
	for (i = 0; i < 3; i++)
		tmnew->cdiv[i] = 256.*pow(tmnew->clf[i], 1./tmnew->mongam);

						/* set input transform */
	tmnew->inppri = tmnew->monpri;
	tmnew->cmat[0][0] = tmnew->cmat[1][1] = tmnew->cmat[2][2] =
			tmnew->inpsf = WHTEFFICACY;
	tmnew->cmat[0][1] = tmnew->cmat[0][2] = tmnew->cmat[1][0] =
	tmnew->cmat[1][2] = tmnew->cmat[2][0] = tmnew->cmat[2][1] = 0.;
	tmnew->hbrmin = 10; tmnew->hbrmax = -10;
	tmnew->histo = NULL;
	tmnew->mbrmin = 10; tmnew->mbrmax = -10;
	tmnew->lumap = NULL;
						/* zero private data */
	for (i = TM_MAXPKG; i--; )
		tmnew->pd[i] = NULL;
						/* make tmnew current */
	tmnew->tmprev = tmTop;
	return(tmTop = tmnew);
}


int
tmSetSpace(pri, sf)		/* set input color space for conversions */
RGBPRIMP	pri;
double	sf;
{
	static char	funcName[] = "tmSetSpace";
	register int	i, j;
						/* error check */
	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (sf <= 1e-12)
		returnErr(TM_E_ILLEGAL);
						/* check if no change */
	if (pri == tmTop->inppri && FEQ(sf, tmTop->inpsf))
		returnOK;
	tmTop->inppri = pri;			/* let's set it */
	tmTop->inpsf = sf;

	if (tmTop->flags & TM_F_BW) {		/* color doesn't matter */
		tmTop->monpri = tmTop->inppri;		/* eliminate xform */
		if (tmTop->inppri == TM_XYZPRIM) {
			tmTop->clf[CIEX] = tmTop->clf[CIEZ] = 0.;
			tmTop->clf[CIEY] = 1.;
		} else {
			comprgb2xyzWBmat(tmTop->cmat, tmTop->monpri);
			tmTop->clf[RED] = tmTop->cmat[1][0];
			tmTop->clf[GRN] = tmTop->cmat[1][1];
			tmTop->clf[BLU] = tmTop->cmat[1][2];
		}
		tmTop->cmat[0][0] = tmTop->cmat[1][1] = tmTop->cmat[2][2] =
				tmTop->inpsf;
		tmTop->cmat[0][1] = tmTop->cmat[0][2] = tmTop->cmat[1][0] =
		tmTop->cmat[1][2] = tmTop->cmat[2][0] = tmTop->cmat[2][1] = 0.;

	} else if (tmTop->inppri == TM_XYZPRIM)	/* input is XYZ */
		compxyz2rgbWBmat(tmTop->cmat, tmTop->monpri);

	else {					/* input is RGB */
		if (tmTop->inppri != tmTop->monpri &&
				PRIMEQ(tmTop->inppri, tmTop->monpri))
			tmTop->inppri = tmTop->monpri;	/* no xform */
		comprgb2rgbWBmat(tmTop->cmat, tmTop->inppri, tmTop->monpri);
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			tmTop->cmat[i][j] *= tmTop->inpsf;
						/* set color divisors */
	for (i = 0; i < 3; i++)
		if (tmTop->clf[i] > .001)
			tmTop->cdiv[i] =
				256.*pow(tmTop->clf[i], 1./tmTop->mongam);
		else
			tmTop->cdiv[i] = 1;
						/* notify packages */
	for (i = tmNumPkgs; i--; )
		if (tmTop->pd[i] != NULL && tmPkg[i]->NewSpace != NULL)
			(*tmPkg[i]->NewSpace)(tmTop);
	returnOK;
}


void
tmClearHisto()			/* clear current histogram */
{
	if (tmTop == NULL || tmTop->histo == NULL)
		return;
	free((MEM_PTR)tmTop->histo);
	tmTop->histo = NULL;
}


int
tmCvColors(ls, cs, scan, len)		/* convert float colors */
TMbright	*ls;
BYTE	*cs;
COLOR	*scan;
int	len;
{
	static char	funcName[] = "tmCvColors";
	static COLOR	csmall = {1e-6, 1e-6, 1e-6};
	COLOR	cmon;
	double	lum, slum;
	register double	d;
	register int	i;

	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | scan == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
	for (i = len; i--; ) {
		if (tmNeedMatrix(tmTop))		/* get monitor RGB */
			colortrans(cmon, tmTop->cmat, scan[i]);
		else {
			cmon[RED] = tmTop->inpsf*scan[i][RED];
			cmon[GRN] = tmTop->inpsf*scan[i][GRN];
			cmon[BLU] = tmTop->inpsf*scan[i][BLU];
		}
							/* world luminance */
		lum =	tmTop->clf[RED]*cmon[RED] +
			tmTop->clf[GRN]*cmon[GRN] +
			tmTop->clf[BLU]*cmon[BLU] ;
							/* check range */
		if (clipgamut(cmon, lum, CGAMUT_LOWER, csmall, cwhite))
			lum =	tmTop->clf[RED]*cmon[RED] +
				tmTop->clf[GRN]*cmon[GRN] +
				tmTop->clf[BLU]*cmon[BLU] ;
		if (lum < MINLUM) {
			ls[i] = MINBRT-1;		/* bogus value */
			lum = MINLUM;
		} else {
			d = TM_BRTSCALE*log(lum);	/* encode it */
			ls[i] = d>0. ? (int)(d+.5) : (int)(d-.5);
		}
		if (cs == TM_NOCHROM)			/* no color? */
			continue;
		if (tmTop->flags & TM_F_MESOPIC && lum < LMESUPPER) {
			slum = scotlum(cmon);		/* mesopic adj. */
			if (lum < LMESLOWER)
				cmon[RED] = cmon[GRN] = cmon[BLU] = slum;
			else {
				d = (lum - LMESLOWER)/(LMESUPPER - LMESLOWER);
				if (tmTop->flags & TM_F_BW)
					cmon[RED] = cmon[GRN] =
							cmon[BLU] = d*lum;
				else
					scalecolor(cmon, d);
				d = (1.-d)*slum;
				cmon[RED] += d;
				cmon[GRN] += d;
				cmon[BLU] += d;
			}
		} else if (tmTop->flags & TM_F_BW) {
			cmon[RED] = cmon[GRN] = cmon[BLU] = lum;
		}
		d = tmTop->clf[RED]*cmon[RED]/lum;
		cs[3*i  ] = d>=.999 ? 255 :
				(int)(256.*pow(d, 1./tmTop->mongam));
		d = tmTop->clf[GRN]*cmon[GRN]/lum;
		cs[3*i+1] = d>=.999 ? 255 :
				(int)(256.*pow(d, 1./tmTop->mongam));
		d = tmTop->clf[BLU]*cmon[BLU]/lum;
		cs[3*i+2] = d>=.999 ? 255 :
				(int)(256.*pow(d, 1./tmTop->mongam));
	}
	returnOK;
}


int
tmCvGrays(ls, scan, len)		/* convert float gray values */
TMbright	*ls;
float	*scan;
int	len;
{
	static char	funcName[] = "tmCvGrays";
	register double	d;
	register int	i;

	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | scan == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
	for (i = len; i--; )
		if (scan[i] <= TM_NOLUM) {
			ls[i] = TM_NOBRT;		/* bogus value */
		} else {
			d = TM_BRTSCALE*log(scan[i]);	/* encode it */
			ls[i] = d>0. ? (int)(d+.5) : (int)(d-.5);
		}
	returnOK;
}


int
tmAddHisto(ls, len, wt)			/* add values to histogram */
register TMbright	*ls;
int	len;
int	wt;
{
	static char	funcName[] = "tmAddHisto";
	int	oldorig=0, oldlen, horig, hlen;
	register int	i, j;

	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (len < 0)
		returnErr(TM_E_ILLEGAL);
	if (len == 0)
		returnOK;
						/* first, grow limits */
	if (tmTop->histo == NULL) {
		for (i = len; i-- && ls[i] < MINBRT; )
			;
		if (i < 0)
			returnOK;
		tmTop->hbrmin = tmTop->hbrmax = ls[i];
		oldlen = 0;
	} else {
		oldorig = (tmTop->hbrmin-MINBRT)/HISTEP;
		oldlen = (tmTop->hbrmax-MINBRT)/HISTEP + 1 - oldorig;
	}
	for (i = len; i--; ) {
		if ((j = ls[i]) < MINBRT)
			continue;
		if (j < tmTop->hbrmin)
			tmTop->hbrmin = j;
		else if (j > tmTop->hbrmax)
			tmTop->hbrmax = j;
	}
	horig = (tmTop->hbrmin-MINBRT)/HISTEP;
	hlen = (tmTop->hbrmax-MINBRT)/HISTEP + 1 - horig;
	if (hlen > oldlen) {			/* (re)allocate histogram */
		register int	*newhist = (int *)calloc(hlen, sizeof(int));
		if (newhist == NULL)
			returnErr(TM_E_NOMEM);
		if (oldlen) {			/* copy and free old */
			for (i = oldlen, j = i+oldorig-horig; i; )
				newhist[--j] = tmTop->histo[--i];
			free((MEM_PTR)tmTop->histo);
		}
		tmTop->histo = newhist;
	}
	if (wt == 0)
		returnOK;
	for (i = len; i--; )			/* add in new counts */
		if (ls[i] >= MINBRT)
			tmTop->histo[ (ls[i]-MINBRT)/HISTEP - horig ] += wt;
	returnOK;
}


static double
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


static int
tmNewMap()
{
	if (tmTop->lumap != NULL && (tmTop->mbrmax - tmTop->mbrmin) !=
					(tmTop->hbrmax - tmTop->hbrmin)) {
		free((MEM_PTR)tmTop->lumap);
		tmTop->lumap = NULL;
	}
	tmTop->mbrmin = tmTop->hbrmin;
	tmTop->mbrmax = tmTop->hbrmax;
	if (tmTop->mbrmin > tmTop->mbrmax)
		return 0;
	if (tmTop->lumap == NULL)
		tmTop->lumap = (unsigned short *)malloc(sizeof(unsigned short)*
					(tmTop->mbrmax-tmTop->mbrmin+1));
	return(tmTop->lumap != NULL);
}


int
tmFixedMapping(expmult, gamval)
double	expmult;
double	gamval;
{
	static char	funcName[] = "tmFixedMapping";
	double		d;
	register int	i;
	
	if (!tmNewMap())
		returnErr(TM_E_NOMEM);
	if (expmult <= .0)
		expmult = 1.;
	if (gamval < MINGAM)
		gamval = tmTop->mongam;
	d = log(expmult/tmTop->inpsf);
	for (i = tmTop->mbrmax-tmTop->mbrmin+1; i--; )
		tmTop->lumap[i] = 256. * exp(
			( d + (tmTop->mbrmin+i)*(1./TM_BRTSCALE) )
			/ gamval );
	returnOK;
}


int
tmComputeMapping(gamval, Lddyn, Ldmax)
double	gamval;
double	Lddyn;
double	Ldmax;
{
	static char	funcName[] = "tmComputeMapping";
	int	*histo;
	float	*cumf;
	int	brt0, histlen, threshold, ceiling, trimmings;
	double	logLddyn, Ldmin, Ldavg, Lwavg, Tr, Lw, Ld;
	int4	histot;
	double	sum;
	register double	d;
	register int	i, j;

	if (tmTop == NULL || tmTop->histo == NULL)
		returnErr(TM_E_TMINVAL);
					/* check arguments */
	if (Lddyn < MINLDDYN) Lddyn = DEFLDDYN;
	if (Ldmax < MINLDMAX) Ldmax = DEFLDMAX;
	if (gamval < MINGAM) gamval = tmTop->mongam;
					/* compute handy values */
	Ldmin = Ldmax/Lddyn;
	logLddyn = log(Lddyn);
	Ldavg = sqrt(Ldmax*Ldmin);
	i = (tmTop->hbrmin-MINBRT)/HISTEP;
	brt0 = MINBRT + HISTEP/2 + i*HISTEP;
	histlen = (tmTop->hbrmax-MINBRT)/HISTEP + 1 - i;
					/* histogram total and mean */
	histot = 0; sum = 0;
	j = brt0 + histlen*HISTEP;
	for (i = histlen; i--; ) {
		histot += tmTop->histo[i];
		sum += (j -= HISTEP) * tmTop->histo[i];
	}
	threshold = histot*.025 + .5;
	if (threshold < 4)
		returnErr(TM_E_TMFAIL);
	Lwavg = tmLuminance( (double)sum / histot );
					/* allocate space for mapping */
	if (!tmNewMap())
		returnErr(TM_E_NOMEM);
					/* use linear tone mapping? */
	if (tmTop->flags & TM_F_LINEAR)
		goto linearmap;
					/* clamp histogram */
	histo = (int *)malloc(histlen*sizeof(int));
	cumf = (float *)malloc((histlen+2)*sizeof(float));
	if (histo == NULL | cumf == NULL)
		returnErr(TM_E_NOMEM);
	cumf[histlen+1] = 1.;		/* guard for assignment code */
	for (i = histlen; i--; )	/* make malleable copy */
		histo[i] = tmTop->histo[i];
	do {				/* iterate to solution */
		sum = 0;		/* cumulative probability */
		for (i = 0; i < histlen; i++) {
			cumf[i] = (double)sum/histot;
			sum += histo[i];
		}
		cumf[histlen] = 1.;
		Tr = histot * (double)(tmTop->hbrmax - tmTop->hbrmin) /
			((double)histlen*TM_BRTSCALE) / logLddyn;
		ceiling = Tr + 1.;
		trimmings = 0;		/* clip to envelope */
		for (i = histlen; i--; ) {
			if (tmTop->flags & TM_F_HCONTR) {
				Lw = tmLuminance(brt0 + i*HISTEP);
				Ld = Ldmin * exp( logLddyn *
					.5*(cumf[i]+cumf[i+1]) );
				ceiling = Tr * (htcontrs(Ld) * Lw) /
					(htcontrs(Lw) * Ld) + 1.;
			}
			if (histo[i] > ceiling) {
				trimmings += histo[i] - ceiling;
				histo[i] = ceiling;
			}
		}
					/* check if we're out of data */
		if ((histot -= trimmings) <= threshold) {
			free((MEM_PTR)histo);
			free((MEM_PTR)cumf);
			goto linearmap;
		}
	} while (trimmings > threshold);
					/* assign tone-mapping */
	for (i = tmTop->mbrmax-tmTop->mbrmin+1; i--; ) {
		j = d = (double)i/(tmTop->mbrmax-tmTop->mbrmin)*histlen;
		d -= (double)j;
		Ld = Ldmin*exp(logLddyn*((1.-d)*cumf[j]+d*cumf[j+1]));
		d = (Ld - Ldmin)/(Ldmax - Ldmin);
		tmTop->lumap[i] = 256.*pow(d, 1./gamval);
	}
	free((MEM_PTR)histo);		/* clean up and return */
	free((MEM_PTR)cumf);
	returnOK;
linearmap:				/* linear tone-mapping */
	if (tmTop->flags & TM_F_HCONTR)
		d = htcontrs(Ldavg) / htcontrs(Lwavg);
	else
		d = Ldavg / Lwavg;
	return(tmFixedMapping(tmTop->inpsf*d/Ldmax, gamval));
}


int
tmMapPixels(ps, ls, cs, len)
register BYTE	*ps;
TMbright	*ls;
register BYTE	*cs;
int	len;
{
	static char	funcName[] = "tmMapPixels";
	register int4	li, pv;

	if (tmTop == NULL || tmTop->lumap == NULL)
		returnErr(TM_E_TMINVAL);
	if (ps == NULL | ls == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
	while (len--) {
		if ((li = *ls++) < tmTop->mbrmin) {
			li = 0;
		} else {
			if (li > tmTop->mbrmax)
				li = tmTop->mbrmax;
			li = tmTop->lumap[li - tmTop->mbrmin];
		}
		if (cs == TM_NOCHROM)
			*ps++ = li>255 ? 255 : li;
		else {
			pv = *cs++ * li / tmTop->cdiv[RED];
			*ps++ = pv>255 ? 255 : pv;
			pv = *cs++ * li / tmTop->cdiv[GRN];
			*ps++ = pv>255 ? 255 : pv;
			pv = *cs++ * li / tmTop->cdiv[BLU];
			*ps++ = pv>255 ? 255 : pv;
		}
	}
	returnOK;
}


struct tmStruct *
tmPop()				/* pop top tone mapping off stack */
{
	register struct tmStruct	*tms;

	if ((tms = tmTop) != NULL)
		tmTop = tms->tmprev;
	return(tms);
}


int
tmPull(tms)			/* pull a tone mapping from stack */
register struct tmStruct	*tms;
{
	register struct tmStruct	*tms2;
					/* special cases first */
	if (tms == NULL | tmTop == NULL)
		return(0);
	if (tms == tmTop) {
		tmTop = tms->tmprev;
		tms->tmprev = NULL;
		return(1);
	}
	for (tms2 = tmTop; tms2->tmprev != NULL; tms2 = tms2->tmprev)
		if (tms == tms2->tmprev) {	/* remove it */
			tms2->tmprev = tms->tmprev;
			tms->tmprev = NULL;
			return(1);
		}
	return(0);			/* not found on stack */
}


struct tmStruct *
tmDup()				/* duplicate top tone mapping */
{
	int	len;
	register int	i;
	register struct tmStruct	*tmnew;

	if (tmTop == NULL)		/* anything to duplicate? */
		return(NULL);
	tmnew = (struct tmStruct *)malloc(sizeof(struct tmStruct));
	if (tmnew == NULL)
		return(NULL);
	*tmnew = *tmTop;		/* copy everything */
	if (tmnew->histo != NULL) {	/* duplicate histogram */
		len = (tmnew->hbrmax-MINBRT)/HISTEP + 1 -
				(tmnew->hbrmin-MINBRT)/HISTEP;
		tmnew->histo = (int *)malloc(len*sizeof(int));
		if (tmnew->histo != NULL)
			for (i = len; i--; )
				tmnew->histo[i] = tmTop->histo[i];
	}
	if (tmnew->lumap != NULL) {	/* duplicate luminance mapping */
		len = tmnew->mbrmax-tmnew->mbrmin+1;
		tmnew->lumap = (unsigned short *)malloc(
						len*sizeof(unsigned short) );
		if (tmnew->lumap != NULL)
			for (i = len; i--; )
				tmnew->lumap[i] = tmTop->lumap[i];
	}
					/* clear package data */
	for (i = tmNumPkgs; i--; )
		tmnew->pd[i] = NULL;
	tmnew->tmprev = tmTop;		/* make copy current */
	return(tmTop = tmnew);
}


int
tmPush(tms)			/* push tone mapping on top of stack */
register struct tmStruct	*tms;
{
	static char	funcName[] = "tmPush";
					/* check validity */
	if (tms == NULL)
		returnErr(TM_E_ILLEGAL);
	if (tms == tmTop)		/* check necessity */
		returnOK;
					/* pull if already in stack */
	(void)tmPull(tms);
					/* push it on top */
	tms->tmprev = tmTop;
	tmTop = tms;
	returnOK;
}


void
tmDone(tms)			/* done with tone mapping -- destroy it */
register struct tmStruct	*tms;
{
	register int	i;
					/* NULL arg. is equiv. to tmTop */
	if (tms == NULL && (tms = tmTop) == NULL)
		return;
					/* take out of stack if present */
	(void)tmPull(tms);
					/* free tables */
	if (tms->histo != NULL)
		free((MEM_PTR)tms->histo);
	if (tms->lumap != NULL)
		free((MEM_PTR)tms->lumap);
					/* free private data */
	for (i = tmNumPkgs; i--; )
		if (tms->pd[i] != NULL)
			(*tmPkg[i]->Free)(tms->pd[i]);
	free((MEM_PTR)tms);		/* free basic structure */
}

/******************** Shared but Private library routines *********************/

BYTE	tmMesofact[BMESUPPER-BMESLOWER];

void
tmMkMesofact()				/* build mesopic lookup factor table */
{
	register int	i;

	if (tmMesofact[BMESUPPER-BMESLOWER-1])
		return;

	for (i = BMESLOWER; i < BMESUPPER; i++)
		tmMesofact[i-BMESLOWER] = 256. *
				(tmLuminance(i) - LMESLOWER) /
				(LMESUPPER - LMESLOWER);
}


int
tmErrorReturn(func, err)		/* error return (with message) */
char	*func;
int	err;
{
	tmLastFunction = func;
	tmLastError = err;
	if (tmTop != NULL && tmTop->flags & TM_F_NOSTDERR)
		return(err);
	fputs(func, stderr);
	fputs(": ", stderr);
	fputs(tmErrorMessage[err], stderr);
	fputs("!\n", stderr);
	return(err);
}
