#ifndef lint
static const char	RCSid[] = "$Id: tonemap.c,v 3.31 2008/07/03 03:30:21 greg Exp $";
#endif
/*
 * Tone mapping functions.
 * See tonemap.h for detailed function descriptions.
 * Added von Kries white-balance calculations 10/01 (GW).
 *
 * Externals declared in tonemap.h
 */

#include "copyright.h"

#include	<stdio.h>
#include	<math.h>
#include	"tmprivat.h"
#include	"tmerrmsg.h"

#define	exp10(x)	exp(M_LN10*(x))

					/* our list of conversion packages */
struct tmPackage	*tmPkg[TM_MAXPKG];
int	tmNumPkgs = 0;			/* number of registered packages */

					/* luminance->brightness lookup */
static TMbright		*tmFloat2BrtLUT = NULL;

#define tmCvLumLUfp(pf)	tmFloat2BrtLUT[*(int32 *)(pf) >> 15]


TMstruct *
tmInit(					/* initialize new tone mapping */
int	flags,
RGBPRIMP	monpri,
double	gamval
)
{
	COLORMAT	cmat;
	TMstruct	*tmnew;
	int	i;
						/* allocate structure */
	tmnew = (TMstruct *)malloc(sizeof(TMstruct));
	if (tmnew == NULL)
		return(NULL);

	tmnew->flags = flags & ~TM_F_UNIMPL;
	if (tmnew->flags & TM_F_BW)
		tmnew->flags &= ~TM_F_MESOPIC;
						/* set monitor transform */
	if (monpri == NULL || monpri == stdprims || tmnew->flags & TM_F_BW) {
		tmnew->monpri = stdprims;
		tmnew->clf[RED] = rgb2xyzmat[1][0];
		tmnew->clf[GRN] = rgb2xyzmat[1][1];
		tmnew->clf[BLU] = rgb2xyzmat[1][2];
	} else {
		comprgb2xyzmat(cmat, tmnew->monpri=monpri);
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
	tmnew->inpdat = NULL;
	tmnew->hbrmin = 10; tmnew->hbrmax = -10;
	tmnew->histo = NULL;
	tmnew->mbrmin = 10; tmnew->mbrmax = -10;
	tmnew->lumap = NULL;
						/* zero private data */
	for (i = TM_MAXPKG; i--; )
		tmnew->pd[i] = NULL;
	tmnew->lastError = TM_E_OK;
	tmnew->lastFunc = "NoErr";
						/* return new TMstruct */
	return(tmnew);
}


int
tmSetSpace(			/* set input color space for conversions */
TMstruct	*tms,
RGBPRIMP	pri,
double	sf,
MEM_PTR	dat
)
{
	static const char funcName[] = "tmSetSpace";
	int	i, j;
						/* error check */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if (sf <= 1e-12)
		returnErr(TM_E_ILLEGAL);
						/* check if no change */
	if (pri == tms->inppri && FEQ(sf, tms->inpsf) && dat == tms->inpdat)
		returnOK;
	tms->inppri = pri;			/* let's set it */
	tms->inpsf = sf;
	tms->inpdat = dat;

	if (tms->flags & TM_F_BW) {		/* color doesn't matter */
		tms->monpri = tms->inppri;		/* eliminate xform */
		if (tms->inppri == TM_XYZPRIM) {
			tms->clf[CIEX] = tms->clf[CIEZ] = 0.;
			tms->clf[CIEY] = 1.;
		} else {
			comprgb2xyzmat(tms->cmat, tms->monpri);
			tms->clf[RED] = tms->cmat[1][0];
			tms->clf[GRN] = tms->cmat[1][1];
			tms->clf[BLU] = tms->cmat[1][2];
		}
		tms->cmat[0][0] = tms->cmat[1][1] = tms->cmat[2][2] =
				tms->inpsf;
		tms->cmat[0][1] = tms->cmat[0][2] = tms->cmat[1][0] =
		tms->cmat[1][2] = tms->cmat[2][0] = tms->cmat[2][1] = 0.;

	} else if (tms->inppri == TM_XYZPRIM)	/* input is XYZ */
		compxyz2rgbWBmat(tms->cmat, tms->monpri);

	else {					/* input is RGB */
		if (tms->inppri != tms->monpri &&
				PRIMEQ(tms->inppri, tms->monpri))
			tms->inppri = tms->monpri;	/* no xform */
		comprgb2rgbWBmat(tms->cmat, tms->inppri, tms->monpri);
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			tms->cmat[i][j] *= tms->inpsf;
						/* set color divisors */
	for (i = 0; i < 3; i++)
		if (tms->clf[i] > .001)
			tms->cdiv[i] =
				256.*pow(tms->clf[i], 1./tms->mongam);
		else
			tms->cdiv[i] = 1;
						/* notify packages */
	for (i = tmNumPkgs; i--; )
		if (tms->pd[i] != NULL && tmPkg[i]->NewSpace != NULL)
			(*tmPkg[i]->NewSpace)(tms);
	returnOK;
}


void
tmClearHisto(				/* clear current histogram */
TMstruct	*tms
)
{
	if (tms == NULL || tms->histo == NULL)
		return;
	free((MEM_PTR)tms->histo);
	tms->histo = NULL;
}


TMbright
tmCvLuminance(				/* convert a single luminance */
double	lum
)
{
	double	d;

#ifdef isfinite
	if (!isfinite(lum) || lum <= TM_NOLUM)
#else
	if (lum <= TM_NOLUM)
#endif
		return(TM_NOBRT);
	d = TM_BRTSCALE*log(lum);
	if (d > 0.)
		return((TMbright)(d+.5));
	return((TMbright)(d-.5));
}


int
tmCvLums(				/* convert luminances using lookup */
TMbright	*ls,
float		*scan,
int		len
)
{
	if (tmFloat2BrtLUT == NULL) {	/* initialize lookup table */
		int32	i;
		tmFloat2BrtLUT = (TMbright *)malloc(sizeof(TMbright)*0x10000);
		if (tmFloat2BrtLUT == NULL)
			return(TM_E_NOMEM);
		for (i = 0; i < 0x10000; i++) {
			int32	l = (i<<1 | 1) << 14;
#ifndef isfinite
			if ((l & 0x7f800000) == 0x7f800000)
				tmFloat2BrtLUT[i] = TM_NOBRT;
			else
#endif
			tmFloat2BrtLUT[i] = tmCvLuminance(*(float *)&l);
		}
	}
	if (len <= 0)
		return(TM_E_OK);
	if ((ls == NULL) | (scan == NULL))
		return(TM_E_ILLEGAL);
	while (len--) {
		if (*scan <= TM_NOLUM) {
			*ls++ = TM_NOBRT;
			++scan;
			continue;
		}
		*ls++ = tmCvLumLUfp(scan++);
	}
	return(TM_E_OK);
}


int
tmCvGrays(				/* convert float gray values */
TMstruct	*tms,
TMbright	*ls,
float		*scan,
int		len
)
{
	static const char funcName[] = "tmCvGrays";
	int	i;

	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (scan == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
	if (tmFloat2BrtLUT == NULL)			/* initialize */
		tmCvLums(NULL, NULL, 0);
	for (i = len; i--; ) {
		float	lum = tms->inpsf * scan[i];
		if (lum <= TM_NOLUM)
			ls[i] = TM_NOBRT;
		else
			ls[i] = tmCvLumLUfp(&lum);
	}
	returnOK;
}


int
tmCvColors(				/* convert float colors */
TMstruct	*tms,
TMbright	*ls,
BYTE	*cs,
COLOR	*scan,
int	len
)
{
	static const char funcName[] = "tmCvColors";
	static BYTE	gamtab[1024];
	static double	curgam = .0;
	COLOR	cmon;
	float	lum, slum, d;
	int	i;

	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (scan == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
	if (tmFloat2BrtLUT == NULL)			/* initialize */
		tmCvLums(NULL, NULL, 0);
	if (cs != TM_NOCHROM && fabs(tms->mongam - curgam) > .02) {
		curgam = tms->mongam;			/* (re)build table */
		for (i = 1024; i--; )
			gamtab[i] = (int)(256.*pow((i+.5)/1024., 1./curgam));
	}
	for (i = len; i--; ) {
		if (tmNeedMatrix(tms)) {		/* get monitor RGB */
			colortrans(cmon, tms->cmat, scan[i]);
		} else {
			cmon[RED] = tms->inpsf*scan[i][RED];
			cmon[GRN] = tms->inpsf*scan[i][GRN];
			cmon[BLU] = tms->inpsf*scan[i][BLU];
		}
#ifdef isfinite
		if (!isfinite(cmon[RED]) || cmon[RED] < .0f) cmon[RED] = .0f;
		if (!isfinite(cmon[GRN]) || cmon[GRN] < .0f) cmon[GRN] = .0f;
		if (!isfinite(cmon[BLU]) || cmon[BLU] < .0f) cmon[BLU] = .0f;
#else
		if (cmon[RED] < .0f) cmon[RED] = .0f;
		if (cmon[GRN] < .0f) cmon[GRN] = .0f;
		if (cmon[BLU] < .0f) cmon[BLU] = .0f;
#endif
							/* world luminance */
		lum =	tms->clf[RED]*cmon[RED] +
			tms->clf[GRN]*cmon[GRN] +
			tms->clf[BLU]*cmon[BLU] ;
		if (lum <= TM_NOLUM) {			/* convert brightness */
			lum = cmon[RED] = cmon[GRN] = cmon[BLU] = TM_NOLUM;
			ls[i] = TM_NOBRT;
		} else
			ls[i] = tmCvLumLUfp(&lum);
		if (cs == TM_NOCHROM)			/* no color? */
			continue;
		if (tms->flags & TM_F_MESOPIC && lum < LMESUPPER) {
			slum = scotlum(cmon);		/* mesopic adj. */
			if (lum < LMESLOWER) {
				cmon[RED] = cmon[GRN] = cmon[BLU] = slum;
			} else {
				d = (lum - LMESLOWER)/(LMESUPPER - LMESLOWER);
				if (tms->flags & TM_F_BW)
					cmon[RED] = cmon[GRN] =
							cmon[BLU] = d*lum;
				else
					scalecolor(cmon, d);
				d = (1.f-d)*slum;
				cmon[RED] += d;
				cmon[GRN] += d;
				cmon[BLU] += d;
			}
		} else if (tms->flags & TM_F_BW) {
			cmon[RED] = cmon[GRN] = cmon[BLU] = lum;
		}
		d = tms->clf[RED]*cmon[RED]/lum;
		cs[3*i  ] = d>=.999f ? 255 : gamtab[(int)(1024.f*d)];
		d = tms->clf[GRN]*cmon[GRN]/lum;
		cs[3*i+1] = d>=.999f ? 255 : gamtab[(int)(1024.f*d)];
		d = tms->clf[BLU]*cmon[BLU]/lum;
		cs[3*i+2] = d>=.999f ? 255 : gamtab[(int)(1024.f*d)];
	}
	returnOK;
}


int
tmAddHisto(				/* add values to histogram */
TMstruct	*tms,
TMbright	*ls,
int	len,
int	wt
)
{
	static const char funcName[] = "tmAddHisto";
	int	oldorig=0, oldlen, horig, hlen;
	int	i, j;

	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if (len < 0)
		returnErr(TM_E_ILLEGAL);
	if (len == 0)
		returnOK;
						/* first, grow limits */
	if (tms->histo == NULL) {
		for (i = len; i-- && ls[i] < MINBRT; )
			;
		if (i < 0)
			returnOK;
		tms->hbrmin = tms->hbrmax = ls[i];
		oldlen = 0;
	} else {
		oldorig = HISTI(tms->hbrmin);
		oldlen = HISTI(tms->hbrmax) + 1 - oldorig;
	}
	for (i = len; i--; ) {
		if ((j = ls[i]) < MINBRT)
			continue;
		if (j < tms->hbrmin)
			tms->hbrmin = j;
		else if (j > tms->hbrmax)
			tms->hbrmax = j;
	}
	horig = HISTI(tms->hbrmin);
	hlen = HISTI(tms->hbrmax) + 1 - horig;
	if (hlen > oldlen) {			/* (re)allocate histogram */
		int	*newhist = (int *)calloc(hlen, sizeof(int));
		if (newhist == NULL)
			returnErr(TM_E_NOMEM);
		if (oldlen) {			/* copy and free old */
			for (i = oldlen, j = i+oldorig-horig; i; )
				newhist[--j] = tms->histo[--i];
			free((MEM_PTR)tms->histo);
		}
		tms->histo = newhist;
	}
	if (wt == 0)
		returnOK;
	for (i = len; i--; )			/* add in new counts */
		if (ls[i] >= MINBRT)
			tms->histo[ HISTI(ls[i]) - horig ] += wt;
	returnOK;
}


static double
htcontrs(		/* human threshold contrast sensitivity, dL(La) */
double	La
)
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
tmNewMap(			/* allocate new tone-mapping array */
TMstruct	*tms
)
{
	if (tms->lumap != NULL && (tms->mbrmax - tms->mbrmin) !=
					(tms->hbrmax - tms->hbrmin)) {
		free((MEM_PTR)tms->lumap);
		tms->lumap = NULL;
	}
	tms->mbrmin = tms->hbrmin;
	tms->mbrmax = tms->hbrmax;
	if (tms->mbrmin > tms->mbrmax)
		return 0;
	if (tms->lumap == NULL)
		tms->lumap = (unsigned short *)malloc(sizeof(unsigned short)*
					(tms->mbrmax-tms->mbrmin+1));
	return(tms->lumap != NULL);
}


int
tmFixedMapping(			/* compute fixed, linear tone-mapping */
TMstruct	*tms,
double	expmult,
double	gamval
)
{
	static const char funcName[] = "tmFixedMapping";
	double		d;
	int	i;
	
	if (!tmNewMap(tms))
		returnErr(TM_E_NOMEM);
	if (expmult <= .0)
		expmult = 1.;
	if (gamval < MINGAM)
		gamval = tms->mongam;
	d = log(expmult/tms->inpsf);
	for (i = tms->mbrmax-tms->mbrmin+1; i--; ) {
		double	val = 256. * exp(
			( d + (tms->mbrmin+i)*(1./TM_BRTSCALE) )
			/ gamval);
		tms->lumap[i] = val >= (double)0xffff ? 0xffff : (int)val;
	}
	returnOK;
}


int
tmComputeMapping(			/* compute histogram tone-mapping */
TMstruct	*tms,
double	gamval,
double	Lddyn,
double	Ldmax
)
{
	static const char funcName[] = "tmComputeMapping";
	int	*histo;
	float	*cumf;
	int	brt0, histlen, threshold, ceiling, trimmings;
	double	logLddyn, Ldmin, Ldavg, Lwavg, Tr, Lw, Ld;
	int32	histot;
	double	sum;
	double	d;
	int	i, j;

	if (tms == NULL || tms->histo == NULL)
		returnErr(TM_E_TMINVAL);
					/* check arguments */
	if (Lddyn < MINLDDYN) Lddyn = DEFLDDYN;
	if (Ldmax < MINLDMAX) Ldmax = DEFLDMAX;
	if (gamval < MINGAM) gamval = tms->mongam;
					/* compute handy values */
	Ldmin = Ldmax/Lddyn;
	logLddyn = log(Lddyn);
	Ldavg = sqrt(Ldmax*Ldmin);
	i = HISTI(tms->hbrmin);
	brt0 = MINBRT + HISTEP/2 + i*HISTEP;
	histlen = HISTI(tms->hbrmax) + 1 - i;
					/* histogram total and mean */
	histot = 0; sum = 0;
	j = brt0 + histlen*HISTEP;
	for (i = histlen; i--; ) {
		histot += tms->histo[i];
		sum += (double)(j -= HISTEP) * tms->histo[i];
	}
	threshold = histot*0.005 + .5;
	if (!histot)
		returnErr(TM_E_TMFAIL);
	Lwavg = tmLuminance( (double)sum / histot );
					/* use linear tone mapping? */
	if (tms->flags & TM_F_LINEAR || threshold < 4 ||
			tms->hbrmax - tms->hbrmin < TM_BRTSCALE*logLddyn)
		goto linearmap;
					/* clamp histogram */
	histo = (int *)malloc(histlen*sizeof(int));
	cumf = (float *)malloc((histlen+2)*sizeof(float));
	if ((histo == NULL) | (cumf == NULL))
		returnErr(TM_E_NOMEM);
	cumf[histlen+1] = 1.;		/* guard for assignment code */
	for (i = histlen; i--; )	/* make malleable copy */
		histo[i] = tms->histo[i];
	do {				/* iterate to solution */
		sum = 0;		/* cumulative probability */
		for (i = 0; i < histlen; i++) {
			cumf[i] = (double)sum/histot;
			sum += histo[i];
		}
		cumf[histlen] = 1.;
		Tr = histot * (double)(tms->hbrmax - tms->hbrmin) /
			((double)histlen*TM_BRTSCALE) / logLddyn;
		ceiling = Tr + 1.;
		trimmings = 0;		/* clip to envelope */
		for (i = histlen; i--; ) {
			if (tms->flags & TM_F_HCONTR) {
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
					/* allocate space for mapping */
	if (!tmNewMap(tms))
		returnErr(TM_E_NOMEM);
					/* assign tone-mapping */
	for (i = tms->mbrmax-tms->mbrmin+1; i--; ) {
		j = d = (double)i/(tms->mbrmax-tms->mbrmin)*histlen;
		d -= (double)j;
		Ld = Ldmin*exp(logLddyn*((1.-d)*cumf[j]+d*cumf[j+1]));
		d = (Ld - Ldmin)/(Ldmax - Ldmin);
		tms->lumap[i] = 256.*pow(d, 1./gamval);
	}
	free((MEM_PTR)histo);		/* clean up and return */
	free((MEM_PTR)cumf);
	returnOK;
linearmap:				/* linear tone-mapping */
	if (tms->flags & TM_F_HCONTR)
		d = htcontrs(Ldavg) / htcontrs(Lwavg);
	else
		d = Ldavg / Lwavg;
	return(tmFixedMapping(tms, tms->inpsf*d/Ldmax, gamval));
}


int
tmMapPixels(			/* apply tone-mapping to pixel(s) */
TMstruct	*tms,
BYTE	*ps,
TMbright	*ls,
BYTE	*cs,
int	len
)
{
	static const char funcName[] = "tmMapPixels";
	int32	li, pv;

	if (tms == NULL || tms->lumap == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ps == NULL) | (ls == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
	while (len--) {
		if ((li = *ls++) < tms->mbrmin) {
			li = 0;
		} else {
			if (li > tms->mbrmax)
				li = tms->mbrmax;
			li = tms->lumap[li - tms->mbrmin];
		}
		if (cs == TM_NOCHROM)
			*ps++ = li>255 ? 255 : li;
		else {
			pv = *cs++ * li / tms->cdiv[RED];
			*ps++ = pv>255 ? 255 : pv;
			pv = *cs++ * li / tms->cdiv[GRN];
			*ps++ = pv>255 ? 255 : pv;
			pv = *cs++ * li / tms->cdiv[BLU];
			*ps++ = pv>255 ? 255 : pv;
		}
	}
	returnOK;
}




TMstruct *
tmDup(				/* duplicate top tone mapping */
TMstruct	*tms
)
{
	int	len;
	int	i;
	TMstruct	*tmnew;

	if (tms == NULL)		/* anything to duplicate? */
		return(NULL);
	tmnew = (TMstruct *)malloc(sizeof(TMstruct));
	if (tmnew == NULL)
		return(NULL);
	*tmnew = *tms;		/* copy everything */
	if (tmnew->histo != NULL) {	/* duplicate histogram */
		len = HISTI(tmnew->hbrmax) + 1 - HISTI(tmnew->hbrmin);
		tmnew->histo = (int *)malloc(len*sizeof(int));
		if (tmnew->histo != NULL)
			for (i = len; i--; )
				tmnew->histo[i] = tms->histo[i];
	}
	if (tmnew->lumap != NULL) {	/* duplicate luminance mapping */
		len = tmnew->mbrmax-tmnew->mbrmin+1;
		tmnew->lumap = (unsigned short *)malloc(
						len*sizeof(unsigned short) );
		if (tmnew->lumap != NULL)
			for (i = len; i--; )
				tmnew->lumap[i] = tms->lumap[i];
	}
					/* clear package data */
	for (i = tmNumPkgs; i--; )
		tmnew->pd[i] = NULL;
					/* return copy */
	return(tmnew);
}


void
tmDone(tms)			/* done with tone mapping -- destroy it */
TMstruct	*tms;
{
	int	i;
					/* NULL arg. is equiv. to tms */
	if (tms == NULL)
		return;
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
	int	i;

	if (tmMesofact[BMESUPPER-BMESLOWER-1])
		return;

	for (i = BMESLOWER; i < BMESUPPER; i++)
		tmMesofact[i-BMESLOWER] = 256. *
				(tmLuminance(i) - LMESLOWER) /
				(LMESUPPER - LMESLOWER);
}


int
tmErrorReturn(				/* error return (with message) */
const char	*func,
TMstruct	*tms,
int	err
)
{
	if (tms != NULL) {
		tms->lastFunc = func;
		tms->lastError = err;
		if (tms->flags & TM_F_NOSTDERR)
			return(err);
	}
	fputs(func, stderr);
	fputs(": ", stderr);
	fputs(tmErrorMessage[err], stderr);
	fputs("!\n", stderr);
	return(err);
}
