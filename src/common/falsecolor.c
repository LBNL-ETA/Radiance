#ifndef lint
static const char RCSid[] = "$Id: falsecolor.c,v 3.9 2008/07/03 03:35:10 greg Exp $";
#endif
/*
 * False color mapping functions.
 * See falsecolor.h for detailed function descriptions.
 *
 * Externals declared in falsecolor.h
 */

#include "copyright.h"

#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#include	"tmprivat.h"
#include	"falsecolor.h"

/* Initialize new false color mapping */
FCstruct *
fcInit(BYTE fcscale[256][3])
{
	FCstruct	*fcs = (FCstruct *)malloc(sizeof(FCstruct));
	
	if (fcs == NULL)
		return(NULL);
	fcs->mbrmin = 10; fcs->mbrmax = -10;
	fcs->lumap = NULL;
	if ((fcs->scale = fcscale) == NULL)
		fcs->scale = fcDefaultScale;
	return(fcs);
}

/* Assign fixed linear false color map */
int
fcFixedLinear(FCstruct *fcs, double Lwmax)
{
	double	mult;
	int	i;

	if ((fcs == NULL) | (Lwmax <= MINLUM))
		return(TM_E_ILLEGAL);
	if (fcs->lumap != NULL)
		free((void *)fcs->lumap);
	fcs->mbrmin = tmCvLuminance(Lwmax/256.);
	fcs->mbrmax = tmCvLuminance(Lwmax);
	fcs->lumap = (BYTE *)malloc(sizeof(BYTE)*(fcs->mbrmax - fcs->mbrmin + 1));
	if (fcs->lumap == NULL)
		return(TM_E_NOMEM);
	mult = 255.999/tmLuminance(fcs->mbrmax);
	for (i = fcs->mbrmin; i <= fcs->mbrmax; i++)
		fcs->lumap[i - fcs->mbrmin] = (int)(mult * tmLuminance(i));
	returnOK;
}

/* Assign fixed logarithmic false color map */
int
fcFixedLog(FCstruct *fcs, double Lwmin, double Lwmax)
{
	int	i;

	if ((fcs == NULL) | (Lwmin <= MINLUM) | (Lwmax <= Lwmin))
		return(TM_E_ILLEGAL);
	if (fcs->lumap != NULL)
		free((void *)fcs->lumap);
	fcs->mbrmin = tmCvLuminance(Lwmin);
	fcs->mbrmax = tmCvLuminance(Lwmax);
	if (fcs->mbrmin >= fcs->mbrmax) {
		fcs->lumap = NULL;
		return(TM_E_ILLEGAL);
	}
	fcs->lumap = (BYTE *)malloc(sizeof(BYTE)*(fcs->mbrmax - fcs->mbrmin + 1));
	if (fcs->lumap == NULL)
		return(TM_E_NOMEM);
	for (i = fcs->mbrmax - fcs->mbrmin; i >= 0; i--)
		fcs->lumap[i] = 256L * i / (fcs->mbrmax - fcs->mbrmin + 1);
	returnOK;
}

/* Compute linear false color map */
int
fcLinearMapping(FCstruct *fcs, TMstruct *tms, double pctile)
{
	int		i, histlen;
	int32		histot, cnt;
	int		brt0;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(0 > pctile) | (pctile >= 50))
		return(TM_E_ILLEGAL);
	i = HISTI(tms->hbrmin);
	brt0 = HISTV(i);
	histlen = HISTI(tms->hbrmax) + 1 - i;
	histot = 0;
	for (i = histlen; i--; )
		histot += tms->histo[i];
	cnt = histot * pctile / 100;
	for (i = histlen; i--; )
		if ((cnt -= tms->histo[i]) < 0)
			break;
	if (i <= 0)
		return(TM_E_TMFAIL);
	return(fcFixedLinear(fcs, tmLuminance(brt0 + i*HISTEP)));
}

/* Compute logarithmic false color map */
int
fcLogMapping(FCstruct *fcs, TMstruct *tms, double pctile)
{
	int		i, histlen;
	int32		histot, cnt;
	int		brt0, wbrmin, wbrmax;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(.0 > pctile) | (pctile >= 50.))
		return(TM_E_ILLEGAL);
	i = HISTI(tms->hbrmin);
	brt0 = HISTV(i);
	histlen = HISTI(tms->hbrmax) + 1 - i;
	histot = 0;
	for (i = histlen; i--; )
		histot += tms->histo[i];
	cnt = histot * pctile * .01;
	for (i = 0; i < histlen; i++)
		if ((cnt -= tms->histo[i]) < 0)
			break;
	if (i >= histlen)
		return(TM_E_TMFAIL);
	wbrmin = brt0 + i*HISTEP;
	cnt = histot * pctile * .01;
	for (i = histlen; i--; )
		if ((cnt -= tms->histo[i]) < 0)
			break;
	wbrmax = brt0 + i*HISTEP;
	if (wbrmax <= wbrmin)
		return(TM_E_TMFAIL);
	return(fcFixedLog(fcs, tmLuminance(wbrmin), tmLuminance(wbrmax)));
}
 
/* Apply false color mapping to pixel values */
int
fcMapPixels(FCstruct *fcs, BYTE *ps, TMbright *ls, int len)
{
	int	li;

	if (fcs == NULL || (fcs->lumap == NULL) | (fcs->scale == NULL))
		return(TM_E_ILLEGAL);
	if ((ps == NULL) | (ls == NULL) | (len < 0))
		return(TM_E_ILLEGAL);
	while (len--) {
		if ((li = *ls++) < fcs->mbrmin)
			li = 0;
		else if (li > fcs->mbrmax)
			li = 255;
		else
			li = fcs->lumap[li - fcs->mbrmin];
		*ps++ = fcs->scale[li][RED];
		*ps++ = fcs->scale[li][GRN];
		*ps++ = fcs->scale[li][BLU];
	}
	returnOK;
}

/* Determine if false color mapping is logarithmic */
int
fcIsLogMap(FCstruct *fcs)
{
	int	miderr;

	if (fcs == NULL || fcs->lumap == NULL)
		return(-1);
	
	miderr = (fcs->mbrmax - fcs->mbrmin)/2;
	miderr = fcs->lumap[miderr] -
			256L * miderr / (fcs->mbrmax - fcs->mbrmin + 1);

	return((-1 <= miderr) & (miderr <= 1));
}

/* Duplicate a false color structure */
FCstruct *
fcDup(FCstruct *fcs)
{
	FCstruct	*fcnew;

	if (fcs == NULL)
		return(NULL);
	fcnew = fcInit(fcs->scale);
	if (fcnew == NULL)
		return(NULL);
	if (fcs->lumap != NULL) {
		fcnew->lumap = (BYTE *)malloc(sizeof(BYTE)*(fcs->mbrmax -
							fcs->mbrmin + 1));
		if (fcnew->lumap == NULL)
			return(fcnew);
		fcnew->mbrmin = fcs->mbrmin; fcnew->mbrmax = fcs->mbrmax;
		memcpy((void *)fcnew->lumap, (void *)fcs->lumap,
				sizeof(BYTE)*(fcs->mbrmax - fcs->mbrmin + 1));
	}
	return(fcnew);
}

/* Free data associated with the given false color mapping structure */
void
fcDone(FCstruct *fcs)
{
	if (fcs == NULL)
		return;
	if (fcs->lumap != NULL)
		free((void *)fcs->lumap);
	free((void *)fcs);
}

BYTE	fcDefaultScale[256][3] = {	/* default false color scale */
	111,8,132,
	108,7,133,
	105,7,134,
	102,6,136,
	98,6,137,
	93,5,139,
	89,4,141,
	84,3,143,
	79,2,145,
	74,1,148,
	68,0,150,
	63,0,153,
	57,0,155,
	52,0,157,
	46,0,160,
	41,0,162,
	36,0,164,
	31,0,166,
	26,0,168,
	22,0,170,
	18,0,172,
	14,2,174,
	11,4,175,
	 8,7,176,
	 7,9,177,
	 6,11,177,
	 5,13,178,
	 4,16,178,
	 3,18,179,
	 2,21,180,
	 1,24,180,
	 1,28,181,
	 0,31,181,
	 0,35,182,
	 0,38,182,
	 0,42,183,
	 0,46,184,
	 0,50,184,
	 0,54,184,
	 0,58,185,
	 0,63,185,
	 0,67,186,
	 0,71,186,
	 0,76,186,
	 0,80,187,
	 0,84,187,
	 0,89,187,
	 0,93,187,
	 1,97,187,
	 1,102,187,
	 1,106,187,
	 2,110,187,
	 2,114,187,
	 3,118,186,
	 3,122,186,
	 4,126,186,
	 4,130,185,
	 4,133,185,
	 5,137,184,
	 5,140,183,
	 6,143,182,
	 6,146,181,
	 6,149,180,
	 7,151,179,
	 7,154,178,
	 7,156,177,
	 8,158,175,
	 8,161,172,
	 9,163,169,
	 9,165,165,
	 9,167,161,
	 9,169,157,
	10,170,153,
	10,172,148,
	10,173,143,
	11,174,138,
	11,174,133,
	11,175,127,
	12,175,122,
	12,176,117,
	13,176,111,
	14,176,106,
	14,176,101,
	15,175,95,
	16,175,90,
	17,175,86,
	18,174,81,
	20,174,77,
	21,173,73,
	22,172,69,
	24,172,66,
	26,171,63,
	28,170,60,
	30,169,58,
	32,168,57,
	34,167,56,
	37,166,55,
	40,165,54,
	42,164,54,
	45,163,54,
	48,162,55,
	52,160,55,
	55,158,56,
	58,157,57,
	62,155,57,
	66,153,59,
	69,152,60,
	73,150,61,
	77,148,63,
	81,146,64,
	84,144,66,
	88,142,67,
	92,139,69,
	96,137,70,
	99,135,72,
	103,133,73,
	107,131,75,
	110,128,76,
	113,126,77,
	117,124,78,
	120,121,79,
	123,119,80,
	126,117,80,
	128,114,81,
	131,112,81,
	133,110,81,
	135,108,80,
	136,106,80,
	137,105,80,
	138,104,79,
	139,102,79,
	140,101,79,
	141,100,78,
	142,98,78,
	143,96,77,
	144,95,76,
	144,93,76,
	145,92,75,
	146,90,74,
	146,89,73,
	147,87,73,
	148,85,72,
	148,84,71,
	149,82,70,
	149,80,69,
	150,79,68,
	150,77,67,
	151,75,66,
	151,73,65,
	151,72,64,
	152,70,63,
	152,68,62,
	153,66,61,
	153,65,60,
	153,63,59,
	154,61,58,
	154,60,57,
	154,58,56,
	154,56,55,
	155,55,54,
	155,53,53,
	155,51,51,
	156,50,50,
	156,48,49,
	156,46,48,
	157,45,47,
	157,43,46,
	157,42,45,
	158,40,44,
	158,39,43,
	158,37,42,
	159,36,41,
	159,34,40,
	159,33,39,
	160,32,38,
	160,31,37,
	161,29,37,
	161,28,36,
	162,27,35,
	162,26,34,
	163,25,33,
	163,24,33,
	164,23,32,
	165,22,31,
	165,21,31,
	168,18,29,
	170,16,28,
	172,13,26,
	175,11,25,
	177,9,24,
	180,7,23,
	183,5,22,
	185,3,21,
	188,2,21,
	191,1,20,
	194,0,19,
	197,0,19,
	199,0,18,
	202,0,17,
	205,0,17,
	207,0,16,
	210,2,16,
	213,3,15,
	215,6,14,
	217,8,13,
	219,11,13,
	220,13,12,
	222,17,11,
	224,20,11,
	226,24,10,
	227,28,9,
	229,32,8,
	231,37,7,
	232,42,6,
	234,47,5,
	236,52,5,
	237,57,4,
	239,63,3,
	240,68,2,
	242,74,2,
	243,79,1,
	245,85,0,
	246,91,0,
	247,96,0,
	248,102,0,
	250,108,0,
	251,113,0,
	252,118,0,
	253,123,0,
	254,128,0,
	254,133,0,
	255,138,0,
	255,143,1,
	255,148,2,
	255,154,3,
	255,159,4,
	255,165,6,
	255,170,7,
	255,176,9,
	255,181,11,
	255,187,13,
	255,192,15,
	255,198,17,
	255,203,20,
	255,208,22,
	255,213,24,
	255,218,26,
	255,223,28,
	255,227,30,
	255,232,32,
	255,236,34,
	254,240,35,
	254,243,37,
	254,246,38,
	254,249,39,
	254,252,40,
};
