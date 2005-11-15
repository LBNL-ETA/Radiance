#ifndef lint
static const char RCSid[] = "$Id: falsecolor.c,v 3.3 2005/11/15 20:04:50 greg Exp $";
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
	int		i;
	int32		histot, cnt;
	TMbright	wbrmin, wbrmax;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(0 > pctile) | (pctile >= 50))
		return(TM_E_ILLEGAL);
	histot = 0;
	for (i = tms->hbrmax - tms->hbrmin; i >= 0; )
		histot += tms->histo[i--];
	cnt = histot * pctile / 100;
	for (i = tms->hbrmax - tms->hbrmin; i >= 0; i--)
		if ((cnt -= tms->histo[i]) < 0)
			break;
	if (i < 0)
		return(TM_E_TMFAIL);
	return(fcFixedLinear(fcs, tmLuminance(tms->hbrmin + i)));
}

/* Compute logarithmic false color map */
int
fcLogMapping(FCstruct *fcs, TMstruct *tms, double pctile)
{
	int		i;
	int32		histot, cnt;
	TMbright	wbrmin, wbrmax;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(0 > pctile) | (pctile >= 50))
		return(TM_E_ILLEGAL);
	histot = 0;
	for (i = tms->hbrmax - tms->hbrmin; i >= 0; )
		histot += tms->histo[i--];
	cnt = histot * pctile / 100;
	for (i = 0; i <= tms->hbrmax - tms->hbrmin; i++)
		if ((cnt -= tms->histo[i]) < 0)
			break;
	if (i >= tms->hbrmax - tms->hbrmin)
		return(TM_E_TMFAIL);
	wbrmin = tms->hbrmin + i;
	cnt = histot * pctile / 100;
	for (i = tms->hbrmax - tms->hbrmin; i >= 0; i--)
		if ((cnt -= tms->histo[i]) < 0)
			break;
	if (i < 0)
		return(TM_E_TMFAIL);
	wbrmax = tms->hbrmin + i;
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
	int	midval;

	if (fcs == NULL || fcs->lumap == NULL)
		return(-1);
	midval = fcs->lumap[(fcs->mbrmax - fcs->mbrmin)/2];
	return((127 <= midval) & (midval <= 129));
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
       119,10,140,
       116,10,142,
       112,10,144,
       109,9,146,
       105,9,148,
       101,8,150,
        97,7,151,
        92,7,153,
        88,6,155,
        83,5,157,
        77,4,158,
        72,2,160,
        67,2,162,
        64,4,164,
        61,5,166,
        58,6,168,
        55,6,170,
        51,7,171,
        47,8,173,
        43,8,175,
        38,9,177,
        33,9,179,
        25,10,180,
        14,10,182,
        10,17,183,
        10,22,183,
         9,27,184,
         9,30,184,
         8,33,185,
         8,36,185,
         7,39,186,
         6,41,187,
         5,43,187,
         4,45,188,
         3,47,188,
         0,50,189,
         0,57,189,
         0,62,189,
         0,67,190,
         0,72,190,
         0,77,190,
         0,81,191,
         0,84,191,
         1,88,191,
         1,92,192,
         1,95,192,
         1,98,192,
         2,102,192,
         3,107,192,
         3,112,192,
         4,116,191,
         5,120,191,
         5,124,191,
         5,128,191,
         6,131,191,
         6,135,190,
         7,138,190,
         7,141,190,
         7,145,189,
         8,147,188,
         8,150,186,
         9,153,185,
         9,155,183,
        10,158,181,
        10,160,180,
        11,163,178,
        11,165,176,
        11,167,175,
        12,170,173,
        12,172,171,
        13,173,168,
        13,174,164,
        14,175,160,
        14,176,156,
        15,177,152,
        15,178,147,
        16,178,143,
        16,179,138,
        17,180,133,
        17,181,127,
        18,182,122,
        18,182,116,
        20,182,112,
        22,181,109,
        24,181,105,
        26,180,102,
        27,180,98,
        29,179,94,
        30,179,90,
        32,178,86,
        33,177,81,
        34,177,77,
        35,176,71,
        37,176,67,
        42,174,67,
        46,173,66,
        49,171,66,
        52,170,65,
        55,169,65,
        58,167,65,
        61,166,64,
        64,164,64,
        66,163,63,
        68,161,63,
        71,160,63,
        75,158,63,
        79,156,64,
        84,154,65,
        88,151,66,
        92,149,68,
        95,147,69,
        99,145,70,
       102,142,71,
       105,140,71,
       109,137,72,
       112,135,73,
       115,132,74,
       118,130,74,
       121,128,74,
       123,126,74,
       126,124,74,
       129,121,75,
       132,119,75,
       134,117,75,
       137,114,75,
       139,112,75,
       142,109,75,
       144,106,75,
       146,104,75,
       147,103,74,
       148,101,73,
       149,100,72,
       150,98,71,
       151,96,70,
       152,95,70,
       152,93,69,
       153,91,68,
       154,90,67,
       155,88,66,
       156,86,65,
       157,84,64,
       158,83,63,
       158,81,62,
       159,79,60,
       160,77,59,
       161,75,58,
       161,73,57,
       162,71,55,
       163,69,54,
       163,67,52,
       164,65,51,
       165,63,49,
       165,61,48,
       166,60,47,
       166,58,46,
       167,56,44,
       167,55,43,
       168,53,42,
       169,51,40,
       169,49,39,
       170,46,37,
       170,44,35,
       171,42,33,
       171,40,32,
       172,39,31,
       172,37,30,
       173,36,29,
       173,34,28,
       174,33,27,
       174,31,25,
       175,29,24,
       175,27,22,
       176,25,21,
       176,23,19,
       177,21,17,
       178,20,16,
       179,19,16,
       179,18,15,
       180,17,14,
       181,16,13,
       182,15,13,
       182,13,12,
       183,12,11,
       184,10,9,
       185,8,8,
       185,5,6,
       187,3,5,
       189,3,5,
       191,3,5,
       193,3,6,
       195,3,6,
       197,2,6,
       198,2,6,
       200,2,6,
       202,2,6,
       204,1,6,
       206,1,6,
       208,0,6,
       209,9,7,
       211,13,7,
       213,15,8,
       214,18,8,
       216,20,9,
       218,21,9,
       219,23,9,
       221,24,10,
       223,26,10,
       224,27,10,
       226,28,11,
       228,33,11,
       229,43,11,
       231,50,11,
       232,56,11,
       234,61,11,
       235,66,11,
       237,71,11,
       238,75,11,
       239,79,11,
       241,83,11,
       242,86,11,
       244,90,11,
       245,96,11,
       246,102,12,
       247,108,12,
       248,114,13,
       249,119,13,
       250,124,13,
       251,129,14,
       252,134,14,
       253,138,15,
       254,143,15,
       255,147,15,
       255,151,16,
       255,157,18,
       255,163,20,
       255,169,22,
       255,174,24,
       255,179,25,
       255,184,26,
       255,189,28,
       255,194,29,
       255,199,30,
       255,203,31,
       255,207,32,
       255,212,34,
       255,216,35,
       255,220,37,
       255,224,38,
       255,227,40,
       255,231,41,
       255,235,42,
       254,238,43,
       254,242,44,
       254,245,46,
       254,249,47,
       254,252,48,
};
