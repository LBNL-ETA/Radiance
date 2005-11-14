#ifndef lint
static const char RCSid[] = "$Id: falsecolor.c,v 3.1 2005/11/14 22:18:18 greg Exp $";
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
	double	v;
	int	i;

	if ((fcs == NULL) | (Lwmax <= MINLUM))
		return(TM_E_ILLEGAL);
	if (fcs->lumap != NULL)
		free((void *)fcs->lumap);
	v = TM_BRTSCALE * log(Lwmax);
	fcs->mbrmax = (int)(v<0. ? v-.5 : v+.5);
	v = TM_BRTSCALE * log(Lwmax/256.);
	fcs->mbrmin = (int)(v<0. ? v-.5 : v+.5);
	fcs->lumap = (BYTE *)malloc(sizeof(BYTE)*(fcs->mbrmax - fcs->mbrmin + 1));
	if (fcs->lumap == NULL)
		return(TM_E_NOMEM);
	v = 255.999/tmLuminance(fcs->mbrmax);
	for (i = fcs->mbrmin; i <= fcs->mbrmax; i++)
		fcs->lumap[i] = (int)(v * tmLuminance(i));
	returnOK;
}

/* Assign fixed logarithmic false color map */
int
fcFixedLog(FCstruct *fcs, double Lwmin, double Lwmax)
{
	double	v;
	int	i;

	if ((fcs == NULL) | (Lwmin <= MINLUM) | (Lwmax <= Lwmin))
		return(TM_E_ILLEGAL);
	if (fcs->lumap != NULL)
		free((void *)fcs->lumap);
	v = TM_BRTSCALE * log(Lwmin);
	fcs->mbrmin = (int)(v<0. ? v-.5 : v+.5);
	v = TM_BRTSCALE * log(Lwmax);
	fcs->mbrmax = (int)(v<0. ? v-.5 : v+.5);
	fcs->lumap = (BYTE *)malloc(sizeof(BYTE)*(fcs->mbrmax - fcs->mbrmin + 1));
	if (fcs->lumap == NULL)
		return(TM_E_NOMEM);
	for (i = fcs->mbrmin; i <= fcs->mbrmax; i++)
		fcs->lumap[i] = 256*(i - fcs->mbrmin) /
					(fcs->mbrmax - fcs->mbrmin + 1);
	returnOK;
}

/* Compute linear false color map */
int
fcLinearMapping(FCstruct *fcs, TMstruct *tms, int pctile)
{
	int		i;
	int32		histot, cnt;
	TMbright	wbrmin, wbrmax;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(0 > pctile) | (pctile >= 100))
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
fcLogMapping(FCstruct *fcs, TMstruct *tms, int pctile)
{
	int		i;
	int32		histot, cnt;
	TMbright	wbrmin, wbrmax;

	if ((fcs == NULL) | (tms == NULL) || (tms->histo == NULL) |
			(0 > pctile) | (pctile >= 100))
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
	return(fcFixedLog(fcs, tmLuminannce(wbrmin), tmLuminance(wbrmax)));
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
        48,0,68, 45,0,70, 42,0,72, 39,0,74, 36,0,76, 33,0,78, 30,0,81,
        27,0,83, 24,0,85, 21,0,87, 18,0,89, 15,0,91, 13,0,94, 12,0,96,
        11,0,99, 9,0,101, 8,0,104, 7,0,106, 6,0,109, 5,0,111, 4,0,114,
         2,0,116, 1,0,119, 0,0,121, 0,0,122, 0,1,123, 0,1,124, 0,2,125,
         0,2,125, 0,3,126, 0,4,127, 0,4,128, 0,5,129, 0,5,129, 0,6,130,
         0,7,131, 0,9,131, 0,11,132, 0,13,133, 0,16,133, 0,18,134,
         0,20,134, 0,22,135, 0,24,135, 0,26,136, 0,29,136, 0,31,137,
         0,34,137, 0,37,136, 0,41,136, 0,45,135, 0,48,135, 0,52,135,
         0,55,134, 0,59,134, 0,62,134, 0,66,133, 0,69,133, 0,73,132,
         0,76,130, 0,79,127, 0,82,125, 0,85,123, 0,88,120, 0,91,118,
         0,94,115, 0,97,113, 0,101,110, 0,104,108, 0,107,106, 0,109,102,
         0,110,97, 0,111,91, 0,112,86, 0,113,81, 0,115,76, 0,116,71,
         0,117,65, 0,118,60, 0,119,55, 0,120,50, 0,122,45, 0,121,42,
         1,120,39, 1,119,36, 1,119,34, 1,118,31, 2,117,28, 2,116,26,
         2,115,23, 2,115,20, 3,114,18, 3,113,15, 3,112,13, 4,110,13,
         5,108,13, 6,106,13, 7,104,12, 9,102,12, 10,100,12, 11,98,12,
        12,97,12, 13,95,12, 14,93,11, 15,91,11, 17,89,12, 19,86,12,
        22,83,12, 24,81,13, 26,78,13, 29,76,14, 31,73,14, 34,70,15,
        36,68,15, 39,65,16, 41,63,16, 44,60,17, 46,58,17, 49,56,17,
        51,54,17, 54,52,17, 57,50,17, 59,48,17, 62,45,17, 64,43,17,
        67,41,17, 70,39,17, 72,37,17, 74,35,17, 75,34,16, 76,33,16,
        77,32,16, 79,31,15, 80,30,15, 81,29,14, 82,28,14, 83,26,13,
        84,25,13, 85,24,13, 87,23,12, 87,22,12, 88,21,11, 89,20,11,
        90,19,10, 91,18,10, 92,17,9, 93,16,9, 94,15,8, 95,14,8,
        95,13,7, 96,12,7, 97,11,7, 98,11,6, 99,10,6, 99,9,5, 100,9,5,
       101,8,5, 102,8,4, 102,7,4, 103,6,4, 104,6,3, 104,5,3, 105,4,2,
       106,4,2, 107,4,2, 107,3,2, 108,3,2, 109,3,2, 109,2,1, 110,2,1,
       111,2,1, 112,1,1, 112,1,1, 113,1,0, 114,1,0, 115,0,0, 116,0,0,
       117,0,0, 118,0,0, 119,0,0, 121,0,0, 122,0,0, 123,0,0, 124,0,0,
       125,0,0, 126,0,0, 128,0,0, 131,0,0, 134,0,0, 137,0,0, 140,0,0,
       144,0,0, 147,0,0, 150,0,0, 153,0,0, 156,0,0, 159,0,0, 162,0,0,
       165,0,0, 168,0,0, 171,0,0, 174,0,0, 177,0,0, 180,1,0, 183,1,0,
       186,1,0, 189,1,0, 192,1,0, 195,2,0, 198,2,0, 201,5,0, 204,7,0,
       207,9,0, 210,11,0, 213,13,0, 216,15,0, 219,17,0, 222,19,0,
       225,21,0, 228,23,0, 230,25,0, 233,29,0, 235,34,0, 237,39,0,
       239,43,0, 241,48,0, 243,52,0, 245,57,0, 247,61,0, 250,66,0,
       252,71,0, 254,75,0, 255,80,0, 255,88,0, 255,95,1, 255,103,1,
       255,110,1, 255,117,1, 255,125,1, 255,132,2, 255,139,2, 255,147,2,
       255,154,2, 255,162,2, 255,169,3, 255,176,3, 255,183,3, 254,190,4,
       254,198,4, 254,205,4, 254,212,4, 253,219,5, 253,226,5, 253,234,5,
       252,241,6, 252,248,6
};
