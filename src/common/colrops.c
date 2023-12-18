#ifndef lint
static const char	RCSid[] = "$Id: colrops.c,v 2.13 2023/12/18 18:54:38 greg Exp $";
#endif
/*
 * Integer operations on COLR scanlines
 */

#include "copyright.h"

#include <stdio.h>
#include <math.h>

#include "rtmisc.h"
#include "color.h"


#define MAXGSHIFT	31		/* maximum shift for gamma table */

static uby8	*g_mant = NULL, *g_nexp = NULL;

static uby8	(*g_bval)[256] = NULL;


int
setcolrcor(			/* set brightness correction */
	double	(*f)(double,double),
	double	a2
)
{
	double	mult;
	int	i, j;
					/* allocate tables */
	if (g_bval == NULL && (g_bval =
			(uby8 (*)[256])bmalloc((MAXGSHIFT+1)*256)) == NULL)
		return(-1);
					/* compute colr -> gamb mapping */
	mult = 1.0/256.0;
	for (i = 0; i <= MAXGSHIFT; i++) {
		for (j = 0; j < 256; j++)
			g_bval[i][j] = 256.0 * (*f)((j+.5)*mult, a2);
		mult *= 0.5;
	}
	return(0);
}


int
setcolrinv(			/* set inverse brightness correction */
	double	(*f)(double,double),
	double	a2
)
{
	double	mult;
	int	i, j;
					/* allocate tables */
	if (g_mant == NULL && (g_mant = (uby8 *)bmalloc(256)) == NULL)
		return(-1);
	if (g_nexp == NULL && (g_nexp = (uby8 *)bmalloc(256)) == NULL)
		return(-1);
					/* compute gamb -> colr mapping */
	i = 0;
	mult = 256.0;
	for (j = 256; j--; ) {
		while ((g_mant[j] = mult * (*f)((j+.5)/256.0, a2)) < 128) {
			i++;
			mult *= 2.0;
		}
		g_nexp[j] = i;
	}
	return(0);
}


int
setcolrgam(				/* set gamma conversion */
	double	g
)
{
	if (setcolrcor(pow, 1.0/g) < 0)
		return(-1);
	return(setcolrinv(pow, g));
}


int
colrs_gambs(			/* convert scanline of colrs to gamma bytes */
	COLR	*scan,
	int	len
)
{
	int	i, expo;

	if (g_bval == NULL)
		return(-1);
	while (len-- > 0) {
		expo = scan[0][EXP] - COLXS;
		if (expo < -MAXGSHIFT) {
			if (expo < -MAXGSHIFT-8) {
				scan[0][RED] =
				scan[0][GRN] =
				scan[0][BLU] = 0;
			} else {
				i = (-MAXGSHIFT-1) - expo;
				scan[0][RED] = 
				g_bval[MAXGSHIFT][((scan[0][RED]>>i)+1)>>1];
				scan[0][GRN] =
				g_bval[MAXGSHIFT][((scan[0][GRN]>>i)+1)>>1];
				scan[0][BLU] =
				g_bval[MAXGSHIFT][((scan[0][BLU]>>i)+1)>>1];
			}
		} else if (expo > 0) {
			if (expo > 8) {
				scan[0][RED] =
				scan[0][GRN] =
				scan[0][BLU] = 255;
			} else {
				i = (scan[0][RED]<<1 | 1) << (expo-1);
				scan[0][RED] = i > 255 ? 255 : g_bval[0][i];
				i = (scan[0][GRN]<<1 | 1) << (expo-1);
				scan[0][GRN] = i > 255 ? 255 : g_bval[0][i];
				i = (scan[0][BLU]<<1 | 1) << (expo-1);
				scan[0][BLU] = i > 255 ? 255 : g_bval[0][i];
			}
		} else {
			scan[0][RED] = g_bval[-expo][scan[0][RED]];
			scan[0][GRN] = g_bval[-expo][scan[0][GRN]];
			scan[0][BLU] = g_bval[-expo][scan[0][BLU]];
		}
		scan[0][EXP] = COLXS;
		scan++;
	}
	return(0);
}


int
gambs_colrs(		/* convert gamma bytes to colr scanline */
	COLR	*scan,
	int	len
)
{
	int	nexpo;

	if ((g_mant == NULL) | (g_nexp == NULL))
		return(-1);
	while (len-- > 0) {
		nexpo = g_nexp[scan[0][RED]];
		if (g_nexp[scan[0][GRN]] < nexpo)
			nexpo = g_nexp[scan[0][GRN]];
		if (g_nexp[scan[0][BLU]] < nexpo)
			nexpo = g_nexp[scan[0][BLU]];
		if (nexpo < g_nexp[scan[0][RED]])
			scan[0][RED] = g_mant[scan[0][RED]]
					>> (g_nexp[scan[0][RED]]-nexpo);
		else
			scan[0][RED] = g_mant[scan[0][RED]];
		if (nexpo < g_nexp[scan[0][GRN]])
			scan[0][GRN] = g_mant[scan[0][GRN]]
					>> (g_nexp[scan[0][GRN]]-nexpo);
		else
			scan[0][GRN] = g_mant[scan[0][GRN]];
		if (nexpo < g_nexp[scan[0][BLU]])
			scan[0][BLU] = g_mant[scan[0][BLU]]
					>> (g_nexp[scan[0][BLU]]-nexpo);
		else
			scan[0][BLU] = g_mant[scan[0][BLU]];
		scan[0][EXP] = COLXS - nexpo;
		scan++;
	}
	return(0);
}


void
shiftcolrs(		/* shift a scanline of colors by 2^adjust */
	COLR	*scan,
	int	len,
	int	adjust
)
{
	int	minexp;

	if (adjust == 0)
		return;
	minexp = adjust < 0 ? -adjust : 0;
	while (len-- > 0) {
		if (scan[0][EXP] <= minexp)
			scan[0][RED] = scan[0][GRN] = scan[0][BLU] =
			scan[0][EXP] = 0;
		else
			scan[0][EXP] += adjust;
		scan++;
	}
}


void
normcolrs(		/* normalize a scanline of colrs */
	COLR  *scan,
	int  len,
	int  adjust
)
{
	int  c;
	int  shift;

	while (len-- > 0) {
		shift = scan[0][EXP] + adjust - COLXS;
		if (shift > 0) {
			if (shift > 8) {
				scan[0][RED] =
				scan[0][GRN] =
				scan[0][BLU] = 255;
			} else {
				shift--;
				c = (scan[0][RED]<<1 | 1) << shift;
				scan[0][RED] = c > 255 ? 255 : c;
				c = (scan[0][GRN]<<1 | 1) << shift;
				scan[0][GRN] = c > 255 ? 255 : c;
				c = (scan[0][BLU]<<1 | 1) << shift;
				scan[0][BLU] = c > 255 ? 255 : c;
			}
		} else if (shift < 0) {
			if (shift < -8) {
				scan[0][RED] =
				scan[0][GRN] =
				scan[0][BLU] = 0;
			} else {
				shift = -1-shift;
				scan[0][RED] = ((scan[0][RED]>>shift)+1)>>1;
				scan[0][GRN] = ((scan[0][GRN]>>shift)+1)>>1;
				scan[0][BLU] = ((scan[0][BLU]>>shift)+1)>>1;
			}
		}
		scan[0][EXP] = COLXS - adjust;
		scan++;
	}
}
