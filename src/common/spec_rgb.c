#ifndef lint
static const char	RCSid[] = "$Id: spec_rgb.c,v 2.11 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Convert colors and spectral ranges.
 * Added von Kries white-balance calculations 10/01 (GW).
 *
 * Externals declared in color.h
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

#include "color.h"
#include <string.h>

#define CEPS	1e-4			/* color epsilon */

#define CEQ(v1,v2)	((v1) <= (v2)+CEPS && (v2) <= (v1)+CEPS)

#define XYEQ(c1,c2)	(CEQ((c1)[CIEX],(c2)[CIEX]) && CEQ((c1)[CIEY],(c2)[CIEY]))


RGBPRIMS  stdprims = STDPRIMS;	/* standard primary chromaticities */

COLOR  cblack = BLKCOLOR;		/* global black color */
COLOR  cwhite = WHTCOLOR;		/* global white color */

float  xyneu[2] = {1./3., 1./3.};	/* neutral xy chromaticities */

/*
 *	The following table contains the CIE tristimulus integrals
 *  for X, Y, and Z.  The table is cumulative, so that
 *  each color coordinate integrates to 1.
 */

#define  STARTWL	380		/* starting wavelength (nanometers) */
#define  INCWL		10		/* wavelength increment */
#define  NINC		40		/* # of values */

static BYTE  chroma[3][NINC] = {
	{							/* X */
		0,   0,   0,   2,   6,   13,  22,  30,  36,  41,
		42,  43,  43,  44,  46,  52,  60,  71,  87,  106,
		128, 153, 178, 200, 219, 233, 243, 249, 252, 254,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255
	}, {							/* Y */
		0,   0,   0,   0,   0,   1,   2,   4,   7,   11,
		17,  24,  34,  48,  64,  84,  105, 127, 148, 169,
		188, 205, 220, 232, 240, 246, 250, 253, 254, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255
	}, {							/* Z */
		0,   0,   2,   10,  32,  66,  118, 153, 191, 220,
		237, 246, 251, 253, 254, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255
	}
};

COLORMAT  xyz2rgbmat = {		/* XYZ to RGB (no white balance) */
	{(CIE_y_g - CIE_y_b - CIE_x_b*CIE_y_g + CIE_y_b*CIE_x_g)/CIE_C_rD,
	 (CIE_x_b - CIE_x_g - CIE_x_b*CIE_y_g + CIE_x_g*CIE_y_b)/CIE_C_rD,
	 (CIE_x_g*CIE_y_b - CIE_x_b*CIE_y_g)/CIE_C_rD},
	{(CIE_y_b - CIE_y_r - CIE_y_b*CIE_x_r + CIE_y_r*CIE_x_b)/CIE_C_gD,
	 (CIE_x_r - CIE_x_b - CIE_x_r*CIE_y_b + CIE_x_b*CIE_y_r)/CIE_C_gD,
	 (CIE_x_b*CIE_y_r - CIE_x_r*CIE_y_b)/CIE_C_gD},
	{(CIE_y_r - CIE_y_g - CIE_y_r*CIE_x_g + CIE_y_g*CIE_x_r)/CIE_C_bD,
	 (CIE_x_g - CIE_x_r - CIE_x_g*CIE_y_r + CIE_x_r*CIE_y_g)/CIE_C_bD,
	 (CIE_x_r*CIE_y_g - CIE_x_g*CIE_y_r)/CIE_C_bD}
};

COLORMAT  rgb2xyzmat = {		/* RGB to XYZ (no white balance) */
	{CIE_x_r*CIE_C_rD/CIE_D,CIE_x_g*CIE_C_gD/CIE_D,CIE_x_b*CIE_C_bD/CIE_D},
	{CIE_y_r*CIE_C_rD/CIE_D,CIE_y_g*CIE_C_gD/CIE_D,CIE_y_b*CIE_C_bD/CIE_D},
	{(1.-CIE_x_r-CIE_y_r)*CIE_C_rD/CIE_D,
	 (1.-CIE_x_g-CIE_y_g)*CIE_C_gD/CIE_D,
	 (1.-CIE_x_b-CIE_y_b)*CIE_C_bD/CIE_D}
};

COLORMAT  vkmat = {		/* Sharp primary matrix */
	{ 1.2694, -0.0988, -0.1706},
	{-0.8364,  1.8006,  0.0357},
	{ 0.0297, -0.0315,  1.0018}
};

COLORMAT  ivkmat = {		/* inverse Sharp primary matrix */
	{ 0.8156,  0.0472,  0.1372},
	{ 0.3791,  0.5769,  0.0440},
	{-0.0123,  0.0167,  0.9955}
};


void
spec_rgb(col, s, e)		/* compute RGB color from spectral range */
COLOR  col;
int  s, e;
{
	COLOR  ciecolor;

	spec_cie(ciecolor, s, e);
	cie_rgb(col, ciecolor);
}


void
spec_cie(col, s, e)		/* compute a color from a spectral range */
COLOR  col;		/* returned color */
int  s, e;		/* starting and ending wavelengths */
{
	register int  i, d, r;
	
	s -= STARTWL;
	if (s < 0)
		s = 0;

	e -= STARTWL;
	if (e <= s) {
		col[CIEX] = col[CIEY] = col[CIEZ] = 0.0;
		return;
	}
	if (e >= INCWL*(NINC - 1))
		e = INCWL*(NINC - 1) - 1;

	d = e / INCWL;			/* interpolate values */
	r = e % INCWL;
	for (i = 0; i < 3; i++)
		col[i] = chroma[i][d]*(INCWL - r) + chroma[i][d + 1]*r;

	d = s / INCWL;
	r = s % INCWL;
	for (i = 0; i < 3; i++)
		col[i] -= chroma[i][d]*(INCWL - r) - chroma[i][d + 1]*r;

	col[CIEX] = (col[CIEX] + 0.5) * (1./(256*INCWL));
	col[CIEY] = (col[CIEY] + 0.5) * (1./(256*INCWL));
	col[CIEZ] = (col[CIEZ] + 0.5) * (1./(256*INCWL));
}


void
cie_rgb(rgb, xyz)		/* convert CIE color to standard RGB */
COLOR	rgb;
COLOR  xyz;
{
	colortrans(rgb, xyz2rgbmat, xyz);
	clipgamut(rgb, xyz[CIEY], CGAMUT_LOWER, cblack, cwhite);
}


int
clipgamut(col, brt, gamut, lower, upper)	/* clip to gamut cube */
COLOR  col;
double  brt;
int  gamut;
COLOR  lower, upper;
{
	int  rflags = 0;
	double  brtmin, brtmax, v, vv;
	COLOR  cgry;
	register int  i;
					/* check for no check */
	if (gamut == 0) return(0);
					/* check brightness limits */
	brtmin = 1./3.*(lower[0]+lower[1]+lower[2]);
	if (gamut & CGAMUT_LOWER && brt < brtmin) {
		copycolor(col, lower);
		return(CGAMUT_LOWER);
	}
	brtmax = 1./3.*(upper[0]+upper[1]+upper[2]);
	if (gamut & CGAMUT_UPPER && brt > brtmax) {
		copycolor(col, upper);
		return(CGAMUT_UPPER);
	}
					/* compute equivalent grey */
	v = (brt - brtmin)/(brtmax - brtmin);
	for (i = 0; i < 3; i++)
		cgry[i] = v*upper[i] + (1.-v)*lower[i];
	vv = 1.;			/* check each limit */
	for (i = 0; i < 3; i++)
		if (gamut & CGAMUT_LOWER && col[i] < lower[i]) {
			v = (lower[i]+CEPS - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_LOWER;
		} else if (gamut & CGAMUT_UPPER && col[i] > upper[i]) {
			v = (upper[i]-CEPS - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_UPPER;
		}
	if (rflags)			/* desaturate to cube face */
		for (i = 0; i < 3; i++)
			col[i] = vv*col[i] + (1.-vv)*cgry[i];
	return(rflags);
}


void
colortrans(c2, mat, c1)		/* convert c1 by mat and put into c2 */
register COLOR  c2;
register COLORMAT  mat;
register COLOR  c1;
{
	COLOR	cout;

	cout[0] = mat[0][0]*c1[0] + mat[0][1]*c1[1] + mat[0][2]*c1[2];
	cout[1] = mat[1][0]*c1[0] + mat[1][1]*c1[1] + mat[1][2]*c1[2];
	cout[2] = mat[2][0]*c1[0] + mat[2][1]*c1[1] + mat[2][2]*c1[2];

	copycolor(c2, cout);
}


void
multcolormat(m3, m2, m1)	/* multiply m1 by m2 and put into m3 */
COLORMAT  m3;			/* m3 can be either m1 or m2 w/o harm */
COLORMAT  m2, m1;
{
	COLORMAT  mt;
	register int  i, j;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			mt[i][j] =	m1[i][0]*m2[0][j] +
					m1[i][1]*m2[1][j] +
					m1[i][2]*m2[2][j] ;
	cpcolormat(m3, mt);
}


void
compxyz2rgbmat(mat, pr)		/* compute conversion from CIE to RGB space */
COLORMAT  mat;
register RGBPRIMS  pr;
{
	double  C_rD, C_gD, C_bD;

	if (pr == stdprims) {	/* can use xyz2rgbmat */
		cpcolormat(mat, xyz2rgbmat);
		return;
	}
	C_rD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) -
			  pr[WHT][CIEY]*(pr[GRN][CIEX] - pr[BLU][CIEX]) +
		  pr[GRN][CIEX]*pr[BLU][CIEY] - pr[BLU][CIEX]*pr[GRN][CIEY] ) ;
	C_gD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) -
			  pr[WHT][CIEY]*(pr[BLU][CIEX] - pr[RED][CIEX]) -
		  pr[RED][CIEX]*pr[BLU][CIEY] + pr[BLU][CIEX]*pr[RED][CIEY] ) ;
	C_bD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) -
			  pr[WHT][CIEY]*(pr[RED][CIEX] - pr[GRN][CIEX]) +
		  pr[RED][CIEX]*pr[GRN][CIEY] - pr[GRN][CIEX]*pr[RED][CIEY] ) ;

	mat[0][0] = (pr[GRN][CIEY] - pr[BLU][CIEY] -
			pr[BLU][CIEX]*pr[GRN][CIEY] +
			pr[BLU][CIEY]*pr[GRN][CIEX])/C_rD ;
	mat[0][1] = (pr[BLU][CIEX] - pr[GRN][CIEX] -
			pr[BLU][CIEX]*pr[GRN][CIEY] +
			pr[GRN][CIEX]*pr[BLU][CIEY])/C_rD ;
	mat[0][2] = (pr[GRN][CIEX]*pr[BLU][CIEY] -
			pr[BLU][CIEX]*pr[GRN][CIEY])/C_rD ;
	mat[1][0] = (pr[BLU][CIEY] - pr[RED][CIEY] -
			pr[BLU][CIEY]*pr[RED][CIEX] +
			pr[RED][CIEY]*pr[BLU][CIEX])/C_gD ;
	mat[1][1] = (pr[RED][CIEX] - pr[BLU][CIEX] -
			pr[RED][CIEX]*pr[BLU][CIEY] +
			pr[BLU][CIEX]*pr[RED][CIEY])/C_gD ;
	mat[1][2] = (pr[BLU][CIEX]*pr[RED][CIEY] -
			pr[RED][CIEX]*pr[BLU][CIEY])/C_gD ;
	mat[2][0] = (pr[RED][CIEY] - pr[GRN][CIEY] -
			pr[RED][CIEY]*pr[GRN][CIEX] +
			pr[GRN][CIEY]*pr[RED][CIEX])/C_bD ;
	mat[2][1] = (pr[GRN][CIEX] - pr[RED][CIEX] -
			pr[GRN][CIEX]*pr[RED][CIEY] +
			pr[RED][CIEX]*pr[GRN][CIEY])/C_bD ;
	mat[2][2] = (pr[RED][CIEX]*pr[GRN][CIEY] -
			pr[GRN][CIEX]*pr[RED][CIEY])/C_bD ;
}


void
comprgb2xyzmat(mat, pr)		/* compute conversion from RGB to CIE space */
COLORMAT  mat;
register RGBPRIMS  pr;
{
	double  C_rD, C_gD, C_bD, D;

	if (pr == stdprims) {	/* can use rgb2xyzmat */
		cpcolormat(mat, rgb2xyzmat);
		return;
	}
	C_rD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) -
			  pr[WHT][CIEY]*(pr[GRN][CIEX] - pr[BLU][CIEX]) +
		  pr[GRN][CIEX]*pr[BLU][CIEY] - pr[BLU][CIEX]*pr[GRN][CIEY] ) ;
	C_gD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) -
			  pr[WHT][CIEY]*(pr[BLU][CIEX] - pr[RED][CIEX]) -
		  pr[RED][CIEX]*pr[BLU][CIEY] + pr[BLU][CIEX]*pr[RED][CIEY] ) ;
	C_bD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) -
			  pr[WHT][CIEY]*(pr[RED][CIEX] - pr[GRN][CIEX]) +
		  pr[RED][CIEX]*pr[GRN][CIEY] - pr[GRN][CIEX]*pr[RED][CIEY] ) ;
	D = pr[RED][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) +
			pr[GRN][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) +
			pr[BLU][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) ;
	mat[0][0] = pr[RED][CIEX]*C_rD/D;
	mat[0][1] = pr[GRN][CIEX]*C_gD/D;
	mat[0][2] = pr[BLU][CIEX]*C_bD/D;
	mat[1][0] = pr[RED][CIEY]*C_rD/D;
	mat[1][1] = pr[GRN][CIEY]*C_gD/D;
	mat[1][2] = pr[BLU][CIEY]*C_bD/D;
	mat[2][0] = (1.-pr[RED][CIEX]-pr[RED][CIEY])*C_rD/D;
	mat[2][1] = (1.-pr[GRN][CIEX]-pr[GRN][CIEY])*C_gD/D;
	mat[2][2] = (1.-pr[BLU][CIEX]-pr[BLU][CIEY])*C_bD/D;
}


void
comprgb2rgbmat(mat, pr1, pr2)	/* compute conversion from RGB1 to RGB2 */
COLORMAT  mat;
RGBPRIMS  pr1, pr2;
{
	COLORMAT  pr1toxyz, xyztopr2;

	if (pr1 == pr2) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return;
	}
	comprgb2xyzmat(pr1toxyz, pr1);
	compxyz2rgbmat(xyztopr2, pr2);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, xyztopr2);
}


void
compxyzWBmat(mat, wht1, wht2)	/* CIE von Kries transform from wht1 to wht2 */
COLORMAT  mat;
float  wht1[2], wht2[2];
{
	COLOR	cw1, cw2;
	if (XYEQ(wht1,wht2)) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return;
	}
	cw1[RED] = wht1[CIEX]/wht1[CIEY];
	cw1[GRN] = 1.;
	cw1[BLU] = (1. - wht1[CIEX] - wht1[CIEY])/wht1[CIEY];
	colortrans(cw1, vkmat, cw1);
	cw2[RED] = wht2[CIEX]/wht2[CIEY];
	cw2[GRN] = 1.;
	cw2[BLU] = (1. - wht2[CIEX] - wht2[CIEY])/wht2[CIEY];
	colortrans(cw2, vkmat, cw2);
	mat[0][0] = cw2[RED]/cw1[RED];
	mat[1][1] = cw2[GRN]/cw1[GRN];
	mat[2][2] = cw2[BLU]/cw1[BLU];
	mat[0][1] = mat[0][2] = mat[1][0] =
	mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
	multcolormat(mat, vkmat, mat);
	multcolormat(mat, mat, ivkmat);
}


void
compxyz2rgbWBmat(mat, pr)	/* von Kries conversion from CIE to RGB space */
COLORMAT  mat;
RGBPRIMS  pr;
{
	COLORMAT	wbmat;

	compxyz2rgbmat(mat, pr);
	if (XYEQ(pr[WHT],xyneu))
		return;
	compxyzWBmat(wbmat, xyneu, pr[WHT]);
	multcolormat(mat, wbmat, mat);
}

void
comprgb2xyzWBmat(mat, pr)	/* von Kries conversion from RGB to CIE space */
COLORMAT  mat;
RGBPRIMS  pr;
{
	COLORMAT	wbmat;
	
	comprgb2xyzmat(mat, pr);
	if (XYEQ(pr[WHT],xyneu))
		return;
	compxyzWBmat(wbmat, pr[WHT], xyneu);
	multcolormat(mat, mat, wbmat);
}

void
comprgb2rgbWBmat(mat, pr1, pr2)	/* von Kries conversion from RGB1 to RGB2 */
COLORMAT  mat;
RGBPRIMS  pr1, pr2;
{
	COLORMAT  pr1toxyz, xyztopr2, wbmat;

	if (pr1 == pr2) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return;
	}
	comprgb2xyzmat(pr1toxyz, pr1);
	compxyzWBmat(wbmat, pr1[WHT], pr2[WHT]);
	compxyz2rgbmat(xyztopr2, pr2);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, wbmat);
	multcolormat(mat, mat, xyztopr2);
}
