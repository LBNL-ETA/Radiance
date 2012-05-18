#ifndef lint
static const char	RCSid[] = "$Id: spec_rgb.c,v 2.24 2012/05/18 22:48:39 greg Exp $";
#endif
/*
 * Convert colors and spectral ranges.
 * Added von Kries white-balance calculations 10/01 (GW).
 *
 * Externals declared in color.h
 */

#include "copyright.h"

#include <stdio.h>
#include <string.h>
#include "color.h"

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

static uby8  chroma[3][NINC] = {
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
spec_rgb(			/* compute RGB color from spectral range */
COLOR  col,
int  s,
int  e
)
{
	COLOR  ciecolor;

	spec_cie(ciecolor, s, e);
	cie_rgb(col, ciecolor);
}


void
spec_cie(			/* compute a color from a spectral range */
COLOR  col,		/* returned color */
int  s,			/* starting and ending wavelengths */
int  e
)
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
cie_rgb(			/* convert CIE color to standard RGB */
COLOR	rgb,
COLOR  xyz
)
{
	colortrans(rgb, xyz2rgbmat, xyz);
	clipgamut(rgb, xyz[CIEY], CGAMUT_LOWER, cblack, cwhite);
}


int
clipgamut(			/* clip to gamut cube */
COLOR  col,
double  brt,
int  gamut,
COLOR  lower,
COLOR  upper
)
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
			v = (lower[i] - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_LOWER;
		} else if (gamut & CGAMUT_UPPER && col[i] > upper[i]) {
			v = (upper[i] - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_UPPER;
		}
	if (rflags)			/* desaturate to cube face */
		for (i = 0; i < 3; i++)
			col[i] = vv*col[i] + (1.-vv)*cgry[i];
	return(rflags);
}


void
colortrans(			/* convert c1 by mat and put into c2 */
register COLOR  c2,
register COLORMAT  mat,
register COLOR  c1
)
{
	COLOR	cout;

	cout[0] = mat[0][0]*c1[0] + mat[0][1]*c1[1] + mat[0][2]*c1[2];
	cout[1] = mat[1][0]*c1[0] + mat[1][1]*c1[1] + mat[1][2]*c1[2];
	cout[2] = mat[2][0]*c1[0] + mat[2][1]*c1[1] + mat[2][2]*c1[2];

	copycolor(c2, cout);
}


void
multcolormat(			/* multiply m1 by m2 and put into m3 */
COLORMAT  m3,			/* m3 can be either m1 or m2 w/o harm */
COLORMAT  m2,
COLORMAT  m1
)
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


int
colorprimsOK(			/* are color primaries reasonable? */
RGBPRIMS  pr
)
{
	int	i, j;
	
	for (i = 0; i < 3; i++) {
		if ((pr[i][CIEX] <= -2.) | (pr[i][CIEY] <= -2.))
			return(0);
		if ((pr[i][CIEX] >= 3.) | (pr[i][CIEY] >= 3.))
			return(0);
		if (pr[i][CIEX] + pr[i][CIEY] <= -2.)
			return(0);
		if (pr[i][CIEX] + pr[i][CIEY] >= 3.)
			return(0);
	}
	if ((pr[3][CIEX] <= 0.) | (pr[3][CIEX] >= 1.) |
			(pr[3][CIEY] <= 0.) | (pr[3][CIEY] >= 1.))
		return(0);
	for (i = 0; i < 4; i++)
		for (j = i+1; j < 4; j++)
			if (CEQ(pr[i][CIEX],pr[j][CIEX]) &&
					CEQ(pr[i][CIEY],pr[j][CIEY]))
				return(0);
	return(1);
}



int
compxyz2rgbmat(			/* compute conversion from CIE to RGB space */
COLORMAT  mat,
register RGBPRIMS  pr
)
{
	double  C_rD, C_gD, C_bD;

	if (pr == stdprims) {	/* can use xyz2rgbmat */
		cpcolormat(mat, xyz2rgbmat);
		return(1);
	}
	if (CEQ(pr[WHT][CIEX],0.) | CEQ(pr[WHT][CIEY],0.))
		return(0);
	C_rD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) -
			  pr[WHT][CIEY]*(pr[GRN][CIEX] - pr[BLU][CIEX]) +
		  pr[GRN][CIEX]*pr[BLU][CIEY] - pr[BLU][CIEX]*pr[GRN][CIEY] ) ;
	if (CEQ(C_rD,0.))
		return(0);
	C_gD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) -
			  pr[WHT][CIEY]*(pr[BLU][CIEX] - pr[RED][CIEX]) -
		  pr[RED][CIEX]*pr[BLU][CIEY] + pr[BLU][CIEX]*pr[RED][CIEY] ) ;
	if (CEQ(C_gD,0.))
		return(0);
	C_bD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) -
			  pr[WHT][CIEY]*(pr[RED][CIEX] - pr[GRN][CIEX]) +
		  pr[RED][CIEX]*pr[GRN][CIEY] - pr[GRN][CIEX]*pr[RED][CIEY] ) ;
	if (CEQ(C_bD,0.))
		return(0);
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
	return(1);
}


int
comprgb2xyzmat(			/* compute conversion from RGB to CIE space */
COLORMAT  mat,
register RGBPRIMS  pr
)
{
	double  C_rD, C_gD, C_bD, D;

	if (pr == stdprims) {	/* can use rgb2xyzmat */
		cpcolormat(mat, rgb2xyzmat);
		return(1);
	}
	if (CEQ(pr[WHT][CIEX],0.) | CEQ(pr[WHT][CIEY],0.))
		return(0);
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
	if (CEQ(D,0.))
		return(0);
	mat[0][0] = pr[RED][CIEX]*C_rD/D;
	mat[0][1] = pr[GRN][CIEX]*C_gD/D;
	mat[0][2] = pr[BLU][CIEX]*C_bD/D;
	mat[1][0] = pr[RED][CIEY]*C_rD/D;
	mat[1][1] = pr[GRN][CIEY]*C_gD/D;
	mat[1][2] = pr[BLU][CIEY]*C_bD/D;
	mat[2][0] = (1.-pr[RED][CIEX]-pr[RED][CIEY])*C_rD/D;
	mat[2][1] = (1.-pr[GRN][CIEX]-pr[GRN][CIEY])*C_gD/D;
	mat[2][2] = (1.-pr[BLU][CIEX]-pr[BLU][CIEY])*C_bD/D;
	return(1);
}


int
comprgb2rgbmat(			/* compute conversion from RGB1 to RGB2 */
COLORMAT  mat,
RGBPRIMS  pr1,
RGBPRIMS  pr2
)
{
	COLORMAT  pr1toxyz, xyztopr2;

	if (pr1 == pr2) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return(1);
	}
	if (!comprgb2xyzmat(pr1toxyz, pr1))
		return(0);
	if (!compxyz2rgbmat(xyztopr2, pr2))
		return(0);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, xyztopr2);
	return(1);
}


int
compxyzWBmat(			/* CIE von Kries transform from wht1 to wht2 */
COLORMAT  mat,
float  wht1[2],
float  wht2[2]
)
{
	COLOR	cw1, cw2;
	if (XYEQ(wht1,wht2)) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return(1);
	}
	if (CEQ(wht1[CIEX],0.) | CEQ(wht1[CIEY],0.))
		return(0);
	cw1[RED] = wht1[CIEX]/wht1[CIEY];
	cw1[GRN] = 1.;
	cw1[BLU] = (1. - wht1[CIEX] - wht1[CIEY])/wht1[CIEY];
	colortrans(cw1, vkmat, cw1);
	if (CEQ(wht2[CIEX],0.) | CEQ(wht2[CIEY],0.))
		return(0);
	cw2[RED] = wht2[CIEX]/wht2[CIEY];
	cw2[GRN] = 1.;
	cw2[BLU] = (1. - wht2[CIEX] - wht2[CIEY])/wht2[CIEY];
	colortrans(cw2, vkmat, cw2);
	if (CEQ(cw1[RED],0.) | CEQ(cw1[GRN],0.) | CEQ(cw1[BLU],0.))
		return(0);
	mat[0][0] = cw2[RED]/cw1[RED];
	mat[1][1] = cw2[GRN]/cw1[GRN];
	mat[2][2] = cw2[BLU]/cw1[BLU];
	mat[0][1] = mat[0][2] = mat[1][0] =
	mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
	multcolormat(mat, vkmat, mat);
	multcolormat(mat, mat, ivkmat);
	return(1);
}


int
compxyz2rgbWBmat(		/* von Kries conversion from CIE to RGB space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	COLORMAT	wbmat;

	if (!compxyz2rgbmat(mat, pr))
		return(0);
	if (XYEQ(pr[WHT],xyneu))
		return(1);
	if (!compxyzWBmat(wbmat, xyneu, pr[WHT]))
		return(0);
	multcolormat(mat, wbmat, mat);
	return(1);
}

int
comprgb2xyzWBmat(		/* von Kries conversion from RGB to CIE space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	COLORMAT	wbmat;
	
	if (!comprgb2xyzmat(mat, pr))
		return(0);
	if (XYEQ(pr[WHT],xyneu))
		return(1);
	if (!compxyzWBmat(wbmat, pr[WHT], xyneu))
		return(0);
	multcolormat(mat, mat, wbmat);
	return(1);
}

int
comprgb2rgbWBmat(		/* von Kries conversion from RGB1 to RGB2 */
COLORMAT  mat,
RGBPRIMS  pr1,
RGBPRIMS  pr2
)
{
	COLORMAT  pr1toxyz, xyztopr2, wbmat;

	if (pr1 == pr2) {
		mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
		mat[0][1] = mat[0][2] = mat[1][0] =
		mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
		return(1);
	}
	if (!comprgb2xyzmat(pr1toxyz, pr1))
		return(0);
	if (!compxyzWBmat(wbmat, pr1[WHT], pr2[WHT]))
		return(0);
	if (!compxyz2rgbmat(xyztopr2, pr2))
		return(0);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, wbmat);
	multcolormat(mat, mat, xyztopr2);
	return(1);
}
