#ifndef lint
static const char RCSid[] = "$Id: ccolor.c,v 3.11 2017/04/13 00:42:01 greg Exp $";
#endif
/*
 * Spectral color handling routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ccolor.h"

#undef frand
#define frand()	(rand()*(1./(RAND_MAX+.5)))

					/* Sharp primary matrix */
float	XYZtoSharp[3][3] = {
	{ 1.2694, -0.0988, -0.1706},
	{-0.8364,  1.8006,  0.0357},
	{ 0.0297, -0.0315,  1.0018}
};
					/* inverse Sharp primary matrix */
float	XYZfromSharp[3][3] = {
	{ 0.8156,  0.0472,  0.1372},
	{ 0.3791,  0.5769,  0.0440},
	{-0.0123,  0.0167,  0.9955}
};

const C_COLOR	c_dfcolor = C_DEFCOLOR;

const C_CHROMA	c_dfchroma = 49750;	/* c_encodeChroma(&c_dfcolor) */

				/* CIE 1931 Standard Observer curves */
const C_COLOR	c_x31 = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY|C_CSEFF,
			{14,42,143,435,1344,2839,3483,3362,2908,1954,956,
			320,49,93,633,1655,2904,4334,5945,7621,9163,10263,
			10622,10026,8544,6424,4479,2835,1649,874,468,227,
			114,58,29,14,7,3,2,1,0}, 106836L, .467, .368, 362.230
			};
const C_COLOR	c_y31 = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY|C_CSEFF,
			{0,1,4,12,40,116,230,380,600,910,1390,2080,3230,
			5030,7100,8620,9540,9950,9950,9520,8700,7570,6310,
			5030,3810,2650,1750,1070,610,320,170,82,41,21,10,
			5,2,1,1,0,0}, 106856L, .398, .542, 493.525
			};
const C_COLOR	c_z31 = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY|C_CSEFF,
			{65,201,679,2074,6456,13856,17471,17721,16692,
			12876,8130,4652,2720,1582,782,422,203,87,39,21,17,
			11,8,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			106770L, .147, .077, 54.363
			};
				/* Derived CIE 1931 Primaries (imaginary) */
static const C_COLOR	cie_xp = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY,
			{-174,-198,-195,-197,-202,-213,-235,-272,-333,
			-444,-688,-1232,-2393,-4497,-6876,-6758,-5256,
			-3100,-815,1320,3200,4782,5998,6861,7408,7754,
			7980,8120,8199,8240,8271,8292,8309,8283,8469,
			8336,8336,8336,8336,8336,8336},
			127424L, 1., .0,
			};
static const C_COLOR	cie_yp = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY,
			{-451,-431,-431,-430,-427,-417,-399,-366,-312,
			-204,57,691,2142,4990,8810,9871,9122,7321,5145,
			3023,1123,-473,-1704,-2572,-3127,-3474,-3704,
			-3846,-3927,-3968,-3999,-4021,-4038,-4012,-4201,
			-4066,-4066,-4066,-4066,-4066,-4066},
			-23035L, .0, 1.,
			};
static const C_COLOR	cie_zp = { 1, NULL, C_CDSPEC|C_CSSPEC|C_CSXY,
			{4051,4054,4052,4053,4054,4056,4059,4064,4071,
			4074,4056,3967,3677,2933,1492,313,-440,-795,
			-904,-918,-898,-884,-869,-863,-855,-855,-851,
			-848,-847,-846,-846,-846,-845,-846,-843,-845,
			-845,-845,-845,-845,-845},
			36057L, .0, .0,
			};


/* convert to sharpened RGB color for low-error operations */
void
c_toSharpRGB(C_COLOR *cin, double cieY, float cout[3])
{
	double	xyz[3];

	c_ccvt(cin, C_CSXY);

	xyz[0] = cin->cx/cin->cy * cieY;
	xyz[1] = cieY;
	xyz[2] = (1. - cin->cx - cin->cy)/cin->cy * cieY;

	cout[0] = XYZtoSharp[0][0]*xyz[0] + XYZtoSharp[0][1]*xyz[1] +
				XYZtoSharp[0][2]*xyz[2];
	cout[1] = XYZtoSharp[1][0]*xyz[0] + XYZtoSharp[1][1]*xyz[1] +
				XYZtoSharp[1][2]*xyz[2];
	cout[2] = XYZtoSharp[2][0]*xyz[0] + XYZtoSharp[2][1]*xyz[1] +
				XYZtoSharp[2][2]*xyz[2];
}

/* convert back from sharpened RGB color */
double
c_fromSharpRGB(float cin[3], C_COLOR *cout)
{
	double	xyz[3], sf;

	xyz[1] = XYZfromSharp[1][0]*cin[0] + XYZfromSharp[1][1]*cin[1] +
				XYZfromSharp[1][2]*cin[2];
	if (xyz[1] <= 1e-6) {
		*cout = c_dfcolor;	/* punting, here... */
		return xyz[1];
	}
	xyz[0] = XYZfromSharp[0][0]*cin[0] + XYZfromSharp[0][1]*cin[1] +
				XYZfromSharp[0][2]*cin[2];
	xyz[2] = XYZfromSharp[2][0]*cin[0] + XYZfromSharp[2][1]*cin[1] +
				XYZfromSharp[2][2]*cin[2];
	
	sf = 1./(xyz[0] + xyz[1] + xyz[2]);

	cout->cx = xyz[0] * sf;
	cout->cy = xyz[1] * sf;
	cout->flags = C_CDXY|C_CSXY;
	
	return(xyz[1]);
}

/* assign arbitrary spectrum and return Y value */
double
c_sset(C_COLOR *clr, double wlmin, double wlmax, const float spec[], int nwl)
{
	double	yval, scale;
	float	va[C_CNSS];
	int	i, pos, n, imax, wl;
	double	wl0, wlstep;
	double	boxpos, boxstep;
					/* check arguments */
	if ((nwl <= 1) | (spec == NULL) | (wlmin >= C_CMAXWL) |
			(wlmax <= C_CMINWL) | (wlmin >= wlmax))
		return(0.);
	wlstep = (wlmax - wlmin)/(nwl-1);
	while (wlmin < C_CMINWL) {
		wlmin += wlstep;
		--nwl; ++spec;
	}
	while (wlmax > C_CMAXWL) {
		wlmax -= wlstep;
		--nwl;
	}
	if ((nwl <= 1) | (wlmin >= wlmax))
		return(0.);
	imax = nwl;			/* box filter if necessary */
	boxpos = 0;
	boxstep = 1;
	if (wlstep < C_CWLI) {
		imax = (wlmax - wlmin)/C_CWLI + 1e-7;
		boxpos = (wlmin - C_CMINWL)/C_CWLI;
		boxstep = wlstep/C_CWLI;
		wlstep = C_CWLI;
	}
	scale = 0.;			/* get values and maximum */
	yval = 0.;
	pos = 0;
	for (i = 0; i < imax; i++) {
		va[i] = 0.; n = 0;
		while (boxpos < i+.5 && pos < nwl) {
			va[i] += spec[pos++];
			n++;
			boxpos += boxstep;
		}
		if (n > 1)
			va[i] /= (double)n;
		if (va[i] > scale)
			scale = va[i];
		else if (va[i] < -scale)
			scale = -va[i];
		yval += va[i] * c_y31.ssamp[i];
	}
	if (scale <= 1e-7)
		return(0.);
	yval /= (double)c_y31.ssum;
	scale = C_CMAXV / scale;
	clr->ssum = 0;			/* convert to our spacing */
	wl0 = wlmin;
	pos = 0;
	for (i = 0, wl = C_CMINWL; i < C_CNSS; i++, wl += C_CWLI)
		if ((wl < wlmin) | (wl > wlmax))
			clr->ssamp[i] = 0;
		else {
			while (wl0 + wlstep < wl+1e-7) {
				wl0 += wlstep;
				pos++;
			}
			if ((wl+1e-7 >= wl0) & (wl-1e-7 <= wl0))
				clr->ssamp[i] = scale*va[pos] + frand();
			else		/* interpolate if necessary */
				clr->ssamp[i] = frand() + scale / wlstep *
						( va[pos]*(wl0+wlstep - wl) +
							va[pos+1]*(wl - wl0) );
			clr->ssum += clr->ssamp[i];
		}
	clr->flags = C_CDSPEC|C_CSSPEC;
	return(yval);
}

/* check if color is grey */
int
c_isgrey(C_COLOR *clr)
{
	if (!(clr->flags & (C_CSXY|C_CSSPEC)))
		return(1);		/* no settings == grey */

	c_ccvt(clr, C_CSXY);

	return(clr->cx >= .323 && clr->cx <= .343 &&
			clr->cy >= .323 && clr->cy <= .343);
}

/* convert color representations */
void
c_ccvt(C_COLOR *clr, int fl)
{
	double	x, y, z;
	int	i;

	fl &= ~clr->flags;			/* ignore what's done */
	if (!fl)				/* everything's done! */
		return;
	if (!(clr->flags & (C_CSXY|C_CSSPEC))) {
		*clr = c_dfcolor;		/* nothing was set! */
		return;
	}
	if (fl & C_CSXY) {		/* cspec -> cxy */
		x = y = z = 0.;
		for (i = 0; i < C_CNSS; i++) {
			x += c_x31.ssamp[i] * clr->ssamp[i];
			y += c_y31.ssamp[i] * clr->ssamp[i];
			z += c_z31.ssamp[i] * clr->ssamp[i];
		}
		x /= (double)c_x31.ssum;
		y /= (double)c_y31.ssum;
		z /= (double)c_z31.ssum;
		z += x + y;
		if (z > 1e-6) {
			clr->cx = x / z;
			clr->cy = y / z;
		} else
			clr->cx = clr->cy = 1./3.;
		clr->flags |= C_CSXY;
	}
	if (fl & C_CSSPEC) {		/* cxy -> cspec */
		x = clr->cx;
		y = clr->cy;
		z = 1. - x - y;
		clr->ssum = 0;
		for (i = 0; i < C_CNSS; i++) {
			clr->ssamp[i] = x*cie_xp.ssamp[i] + y*cie_yp.ssamp[i]
					+ z*cie_zp.ssamp[i] + frand();
			if (clr->ssamp[i] < 0)		/* out of gamut! */
				clr->ssamp[i] = 0;
			else
				clr->ssum += clr->ssamp[i];
		}
		clr->flags |= C_CSSPEC;
	}
	if (fl & C_CSEFF) {		/* compute efficacy */
		if (clr->flags & C_CSSPEC) {		/* from spectrum */
			y = 0.;
			for (i = 0; i < C_CNSS; i++)
				y += c_y31.ssamp[i] * clr->ssamp[i];
			clr->eff = C_CLPWM * y / (clr->ssum + 0.0001);
		} else /* clr->flags & C_CSXY */ {	/* from (x,y) */
			clr->eff = clr->cx*c_x31.eff + clr->cy*c_y31.eff +
					(1. - clr->cx - clr->cy)*c_z31.eff;
		}
		clr->flags |= C_CSEFF;
	}
}

/* mix two colors according to weights given -- cres may be c1 or c2 */
void
c_cmix(C_COLOR *cres, double w1, C_COLOR *c1, double w2, C_COLOR *c2)
{
	double	scale;
	int	i;

	if (w1 == 0) {
		*cres = *c2;
		return;
	}
	if (w2 == 0) {
		*cres = *c1;
		return;
	}
	if ((c1->flags|c2->flags) & C_CDSPEC) {		/* spectral mixing */
		float	cmix[C_CNSS];

		c_ccvt(c1, C_CSSPEC|C_CSEFF);
		c_ccvt(c2, C_CSSPEC|C_CSEFF);
		if (c1->ssum*c2->ssum == 0) {
			*cres = c_dfcolor;
			return;
		}
		w1 /= c1->eff*c1->ssum;
		w2 /= c2->eff*c2->ssum;
		scale = 0.;
		for (i = 0; i < C_CNSS; i++) {
			cmix[i] = w1*c1->ssamp[i] + w2*c2->ssamp[i];
			if (cmix[i] > scale)
				scale = cmix[i];
			else if (cmix[i] < -scale)
				scale = -cmix[i];
		}
		scale = C_CMAXV / scale;
		cres->ssum = 0;
		for (i = 0; i < C_CNSS; i++)
			cres->ssum += cres->ssamp[i] = scale*cmix[i] + frand();
		cres->flags = C_CDSPEC|C_CSSPEC;
	} else {					/* CIE xy mixing */
		c_ccvt(c1, C_CSXY);
		c_ccvt(c2, C_CSXY);
		w1 *= c2->cy;
		w2 *= c1->cy;
		scale = w1 + w2;
		if (scale == 0.) {
			*cres = c_dfcolor;
			return;
		}
		scale = 1. / scale;
		cres->cx = (w1*c1->cx + w2*c2->cx) * scale;
		cres->cy = (w1*c1->cy + w2*c2->cy) * scale;
		cres->flags = C_CDXY|C_CSXY;
	}
}

/* multiply two colors -- cres may be c1 or c2 */
double
c_cmult(C_COLOR *cres, C_COLOR *c1, double y1, C_COLOR *c2, double y2)
{
	double	yres;
	int	i;

	if ((c1->flags|c2->flags) & C_CDSPEC) {
		long	cmix[C_CNSS], cmax = 0;		/* spectral multiply */

		c_ccvt(c1, C_CSSPEC|C_CSEFF);
		c_ccvt(c2, C_CSSPEC|C_CSEFF);
		for (i = 0; i < C_CNSS; i++) {
			cmix[i] = c1->ssamp[i] * c2->ssamp[i];
			if (cmix[i] > cmax)
				cmax = cmix[i];
			else if (cmix[i] < -cmax)
				cmax = -cmix[i];
		}
		cmax /= C_CMAXV;
		if (!cmax) {
			*cres = c_dfcolor;
			return(0.);
		}
		cres->ssum = 0;
		for (i = 0; i < C_CNSS; i++)
			cres->ssum += cres->ssamp[i] = cmix[i] / cmax;
		cres->flags = C_CDSPEC|C_CSSPEC;

		c_ccvt(cres, C_CSEFF);			/* below is correct */
		yres = y1 * y2 * c_y31.ssum * C_CLPWM /
			(c1->eff*c1->ssum * c2->eff*c2->ssum) *
			cres->eff*( cres->ssum*(double)cmax +
						C_CNSS/2.0*(cmax-1) );
	} else {
		float	rgb1[3], rgb2[3];		/* CIE xy multiply */

		c_toSharpRGB(c1, y1, rgb1);
		c_toSharpRGB(c2, y2, rgb2);
		
		rgb2[0] *= rgb1[0]; rgb2[1] *= rgb1[1]; rgb2[2] *= rgb1[2];

		yres = c_fromSharpRGB(rgb2, cres);
	}
	return(yres);
}

#define	C1		3.741832e-16	/* W-m^2 */
#define C2		1.4388e-2	/* m-K */

#define bbsp(l,t)	(C1/((l)*(l)*(l)*(l)*(l)*(exp(C2/((t)*(l)))-1.)))
#define bblm(t)		(C2*0.2/(t))

/* set black body spectrum */
int
c_bbtemp(C_COLOR *clr, double tk)
{
	double	sf, wl;
	int	i;

	if (tk < 1000)
		return(0);
	wl = bblm(tk);			/* scalefactor based on peak */
	if (wl < C_CMINWL*1e-9)
		wl = C_CMINWL*1e-9;
	else if (wl > C_CMAXWL*1e-9)
		wl = C_CMAXWL*1e-9;
	sf = C_CMAXV/bbsp(wl,tk);
	clr->ssum = 0;
	for (i = 0; i < C_CNSS; i++) {
		wl = (C_CMINWL + i*C_CWLI)*1e-9;
		clr->ssum += clr->ssamp[i] = sf*bbsp(wl,tk) + frand();
	}
	clr->flags = C_CDSPEC|C_CSSPEC;
	return(1);
}

#undef	C1
#undef	C2
#undef	bbsp
#undef	bblm

#define UV_NORMF	410.

/* encode (x,y) chromaticity */
C_CHROMA
c_encodeChroma(C_COLOR *clr)
{
	double	df;
	int	ub, vb;

	c_ccvt(clr, C_CSXY);
	df = UV_NORMF/(-2.*clr->cx + 12.*clr->cy + 3.);
	ub = 4.*clr->cx*df + frand();
	if (ub > 0xff) ub = 0xff;
	else ub *= (ub > 0);
	vb = 9.*clr->cy*df + frand();
	if (vb > 0xff) vb = 0xff;
	else vb *= (vb > 0);

	return(vb<<8 | ub);
}

/* decode (x,y) chromaticity */
void
c_decodeChroma(C_COLOR *cres, C_CHROMA ccode)
{
	double	up = (ccode & 0xff)*(1./UV_NORMF);
	double	vp = (ccode>>8 & 0xff)*(1./UV_NORMF);
	double	df = 1./(6.*up - 16.*vp + 12.);

	cres->cx = 9.*up * df;
	cres->cy = 4.*vp * df;
	cres->flags = C_CDXY|C_CSXY;
}

#undef	UV_NORMF
