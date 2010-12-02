#ifndef lint
static const char	RCSid[] = "$Id: tmapluv.c,v 3.15 2010/12/02 18:10:40 greg Exp $";
#endif
/*
 * Routines for tone-mapping LogLuv encoded pixels.
 *
 * Externals declared in tmaptiff.h
 */

#include "copyright.h"

#define LOGLUV_PUBLIC		1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tiffio.h"
#include "tmprivat.h"
#include "tmaptiff.h"

#define uvflgop(p,uv,op)	((p)->rgbflg[(uv)>>5] op (1L<<((uv)&0x1f)))
#define isuvset(p,uv)		uvflgop(p,uv,&)
#define setuv(p,uv)		uvflgop(p,uv,|=)
#define clruv(p,uv)		uvflgop(p,uv,&=~)
#define clruvall(p)		memset((MEM_PTR)(p)->rgbflg,'\0',sizeof((p)->rgbflg))

#ifdef NOPROTO
static MEM_PTR	luv32Init();
static void	luv32NewSpace();
static MEM_PTR	luv24Init();
static void	luv24NewSpace();
#else
static MEM_PTR	luv32Init(TMstruct *);
static void	luv32NewSpace(TMstruct *);
static MEM_PTR	luv24Init(TMstruct *);
static void	luv24NewSpace(TMstruct *);
#endif

typedef struct {
	int	offset;			/* computed luminance offset */
	BYTE	rgbval[1<<16][3];	/* computed RGB value for given uv */
	uint32	rgbflg[1<<(16-5)];	/* flags for computed values */
} LUV32DATA;			/* LogLuv 32-bit conversion data */

#define UVNEU			((int)(UVSCALE*U_NEU)<<8 \
					| (int)(UVSCALE*V_NEU))

static struct tmPackage  luv32Pkg = {	/* 32-bit package functions */
	luv32Init, luv32NewSpace, free
};
static int	luv32Reg = -1;		/* 32-bit package reg. number */

typedef struct {
	int	offset;			/* computed luminance offset */
	BYTE	rgbval[1<<14][3];	/* computed rgb value for uv index */
	uint32	rgbflg[1<<(14-5)];	/* flags for computed values */
} LUV24DATA;			/* LogLuv 24-bit conversion data */

static struct tmPackage  luv24Pkg = {	/* 24-bit package functions */
	luv24Init, luv24NewSpace, free
};
static int	luv24Reg = -1;		/* 24-bit package reg. number */

static int	uv14neu = -1;		/* neutral index for 14-bit (u',v') */


static void
uv2rgb(rgb, tms, uvp)			/* compute RGB from uv coordinate */
BYTE	rgb[3];
register TMstruct	*tms;
double	uvp[2];
{	/* Should check that tms->inppri==TM_XYZPRIM beforehand... */
	double	d, x, y;
	COLOR	XYZ, RGB;
					/* convert to XYZ */
	d = 1./(6.*uvp[0] - 16.*uvp[1] + 12.);
	x = 9.*uvp[0] * d;
	y = 4.*uvp[1] * d;
	XYZ[CIEY] = 1./tms->inpsf;
	XYZ[CIEX] = x/y * XYZ[CIEY];
	XYZ[CIEZ] = (1.-x-y)/y * XYZ[CIEY];
					/* convert to RGB and clip */
	colortrans(RGB, tms->cmat, XYZ);
	clipgamut(RGB, 1., CGAMUT_LOWER, cblack, cwhite);
					/* perform final scaling & gamma */
	d = tms->clf[RED] * RGB[RED];
	rgb[RED] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tms->mongam));
	d = tms->clf[GRN] * RGB[GRN];
	rgb[GRN] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tms->mongam));
	d = tms->clf[BLU] * RGB[BLU];
	rgb[BLU] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tms->mongam));
}


static TMbright
compmeshift(li, uvp)			/* compute mesopic color shift */
TMbright	li;		/* encoded world luminance */
double	uvp[2];			/* world (u',v') -> returned desaturated */
{
	double	scotrat;
	register double	d;

	if (li >= BMESUPPER)
		return(li);
	scotrat = (.767676768 - 1.02356902*uvp[1])/uvp[0] - .343434343;
	if (li <= BMESLOWER) {
		d = 0.;
		uvp[0] = U_NEU; uvp[1] = V_NEU;
	} else {
		d = (tmMesofact[li-BMESLOWER] + .5) * (1./256.);
		uvp[0] = d*uvp[0] + (1.-d)*U_NEU;
		uvp[1] = d*uvp[1] + (1.-d)*V_NEU;
	}
	/*
	d = li + (double)TM_BRTSCALE*log(d + (1.-d)*scotrat);
	*/
	d = d + (1.-d)*scotrat;
	d -= 1.;			/* Taylor expansion of log(x) about 1 */
	d = d*(1. + d*(-.5 + d*(1./3. + d*-.125)));
	d = li + (double)TM_BRTSCALE*d;
	return((TMbright)(d>0. ? d+.5 : d-.5));
}


int
tmCvLuv32(				/* convert raw 32-bit LogLuv values */
TMstruct	*tms,
TMbright	*ls,
BYTE	*cs,
uint32	*luvs,
int	len
)
{
	static const char	funcName[] = "tmCvLuv32";
	double	uvp[2];
	register LUV32DATA	*ld;
	register int	i, j;
					/* check arguments */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (luvs == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
					/* check package registration */
	if (luv32Reg < 0) {
		if ((luv32Reg = tmRegPkg(&luv32Pkg)) < 0)
			returnErr(TM_E_CODERR1);
		tmMkMesofact();
	}
					/* get package data */
	if ((ld = (LUV32DATA *)tmPkgData(tms,luv32Reg)) == NULL)
		returnErr(TM_E_NOMEM);
					/* convert each pixel */
	for (i = len; i--; ) {
		j = luvs[i] >> 16;		/* get luminance */
		if (j & 0x8000)			/* negative luminance */
			ls[i] = TM_NOBRT;	/* assign bogus value */
		else				/* else convert to lnL */
			ls[i] = (BRT2SCALE(j) >> 8) - ld->offset;
		if (cs == TM_NOCHROM)		/* no color? */
			continue;
						/* get chrominance */
		if (tms->flags & TM_F_MESOPIC && ls[i] < BMESUPPER) {
			uvp[0] = 1./UVSCALE*((luvs[i]>>8 & 0xff) + .5);
			uvp[1] = 1./UVSCALE*((luvs[i] & 0xff) + .5);
			ls[i] = compmeshift(ls[i], uvp);
			j = tms->flags&TM_F_BW || ls[i]<BMESLOWER
					? UVNEU
					: (int)(uvp[0]*UVSCALE)<<8
						| (int)(uvp[1]*UVSCALE);
		} else {
			j = tms->flags&TM_F_BW ? UVNEU : luvs[i]&0xffff;
		}
		if (!isuvset(ld, j)) {
			uvp[0] = 1./UVSCALE*((j>>8) + .5);
			uvp[1] = 1./UVSCALE*((j & 0xff) + .5);
			uv2rgb(ld->rgbval[j], tms, uvp);
			setuv(ld, j);
		}
		cs[3*i  ] = ld->rgbval[j][RED];
		cs[3*i+1] = ld->rgbval[j][GRN];
		cs[3*i+2] = ld->rgbval[j][BLU];
	}
	returnOK;
}


int
tmCvLuv24(			/* convert raw 24-bit LogLuv values */
TMstruct	*tms,
TMbright	*ls,
BYTE	*cs,
uint32	*luvs,
int	len
)
{
	char	funcName[] = "tmCvLuv24";
	double	uvp[2];
	register LUV24DATA	*ld;
	register int	i, j;
					/* check arguments */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (luvs == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
					/* check package registration */
	if (luv24Reg < 0) {
		if ((luv24Reg = tmRegPkg(&luv24Pkg)) < 0)
			returnErr(TM_E_CODERR1);
		tmMkMesofact();
	}
					/* get package data */
	if ((ld = (LUV24DATA *)tmPkgData(tms,luv24Reg)) == NULL)
		returnErr(TM_E_NOMEM);
					/* convert each pixel */
	for (i = len; i--; ) {
		j = luvs[i] >> 14;		/* get luminance */
		ls[i] = (BRT2SCALE(j) >> 6) - ld->offset;
		if (cs == TM_NOCHROM)		/* no color? */
			continue;
						/* get chrominance */
		if (tms->flags & TM_F_MESOPIC && ls[i] < BMESUPPER) {
			if (uv_decode(&uvp[0], &uvp[1], luvs[i]&0x3fff) < 0) {
				uvp[0] = U_NEU;		/* should barf? */
				uvp[1] = V_NEU;
			}
			ls[i] = compmeshift(ls[i], uvp);
			if (tms->flags&TM_F_BW || ls[i]<BMESLOWER
					|| (j = uv_encode(uvp[0], uvp[1],
						SGILOGENCODE_NODITHER)) < 0)
				j = uv14neu;
		} else {
			j = tms->flags&TM_F_BW ? uv14neu :
					(int)(luvs[i]&0x3fff);
		}
		if (!isuvset(ld, j)) {
			if (uv_decode(&uvp[0], &uvp[1], j) < 0) {
				uvp[0] = U_NEU; uvp[1] = V_NEU;
			}
			uv2rgb(ld->rgbval[j], tms, uvp);
			setuv(ld, j);
		}
		cs[3*i  ] = ld->rgbval[j][RED];
		cs[3*i+1] = ld->rgbval[j][GRN];
		cs[3*i+2] = ld->rgbval[j][BLU];
	}
	returnOK;
}


int
tmCvL16(				/* convert 16-bit LogL values */
TMstruct	*tms,
TMbright	*ls,
uint16	*l16s,
int	len
)
{
	static const char	funcName[] = "tmCvL16";
	static double	lastsf;
	static int	offset;
	register int	i;
					/* check arguments */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (l16s == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
					/* check scaling offset */
	if (!FEQ(tms->inpsf, lastsf)) {
		offset = BRT2SCALE(64) - tmCvLuminance(tms->inpsf);
		lastsf = tms->inpsf;
	}
					/* convert each pixel */
	for (i = len; i--; ) {
		if (l16s[i] & 0x8000)		/* negative luminance */
			ls[i] = TM_NOBRT;	/* assign bogus value */
		else				/* else convert to lnL */
			ls[i] = (BRT2SCALE(l16s[i]) >> 8) - offset;
	}
	returnOK;
}


static void
luv32NewSpace(tms)		/* initialize 32-bit LogLuv color space */
TMstruct	*tms;
{
	register LUV32DATA	*ld;

	if (tms->inppri != TM_XYZPRIM) {		/* panic time! */
		fputs("Improper input color space in luv32NewSpace!\n", stderr);
		exit(1);
	}
	ld = (LUV32DATA *)tms->pd[luv32Reg];
	ld->offset = BRT2SCALE(64) - tmCvLuminance(tms->inpsf);
	clruvall(ld);
}


static MEM_PTR
luv32Init(tms)			/* allocate data for 32-bit LogLuv decoder */
TMstruct	*tms;
{
	register LUV32DATA	*ld;

	ld = (LUV32DATA *)malloc(sizeof(LUV32DATA));
	if (ld == NULL)
		return(NULL);
	tms->pd[luv32Reg] = (MEM_PTR)ld;
	luv32NewSpace(tms);
	return((MEM_PTR)ld);
}


static void
luv24NewSpace(tms)		/* initialize 24-bit LogLuv color space */
TMstruct *tms;
{
	register LUV24DATA	*ld;

	if (tms->inppri != TM_XYZPRIM) {		/* panic time! */
		fputs("Improper input color space in luv24NewSpace!\n", stderr);
		exit(1);
	}
	ld = (LUV24DATA *)tms->pd[luv24Reg];
	ld->offset = BRT2SCALE(12) - tmCvLuminance(tms->inpsf);
	clruvall(ld);
}


static MEM_PTR
luv24Init(tms)			/* allocate data for 24-bit LogLuv decoder */
TMstruct	*tms;
{
	register LUV24DATA	*ld;

	ld = (LUV24DATA *)malloc(sizeof(LUV24DATA));
	if (ld == NULL)
		return(NULL);
	tms->pd[luv24Reg] = (MEM_PTR)ld;
	if (uv14neu < 0)		/* initialize neutral color index */
		uv14neu = uv_encode(U_NEU, V_NEU, SGILOGENCODE_NODITHER);
	luv24NewSpace(tms);
	return((MEM_PTR)ld);
}
