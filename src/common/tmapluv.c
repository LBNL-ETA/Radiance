#ifndef lint
static const char	RCSid[] = "$Id: tmapluv.c,v 3.4 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Routines for tone-mapping LogLuv encoded pixels.
 *
 * Externals declared in tmaptiff.h
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

#define LOGLUV_PUBLIC		1

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tiffio.h"
#include "tmprivat.h"
#include "tmaptiff.h"

#ifndef BSD
#define bzero(d,n)		(void)memset(d,0,n)
#endif
#define uvflgop(p,uv,op)	((p)->rgbflg[(uv)>>5] op (1L<<((uv)&0x1f)))
#define isuvset(p,uv)		uvflgop(p,uv,&)
#define setuv(p,uv)		uvflgop(p,uv,|=)
#define clruv(p,uv)		uvflgop(p,uv,&=~)
#define clruvall(p)		bzero((MEM_PTR)(p)->rgbflg,sizeof((p)->rgbflg))

#ifdef NOPROTO
static MEM_PTR	luv32Init();
static void	luv32NewSpace();
static MEM_PTR	luv24Init();
static void	luv24NewSpace();
#else
static MEM_PTR	luv32Init(struct tmStruct *);
static void	luv32NewSpace(struct tmStruct *);
static MEM_PTR	luv24Init(struct tmStruct *);
static void	luv24NewSpace(struct tmStruct *);
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
uv2rgb(rgb, tm, uvp)			/* compute RGB from uv coordinate */
BYTE	rgb[3];
register struct tmStruct	*tm;
double	uvp[2];
{	/* Should check that tm->inppri==TM_XYZPRIM beforehand... */
	double	d, x, y;
	COLOR	XYZ, RGB;
					/* convert to XYZ */
	d = 1./(6.*uvp[0] - 16.*uvp[1] + 12.);
	x = 9.*uvp[0] * d;
	y = 4.*uvp[1] * d;
	XYZ[CIEY] = 1./tm->inpsf;
	XYZ[CIEX] = x/y * XYZ[CIEY];
	XYZ[CIEZ] = (1.-x-y)/y * XYZ[CIEY];
					/* convert to RGB and clip */
	colortrans(RGB, tm->cmat, XYZ);
	clipgamut(RGB, 1., CGAMUT_LOWER, cblack, cwhite);
					/* perform final scaling & gamma */
	d = tm->clf[RED] * RGB[RED];
	rgb[RED] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tm->mongam));
	d = tm->clf[GRN] * RGB[GRN];
	rgb[GRN] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tm->mongam));
	d = tm->clf[BLU] * RGB[BLU];
	rgb[BLU] = d>=.999 ? 255 : (int)(256.*pow(d, 1./tm->mongam));
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
tmCvLuv32(ls, cs, luvs, len)		/* convert raw 32-bit LogLuv values */
TMbright	*ls;
BYTE	*cs;
uint32	*luvs;
int	len;
{
	static char	funcName[] = "tmCvLuv32";
	double	uvp[2];
	register LUV32DATA	*ld;
	register int	i, j;
					/* check arguments */
	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | luvs == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
					/* check package registration */
	if (luv32Reg < 0) {
		if ((luv32Reg = tmRegPkg(&luv32Pkg)) < 0)
			returnErr(TM_E_CODERR1);
		tmMkMesofact();
	}
					/* get package data */
	if ((ld = (LUV32DATA *)tmPkgData(tmTop,luv32Reg)) == NULL)
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
		if (tmTop->flags & TM_F_MESOPIC && ls[i] < BMESUPPER) {
			uvp[0] = 1./UVSCALE*((luvs[i]>>8 & 0xff) + .5);
			uvp[1] = 1./UVSCALE*((luvs[i] & 0xff) + .5);
			ls[i] = compmeshift(ls[i], uvp);
			j = tmTop->flags&TM_F_BW || ls[i]<BMESLOWER
					? UVNEU
					: (int)(uvp[0]*UVSCALE)<<8
						| (int)(uvp[1]*UVSCALE);
		} else {
			j = tmTop->flags&TM_F_BW ? UVNEU : luvs[i]&0xffff;
		}
		if (!isuvset(ld, j)) {
			uvp[0] = 1./UVSCALE*((j>>8) + .5);
			uvp[1] = 1./UVSCALE*((j & 0xff) + .5);
			uv2rgb(ld->rgbval[j], tmTop, uvp);
			setuv(ld, j);
		}
		cs[3*i  ] = ld->rgbval[j][RED];
		cs[3*i+1] = ld->rgbval[j][GRN];
		cs[3*i+2] = ld->rgbval[j][BLU];
	}
	returnOK;
}


int
tmCvLuv24(ls, cs, luvs, len)		/* convert raw 24-bit LogLuv values */
TMbright	*ls;
BYTE	*cs;
uint32	*luvs;
int	len;
{
	char	funcName[] = "tmCvLuv24";
	double	uvp[2];
	register LUV24DATA	*ld;
	register int	i, j;
					/* check arguments */
	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | luvs == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
					/* check package registration */
	if (luv24Reg < 0) {
		if ((luv24Reg = tmRegPkg(&luv24Pkg)) < 0)
			returnErr(TM_E_CODERR1);
		tmMkMesofact();
	}
					/* get package data */
	if ((ld = (LUV24DATA *)tmPkgData(tmTop,luv24Reg)) == NULL)
		returnErr(TM_E_NOMEM);
					/* convert each pixel */
	for (i = len; i--; ) {
		j = luvs[i] >> 14;		/* get luminance */
		ls[i] = (BRT2SCALE(j) >> 6) - ld->offset;
		if (cs == TM_NOCHROM)		/* no color? */
			continue;
						/* get chrominance */
		if (tmTop->flags & TM_F_MESOPIC && ls[i] < BMESUPPER) {
			if (uv_decode(uvp, uvp+1, luvs[i]&0x3fff) < 0) {
				uvp[0] = U_NEU;		/* should barf? */
				uvp[1] = V_NEU;
			}
			ls[i] = compmeshift(ls[i], uvp);
			if (tmTop->flags&TM_F_BW || ls[i]<BMESLOWER
					|| (j = uv_encode(uvp[0], uvp[1],
						SGILOGENCODE_NODITHER)) < 0)
				j = uv14neu;
		} else {
			j = tmTop->flags&TM_F_BW ? uv14neu : luvs[i]&0x3fff;
		}
		if (!isuvset(ld, j)) {
			if (uv_decode(&uvp[0], &uvp[1], j) < 0) {
				uvp[0] = U_NEU; uvp[1] = V_NEU;
			}
			uv2rgb(ld->rgbval[j], tmTop, uvp);
			setuv(ld, j);
		}
		cs[3*i  ] = ld->rgbval[j][RED];
		cs[3*i+1] = ld->rgbval[j][GRN];
		cs[3*i+2] = ld->rgbval[j][BLU];
	}
	returnOK;
}


int
tmCvL16(ls, l16s, len)			/* convert 16-bit LogL values */
TMbright	*ls;
uint16	*l16s;
int	len;
{
	static char	funcName[] = "tmCvL16";
	static double	lastsf;
	static int	offset;
	register int	i;
					/* check arguments */
	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | l16s == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
					/* check scaling offset */
	if (!FEQ(tmTop->inpsf, lastsf)) {
		offset = BRT2SCALE(64);
		if (tmTop->inpsf > 1.0001)
			offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)+.5);
		else if (tmTop->inpsf < 0.9999)
			offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)-.5);
		lastsf = tmTop->inpsf;
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
luv32NewSpace(tm)		/* initialize 32-bit LogLuv color space */
struct tmStruct	*tm;
{
	register LUV32DATA	*ld;

	if (tm->inppri != TM_XYZPRIM) {		/* panic time! */
		fputs("Improper input color space in luv32NewSpace!\n", stderr);
		exit(1);
	}
	ld = (LUV32DATA *)tm->pd[luv32Reg];
	ld->offset = BRT2SCALE(64);
	if (tm->inpsf > 1.0001)
		ld->offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)+.5);
	else if (tm->inpsf < 0.9999)
		ld->offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)-.5);
	clruvall(ld);
}


static MEM_PTR
luv32Init(tm)			/* allocate data for 32-bit LogLuv decoder */
struct tmStruct	*tm;
{
	register LUV32DATA	*ld;

	ld = (LUV32DATA *)malloc(sizeof(LUV32DATA));
	if (ld == NULL)
		return(NULL);
	tm->pd[luv32Reg] = (MEM_PTR)ld;
	luv32NewSpace(tm);
	return((MEM_PTR)ld);
}


static void
luv24NewSpace(tm)		/* initialize 24-bit LogLuv color space */
struct tmStruct	*tm;
{
	register LUV24DATA	*ld;

	if (tm->inppri != TM_XYZPRIM) {		/* panic time! */
		fputs("Improper input color space in luv24NewSpace!\n", stderr);
		exit(1);
	}
	ld = (LUV24DATA *)tm->pd[luv24Reg];
	ld->offset = BRT2SCALE(12);
	if (tm->inpsf > 1.0001)
		ld->offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)+.5);
	else if (tm->inpsf < 0.9999)
		ld->offset -= (int)(TM_BRTSCALE*log(tmTop->inpsf)-.5);
	clruvall(ld);
}


static MEM_PTR
luv24Init(tm)			/* allocate data for 24-bit LogLuv decoder */
struct tmStruct	*tm;
{
	register LUV24DATA	*ld;

	ld = (LUV24DATA *)malloc(sizeof(LUV24DATA));
	if (ld == NULL)
		return(NULL);
	tm->pd[luv24Reg] = (MEM_PTR)ld;
	if (uv14neu < 0)		/* initialize neutral color index */
		uv14neu = uv_encode(U_NEU, V_NEU, SGILOGENCODE_NODITHER);
	luv24NewSpace(tm);
	return((MEM_PTR)ld);
}
