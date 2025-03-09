/* RCSid $Id: color.h,v 2.52 2025/03/09 19:11:51 greg Exp $ */
/*
 *  color.h - header for routines using pixel color and spectral values
 *		(Notes by Randolph Fritz)
 *
 *  COLOR REPRESENTATION OVERVIEW
 *  =============================
 *  Internally, Radiance represents light in multiple spectral
 *  bands. Four spectral models are used: monochrome, RGB, XYZ,
 *  and multiband.
 *
 *  Units
 *  -----
 *  Radiance -- W/sr/m2
 *  Irradiance -- W/m2
 *  Luminance -- lm/sr/m2
 *  Illuminance -- lm/m2
 *
 *  Colors
 *  ------
 *  In the monochrome, RGB, and multiband formats, units of
 *  radiance and irradiance are used. In the XYZ format, units
 *  of luminance and illuminance are used, or sometimes an
 *  intermediate form, where the units are multiplied by the
 *  constant WHTEFFICACY, the scotopic luminous efficacy of
 *  white light. WHTEFFICACY is 179, an approximation of the
 *  luminous efficacy of the equal energy spectrum.
 *
 *  In the multiband format, up to MAXCSAMP (by default 24)
 *  spectral bands may be used. 24 was chosen after testing,
 *  which found virtually no benefit for more than 18 bands. 24
 *  was chosen to allow for additional infrared or ultraviolet
 *  bands. If even that is not enough, MAXCSAMP can be increased
 *  at compile time.
 *
 *  Numbers
 *  -------
 *  Calculations are done using 32-bit floating point numbers.
 *  Values are stored in the compressed real pixels[1] format,
 *  which allocates one byte per band, plus one for an exponent.
 *
 *  Real Pixel Format
 *  -----------------
 *  The real pixel format gives at most eight bits of floating
 *  point precision to each spectral band; the brightest bands
 *  have full precision, darker ones less.
 *
 *  In the real pixel format, each spectral band is allotted a
 *  one byte mantissa and a common single-byte exponent is used.
 *  The exponent has the range [-128,127] and each mantissa is
 *  assumed to have a binary point at the left, so they have the
 *  range [0,255/256]. In addition, the mantissas are
 *  normalized, so that at least one mantissa always is in the
 *  range [128/256, 255/256] -- one mantissa will always have
 *  its high-order bit set.
 *
 *  References
 *  ----------
 *  [1] Ward, Greg. "Real Pixels." In Graphics Gems II, edited
 *  by Arvo, James, 80--83. Graphics Gems Series. Boston:
 *  Academic Press, 1991.
 *
 *
 *  IMPLEMENTATION DETAILS
 *  ======================
 *  Two color representations are used, one for calculation and
 *  another for storage.  Calculation is done with an array of 32-bit
 *  floats for speed.  Stored color values use single byte mantissas
 *  and a common exponent.  By convention, types containing 32-bit
 *  floats are denoted by COLOR and compressed types are denoted by
 *  COLR.  Tristimulus -- RGB or XYZ -- values can be stored in a
 *  single 32-bit word.
 *
 *  Spectral colors have between 3 and MAXCSAMP samples, and cover
 *  wavelengths from WLPART[0] to WLPART[3] (max to min nanometers).
 *  Wavelengths WLPART[1] and WLPART[2] mark the separation
 *  of red/green and green/blue intervals, respectively.
 *  The wavelength range may go well beyond visible in either
 *  direction, but some routines will average over each interval
 *  and designate the means as R, G, and B, regardless.
 *  The number of samples is set in CNDX[3], and CNDX[0,1,2]
 *  give peak wavelength indices corresponding to stdprims.
 *  For best accuracy, internal calculations are promoted to
 *  the current number of samples and final conversions
 *  to tristimulus should use scolor_rgb() or scolor_cie().
 *
 *  A Radiance file format is provided for spectral pictures, and
 *  spectral colors must be converted by caller if the sampling
 *  doesn't match.
 */
#ifndef _RAD_COLOR_H_
#define _RAD_COLOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAXCSAMP
#define MAXCSAMP	24	/* maximum # spectral samples */
#endif

/* Subscripts for tristimulus colors */
#define  RED		0
#define  GRN		1
#define  BLU		2
#define  CIEX		0	/* or, if input is XYZ... */
#define  CIEY		1
#define  CIEZ		2
#define  EXP		3	/* exponent same for either format */
#define  COLXS		128	/* excess used for exponent */
#define  WHT		3	/* used for RGBPRIMS type */

#undef uby8
#define uby8  unsigned char	/* 8-bit unsigned integer */

typedef uby8 COLRV;
typedef COLRV COLR[4];		/* red, green, blue (or X,Y,Z), exponent */
typedef COLRV SCOLR[MAXCSAMP+1];	/* spectral color, common exponent */

typedef float COLORV;
typedef COLORV  COLOR[3];	/* red, green, blue (or X,Y,Z) */
typedef COLORV	SCOLOR[MAXCSAMP];	/* spectral color */

typedef float  RGBPRIMS[4][2];	/* (x,y) chromaticities for RGBW */
typedef float  (*RGBPRIMP)[2];	/* pointer to RGBPRIMS array */

typedef float  COLORMAT[3][3];	/* color coordinate conversion matrix */

#define  copycolr(c1,c2)	(c1[0]=c2[0],c1[1]=c2[1], \
				c1[2]=c2[2],c1[3]=c2[3])

#define  colval(col,pri)	(col)[pri]

#define  setcolor(col,r,g,b)	((col)[RED]=(r),(col)[GRN]=(g),(col)[BLU]=(b))

#define  scalecolor(col,sf)	((col)[0]*=(sf),(col)[1]*=(sf),(col)[2]*=(sf))

#define  opcolor(c1,op,c2)	((c1)[0]op(c2)[0],(c1)[1]op(c2)[1],(c1)[2]op(c2)[2])

#define  copycolor(c1,c2)	opcolor(c1,=,c2)

#define  addcolor(c1,c2)	opcolor(c1,+=,c2)

#define  multcolor(c1,c2)	opcolor(c1,*=,c2)

#define NCSAMP	CNDX[EXP]	/* current number of spectral samples */
#define LSCOLR	(NCSAMP+1)

#define  copyscolr(sc1,sc2)	memcpy(sc1,sc2,LSCOLR)

#define  scolval(sc,pri)	(sc)[CNDX[pri]]

#define  cx2real(m,e)		(((m)+.5f)*cxponent[e])

#define  copyscolor(sc1,sc2)	memcpy(sc1,sc2,sizeof(COLORV)*NCSAMP)

#define  scalescolor(sc,sf)	{const float _f=sf; int _i=NCSAMP; \
					while (_i--) (sc)[_i] *= _f;}

				/* faster, use principal colors for RGB */
#define pcolor_color(col,scol)	setcolor(col,scolval(scol,RED),\
					scolval(scol,GRN),scolval(scol,BLU))
#define pcolor_colr(clr,scol)	setcolr(clr,scolval(scol,RED),\
					scolval(scol,GRN),scolval(scol,BLU))

#define  sopscolor(sc1,op,sc2)	{int _i=NCSAMP; while (_i--) (sc1)[_i] op (sc2)[_i];}

#define  saddscolor(sc1,sc2)	sopscolor(sc1,+=,sc2)

#define  smultscolor(sc1,sc2)	sopscolor(sc1,*=,sc2)

#define  scolrblack(c)		memset(c,0,LSCOLR)

#define  scolorblack(c)		memset(c,0,sizeof(COLORV)*NCSAMP)

#define scolor_color(col,scol)	scolor2color(col,scol,NCSAMP,WLPART)
#define scolor_colr(clr,scol)	scolor2colr(clr,scol,NCSAMP,WLPART)
#define scolor_scolr(sclr,scol)	scolor2scolr(sclr,scol,NCSAMP)
#define scolr_scolor(scol,sclr)	scolr2scolor(scol,sclr,NCSAMP)
#define scolr_color(col,sclr)	scolr2color(col,sclr,NCSAMP,WLPART)

#define  sopcolor(sc1,op,c2)	{SCOLOR _sct;\
				  setscolor(_sct,(c2)[RED],(c2)[GRN],(c2)[BLU]);\
					sopscolor(sc1,op,_sct);}

#define  saddcolor(sc1,c2)	sopcolor(sc1,+=,c2)

#define  smultcolor(sc1,c2)	sopcolor(sc1,*=,c2)

#define  opscolor(c1,op,sc2)	{COLOR _ct; scolor_color(_ct,sc2);\
					opcolor(c1,op,_ct);}

#define  addscolor(c1,sc2)	opscolor(c1,+=,sc2)

#define  multscolor(c1,sc2)	opscolor(c1,*=,sc2)

#if defined(NTSC_RGB)
#define  CIE_x_r		0.670		/* standard NTSC primaries */
#define  CIE_y_r		0.330
#define  CIE_x_g		0.210
#define  CIE_y_g		0.710
#define  CIE_x_b		0.140
#define  CIE_y_b		0.080
#define  CIE_x_w		(1./3.)		/* use EE white */
#define  CIE_y_w		(1./3.)
#elif defined(SHARP_RGB)
#define  CIE_x_r		0.6898		/* "sharp" RGB primaries */
#define  CIE_y_r		0.3206
#define  CIE_x_g		0.0736
#define  CIE_y_g		0.9003
#define  CIE_x_b		0.1166
#define  CIE_y_b		0.0374
#define  CIE_x_w		(1./3.)		/* use EE white */
#define  CIE_y_w		(1./3.)
#else
#define  CIE_x_r		0.640		/* nominal CRT primaries */
#define  CIE_y_r		0.330
#define  CIE_x_g		0.290
#define  CIE_y_g		0.600
#define  CIE_x_b		0.150
#define  CIE_y_b		0.060
#define  CIE_x_w		(1./3.)		/* use EE white */
#define  CIE_y_w		(1./3.)
#endif

#define  STDPRIMS	{{CIE_x_r,CIE_y_r},{CIE_x_g,CIE_y_g}, \
				{CIE_x_b,CIE_y_b},{CIE_x_w,CIE_y_w}}

#define CIE_D		(	CIE_x_r*(CIE_y_g - CIE_y_b) + \
				CIE_x_g*(CIE_y_b - CIE_y_r) + \
				CIE_x_b*(CIE_y_r - CIE_y_g)	)
#define CIE_C_rD	( (1./CIE_y_w) * \
				( CIE_x_w*(CIE_y_g - CIE_y_b) - \
				  CIE_y_w*(CIE_x_g - CIE_x_b) + \
				  CIE_x_g*CIE_y_b - CIE_x_b*CIE_y_g	) )
#define CIE_C_gD	( (1./CIE_y_w) * \
				( CIE_x_w*(CIE_y_b - CIE_y_r) - \
				  CIE_y_w*(CIE_x_b - CIE_x_r) - \
				  CIE_x_r*CIE_y_b + CIE_x_b*CIE_y_r	) )
#define CIE_C_bD	( (1./CIE_y_w) * \
				( CIE_x_w*(CIE_y_r - CIE_y_g) - \
				  CIE_y_w*(CIE_x_r - CIE_x_g) + \
				  CIE_x_r*CIE_y_g - CIE_x_g*CIE_y_r	) )

#define CIE_rf		(CIE_y_r*CIE_C_rD/CIE_D)
#define CIE_gf		(CIE_y_g*CIE_C_gD/CIE_D)
#define CIE_bf		(CIE_y_b*CIE_C_bD/CIE_D)

/* Default CIE_rf=.265074126, CIE_gf=.670114631 and CIE_bf=.064811243 */

/***** The following definitions are not for XYZ colors *****/

#define  bright(col)	(CIE_rf*(col)[RED]+CIE_gf*(col)[GRN]+CIE_bf*(col)[BLU])
#define  pbright(col)	(CIE_rf*scolval(col,RED) + CIE_gf*scolval(col,GRN) + \
				CIE_bf*scolval(col,BLU))
#define  normbright(c)	( ( (long)(CIE_rf*256.+.5)*(c)[RED] + \
			    (long)(CIE_gf*256.+.5)*(c)[GRN] + \
			    (long)(CIE_bf*256.+.5)*(c)[BLU] ) >> 8 )
#define  normpbright(c)	( ( (long)(CIE_rf*256.+.5)*(c)[CNDX[RED]] + \
			    (long)(CIE_gf*256.+.5)*(c)[CNDX[GRN]] + \
			    (long)(CIE_bf*256.+.5)*(c)[CNDX[BLU]] ) >> 8 )

				/* luminous efficacies over visible spectrum */
#define  MAXEFFICACY		683.		/* defined maximum at 550 nm */
#define  WHTEFFICACY		179.		/* equal energy white 380-780nm */
#define  D65EFFICACY		203.		/* standard illuminant D65 */
#define  INCEFFICACY		160.		/* illuminant A (incand.) */
#define  SUNEFFICACY		208.		/* illuminant B (solar dir.) */
#define  SKYEFFICACY		D65EFFICACY	/* skylight (should be 110) */
#define  DAYEFFICACY		D65EFFICACY	/* combined sky and solar */
#define  WHTSCOTOPIC		412.		/* scotopic EE white 380-780nm */
#define  WHTMELANOPIC		179.		/* melanopic EE white 380-780nm */

#define  luminance(col)		(WHTEFFICACY * bright(col))
#define  pluminance(scol)	(WHTEFFICACY * pbright(scol))

#define  scolor_rgb(col,scol)	scolor2rgb(col,scol,NCSAMP,WLPART)

/***** ...end of stuff specific to RGB colors *****/

#define  scolor_cie(col,scol)	scolor2cie(col,scol,NCSAMP,WLPART)

#define  sluminance(scol)	(WHTEFFICACY * scolor_photopic(scol))

#define  intens(col)		( (col)[0] > (col)[1] \
				? (col)[0] > (col)[2] ? (col)[0] : (col)[2] \
				: (col)[1] > (col)[2] ? (col)[1] : (col)[2] )

#define  colrval(c,p)		cx2real((c)[p],(c)[EXP])

#define  scolrval(c,p)		cx2real((c)[CNDX[p]],(c)[CNDX[EXP]])

#define  WHTCOLOR		{1.0,1.0,1.0}
#define  BLKCOLOR		{0.0,0.0,0.0}
#define  WHTCOLR		{128,128,128,COLXS+1}
#define  BLKCOLR		{0,0,0,0}

				/* picture format identifier */
#define  COLRFMT		"32-bit_rle_rgbe"
#define  CIEFMT			"32-bit_rle_xyze"
#define  PICFMT			"32-bit_rle_???e"	/* matches either */
#define  SPECFMT		"Radiance_spectra"	/* spectral data w/ exponent */

				/* Number of spectral components */
#define  NCOMPSTR		"NCOMP="
#define  LNCOMPSTR		6
#define  isncomp(hl)		!strncmp(hl,NCOMPSTR,LNCOMPSTR)
#define  ncompval(hl)		atoi((hl)+LNCOMPSTR)
#define  fputncomp(nc,fp)	fprintf(fp,"%s%d\n",NCOMPSTR,nc)

				/* 4 wavelength partitions for (IR+)R,G,B(+UV) */
#define  WLSPLTSTR		"WAVELENGTH_SPLITS="
#define  LWLSPLTSTR		18
#define  iswlsplit(hl)		!strncmp(hl,WLSPLTSTR,LWLSPLTSTR)
#define  wlsplitval(w,hl)	(sscanf((hl)+LWLSPLTSTR,"%f %f %f %f",\
					&(w)[0],&(w)[1],&(w)[2],&(w)[3]) == 4)
#define  fputwlsplit(w,fp)	fprintf(fp,"%s %g %g %g %g\n",WLSPLTSTR,\
					(w)[0],(w)[1],(w)[2],(w)[3])

				/* macros for exposures */
#define  EXPOSSTR		"EXPOSURE="
#define  LEXPOSSTR		9
#define  isexpos(hl)		!strncmp(hl,EXPOSSTR,LEXPOSSTR)
#define  exposval(hl)		atof((hl)+LEXPOSSTR)
#define  fputexpos(ex,fp)	fprintf(fp,"%s%.4e\n",EXPOSSTR,ex)

				/* macros for pixel aspect ratios */
#define  ASPECTSTR		"PIXASPECT="
#define  LASPECTSTR		10
#define  isaspect(hl)		!strncmp(hl,ASPECTSTR,LASPECTSTR)
#define  aspectval(hl)		atof((hl)+LASPECTSTR)
#define  fputaspect(pa,fp)	fprintf(fp,"%s%f\n",ASPECTSTR,pa)

				/* macros for primary specifications */
#define  PRIMARYSTR		"PRIMARIES="
#define  LPRIMARYSTR		10
#define  isprims(hl)		!strncmp(hl,PRIMARYSTR,LPRIMARYSTR)
#define  primsval(p,hl)		(sscanf((hl)+LPRIMARYSTR, \
					"%f %f %f %f %f %f %f %f", \
					&(p)[RED][CIEX],&(p)[RED][CIEY], \
					&(p)[GRN][CIEX],&(p)[GRN][CIEY], \
					&(p)[BLU][CIEX],&(p)[BLU][CIEY], \
					&(p)[WHT][CIEX],&(p)[WHT][CIEY]) == 8)
#define  fputprims(p,fp)	fprintf(fp, \
				"%s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",\
					PRIMARYSTR, \
					(p)[RED][CIEX],(p)[RED][CIEY], \
					(p)[GRN][CIEX],(p)[GRN][CIEY], \
					(p)[BLU][CIEX],(p)[BLU][CIEY], \
					(p)[WHT][CIEX],(p)[WHT][CIEY])

				/* macros for color correction */
#define  COLCORSTR		"COLORCORR="
#define  LCOLCORSTR		10
#define  iscolcor(hl)		!strncmp(hl,COLCORSTR,LCOLCORSTR)
#define  colcorval(cc,hl)	(sscanf((hl)+LCOLCORSTR,"%f %f %f", \
					&(cc)[RED],&(cc)[GRN],&(cc)[BLU]) == 3)
#define  fputcolcor(cc,fp)	fprintf(fp,"%s %f %f %f\n",COLCORSTR, \
					(cc)[RED],(cc)[GRN],(cc)[BLU])

/*
 * Conversions to and from XYZ space generally don't apply WHTEFFICACY.
 * If you need Y to be luminance (cd/m^2), this must be applied when
 * converting from radiance (watts/sr/m^2).
 */

extern const float  cxponent[256];	/* exponent look-up */
extern int  CNDX[4];		/* RGBE indices for SCOLOR, SCOLR */
extern float  WLPART[4];	/* RGB wavelength limits+partitions (nm) */
extern RGBPRIMS  stdprims;	/* standard primary chromaticities */
extern RGBPRIMS  xyzprims;	/* to indicate XYZ input or output */
extern COLORMAT  rgb2xyzmat;	/* RGB to XYZ conversion matrix */
extern COLORMAT  xyz2rgbmat;	/* XYZ to RGB conversion matrix */
extern const COLOR  cblack, cwhite;	/* black (0,0,0) and white (1,1,1) */
extern const SCOLOR scblack;	/* black spectral color (all 0's) */

#define  CGAMUT_LOWER		01
#define  CGAMUT_UPPER		02
#define  CGAMUT			(CGAMUT_LOWER|CGAMUT_UPPER)

#define  rgb_cie(xyz,rgb)	colortrans(xyz,rgb2xyzmat,rgb)

#define  cpcolormat(md,ms)	memcpy((void *)md,(void *)ms,sizeof(COLORMAT))

#define  colr_color(col,clr)	((col)[RED]=colrval(clr,RED),\
				(col)[GRN]=colrval(clr,GRN),\
				(col)[BLU]=colrval(clr,BLU))

					/* defined in color.c */
extern void	*tempbuffer(size_t len);
				/* in cn[3]=nsamps, wlpt[0],wlpt[3]=extrema */
extern int	setspectrsamp(int cn[4], float wlpt[4]);
extern int	fwritecolrs(COLR *scanline, int len, FILE *fp);
extern int	fwritescan(COLOR *scanline, int len, FILE *fp);
extern int	fwritescolrs(const COLRV *sscanline, int nc, int len, FILE *fp);
extern int	fwritesscan(const COLORV *sscanline, int nc, int len, FILE *fp);
extern int	freadcolrs(COLR *scanline, int len, FILE *fp);
extern int	freadscan(COLOR *scanline, int len, FILE *fp);
extern int	freadscolrs(COLRV *scanline, int nc, int len, FILE *fp);
extern int	freadsscan(COLORV *sscanline, int nc, int len, FILE *fp);
extern int	fread2colrs(COLR *scanline, int len, FILE *fp,
				int nc, const float wlpt[4]);
extern int	fread2scan(COLOR *scanline, int len, FILE *fp,
				int nc, const float wlpt[4]);
				/* spectrum conversion, zero-fill past ends */
extern void	convertscolor(SCOLOR dst, int dnc, double dwl0, double dwl1,
			const COLORV src[], int snc, double swl0, double swl1);
				/* the following use avg spectral ranges */
				/* compare scolor_rgb() and scolor_cie() */
				/* also, pcolor_color() and pcolor_colr() */
extern void	setscolor(SCOLOR scol, double r, double g, double b);
extern void	scolor2color(COLOR col, const SCOLOR scol, int ncs, const float wlpt[4]);
extern void	scolor2colr(COLR clr, const SCOLOR scol, int ncs, const float wlpt[4]);
extern void	scolr2colr(COLR clr, const SCOLR sclr, int ncs, const float wlpt[4]);
extern void	scolor2scolr(SCOLR sclr, const SCOLOR scol, int ncs);
extern void	scolr2scolor(SCOLOR scol, const SCOLR sclr, int ncs);
extern void	scolr2color(COLOR col, const SCOLR sclr, int ncs, const float wlpt[4]);
extern void	setcolr(COLR clr, double r, double g, double b);
extern void	setscolr(SCOLR sclr, double r, double g, double b);
extern double	scolor_mean(const SCOLOR scol);
extern double	sintens(const SCOLOR scol);
extern int	bigdiff(const COLOR c1, const COLOR c2, double md);
extern int	sbigsdiff(const SCOLOR sc1, const SCOLOR sc2, double md);
					/* defined in spec_rgb.c */
extern void	scolor_out(COLORV *cout, RGBPRIMS pr, const SCOLOR cres);
extern void	spec_cie(COLOR col, int s, int e);
extern void	spec_rgb(COLOR col, int s, int e);
extern void	cie_rgb(COLOR rgb, const COLOR xyz);
extern int	clipgamut(COLOR col, double brt, int gamut,
				const COLOR lower, const COLOR upper);
extern void	colortrans(COLOR c2, const COLORMAT mat, const COLOR c1);
extern void	multcolormat(COLORMAT m3, const COLORMAT m2,
					const COLORMAT m1);
extern int	colorprimsOK(RGBPRIMS pr);
extern int	compxyz2rgbmat(COLORMAT mat, RGBPRIMS pr);
extern int	comprgb2xyzmat(COLORMAT mat, RGBPRIMS pr);
extern int	comprgb2rgbmat(COLORMAT mat, RGBPRIMS pr1, RGBPRIMS pr2);
extern int	compxyzWBmat(COLORMAT mat, const float wht1[2], const float wht2[2]);
extern int	compxyz2rgbWBmat(COLORMAT mat, RGBPRIMS pr);
extern int	comprgb2xyzWBmat(COLORMAT mat, RGBPRIMS pr);
extern int	comprgb2rgbWBmat(COLORMAT mat, RGBPRIMS pr1, RGBPRIMS pr2);
					/* any uniform spectrum to working */
extern void	convertscolorcol(SCOLOR rcol, const COLORV src[], int snc,
					double swl0, double swl1);
					/* most accurate spectral->tristim */
extern void	scolor2cie(COLOR col, const SCOLOR scol, int ncs, const float wlpt[4]);
extern void	scolor2rgb(COLOR col, const SCOLOR scol, int ncs, const float wlpt[4]);
extern double	scolor2photopic(const SCOLOR scol, int ncs, const float wlpt[4]);
extern double	scolor2scotopic(const SCOLOR scol, int ncs, const float wlpt[4]);
extern double	scolor2melanopic(const SCOLOR scol, int ncs, const float wlpt[4]);
extern double	scolor_photopic(const SCOLOR scol);
extern double	scolor_scotopic(const SCOLOR scol);
extern double	scolor_melanopic(const SCOLOR scol);
					/* defined in colrops.c */
extern int	setcolrcor(double (*f)(double, double), double a2);
extern int	setcolrinv(double (*f)(double, double), double a2);
extern int	setcolrgam(double g);
extern int	colrs_gambs(COLR *scan, int len);
extern int	gambs_colrs(COLR *scan, int len);
extern void	shiftcolrs(COLR *scan, int len, int adjust);
extern void	normcolrs(COLR *scan, int len, int adjust);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_COLOR_H_ */

