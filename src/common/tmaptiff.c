#ifndef lint
static const char	RCSid[] = "$Id: tmaptiff.c,v 3.11 2011/06/28 21:11:04 greg Exp $";
#endif
/*
 * Perform tone mapping on TIFF input.
 *
 * Externals declared in tmaptiff.h
 */

#include "copyright.h"

#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
#include "tmprivat.h"
#include "tmaptiff.h"

					/* input cases we handle */
#define TC_LOGLUV32	1
#define TC_LOGLUV24	2
#define TC_LOGL16	3
#define TC_GRYFLOAT	4
#define TC_RGBFLOAT	5
#define TC_GRYSHORT	6
#define TC_RGBSHORT	7

/* figure out what kind of TIFF we have and if we can tone-map it */
static int
getTIFFtype(TIFF *tif)
{
	uint16	comp, phot, pconf;
	uint16	samp_fmt, bits_samp;
	
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &pconf);
	if (pconf != PLANARCONFIG_CONTIG)
		return(0);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &phot);
	TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &comp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &samp_fmt);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_samp);
	switch (phot) {
	case PHOTOMETRIC_LOGLUV:
		if (comp == COMPRESSION_SGILOG)
			return(TC_LOGLUV32);
		if (comp == COMPRESSION_SGILOG24)
			return(TC_LOGLUV24);
		return(0);
	case PHOTOMETRIC_LOGL:
		if (comp == COMPRESSION_SGILOG)
			return(TC_LOGL16);
		return(0);
	case PHOTOMETRIC_MINISBLACK:
		if (samp_fmt == SAMPLEFORMAT_UINT) {
			if (bits_samp == 16)
				return(TC_GRYSHORT);
			return(0);
		}
		if (samp_fmt == SAMPLEFORMAT_IEEEFP) {
			if (bits_samp == 8*sizeof(float))
				return(TC_GRYFLOAT);
			return(0);
		}
		return(0);
	case PHOTOMETRIC_RGB:
		if (samp_fmt == SAMPLEFORMAT_UINT) {
			if (bits_samp == 16)
				return(TC_RGBSHORT);
			return(0);
		}
		if (samp_fmt == SAMPLEFORMAT_IEEEFP) {
			if (bits_samp == 8*sizeof(float))
				return(TC_RGBFLOAT);
			return(0);
		}
		return(0);
	}
	return(0);
}

/* load and convert TIFF */
int
tmLoadTIFF(TMstruct *tms, TMbright **lpp, uby8 **cpp,
		int *xp, int *yp, char *fname, TIFF *tp)
{
	char	*funcName = fname==NULL ? "tmLoadTIFF" : fname;
	RGBPRIMP	inppri = tms->monpri;
	RGBPRIMS	myprims;
	float	*fa;
	TIFF	*tif;
	int	err;
	union {uint16 *w; uint32 *l; float *f; MEM_PTR p;} sl;
	uint32	width, height;
	int	tcase;
	double	stonits;
	int	y;
					/* check arguments */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((lpp == NULL) | (xp == NULL) | (yp == NULL) |
			((fname == NULL) & (tp == NULL)))
		returnErr(TM_E_ILLEGAL);
					/* check/get TIFF tags */
	sl.p = NULL; *lpp = NULL;
	if (cpp != TM_NOCHROMP) *cpp = TM_NOCHROM;
	err = TM_E_BADFILE;
	if ((tif = tp) == NULL && (tif = TIFFOpen(fname, "r")) == NULL)
		returnErr(TM_E_BADFILE);
	tcase = getTIFFtype(tif);
	if (!tcase)
		goto done;
	if (!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ||
			!TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height))
		goto done;
	*xp = width; *yp = height;
	if (TIFFGetField(tif, TIFFTAG_PRIMARYCHROMATICITIES, &fa)) {
		myprims[RED][CIEX] = fa[0];
		myprims[RED][CIEY] = fa[1];
		myprims[GRN][CIEX] = fa[2];
		myprims[GRN][CIEY] = fa[3];
		myprims[BLU][CIEX] = fa[4];
		myprims[BLU][CIEY] = fa[5];
		myprims[WHT][CIEX] = 1./3.;
		myprims[WHT][CIEY] = 1./3.;
		if (TIFFGetField(tif, TIFFTAG_WHITEPOINT, &fa)) {
			myprims[WHT][CIEX] = fa[0];
			myprims[WHT][CIEY] = fa[1];
		}
		inppri = myprims;
	}
	if (!TIFFGetField(tif, TIFFTAG_STONITS, &stonits))
		stonits = 1.;
	switch (tcase) {		/* set up conversion */
	case TC_LOGLUV32:
	case TC_LOGLUV24:
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_RAW);
		sl.l = (uint32 *)malloc(width*sizeof(uint32));
		tmSetSpace(tms, TM_XYZPRIM, stonits, NULL);
		break;
	case TC_LOGL16:
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_16BIT);
		sl.w = (uint16 *)malloc(width*sizeof(uint16));
		tmSetSpace(tms, tms->monpri, stonits, NULL);
		break;
	case TC_RGBFLOAT:
		sl.f = (float *)malloc(width*3*sizeof(float));
		tmSetSpace(tms, inppri, stonits, NULL);
		break;
	case TC_GRYFLOAT:
		sl.f = (float *)malloc(width*sizeof(float));
		tmSetSpace(tms, tms->monpri, stonits, NULL);
		break;
	case TC_RGBSHORT:
		sl.w = (uint16 *)malloc(width*3*sizeof(uint16));
		tmSetSpace(tms, inppri, stonits, NULL);
		break;
	case TC_GRYSHORT:
		sl.w = (uint16 *)malloc(width*sizeof(uint16));
		tmSetSpace(tms, tms->monpri, stonits, NULL);
		break;
	default:
		err = TM_E_CODERR1;
		goto done;
	}
	*lpp = (TMbright *)malloc(width*height*sizeof(TMbright));
	if ((sl.p == NULL) | (*lpp == NULL)) {
		err = TM_E_NOMEM;
		goto done;
	}
	switch (tcase) {		/* allocate color if needed */
	case TC_LOGLUV32:
	case TC_LOGLUV24:
	case TC_RGBFLOAT:
	case TC_RGBSHORT:
		if (cpp == TM_NOCHROMP)
			break;
		*cpp = (uby8 *)malloc(width*height*3*sizeof(uby8));
		if (*cpp == NULL) {
			err = TM_E_NOMEM;
			goto done;
		}
		break;
	}
					/* read and convert each scanline */
	for (y = 0; y < height; y++) {
		if (TIFFReadScanline(tif, sl.p, y, 0) < 0) {
			err = TM_E_BADFILE;
			break;
		}
		switch (tcase) {
		case TC_LOGLUV32:
			err = tmCvLuv32(tms, *lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					sl.l, width);
			break;
		case TC_LOGLUV24:
			err = tmCvLuv24(tms, *lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					sl.l, width);
			break;
		case TC_LOGL16:
			err = tmCvL16(tms, *lpp + y*width, sl.w, width);
			break;
		case TC_RGBFLOAT:
			err = tmCvColors(tms, *lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					(COLOR *)sl.f, width);
			break;
		case TC_GRYFLOAT:
			err = tmCvGrays(tms, *lpp + y*width, sl.f, width);
			break;
		case TC_RGBSHORT:
			err = tmCvRGB48(tms, *lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					(uint16 (*)[3])sl.w, width, DEFGAM);
			break;
		case TC_GRYSHORT:
			err = tmCvGray16(tms, *lpp + y*width,
					sl.w, width, DEFGAM);
			break;
		default:
			err = TM_E_CODERR1;
			break;
		}
		if (err != TM_E_OK)
			break;
	}
done:					/* clean up */
	if (tp == NULL)
		TIFFClose(tif);
	if (sl.p != NULL)
		free(sl.p);
	if (err != TM_E_OK) {		/* free buffers on error */
		if (*lpp != NULL)
			free((MEM_PTR)*lpp);
		*lpp = NULL;
		if (cpp != TM_NOCHROMP) {
			if (*cpp != TM_NOCHROM)
				free((MEM_PTR)*cpp);
			*cpp = NULL;
		}
		*xp = *yp = 0;
		returnErr(err);
	}
	returnOK;
}


/*
 * Load and tone-map a SGILOG TIFF.
 * Beware of greyscale input -- you must check the PHOTOMETRIC tag to
 * determine that the returned array contains only grey values, not RGB.
 * As in tmMapPicture(), grey values are also returned if flags&TM_F_BW.
 */
int
tmMapTIFF(uby8 **psp, int *xp, int *yp, int flags, RGBPRIMP monpri,
	double gamval, double Lddyn, double Ldmax, char *fname, TIFF *tp)
{
	char	*funcName = fname==NULL ? "tmMapTIFF" : fname;
	TMstruct	*tms = NULL;
	TMbright	*lp;
	uby8	*cp;
	int	err;
					/* check arguments */
	if ((psp == NULL) | (xp == NULL) | (yp == NULL) | (monpri == NULL) |
			((fname == NULL) & (tp == NULL)))
		returnErr(TM_E_ILLEGAL);
	if (gamval < MINGAM) gamval = DEFGAM;
	if (Lddyn < MINLDDYN) Lddyn = DEFLDDYN;
	if (Ldmax < MINLDMAX) Ldmax = DEFLDMAX;
	if (flags & TM_F_BW) monpri = stdprims;
					/* initialize tone mapping */
	if ((tms = tmInit(flags, monpri, gamval)) == NULL)
		returnErr(TM_E_NOMEM);
					/* load and convert TIFF */
	cp = TM_NOCHROM;
	err = tmLoadTIFF(tms, &lp, flags&TM_F_BW ? TM_NOCHROMP : &cp,
			xp, yp, fname, tp);
	if (err != TM_E_OK) {
		tmDone(tms);
		return(err);
	}
	if (cp == TM_NOCHROM) {
		*psp = (uby8 *)malloc(*xp * *yp * sizeof(uby8));
		if (*psp == NULL) {
			free((MEM_PTR)lp);
			tmDone(tms);
			returnErr(TM_E_NOMEM);
		}
	} else
		*psp = cp;
					/* compute color mapping */
	err = tmAddHisto(tms, lp, *xp * *yp, 1);
	if (err != TM_E_OK)
		goto done;
	err = tmComputeMapping(tms, gamval, Lddyn, Ldmax);
	if (err != TM_E_OK)
		goto done;
					/* map pixels */
	err = tmMapPixels(tms, *psp, lp, cp, *xp * *yp);

done:					/* clean up */
	free((MEM_PTR)lp);
	tmDone(tms);
	if (err != TM_E_OK) {		/* free memory on error */
		free((MEM_PTR)*psp);
		*psp = NULL;
		*xp = *yp = 0;
		returnErr(err);
	}
	returnOK;
}
