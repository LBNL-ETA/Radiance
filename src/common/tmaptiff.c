#ifndef lint
static const char	RCSid[] = "$Id: tmaptiff.c,v 3.2 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Perform tone mapping on TIFF input.
 *
 * Externals declared in tmaptiff.h
 */

#include "copyright.h"

#include <stdio.h>
#include "tiffio.h"
#include "tmprivat.h"
#include "tmaptiff.h"


int
tmLoadTIFF(lpp, cpp, xp, yp, fname, tp)	/* load and convert TIFF */
TMbright	**lpp;
BYTE	**cpp;
int	*xp, *yp;
char	*fname;
TIFF	*tp;
{
	char	*funcName = fname==NULL ? "tmLoadTIFF" : fname;
	TIFF	*tif;
	int	err;
	union {uint16 *w; uint32 *l; MEM_PTR p;} sl;
	uint16	comp, phot, pconf;
	uint32	width, height;
	double	stonits;
	int	y;
					/* check arguments */
	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (lpp == NULL | xp == NULL | yp == NULL |
			(fname == NULL & tp == NULL))
		returnErr(TM_E_ILLEGAL);
					/* check/get TIFF tags */
	sl.p = NULL; *lpp = NULL;
	if (cpp != TM_NOCHROMP) *cpp = TM_NOCHROM;
	err = TM_E_BADFILE;
	if ((tif = tp) == NULL && (tif = TIFFOpen(fname, "r")) == NULL)
		returnErr(TM_E_BADFILE);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &phot);
	TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &comp);
	if (!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ||
			!TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height))
		goto done;
	*xp = width; *yp = height;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &pconf);
	if (pconf != PLANARCONFIG_CONTIG)
		goto done;
	if (!TIFFGetField(tif, TIFFTAG_STONITS, &stonits))
		stonits = 1.;
	if (phot == PHOTOMETRIC_LOGLUV) {
		if (comp != COMPRESSION_SGILOG && comp != COMPRESSION_SGILOG24)
			goto done;
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_RAW);
		sl.l = (uint32 *)malloc(width*sizeof(uint32));
		if (cpp != TM_NOCHROMP) {
			*cpp = (BYTE *)malloc(width*height*3*sizeof(BYTE));
			if (*cpp == NULL) {
				err = TM_E_NOMEM;
				goto done;
			}
		}
	} else if (phot == PHOTOMETRIC_LOGL) {
		if (comp != COMPRESSION_SGILOG)
			goto done;
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_16BIT);
		sl.w = (uint16 *)malloc(width*sizeof(uint16));
	} else
		goto done;
	*lpp = (TMbright *)malloc(width*height*sizeof(TMbright));
	if (sl.p == NULL | *lpp == NULL) {
		err = TM_E_NOMEM;
		goto done;
	}
					/* set input color space */
	tmSetSpace(TM_XYZPRIM, stonits);
					/* read and convert each scanline */
	for (y = 0; y < height; y++) {
		if (TIFFReadScanline(tif, sl.p, y, 0) < 0) {
			err = TM_E_BADFILE;
			break;
		}
		if (phot == PHOTOMETRIC_LOGL)
			err = tmCvL16(*lpp + y*width, sl.w, width);
		else if (comp == COMPRESSION_SGILOG24)
			err = tmCvLuv24(*lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					sl.l, width);
		else
			err = tmCvLuv32(*lpp + y*width,
				cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp+y*3*width,
					sl.l, width);
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
tmMapTIFF(psp, xp, yp, flags, monpri, gamval, Lddyn, Ldmax, fname, tp)
BYTE	**psp;
int	*xp, *yp;
int	flags;
RGBPRIMP	monpri;
double	gamval, Lddyn, Ldmax;
char	*fname;
TIFF	*tp;
{
	char	*funcName = fname==NULL ? "tmMapTIFF" : fname;
	TMbright	*lp;
	BYTE	*cp;
	int	err;
					/* check arguments */
	if (psp == NULL | xp == NULL | yp == NULL | monpri == NULL |
			(fname == NULL & tp == NULL))
		returnErr(TM_E_ILLEGAL);
	if (gamval < MINGAM) gamval = DEFGAM;
	if (Lddyn < MINLDDYN) Lddyn = DEFLDDYN;
	if (Ldmax < MINLDMAX) Ldmax = DEFLDMAX;
	if (flags & TM_F_BW) monpri = stdprims;
					/* initialize tone mapping */
	if (tmInit(flags, monpri, gamval) == NULL)
		returnErr(TM_E_NOMEM);
					/* load and convert TIFF */
	cp = TM_NOCHROM;
	err = tmLoadTIFF(&lp, flags&TM_F_BW ? TM_NOCHROMP : &cp,
			xp, yp, fname, tp);
	if (err != TM_E_OK) {
		tmDone(NULL);
		return(err);
	}
	if (cp == TM_NOCHROM) {
		*psp = (BYTE *)malloc(*xp * *yp * sizeof(BYTE));
		if (*psp == NULL) {
			free((MEM_PTR)lp);
			tmDone(NULL);
			returnErr(TM_E_NOMEM);
		}
	} else
		*psp = cp;
					/* compute color mapping */
	err = tmAddHisto(lp, *xp * *yp, 1);
	if (err != TM_E_OK)
		goto done;
	err = tmComputeMapping(gamval, Lddyn, Ldmax);
	if (err != TM_E_OK)
		goto done;
					/* map pixels */
	err = tmMapPixels(*psp, lp, cp, *xp * *yp);

done:					/* clean up */
	free((MEM_PTR)lp);
	tmDone(NULL);
	if (err != TM_E_OK) {		/* free memory on error */
		free((MEM_PTR)*psp);
		*psp = NULL;
		*xp = *yp = 0;
		returnErr(err);
	}
	returnOK;
}
