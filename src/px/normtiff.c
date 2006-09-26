#ifndef lint
static const char	RCSid[] = "$Id: normtiff.c,v 3.10 2006/09/26 12:26:02 greg Exp $";
#endif
/*
 * Tone map SGILOG TIFF or Radiance picture and output 24-bit RGB TIFF
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "tiffio.h"
#include "color.h"
#include "tonemap.h"
#include "tmaptiff.h"
#include "resolu.h"


TIFF	*tifout;			/* TIFF output */
int	flags = TM_F_CAMERA;		/* tone-mapping flags */
RGBPRIMP	rgbp = stdprims;	/* display primaries */
RGBPRIMS	myprims;		/* overriding display primaries */
double	ldmax = 100.;			/* maximum display luminance */
double	lddyn = 32.;			/* display dynamic range */
double	gamv = 2.2;			/* display gamma value */

short	ortab[8] = {		/* orientation conversion table */
	YMAJOR|YDECR,
	YMAJOR|YDECR|XDECR,
	YMAJOR|XDECR,
	YMAJOR,
	YDECR,
	XDECR|YDECR,
	XDECR,
	0
};

typedef struct {
	FILE	*fp;		/* file pointer */
	char	fmt[32];	/* picture format */
	double	pa;		/* pixel aspect ratio */
	RESOLU	rs;		/* picture resolution */
} PICTURE;

uint16	comp = COMPRESSION_NONE;	/* TIFF compression mode */

#define closepicture(p)		(fclose((p)->fp),free((void *)(p)))

static gethfunc headline;

static int headline(char *s, void *pp);
static PICTURE *openpicture(char *fname);
static int tmap_picture(char *fname, PICTURE *pp);
static int tmap_tiff(char *fname, TIFF *tp);
static int putimage(uint16 or, uint32 xs, uint32 ys, float xr, float yr,
	uint16 ru, BYTE	*pd);


int
main(
	int	argc,
	char	*argv[]
)
{
	PICTURE	*pin = NULL;
	TIFF	*tin = NULL;
	int	i, rval;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'h':			/* human observer settings */
			flags = TM_F_HUMAN;
			break;
		case 's':			/* toggle human contrast */
			flags ^= TM_F_HCONTR;
			break;
		case 'c':			/* toggle mesopic sens. */
			flags ^= TM_F_MESOPIC;
			break;
		case 'l':			/* toggle linear mapping */
			flags ^= TM_F_LINEAR;
			break;
		case 'b':			/* toggle greyscale output */
			flags ^= TM_F_BW;
			break;
		case 'g':			/* set display gamma */
			if (argc-i < 2) goto userr;
			gamv = atof(argv[++i]);
			break;
		case 'u':			/* set display maximum */
			if (argc-i < 2) goto userr;
			ldmax = atof(argv[++i]);
			break;
		case 'd':			/* set display dynamic range */
			if (argc-i < 2) goto userr;
			lddyn = atof(argv[++i]);
			break;
		case 'z':			/* LZW compression */
			comp = COMPRESSION_LZW;
			break;
		case 'p':			/* set display primaries */
			if (argc-i < 9) goto userr;
			myprims[RED][CIEX] = atof(argv[++i]);
			myprims[RED][CIEY] = atof(argv[++i]);
			myprims[GRN][CIEX] = atof(argv[++i]);
			myprims[GRN][CIEY] = atof(argv[++i]);
			myprims[BLU][CIEX] = atof(argv[++i]);
			myprims[BLU][CIEY] = atof(argv[++i]);
			myprims[WHT][CIEX] = atof(argv[++i]);
			myprims[WHT][CIEY] = atof(argv[++i]);
			rgbp = myprims;
			break;
		default:
			goto userr;
		}
	if (argc-i < 2) goto userr;
	if ((pin = openpicture(argv[i])) == NULL &&
			(tin = TIFFOpen(argv[i], "r")) == NULL) {
		fprintf(stderr, "%s: cannot open or interpret file \"%s\"\n",
				argv[0], argv[i]);
		exit(1);
	}
	if ((tifout = TIFFOpen(argv[i+1], "w")) == NULL) {
		fprintf(stderr, "%s: cannot open output TIFF \"%s\"\n",
				argv[0], argv[i+1]);
		exit(1);
	}
	if (pin != NULL) {
		rval = tmap_picture(argv[i], pin);
		closepicture(pin);
	} else {
		rval = tmap_tiff(argv[i], tin);
		TIFFClose(tin);
	}
	TIFFClose(tifout);
	exit(rval==0 ? 0 : 1);
userr:
	fprintf(stderr,
"Usage: %s [-h][-s][-c][-l][-b][-g gv][-d ld][-u lm][-z][-p xr yr xg yg xb yb xw yw] input.{tif|pic} output.tif\n",
			argv[0]);
	exit(1);
}


static int
headline(				/* process line from header */
	char	*s,
	void *pp
)
{
	register char	*cp;

	for (cp = s; *cp; cp++)
		if (*cp & 0x80)
			return(-1);	/* non-ascii in header */
	if (isaspect(s))
		((PICTURE *)pp)->pa *= aspectval(s);
	else
		formatval(((PICTURE *)pp)->fmt, s);
	return(0);
}


static PICTURE *
openpicture(			/* open/check Radiance picture file */
	char	*fname
)
{
	FILE	*fp;
	register PICTURE	*pp;
	register char	*cp;
					/* check filename suffix */
	if (fname == NULL) return(NULL);
	for (cp = fname; *cp; cp++)
		;
	while (cp > fname && cp[-1] != '.')
		if (*--cp == '/') {
			cp = fname;
			break;
		}
	if (cp > fname && !strncmp(cp, "tif", 3))
		return(NULL);		/* assume it's a TIFF */
					/* else try opening it */
	if ((fp = fopen(fname, "r")) == NULL)
		return(NULL);
					/* allocate struct */
	if ((pp = (PICTURE *)malloc(sizeof(PICTURE))) == NULL)
		return(NULL);		/* serious error -- should exit? */
	pp->fp = fp; pp->fmt[0] = '\0'; pp->pa = 1.;
					/* load header */
	if (getheader(fp, headline, pp) < 0) {
		closepicture(pp);
		return(NULL);
	}
	if (!pp->fmt[0])		/* assume RGBE if unspecified */
		strcpy(pp->fmt, COLRFMT);
	if (!globmatch(PICFMT, pp->fmt) || !fgetsresolu(&pp->rs, fp)) {
		closepicture(pp);	/* failed test -- close file */
		return(NULL);
	}
	rewind(fp);			/* passed test -- rewind file */
	return(pp);
}


static int
tmap_picture(			/* tone map Radiance picture */
	char	*fname,
	register PICTURE	*pp
)
{
	uint16	orient;
	double	paspect = (pp->rs.rt & YMAJOR) ? pp->pa : 1./pp->pa;
	int	xsiz, ysiz;
	BYTE	*pix;
					/* read and tone map picture */
	if (tmMapPicture(&pix, &xsiz, &ysiz, flags,
			rgbp, gamv, lddyn, ldmax, fname, pp->fp) != TM_E_OK)
		return(-1);
					/* figure out TIFF orientation */
	for (orient = 8; --orient; )
		if (ortab[orient] == pp->rs.rt)
			break;
	orient++;
					/* put out our image */
	if (putimage(orient, (uint32)xsiz, (uint32)ysiz,
			72., 72./paspect, 2, pix) != 0)
		return(-1);
					/* free data and we're done */
	free((void *)pix);
	return(0);
}


static int
tmap_tiff(			/* tone map SGILOG TIFF */
	char	*fname,
	TIFF	*tp
)
{
	float	xres, yres;
	uint16	orient, resunit, phot;
	int	xsiz, ysiz;
	BYTE	*pix;
					/* check to make sure it's SGILOG */
	TIFFGetFieldDefaulted(tp, TIFFTAG_PHOTOMETRIC, &phot);
	if ((phot == PHOTOMETRIC_LOGL) | (phot == PHOTOMETRIC_MINISBLACK))
		flags |= TM_F_BW;
					/* read and tone map TIFF */
	if (tmMapTIFF(&pix, &xsiz, &ysiz, flags,
			rgbp, gamv, lddyn, ldmax, fname, tp) != TM_E_OK)
		return(-1);
					/* get relevant tags */
	TIFFGetFieldDefaulted(tp, TIFFTAG_RESOLUTIONUNIT, &resunit);
	TIFFGetFieldDefaulted(tp, TIFFTAG_XRESOLUTION, &xres);
	TIFFGetFieldDefaulted(tp, TIFFTAG_YRESOLUTION, &yres);
	TIFFGetFieldDefaulted(tp, TIFFTAG_ORIENTATION, &orient);
					/* put out our image */
	if (putimage(orient, (uint32)xsiz, (uint32)ysiz,
			xres, yres, resunit, pix) != 0)
		return(-1);
					/* free data and we're done */
	free((void *)pix);
	return(0);
}


static int
putimage(	/* write out our image */
	uint16	or,
	uint32	xs,
	uint32	ys,
	float	xr,
	float	yr,
	uint16	ru,
	BYTE	*pd
)
{
	register int	y;
	uint32	rowsperstrip;

	TIFFSetField(tifout, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	if (flags & TM_F_BW) {
		TIFFSetField(tifout, TIFFTAG_PHOTOMETRIC,
				PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tifout, TIFFTAG_SAMPLESPERPIXEL, 1);
	} else {
		TIFFSetField(tifout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tifout, TIFFTAG_SAMPLESPERPIXEL, 3);
	}
	if (rgbp != stdprims) {
		TIFFSetField(tifout, TIFFTAG_PRIMARYCHROMATICITIES,
				(float *)rgbp);
		TIFFSetField(tifout, TIFFTAG_WHITEPOINT, (float *)rgbp[WHT]);
	}
	TIFFSetField(tifout, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tifout, TIFFTAG_IMAGEWIDTH, xs);
	TIFFSetField(tifout, TIFFTAG_IMAGELENGTH, ys);
	TIFFSetField(tifout, TIFFTAG_RESOLUTIONUNIT, ru);
	TIFFSetField(tifout, TIFFTAG_COMPRESSION, comp);
	TIFFSetField(tifout, TIFFTAG_XRESOLUTION, xr);
	TIFFSetField(tifout, TIFFTAG_YRESOLUTION, yr);
	TIFFSetField(tifout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tifout, TIFFTAG_ORIENTATION, or);
					/* compute good strip size */
	rowsperstrip = 8192/TIFFScanlineSize(tifout);
	if (rowsperstrip < 1) rowsperstrip = 1;
	TIFFSetField(tifout, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
					/* write out scanlines */
	if (flags & TM_F_BW) {
		for (y = 0; y < ys; y++)
			if (TIFFWriteScanline(tifout, pd + y*xs, y, 0) < 0)
				goto writerr;
	} else {
		for (y = 0; y < ys; y++)
			if (TIFFWriteScanline(tifout, pd + y*3*xs, y, 0) < 0)
				goto writerr;
	}
	return(0);			/* all done! */
writerr:
	fputs("Error writing TIFF output\n", stderr);
	return(-1);
}
