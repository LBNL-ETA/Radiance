#ifndef lint
static const char RCSid[] = "$Id: dctimestep.c,v 2.53 2025/03/27 16:34:23 greg Exp $";
#endif
/*
 * Compute time-step result using Daylight Coefficient method.
 *
 *	G. Ward
 */

#include <ctype.h>
#include "platform.h"
#include "standard.h"
#include "cmatrix.h"
#include "platform.h"
#include "resolu.h"

char	*progname;			/* global argv[0] */

/* Sum together a set of images and write result to fout */
static int
sum_images(const char *fspec, const CMATRIX *cv, FILE *fout)
{
	static int	runcnt = 0;
	int	myDT = DTfromHeader;
	COLR	*scanline = NULL;
	CMATRIX	*pmat = NULL;
	int	myXR=0, myYR=0;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = cv->nrows; i-- > 0; ) {
		const int	r = runcnt&1 ? i : cv->nrows-1 - i;
		const COLORV	*scv = cv_lval(cv,r);
		int		flat_file = 0;
		char		fname[1024];
		FILE		*fp;
		long		data_start;
		int		dt, xr, yr;
		COLORV		*psp;
		char		*err;
							/* check for zero */
		if ((scv[RED] == 0) & (scv[GRN] == 0) & (scv[BLU] == 0) &&
				(myDT != DTfromHeader) | (i > 0))
			continue;
							/* open next picture */
		sprintf(fname, fspec, r);
		if ((fp = fopen(fname, "rb")) == NULL) {
			sprintf(errmsg, "cannot open picture '%s'", fname);
			error(SYSTEM, errmsg);
		}
#ifdef getc_unlocked
		flockfile(fp);
#endif
		dt = DTfromHeader;
		if ((err = cm_getheader(&dt, NULL, NULL, NULL, NULL, fp)) != NULL)
			error(USER, err);
		if ((dt != DTrgbe) & (dt != DTxyze) ||
				!fscnresolu(&xr, &yr, fp)) {
			sprintf(errmsg, "file '%s' not a picture", fname);
			error(USER, errmsg);
		}
		if (myDT == DTfromHeader) {		/* on first one */
			myDT = dt;
			myXR = xr; myYR = yr;
			scanline = (COLR *)malloc(sizeof(COLR)*myXR);
			if (scanline == NULL)
				error(SYSTEM, "out of memory in sum_images()");
			pmat = cm_alloc(myYR, myXR);
			memset(pmat->cmem, 0, sizeof(COLOR)*myXR*myYR);
							/* finish header */
			fputformat(cm_fmt_id[myDT], fout);
			fputc('\n', fout);
			fflush(fout);
		} else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
					fname);
			error(USER, errmsg);
		}
							/* flat file check */
		if ((data_start = ftell(fp)) > 0 && fseek(fp, 0L, SEEK_END) == 0) {
			flat_file = (ftell(fp) >= data_start + sizeof(COLR)*xr*yr);
			if (fseek(fp, data_start, SEEK_SET) < 0) {
				sprintf(errmsg, "cannot seek on picture '%s'", fname);
				error(SYSTEM, errmsg);
			}
		}
		psp = pmat->cmem;
		for (y = 0; y < yr; y++) {		/* read it in */
			COLOR	col;
			int	x;
			if (flat_file ? getbinary(scanline, sizeof(COLR), xr, fp) != xr :
					freadcolrs(scanline, xr, fp) < 0) {
				sprintf(errmsg, "error reading picture '%s'",
						fname);
				error(SYSTEM, errmsg);
			}
							/* sum in scanline */
			for (x = 0; x < xr; x++, psp += 3) {
				if (!scanline[x][EXP])
					continue;	/* skip zeroes */
				colr_color(col, scanline[x]);
				multcolor(col, scv);
				addcolor(psp, col);
			}
		}
		fclose(fp);				/* done this picture */
	}
	free(scanline);
	i = cm_write(pmat, myDT, fout);			/* write picture */
	cm_free(pmat);					/* free data */
	++runcnt;					/* for zig-zagging */
	return(i);
}

/* adjust matrix dimensions according to user size(s) */
static int
alt_dim(CMATRIX *cm, int nr, int nc)
{
	if ((nr <= 0) & (nc <= 0))
		return(0);
	if ((nr == cm->nrows) & (nc == cm->ncols))
		return(0);
	if (nr > 0) {
		if (nc <= 0)
			nc = cm->nrows*cm->ncols/nr;
		if (nr*nc != cm->nrows*cm->ncols) {
			fprintf(stderr, "Bad dimensions: %dx%d != %dx%d\n",
					nr, nc, cm->nrows, cm->ncols);
			return(-1);
		}
	} else /* nc > 0 */ {
		nr = cm->nrows*cm->ncols/nc;
		if (nc*nr != cm->nrows*cm->ncols) {
			fprintf(stderr, "Bad dimensions: %d does not divide %dx%d evenly\n",
					nc, cm->nrows, cm->ncols);
			return(-1);
		}
	}
	cm->nrows = nr;
	cm->ncols = nc;
	return(1);
}

/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char *s)
{
	if (s == NULL)
		return(0);

	while (*s) {
		while (*s != '%')
			if (!*s++)
				return(0);
		if (*++s == '%') {		/* ignore "%%" */
			++s;
			continue;
		}
		while (isdigit(*s))		/* field length */
			++s;
						/* field we'll use? */
		if ((*s == 'd') | (*s == 'i') | (*s == 'o') |
					(*s == 'x') | (*s == 'X'))
			return(1);
	}
	return(0);				/* didn't find one */
}

int
main(int argc, char *argv[])
{
	int		skyfmt = DTfromHeader;
	int		outfmt = DTascii;
	int		headout = 1;
	int		nsteps = 0;
	char		*ofspec = NULL;
	FILE		*ofp = stdout;
	int		xres=0, yres=0;
	CMATRIX		*cmtx;		/* component vector/matrix result */
	char		fnbuf[256];
	int		a, i;

	progname = argv[0];
					/* get options */
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'n':
			nsteps = atoi(argv[++a]);
			if (nsteps < 0)
				goto userr;
			skyfmt = nsteps ? DTascii : DTfromHeader;
			break;
		case 'h':
			headout = !headout;
			break;
		case 'i':
			switch (argv[a][2]) {
			case 'f':
				skyfmt = DTfloat;
				break;
			case 'd':
				skyfmt = DTdouble;
				break;
			case 'a':
				skyfmt = DTascii;
				break;
			default:
				goto userr;
			}
			break;
		case 'o':
			switch (argv[a][2]) {
			case '\0':	/* output specification (not format) */
				ofspec = argv[++a];
				break;
			case 'f':
				outfmt = DTfloat;
				break;
			case 'd':
				outfmt = DTdouble;
				break;
			case 'a':
				outfmt = DTascii;
				break;
			case 'c':
				outfmt = DTrgbe;
				break;
			default:
				goto userr;
			}
			break;
		case 'x':
			xres = atoi(argv[++a]);
			break;
		case 'y':
			yres = atoi(argv[++a]);
			break;
		default:
			goto userr;
		}
	if ((argc-a < 1) | (argc-a > 4))
		goto userr;

	if (argc-a > 2) {			/* VTDs expression */
		CMATRIX		*smtx, *Dmat, *Tmat, *imtx;
		const char	*ccp;
						/* get sky vector/matrix */
		smtx = cm_load(argv[a+3], 0, nsteps, skyfmt);
		nsteps = smtx->ncols;
						/* load BSDF */
		if (argv[a+1][0] != '!' &&
				(ccp = strrchr(argv[a+1], '.')) > argv[a+1] &&
				!strcasecmp(ccp+1, "XML"))
			Tmat = cm_loadBTDF(argv[a+1]);
		else
			Tmat = cm_load(argv[a+1], 0, 0, DTfromHeader);
						/* load Daylight matrix */
		Dmat = cm_load(argv[a+2], Tmat->ncols,
					smtx->nrows, DTfromHeader);
						/* multiply vector through */
		imtx = cm_multiply(Dmat, smtx);
		cm_free(Dmat); cm_free(smtx);
		cmtx = cm_multiply(Tmat, imtx);
		cm_free(Tmat); 
		cm_free(imtx);
	} else {				/* sky vector/matrix only */
		cmtx = cm_load(argv[a+1], 0, nsteps, skyfmt);
		nsteps = cmtx->ncols;
	}
						/* prepare output stream */
	if ((ofspec != NULL) & (nsteps == 1) && hasNumberFormat(ofspec)) {
		sprintf(fnbuf, ofspec, 0);
		ofspec = fnbuf;
	}
	if (ofspec != NULL && !hasNumberFormat(ofspec)) {
		if ((ofp = fopen(ofspec, "w")) == NULL) {
			fprintf(stderr, "%s: cannot open '%s' for output\n",
					progname, ofspec);
			return(1);
		}
		ofspec = NULL;			/* only need to open once */
	}
	if (hasNumberFormat(argv[a])) {		/* loading image vector(s) */
		if (outfmt != DTrgbe) {
			error(WARNING, "changing output type to -oc");
			outfmt = DTrgbe;
		}
		if (ofspec == NULL) {
			SET_FILE_BINARY(ofp);
			newheader("RADIANCE", ofp);
			printargs(argc, argv, ofp);
			fputnow(ofp);
		}
		if (nsteps > 1)			/* multiple output frames? */
			for (i = 0; i < nsteps; i++) {
				CMATRIX	*cvec = cm_column(cmtx, i);
				if (ofspec != NULL) {
					sprintf(fnbuf, ofspec, i);
					if ((ofp = fopen(fnbuf, "wb")) == NULL) {
						fprintf(stderr,
					"%s: cannot open '%s' for output\n",
							progname, fnbuf);
						return(1);
					}
					newheader("RADIANCE", ofp);
					printargs(argc, argv, ofp);
					fputnow(ofp);
				}
				fprintf(ofp, "FRAME=%d\n", i);
				if (!sum_images(argv[a], cvec, ofp))
					return(1);
				if (ofspec != NULL) {
					if (fclose(ofp) == EOF) {
						fprintf(stderr,
						"%s: error writing to '%s'\n",
							progname, fnbuf);
						return(1);
					}
					ofp = stdout;
				}
				cm_free(cvec);
			}
		else if (!sum_images(argv[a], cmtx, ofp))
			return(1);
	} else {				/* loading view matrix */
		CMATRIX	*Vmat = cm_load(argv[a], 0, cmtx->nrows, DTfromHeader);
		CMATRIX	*rmtx = cm_multiply(Vmat, cmtx);
		cm_free(Vmat);
		if (ofspec != NULL) {		/* multiple vector files? */
			const char	*wtype = (outfmt==DTascii) ? "w" : "wb";
			for (i = 0; i < nsteps; i++) {
				CMATRIX	*rvec = cm_column(rmtx, i);
				if (alt_dim(rvec, yres, xres) < 0)
					return(1);
				sprintf(fnbuf, ofspec, i);
				if ((ofp = fopen(fnbuf, wtype)) == NULL) {
					fprintf(stderr,
					"%s: cannot open '%s' for output\n",
							progname, fnbuf);
					return(1);
				}
#ifdef getc_unlocked
				flockfile(ofp);
#endif
				if (headout) {	/* header output */
					newheader("RADIANCE", ofp);
					printargs(argc, argv, ofp);
					fputnow(ofp);
					fprintf(ofp, "FRAME=%d\n", i);
					if ((outfmt != DTrgbe) & (outfmt != DTxyze)) {
						fprintf(ofp, "NROWS=%d\n", rvec->nrows);
						fprintf(ofp, "NCOLS=%d\n", rvec->ncols);
						fputncomp(3, ofp);
					}
					if ((outfmt == DTfloat) | (outfmt == DTdouble))
						fputendian(ofp);
					fputformat(cm_fmt_id[outfmt], ofp);
					fputc('\n', ofp);
				}
				cm_write(rvec, outfmt, ofp);
				if (fclose(ofp) == EOF) {
					fprintf(stderr,
						"%s: error writing to '%s'\n",
							progname, fnbuf);
					return(1);
				}
				ofp = stdout;
				cm_free(rvec);
			}
		} else {
#ifdef getc_unlocked
			flockfile(ofp);
#endif
			if (outfmt != DTascii)
				SET_FILE_BINARY(ofp);
			if (alt_dim(rmtx, yres, xres) < 0)
				return(1);
			if (headout) {		/* header output */
				newheader("RADIANCE", ofp);
				printargs(argc, argv, ofp);
				fputnow(ofp);
				if ((outfmt != DTrgbe) & (outfmt != DTxyze)) {
					fprintf(ofp, "NROWS=%d\n", rmtx->nrows);
					fprintf(ofp, "NCOLS=%d\n", rmtx->ncols);
					fputncomp(3, ofp);
				}
				if ((outfmt == DTfloat) | (outfmt == DTdouble))
					fputendian(ofp);
				fputformat(cm_fmt_id[outfmt], ofp);
				fputc('\n', ofp);
			}
			cm_write(rmtx, outfmt, ofp);
		}
		cm_free(rmtx);
	}
	if (fflush(ofp) == EOF) {		/* final clean-up */
		fprintf(stderr, "%s: write error on output\n", progname);
		return(1);
	}
	cm_free(cmtx);
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] DCspec [skyf]\n",
				progname);
	fprintf(stderr, "   or: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] Vspec Tbsdf Dmat.dat [skyf]\n",
				progname);
	return(1);
}
