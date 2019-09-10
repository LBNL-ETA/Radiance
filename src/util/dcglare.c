/*
 * Compute time-step glare using imageless DGP calculation method.
 *
 *	N. Jones
 */

/*
 * Copyright (c) 2017-2019 Nathaniel Jones
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <ctype.h>
#include "platform.h"
#include "standard.h"
#include "cmatrix.h"
#include "resolu.h"
#include "cmglare.h"

char	*progname;			/* global argv[0] */

/* Sum together a set of images and write result to fout */
static int
sum_images(const char *fspec, const CMATRIX *cv, FILE *fout)
{
	int	myDT = DTfromHeader;
	COLR	*scanline = NULL;
	CMATRIX	*pmat = NULL;
	int	myXR=0, myYR=0;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = 0; i < cv->nrows; i++) {
		const COLORV	*scv = cv_lval(cv,i);
		int		flat_file = 0;
		char		fname[1024];
		FILE		*fp;
		long		data_start;
		int		dt, xr, yr;
		COLORV		*psp;
		char		*err;
							/* check for zero */
		if ((scv[RED] == 0) & (scv[GRN] == 0) & (scv[BLU] == 0) &&
				(myDT != DTfromHeader) | (i < cv->nrows-1))
			continue;
							/* open next picture */
		sprintf(fname, fspec, i);
		if ((fp = fopen(fname, "rb")) == NULL) {
			sprintf(errmsg, "cannot open picture '%s'", fname);
			error(SYSTEM, errmsg);
		}
		dt = DTfromHeader;
		if ((err = cm_getheader(&dt, NULL, NULL, NULL, fp)) != NULL)
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
			fputformat((char *)cm_fmt_id[myDT], fout);
			fputc('\n', fout);
			fflush(fout);
		} else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
					fname);
			error(USER, errmsg);
		}
							/* flat file check */
		if ((data_start = ftell(fp)) > 0 && fseek(fp, 0L, SEEK_END) == 0) {
			flat_file = (ftell(fp) == data_start + sizeof(COLR)*xr*yr);
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
	return(i);
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
	CMATRIX		*cmtx;		/* component vector/matrix result */
	char		fnbuf[256];
	int		a, i;
#ifdef DC_GLARE
	char	*direct_path = NULL;
	char	*schedule_path = NULL;
	int		*occupancy = NULL;
	int		start_hour = 0;
	int		end_hour = 24;
	double	dgp_limit = -1;
	double	dgp_threshold = 2000;
	char	*view_path = NULL;
	FILE	*fp;
	float	*dgp_values = NULL;
	FVECT	vdir, vup;
	FVECT	*views = NULL;
	int		viewfmt = DTascii;

	vdir[0] = vdir[1] = vdir[2] = vup[0] = vup[1] = 0;
	vup[2] = 1;

	clock_t timer = clock();
#endif /* DC_GLARE */

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
#ifndef DC_GLARE
			case '\0':	/* output specification (not format) */
				ofspec = argv[++a];
				break;
#endif /* DC_GLARE */
			case 'f':
				outfmt = DTfloat;
				break;
			case 'd':
				outfmt = DTdouble;
				break;
			case 'a':
				outfmt = DTascii;
				break;
			default:
				goto userr;
			}
			break;
#ifdef DC_GLARE
		case 's':
			switch (argv[a][2]) {
			case 'f':	/* occupancy schedule file */
				schedule_path = argv[++a];
				break;
			case 's':	/* occupancy start hour */
				start_hour = atoi(argv[++a]);
				break;
			case 'e':	/* occupancy end hour */
				end_hour = atoi(argv[++a]);
				break;
			default:
				goto userr;
			}
			break;
		case 'l':	/* perceptible glare threshold */
			dgp_limit = atof(argv[++a]);
			break;
		case 'b':	/* luminance threshold */
			dgp_threshold = atof(argv[++a]);
			break;
		case 'v':
			switch (argv[a][2]) {
			case 'd':	/* forward */
				//check(3, "fff");
				vdir[0] = atof(argv[++a]);
				vdir[1] = atof(argv[++a]);
				vdir[2] = atof(argv[++a]);
				if (normalize(vdir) == 0.0) goto userr;
				break;
			case 'u':	/* up */
				//check(3, "fff");
				vup[0] = atof(argv[++a]);
				vup[1] = atof(argv[++a]);
				vup[2] = atof(argv[++a]);
				if (normalize(vup) == 0.0) goto userr;
				break;
			case 'f':	/* view directions file */
				//check(2, "s");
				view_path = argv[++a];
				break;
			case 'i':
				switch (argv[a][3]) {
				case 'f':
					viewfmt = DTfloat;
					break;
				case 'd':
					viewfmt = DTdouble;
					break;
				case 'a':
					viewfmt = DTascii;
					break;
				default:
					goto userr;
				}
			default:
				goto userr;
			}
			break;
#endif /* DC_GLARE */
		default:
			goto userr;
		}
#ifdef DC_GLARE
	if ((argc-a < 2) | (argc-a > 5))
		goto userr;
	/* single bounce daylight coefficients file */
	direct_path = argv[++a];
#else
	if ((argc-a < 1) | (argc-a > 4))
		goto userr;
#endif /* DC_GLARE */

	if (argc-a > 2) {			/* VTDs expression */
		CMATRIX		*smtx, *Dmat, *Tmat, *imtx;
		const char	*ccp;
						/* get sky vector/matrix */
		smtx = cm_load(argv[a+3], 0, nsteps, skyfmt);
		nsteps = smtx->ncols;
						/* load BSDF */
		if (argv[a+1][0] != '!' &&
				(ccp = strrchr(argv[a+1], '.')) != NULL &&
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
		//TIMER(timer, "read args");
		cmtx = cm_load(argv[a+1], 0, nsteps, skyfmt);
		nsteps = cmtx->ncols;
		//TIMER(timer, "load sky matrix");
	}
						/* prepare output stream */
	if ((ofspec != NULL) & (nsteps == 1) && hasNumberFormat(ofspec)) {
		sprintf(fnbuf, ofspec, 1);
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
	if (hasNumberFormat(argv[a])) {		/* generating image(s) */
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
	} else {				/* generating vector/matrix */
		CMATRIX	*Vmat = cm_load(argv[a], 0, cmtx->nrows, DTfromHeader);
		//TIMER(timer, "load view matrix");
		CMATRIX	*rmtx = cm_multiply(Vmat, cmtx);
		//TIMER(timer, "multiply");
		cm_free(Vmat);
#ifdef DC_GLARE
		if (direct_path) { /* Do glare autonomy calculation */
			/* Load occupancy schedule */
			occupancy = (int*)malloc(nsteps * sizeof(int));
			if (!occupancy) {
				fprintf(stderr,
					"%s: out of memory for schedule\n",
					progname);
				return(1);
			}
			if (schedule_path) {
				if ((fp = fopen(schedule_path, "r")) == NULL) {
					fprintf(stderr,
						"%s: cannot open input file \"%s\"\n",
						progname, schedule_path);
					return(1);
				}
				if (cm_load_schedule(nsteps, occupancy, fp) < 0) return(1);
			}
			else {
				for (i = 0; i < nsteps; i++) {
					/* Assume hourly spacing */
					int hour = i % 24;
					occupancy[i] = ((hour >= start_hour) & (hour < end_hour));
				}
			}
			//TIMER(timer, "load occupancy schedule");

			/* Load view directions */
			if (view_path == NULL) {
				if (vdir[0] == 0.0f && vdir[1] == 0.0f && vdir[2] == 0.0f) {
					fprintf(stderr,
						"%s: missing view direction\n",
						progname);
					return(1);
				}
			}
			else {
				if ((fp = fopen(view_path, "r")) == NULL) {
					fprintf(stderr,
						"%s: cannot open input file \"%s\"\n",
						progname, view_path);
					return(1);
				}
				if (viewfmt != DTascii)
					SET_FILE_BINARY(fp);
				views = cm_load_views(rmtx->nrows, viewfmt, fp);
				if (!views) return(1);
				//TIMER(timer, "load views");
			}

			/* Calculate glare values */
			Vmat = cm_load(direct_path, 0, cmtx->nrows, DTfromHeader);
			//TIMER(timer, "load direct matrix");
			dgp_values = cm_glare(Vmat, rmtx, cmtx, occupancy, dgp_limit, dgp_threshold, views, vdir, vup);
			//TIMER(timer, "calculate dgp");
			free(views);
			cm_free(Vmat);

			/* Check successful calculation */
			if (!dgp_values) return(1);
		}
#endif /* DC_GLARE */
		if (ofspec != NULL) {		/* multiple vector files? */
			const char	*wtype = (outfmt==DTascii) ? "w" : "wb";
			for (i = 0; i < nsteps; i++) {
				CMATRIX	*rvec = cm_column(rmtx, i);
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
					fprintf(ofp, "NROWS=%d\n", rvec->nrows);
					fputs("NCOLS=1\nNCOMP=3\n", ofp);
					if ((outfmt == 'f') | (outfmt == 'd'))
						fputendian(ofp);
					fputformat((char *)cm_fmt_id[outfmt], ofp);
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
			if (headout) {		/* header output */
				newheader("RADIANCE", ofp);
				printargs(argc, argv, ofp);
				fputnow(ofp);
				fprintf(ofp, "NROWS=%d\n", rmtx->nrows);
#ifdef DC_GLARE
				fprintf(ofp, "NCOLS=%d\n", (!dgp_values || dgp_limit < 0) ? rmtx->ncols : 1);
				fprintf(ofp, "NCOMP=%d\n", dgp_values ? 1 : 3);
#else
				fprintf(ofp, "NCOLS=%d\n", rmtx->ncols);
				fputs("NCOMP=3\n", ofp);
#endif /* DC_GLARE */
				if ((outfmt == 'f') | (outfmt == 'd'))
					fputendian(ofp);
				fputformat((char *)cm_fmt_id[outfmt], ofp);
				fputc('\n', ofp);
			}
#ifdef DC_GLARE
			if (dgp_values) { /* Write glare autonomy */
				cm_write_glare(dgp_values, rmtx->nrows, dgp_limit < 0 ? rmtx->ncols : 1, outfmt, ofp);
				free(dgp_values);
				//TIMER(timer, "write");
			} 
			else
#endif /* DC_GLARE */
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
#ifdef DC_GLARE
	fprintf(stderr, "Usage: %s [-n nsteps][-i{f|d|h}][-o{f|d}][-l limit][-b threshold][{-sf occupancy|-ss start -se end}]{-vf views [-vi{f|d}]|-vd x y z}[-vu x y z] DCdirect DCtotal [skyf]\n",
		progname);
	fprintf(stderr, "   or: %s [-n nsteps][-i{f|d|h}][-o{f|d}][-l limit][-b threshold][{-sf occupancy|-ss start -se end}]{-vf views [-vi{f|d}]|-vd x y z}[-vu x y z] DCdirect Vspec Tbsdf Dmat.dat [skyf]\n",
		progname);
#else
	fprintf(stderr, "Usage: %s [-n nsteps][-o ospec][-i{f|d|h}][-o{f|d}] DCspec [skyf]\n",
				progname);
	fprintf(stderr, "   or: %s [-n nsteps][-o ospec][-i{f|d|h}][-o{f|d}] Vspec Tbsdf Dmat.dat [skyf]\n",
				progname);
#endif /* DC_GLARE */
	return(1);
}
