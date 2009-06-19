#ifndef lint
static const char RCSid[] = "$Id: dctimestep.c,v 2.3 2009/06/19 14:21:42 greg Exp $";
#endif
/*
 * Compute time-step result using Daylight Coefficient method.
 *
 *	G. Ward
 */

#include <ctype.h>
#include "standard.h"
#include "platform.h"
#include "color.h"
#include "resolu.h"
#include "bsdf.h"

char	*progname;			/* global argv[0] */

/* Data types for file loading */
enum {DTfromHeader, DTascii, DTfloat, DTdouble, DTrgbe, DTxyze};

/* A color coefficient matrix -- vectors have ncols==1 */
typedef struct {
	int	nrows, ncols;
	COLORV	cmem[3];		/* extends struct */
} CMATRIX;

#define COLSPEC	(sizeof(COLORV)==sizeof(float) ? "%f %f %f" : "%lf %lf %lf")

#define cm_lval(cm,r,c)	((cm)->cmem + 3*((r)*(cm)->ncols + (c)))

#define cv_lval(cm,i)	((cm)->cmem + 3*(i))

/* Allocate a color coefficient matrix */
static CMATRIX *
cm_alloc(int nrows, int ncols)
{
	CMATRIX	*cm;

	if ((nrows <= 0) | (ncols <= 0))
		return(NULL);
	cm = (CMATRIX *)malloc(sizeof(CMATRIX) +
				3*sizeof(COLORV)*(nrows*ncols - 1));
	if (cm == NULL)
		error(SYSTEM, "out of memory in cm_alloc()");
	cm->nrows = nrows;
	cm->ncols = ncols;
	return(cm);
}

#define cm_free(cm)	free(cm)

/* Resize color coefficient matrix */
static CMATRIX *
cm_resize(CMATRIX *cm, int nrows)
{
	if (nrows == cm->nrows)
		return(cm);
	if (nrows <= 0) {
		cm_free(cm);
		return(NULL);
	}
	cm = (CMATRIX *)realloc(cm, sizeof(CMATRIX) +
			3*sizeof(COLORV)*(nrows*cm->ncols - 1));
	if (cm == NULL)
		error(SYSTEM, "out of memory in cm_resize()");
	cm->nrows = nrows;
	return(cm);
}

/* Load header to obtain data type */
static int
getDT(char *s, void *p)
{
	char	fmt[32];
	
	if (formatval(fmt, s)) {
		if (!strcmp(fmt, "ascii"))
			*((int *)p) = DTascii;
		else if (!strcmp(fmt, "float"))
			*((int *)p) = DTfloat;
		else if (!strcmp(fmt, "double"))
			*((int *)p) = DTdouble;
		else if (!strcmp(fmt, COLRFMT))
			*((int *)p) = DTrgbe;
		else if (!strcmp(fmt, CIEFMT))
			*((int *)p) = DTxyze;
	}
	return(0);
}

static int
getDTfromHeader(FILE *fp)
{
	int	dt = DTfromHeader;
	
	if (getheader(fp, getDT, &dt) < 0)
		error(SYSTEM, "header read error");
	if (dt == DTfromHeader)
		error(USER, "missing data format in header");
	return(dt);
}

/* Allocate and load a matrix from the given file (or stdin if NULL) */
static CMATRIX *
cm_load(const char *fname, int nrows, int ncols, int dtype)
{
	CMATRIX	*cm;
	FILE	*fp = stdin;

	if (fname == NULL)
		fname = "<stdin>";
	else if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open file '%s'", fname);
		error(SYSTEM, errmsg);
	}
	if (dtype != DTascii)
		SET_FILE_BINARY(fp);
	if (dtype == DTfromHeader)
		dtype = getDTfromHeader(fp);
	switch (dtype) {
	case DTascii:
	case DTfloat:
	case DTdouble:
		break;
	default:
		error(USER, "unexpected data type in cm_load()");
	}
	if (nrows <= 0) {			/* don't know length? */
		int	guessrows = 147;	/* usually big enough */
		if ((dtype != DTascii) & (fp != stdin)) {
			long	startpos = ftell(fp);
			if (fseek(fp, 0L, SEEK_END) == 0) {
				long	endpos = ftell(fp);
				long	elemsiz = 3*(dtype==DTfloat ?
					    sizeof(float) : sizeof(double));

				if ((endpos - startpos) % (ncols*elemsiz)) {
					sprintf(errmsg,
					"improper length for binary file '%s'",
							fname);
					error(USER, errmsg);
				}
				guessrows = (endpos - startpos)/(ncols*elemsiz);
				if (fseek(fp, startpos, SEEK_SET) < 0) {
					sprintf(errmsg,
						"fseek() error on file '%s'",
							fname);
					error(SYSTEM, errmsg);
				}
				nrows = guessrows;	/* we're confident */
			}
		}
		cm = cm_alloc(guessrows, ncols);
	} else
		cm = cm_alloc(nrows, ncols);
	if (cm == NULL)
		return(NULL);
	if (dtype == DTascii) {				/* read text file */
		int	maxrow = (nrows > 0 ? nrows : 32000);
		int	r, c;
		for (r = 0; r < maxrow; r++) {
		    if (r >= cm->nrows)		/* need more space? */
			cm = cm_resize(cm, 2*cm->nrows);
		    for (c = 0; c < ncols; c++) {
		        COLORV	*cv = cm_lval(cm,r,c);
			if (fscanf(fp, COLSPEC, cv, cv+1, cv+2) != 3)
				if ((nrows <= 0) & (r > 0) & (c == 0)) {
					cm = cm_resize(cm, maxrow=r);
					break;
				} else
					goto EOFerror;
		    }
		}
		while ((c = getc(fp)) != EOF)
			if (!isspace(c)) {
				sprintf(errmsg,
				"unexpected data at end of ascii file %s",
						fname);
				error(WARNING, errmsg);
				break;
			}
	} else {					/* read binary file */
		if (sizeof(COLORV) == (dtype==DTfloat ? sizeof(float) :
							sizeof(double))) {
			int	nread = 0;
			do {				/* read all we can */
				nread += fread(cm->cmem + 3*nread,
						3*sizeof(COLORV),
						cm->nrows*cm->ncols - nread,
						fp);
				if (nrows <= 0) {	/* unknown length */
					if (nread == cm->nrows*cm->ncols)
							/* need more space? */
						cm = cm_resize(cm, 2*cm->nrows);
					else if (nread % cm->ncols == 0)
							/* seem to be  done */
						cm = cm_resize(cm, nread/cm->ncols);
					else		/* ended mid-row */
						goto EOFerror;
				} else if (nread < cm->nrows*cm->ncols)
					goto EOFerror;
			} while (nread < cm->nrows*cm->ncols);

		} else if (dtype == DTdouble) {
			double	dc[3];			/* load from double */
			COLORV	*cvp = cm->cmem;
			int	n = nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (fread(dc, sizeof(double), 3, fp) != 3)
					goto EOFerror;
				copycolor(cvp, dc);
				cvp += 3;
			}
		} else /* dtype == DTfloat */ {
			float	fc[3];			/* load from float */
			COLORV	*cvp = cm->cmem;
			int	n = nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (fread(fc, sizeof(float), 3, fp) != 3)
					goto EOFerror;
				copycolor(cvp, fc);
				cvp += 3;
			}
		}
		if (getc(fp) != EOF) {
				sprintf(errmsg,
				"unexpected data at end of binary file %s",
						fname);
				error(WARNING, errmsg);
		}
	}
	if (fp != stdin)
		fclose(fp);
	return(cm);
EOFerror:
	sprintf(errmsg, "unexpected EOF reading %s",
			fname);
	error(USER, errmsg);
not_handled:
	error(INTERNAL, "unhandled data size or length in cm_load()");
	return(NULL);	/* gratis return */
}

/* Multiply two matrices (or a matrix and a vector) and allocate the result*/
static CMATRIX *
cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2)
{
	CMATRIX	*cmr;
	int	dr, dc, i;

	if ((cm1->ncols <= 0) | (cm1->ncols != cm2->nrows))
		error(INTERNAL, "matrix dimension mismatch in cm_multiply()");
	cmr = cm_alloc(cm1->nrows, cm2->ncols);
	if (cmr == NULL)
		return(NULL);
	for (dr = 0; dr < cmr->nrows; dr++)
	    for (dc = 0; dc < cmr->ncols; dc++) {
		COLORV	*dp = cm_lval(cmr,dr,dc);
		dp[0] = dp[1] = dp[2] = 0;
		for (i = 0; i < cm1->ncols; i++) {
		    const COLORV	*cp1 = cm_lval(cm1,dr,i);
		    const COLORV	*cp2 = cm_lval(cm2,i,dc);
		    dp[0] += cp1[0] * cp2[0];
		    dp[1] += cp1[1] * cp2[1];
		    dp[2] += cp1[2] * cp2[2];
		}
	    }
	return(cmr);
}

/* print out matrix as ASCII text -- no header */
static void
cm_print(const CMATRIX *cm, FILE *fp)
{
	int		r, c;
	const COLORV	*mp = cm->cmem;
	
	for (r = 0; r < cm->nrows; r++) {
		for (c = 0; c < cm->ncols; c++, mp += 3)
			fprintf(fp, "\t%.6e %.6e %.6e", mp[0], mp[1], mp[2]);
		fputc('\n', fp);
	}
}

/* convert a BSDF to our matrix representation */
static CMATRIX *
cm_bsdf(const struct BSDF_data *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->nout, bsdf->ninc);
	COLORV	*mp = cm->cmem;
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
	
	for (r = 0; r < cm->nrows; r++)
		for (c = 0; c < cm->ncols; c++, mp += 3) {
			float	f = BSDF_value(bsdf,c,r);
			float	dom = getBSDF_incohm(bsdf,c);
			FVECT	v;
			
			if (f <= .0) {
				nneg += (f < -FTINY);
				continue;
			}
			if (dom <= .0) {
				nbadohm++;
				continue;
			}
			if (!getBSDF_incvec(v,bsdf,c) || v[2] > FTINY)
				error(USER, "illegal incoming BTDF direction");
				
			mp[0] = mp[1] = mp[2] = f * dom * -v[2];
		}
	if (nbadohm | nneg) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles");
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Sum together a set of images and write result to stdout */
static int
sum_images(const char *fspec, const CMATRIX *cv)
{
	int	myDT = DTfromHeader;
	CMATRIX	*pmat;
	COLOR	*scanline;
	int	myXR, myYR;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = 0; i < cv->nrows; i++) {
		const COLORV	*scv = cv_lval(cv,i);
		char		fname[1024];
		FILE		*fp;
		int		dt, xr, yr;
		COLORV		*psp;
							/* open next picture */
		sprintf(fname, fspec, i);
		if ((fp = fopen(fname, "r")) == NULL) {
			sprintf(errmsg, "cannot open picture '%s'", fname);
			error(SYSTEM, errmsg);
		}
		SET_FILE_BINARY(fp);
		dt = getDTfromHeader(fp);
		if ((dt != DTrgbe) & (dt != DTxyze) ||
				!fscnresolu(&xr, &yr, fp)) {
			sprintf(errmsg, "file '%s' not a picture", fname);
			error(USER, errmsg);
		}
		if (myDT == DTfromHeader) {		/* on first one */
			myDT = dt;
			myXR = xr; myYR = yr;
			scanline = (COLOR *)malloc(sizeof(COLOR)*myXR);
			if (scanline == NULL)
				error(SYSTEM, "out of memory in sum_images()");
			pmat = cm_alloc(myYR, myXR);
			memset(pmat->cmem, 0, sizeof(COLOR)*myXR*myYR);
							/* finish header */
			fputformat(myDT==DTrgbe ? COLRFMT : CIEFMT, stdout);
			fputc('\n', stdout);
			fprtresolu(myXR, myYR, stdout);
			fflush(stdout);
		} else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
					fname);
			error(USER, errmsg);
		}
		psp = pmat->cmem;
		for (y = 0; y < yr; y++) {		/* read it in */
			int	x;
			if (freadscan(scanline, xr, fp) < 0) {
				sprintf(errmsg, "error reading picture '%s'",
						fname);
				error(SYSTEM, errmsg);
			}
							/* sum in scanline */
			for (x = 0; x < xr; x++, psp += 3) {
				multcolor(scanline[x], scv);
				addcolor(psp, scanline[x]);
			}
		}
		fclose(fp);				/* done this picture */
	}
	free(scanline);
							/* write scanlines */
	for (y = 0; y < myYR; y++)
		if (fwritescan((COLOR *)cm_lval(pmat, y, 0), myXR, stdout) < 0)
			return(0);
	cm_free(pmat);					/* all done */
	return(fflush(stdout) == 0);
}

/* check to see if a string contains a %d specification */
int
hasDecimalSpec(const char *s)
{
	while (*s && *s != '%')
		s++;
	if (!*s)
		return(0);
	do
		++s;
	while (isdigit(*s));

	return(*s == 'd');
}

int
main(int argc, char *argv[])
{
	CMATRIX			*tvec, *Dmat, *Tmat, *ivec, *cvec;
	struct BSDF_data	*btdf;

	progname = argv[0];

	if ((argc < 4) | (argc > 5)) {
		fprintf(stderr, "Usage: %s Vspec Tbsdf.xml Dmat.dat [tregvec]\n",
				progname);
		return(1);
	}
	tvec = cm_load(argv[4], 0, 1, DTascii);	/* argv[4]==NULL iff argc==4 */
	Dmat = cm_load(argv[3], 0, tvec->nrows, DTfromHeader);
	btdf = load_BSDF(argv[2]);
	if (btdf == NULL)
		return(1);
	if (btdf->ninc != Dmat->nrows) {
		sprintf(errmsg, "Incoming BTDF dir (%d) mismatch to D (%d)",
				btdf->ninc, Dmat->nrows);
		error(USER, errmsg);
	}
						/* multiply vector through */
	ivec = cm_multiply(Dmat, tvec);
	cm_free(Dmat); cm_free(tvec);
	Tmat = cm_bsdf(btdf);			/* convert BTDF to matrix */
	free_BSDF(btdf);
	cvec = cm_multiply(Tmat, ivec);		/* cvec = component vector */
	cm_free(Tmat); cm_free(ivec);
	if (hasDecimalSpec(argv[1])) {		/* generating image */
		SET_FILE_BINARY(stdout);
		newheader("RADIANCE", stdout);
		printargs(argc, argv, stdout);
		fputnow(stdout);
		if (!sum_images(argv[1], cvec))
			return(1);
	} else {				/* generating vector */
		CMATRIX	*Vmat = cm_load(argv[1], 0, cvec->nrows, DTfromHeader);
		CMATRIX	*rvec = cm_multiply(Vmat, cvec);
		cm_free(Vmat);
		cm_print(rvec, stdout);
		cm_free(rvec);
	}
	cm_free(cvec);				/* final clean-up */
	return(0);
}
