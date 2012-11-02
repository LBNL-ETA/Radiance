#ifndef lint
static const char RCSid[] = "$Id: dctimestep.c,v 2.23 2012/11/02 00:26:55 greg Exp $";
#endif
/*
 * Compute time-step result using Daylight Coefficient method.
 *
 *	G. Ward
 */

#include <ctype.h>
#include "standard.h"
#include "platform.h"
#include "paths.h"
#include "color.h"
#include "resolu.h"
#include "bsdf.h"
#include "bsdf_m.h"

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
		error(USER, "attempt to create empty matrix");
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

	if (ncols <= 0)
		error(USER, "Non-positive number of columns");
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
		    if (r >= cm->nrows)			/* need more space? */
			cm = cm_resize(cm, 2*cm->nrows);
		    for (c = 0; c < ncols; c++) {
		        COLORV	*cv = cm_lval(cm,r,c);
			if (fscanf(fp, COLSPEC, cv, cv+1, cv+2) != 3)
				if ((nrows <= 0) & (r > 0) & !c) {
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
					else if (nread && !(nread % cm->ncols))
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

/* Scale a matrix by a single value */
static CMATRIX *
cm_scale(const CMATRIX *cm1, const COLOR sca)
{
	CMATRIX	*cmr;
	int	dr, dc;

	cmr = cm_alloc(cm1->nrows, cm1->ncols);
	if (cmr == NULL)
		return(NULL);
	for (dr = 0; dr < cmr->nrows; dr++)
	    for (dc = 0; dc < cmr->ncols; dc++) {
	        const COLORV	*sp = cm_lval(cm1,dr,dc);
		COLORV		*dp = cm_lval(cmr,dr,dc);
		dp[0] = sp[0] * sca[0];
		dp[1] = sp[1] * sca[1];
		dp[2] = sp[2] * sca[2];
	    }
	return(cmr);
}

/* Multiply two matrices (or a matrix and a vector) and allocate the result */
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

/* Convert a BSDF to our matrix representation */
static CMATRIX *
cm_bsdf(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->nout, bsdf->ninc);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const double	dom = mBSDF_incohm(bsdf,c);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			float	f = mBSDF_value(bsdf,c,r);
			COLORV	*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Convert between input and output indices for reciprocity */
static int
recip_out_from_in(const SDMat *bsdf, int in_recip)
{
	FVECT	v;

	if (!mBSDF_incvec(v, bsdf, in_recip+.5))
		return(in_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_outndx(bsdf, v));
}

/* Convert between output and input indices for reciprocity */
static int
recip_in_from_out(const SDMat *bsdf, int out_recip)
{
	FVECT	v;

	if (!mBSDF_outvec(v, bsdf, out_recip+.5))
		return(out_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_incndx(bsdf, v));
}

/* Convert a BSDF to our matrix representation, applying reciprocity */
static CMATRIX *
cm_bsdf_recip(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->ninc, bsdf->nout);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const int	ro = recip_out_from_in(bsdf,c);
		const double	dom = mBSDF_outohm(bsdf,ro);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			const int	ri = recip_in_from_out(bsdf,r);
			float		f = mBSDF_value(bsdf,ri,ro);
			COLORV		*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Load and convert a matrix BSDF from the given XML file */
static CMATRIX *
cm_loadBSDF(char *fname, COLOR cLamb)
{
	CMATRIX		*Tmat;
	char		*fpath;
	int		recip;
	SDError		ec;
	SDData		myBSDF;
	SDSpectralDF	*tdf;
	COLOR		bsdfLamb, specCol;
					/* find path to BSDF file */
	fpath = getpath(fname, getrlibpath(), R_OK);
	if (fpath == NULL) {
		sprintf(errmsg, "cannot find BSDF file '%s'", fname);
		error(USER, errmsg);
	}
	SDclearBSDF(&myBSDF, fname);	/* load XML and check type */
	ec = SDloadFile(&myBSDF, fpath);
	if (ec)
		error(USER, transSDError(ec));
	ccy2rgb(&myBSDF.tLamb.spec, myBSDF.tLamb.cieY/PI, bsdfLamb);
	recip = (myBSDF.tb == NULL);
	tdf = recip ? myBSDF.tf : myBSDF.tb;
	if (tdf == NULL) {		/* no non-Lambertian transmission? */
		if (cLamb != NULL)
			copycolor(cLamb, bsdfLamb);
		SDfreeBSDF(&myBSDF);
		return(NULL);
	}
	if (tdf->ncomp != 1 || tdf->comp[0].func != &SDhandleMtx) {
		sprintf(errmsg, "unsupported BSDF '%s'", fpath);
		error(USER, errmsg);
	}
					/* convert BTDF to matrix */
	ccy2rgb(&tdf->comp[0].cspec[0], 1., specCol);
	Tmat = recip ? cm_bsdf_recip(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist)
			: cm_bsdf(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist);
	if (cLamb != NULL)		/* Lambertian is included */
		setcolor(cLamb, .0, .0, .0);
					/* free BSDF and return */
	SDfreeBSDF(&myBSDF);
	return(Tmat);
}

/* Sum together a set of images and write result to stdout */
static int
sum_images(const char *fspec, const CMATRIX *cv, FILE *fout)
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
							/* check for zero */
		if ((scv[RED] == 0) & (scv[GRN] == 0) & (scv[BLU] == 0) &&
				(myDT != DTfromHeader) | (i < cv->nrows-1))
			continue;
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
			fputformat(myDT==DTrgbe ? COLRFMT : CIEFMT, fout);
			fputc('\n', fout);
			fprtresolu(myXR, myYR, fout);
			fflush(fout);
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
		if (fwritescan((COLOR *)cm_lval(pmat, y, 0), myXR, fout) < 0)
			return(0);
	cm_free(pmat);					/* all done */
	return(fflush(fout) == 0);
}

/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char *s)
{
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
	CMATRIX			*cvec;		/* component vector result */

	progname = argv[0];

	if ((argc < 2) | (argc > 5)) {
		fprintf(stderr, "Usage: %s DCspec [tregvec]\n", progname);
		fprintf(stderr, "   or: %s Vspec Tbsdf.xml Dmat.dat [tregvec]\n",
					progname);
		return(1);
	}

	if (argc > 3) {				/* VTDs expression */
		CMATRIX	*svec, *Dmat, *Tmat, *ivec;
		COLOR	tLamb;
						/* get sky vector */
		svec = cm_load(argv[4], 0, 1, DTascii);
						/* load BSDF */
		Tmat = cm_loadBSDF(argv[2], tLamb);
						/* load Daylight matrix */
		Dmat = cm_load(argv[3], Tmat==NULL ? 0 : Tmat->ncols,
					svec->nrows, DTfromHeader);
						/* multiply vector through */
		ivec = cm_multiply(Dmat, svec);
		cm_free(Dmat); cm_free(svec);
		if (Tmat == NULL) {		/* diffuse only */
			cvec = cm_scale(ivec, tLamb);
		} else {			/* else apply BTDF matrix */
			cvec = cm_multiply(Tmat, ivec);
			cm_free(Tmat); 
		}
		cm_free(ivec);
	} else {				/* else just use sky vector */
		cvec = cm_load(argv[2], 0, 1, DTascii);
	}

	if (hasNumberFormat(argv[1])) {		/* generating image */
		SET_FILE_BINARY(stdout);
		newheader("RADIANCE", stdout);
		printargs(argc, argv, stdout);
		fputnow(stdout);
		if (!sum_images(argv[1], cvec, stdout))
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
