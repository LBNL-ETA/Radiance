#ifndef lint
static const char RCSid[] = "$Id: cmatrix.c,v 2.15 2016/02/17 23:26:06 greg Exp $";
#endif
/*
 * Color matrix routines.
 *
 *	G. Ward
 */

#include <ctype.h>
#include "standard.h"
#include "cmatrix.h"
#include "platform.h"
#include "paths.h"
#include "resolu.h"

const char	*cm_fmt_id[] = {
			"unknown", "ascii", COLRFMT, CIEFMT,
			"float", "double"
		};

const int	cm_elem_size[] = {
			0, 0, 3*sizeof(float), 3*sizeof(double), 4, 4
		};

/* Allocate a color coefficient matrix */
CMATRIX *
cm_alloc(int nrows, int ncols)
{
	CMATRIX	*cm;

	if ((nrows <= 0) | (ncols <= 0))
		error(USER, "attempt to create empty matrix");
	cm = (CMATRIX *)malloc(sizeof(CMATRIX) +
				sizeof(COLOR)*(nrows*ncols - 1));
	if (cm == NULL)
		error(SYSTEM, "out of memory in cm_alloc()");
	cm->nrows = nrows;
	cm->ncols = ncols;
	return(cm);
}

static void
adjacent_ra_sizes(size_t bounds[2], size_t target)
{
	bounds[0] = 0; bounds[1] = 2048;
	while (bounds[1] < target) {
		bounds[0] = bounds[1];
		bounds[1] += bounds[1]>>1;
	}
}

/* Resize color coefficient matrix */
CMATRIX *
cm_resize(CMATRIX *cm, int nrows)
{
	size_t	old_size, new_size, ra_bounds[2];

	if (nrows == cm->nrows)
		return(cm);
	if (nrows <= 0) {
		cm_free(cm);
		return(NULL);
	}
	old_size = sizeof(CMATRIX) + sizeof(COLOR)*(cm->nrows*cm->ncols - 1);
	adjacent_ra_sizes(ra_bounds, old_size);
	new_size = sizeof(CMATRIX) + sizeof(COLOR)*(nrows*cm->ncols - 1);
	if (nrows < cm->nrows ? new_size <= ra_bounds[0] :
				new_size > ra_bounds[1]) {
		adjacent_ra_sizes(ra_bounds, new_size);
		cm = (CMATRIX *)realloc(cm, ra_bounds[1]);
		if (cm == NULL)
			error(SYSTEM, "out of memory in cm_resize()");
	}
	cm->nrows = nrows;
	return(cm);
}

typedef struct {
	int	dtype;		/* data type */
	int	nrows, ncols;	/* matrix size */
	char	*err;		/* error message */
} CMINFO;		/* header info record */

static int
get_cminfo(char *s, void *p)
{
	CMINFO	*ip = (CMINFO *)p;
	char	fmt[32];
	int	i;

	if (!strncmp(s, "NCOMP=", 6) && atoi(s+6) != 3) {
		ip->err = "unexpected # components (must be 3)";
		return(-1);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		ip->nrows = atoi(s+6);
		return(0);
	}
	if (!strncmp(s, "NCOLS=", 6)) {
		ip->ncols = atoi(s+6);
		return(0);
	}
	if (!formatval(fmt, s))
		return(0);
	for (i = 1; i < DTend; i++)
		if (!strcmp(fmt, cm_fmt_id[i]))
			ip->dtype = i;
	return(0);
}

/* Load header to obtain/check data type and number of columns */
char *
cm_getheader(int *dt, int *nr, int *nc, FILE *fp)
{
	CMINFO	cmi;
						/* read header */
	cmi.dtype = DTfromHeader;
	cmi.nrows = cmi.ncols = 0;
	cmi.err = "unexpected EOF in header";
	if (getheader(fp, get_cminfo, &cmi) < 0)
		return(cmi.err);
	if (dt != NULL) {			/* get/check data type? */
		if (cmi.dtype == DTfromHeader) {
			if (*dt == DTfromHeader)
				return("missing/unknown data format in header");
		} else if (*dt == DTfromHeader)
			*dt = cmi.dtype;
		else if (*dt != cmi.dtype)
			return("unexpected data format in header");
	}
	if (nr != NULL) {			/* get/check #rows? */
		if (*nr <= 0)
			*nr = cmi.nrows;
		else if ((cmi.nrows > 0) & (*nr != cmi.nrows))
			return("unexpected row count in header");
	}
	if (nc != NULL) {			/* get/check #columns? */
		if (*nc <= 0)
			*nc = cmi.ncols;
		else if ((cmi.ncols > 0) & (*nc != cmi.ncols))
			return("unexpected column count in header");
	}
	return(NULL);
}

/* Allocate and load a matrix from the given input (or stdin if NULL) */
CMATRIX *
cm_load(const char *inspec, int nrows, int ncols, int dtype)
{
	const int	ROWINC = 2048;
	FILE		*fp = stdin;
	CMATRIX		*cm;

	if (inspec == NULL)
		inspec = "<stdin>";
	else if (inspec[0] == '!') {
		fp = popen(inspec+1, "r");
		if (fp == NULL) {
			sprintf(errmsg, "cannot start command '%s'", inspec);
			error(SYSTEM, errmsg);
		}
	} else if ((fp = fopen(inspec, "r")) == NULL) {
		sprintf(errmsg, "cannot open file '%s'", inspec);
		error(SYSTEM, errmsg);
	}
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (dtype != DTascii)
		SET_FILE_BINARY(fp);		/* doesn't really work */
	if (!dtype | !ncols) {			/* expecting header? */
		char	*err = cm_getheader(&dtype, &nrows, &ncols, fp);
		if (err != NULL)
			error(USER, err);
		if (ncols <= 0)
			error(USER, "unspecified number of columns");
	}
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
		if ((dtype != DTascii) & (fp != stdin) & (inspec[0] != '!')) {
			long	startpos = ftell(fp);
			if (fseek(fp, 0L, SEEK_END) == 0) {
				long	endpos = ftell(fp);
				long	elemsiz = 3*(dtype==DTfloat ?
					    sizeof(float) : sizeof(double));

				if ((endpos - startpos) % (ncols*elemsiz)) {
					sprintf(errmsg,
					"improper length for binary file '%s'",
							inspec);
					error(USER, errmsg);
				}
				guessrows = (endpos - startpos)/(ncols*elemsiz);
				if (fseek(fp, startpos, SEEK_SET) < 0) {
					sprintf(errmsg,
						"fseek() error on file '%s'",
							inspec);
					error(SYSTEM, errmsg);
				}
				nrows = guessrows;	/* we're confident */
			}
		}
		cm = cm_alloc(guessrows, ncols);
	} else
		cm = cm_alloc(nrows, ncols);
	if (cm == NULL)					/* XXX never happens */
		return(NULL);
	if (dtype == DTascii) {				/* read text file */
		int	maxrow = (nrows > 0 ? nrows : 32000);
		int	r, c;
		for (r = 0; r < maxrow; r++) {
		    if (r >= cm->nrows)			/* need more space? */
			cm = cm_resize(cm, cm->nrows+ROWINC);
		    for (c = 0; c < ncols; c++) {
		        COLORV	*cv = cm_lval(cm,r,c);
			if (fscanf(fp, COLSPEC, cv, cv+1, cv+2) != 3) {
				if ((nrows <= 0) & (r > 0) & !c) {
					cm = cm_resize(cm, maxrow=r);
					break;
				} else
					goto EOFerror;
			}
		    }
		}
		while ((c = getc(fp)) != EOF)
			if (!isspace(c)) {
				sprintf(errmsg,
				"unexpected data at end of ascii input '%s'",
						inspec);
				error(WARNING, errmsg);
				break;
			}
	} else {					/* read binary file */
		if (sizeof(COLOR) == cm_elem_size[dtype]) {
			int	nread = 0;
			do {				/* read all we can */
				nread += fread(cm->cmem + 3*nread,
						sizeof(COLOR),
						cm->nrows*cm->ncols - nread,
						fp);
				if (nrows <= 0) {	/* unknown length */
					if (nread == cm->nrows*cm->ncols)
							/* need more space? */
						cm = cm_resize(cm, cm->nrows+ROWINC);
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
		if (fgetc(fp) != EOF) {
				sprintf(errmsg,
				"unexpected data at end of binary input '%s'",
						inspec);
				error(WARNING, errmsg);
		}
	}
	if (fp != stdin) {
		if (inspec[0] != '!')
			fclose(fp);
		else if (pclose(fp)) {
			sprintf(errmsg, "error running command '%s'", inspec);
			error(WARNING, errmsg);
		}
	}
#ifdef getc_unlocked
	else
		funlockfile(fp);
#endif
	return(cm);
EOFerror:
	sprintf(errmsg, "unexpected EOF reading %s", inspec);
	error(USER, errmsg);
not_handled:
	error(INTERNAL, "unhandled data size or length in cm_load()");
	return(NULL);	/* gratis return */
}

/* Extract a column vector from a matrix */
CMATRIX *
cm_column(const CMATRIX *cm, int c)
{
	CMATRIX	*cvr;
	int	dr;

	if ((c < 0) | (c >= cm->ncols))
		error(INTERNAL, "column requested outside matrix");
	cvr = cm_alloc(cm->nrows, 1);
	if (cvr == NULL)
		return(NULL);
	for (dr = 0; dr < cm->nrows; dr++) {
		const COLORV	*sp = cm_lval(cm,dr,c);
		COLORV		*dp = cv_lval(cvr,dr);
		dp[0] = sp[0];
		dp[1] = sp[1];
		dp[2] = sp[2];
	}
	return(cvr);
}

/* Scale a matrix by a single value */
CMATRIX *
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
CMATRIX *
cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2)
{
	char	*rowcheck=NULL, *colcheck=NULL;
	CMATRIX	*cmr;
	int	dr, dc, i;

	if ((cm1->ncols <= 0) | (cm1->ncols != cm2->nrows))
		error(INTERNAL, "matrix dimension mismatch in cm_multiply()");
	cmr = cm_alloc(cm1->nrows, cm2->ncols);
	if (cmr == NULL)
		return(NULL);
				/* optimization: check for zero rows & cols */
	if (((cm1->nrows > 5) | (cm2->ncols > 5)) & (cm1->ncols > 5)) {
		static const COLOR	czero;
		rowcheck = (char *)calloc(cmr->nrows, 1);
		for (dr = cm1->nrows*(rowcheck != NULL); dr--; )
		    for (dc = cm1->ncols; dc--; )
			if (memcmp(cm_lval(cm1,dr,dc), czero, sizeof(COLOR))) {
				rowcheck[dr] = 1;
				break;
			}
		colcheck = (char *)calloc(cmr->ncols, 1);
		for (dc = cm2->ncols*(colcheck != NULL); dc--; )
		    for (dr = cm2->nrows; dr--; )
			if (memcmp(cm_lval(cm2,dr,dc), czero, sizeof(COLOR))) {
				colcheck[dc] = 1;
				break;
			}
	}
	for (dr = 0; dr < cmr->nrows; dr++)
	    for (dc = 0; dc < cmr->ncols; dc++) {
		COLORV	*dp = cm_lval(cmr,dr,dc);
		double	res[3];
		dp[0] = dp[1] = dp[2] = 0;
		if (rowcheck != NULL && !rowcheck[dr])
			continue;
		if (colcheck != NULL && !colcheck[dc])
			continue;
		res[0] = res[1] = res[2] = 0;
		for (i = 0; i < cm1->ncols; i++) {
		    const COLORV	*cp1 = cm_lval(cm1,dr,i);
		    const COLORV	*cp2 = cm_lval(cm2,i,dc);
		    res[0] += cp1[0] * cp2[0];
		    res[1] += cp1[1] * cp2[1];
		    res[2] += cp1[2] * cp2[2];
		}
		copycolor(dp, res);
	    }
	if (rowcheck != NULL) free(rowcheck);
	if (colcheck != NULL) free(colcheck);
	return(cmr);
}

/* write out matrix to file (precede by resolution string if picture) */
int
cm_write(const CMATRIX *cm, int dtype, FILE *fp)
{
	static const char	tabEOL[2] = {'\t','\n'};
	const COLORV		*mp = cm->cmem;
	int			r, c;

	switch (dtype) {
	case DTascii:
		for (r = 0; r < cm->nrows; r++)
			for (c = 0; c < cm->ncols; c++, mp += 3)
				fprintf(fp, "%.6e %.6e %.6e%c",
						mp[0], mp[1], mp[2],
						tabEOL[c >= cm->ncols-1]);
		break;
	case DTfloat:
	case DTdouble:
		if (sizeof(COLOR) == cm_elem_size[dtype]) {
			r = cm->ncols*cm->nrows;
			while (r > 0) {
				c = fwrite(mp, sizeof(COLOR), r, fp);
				if (c <= 0)
					return(0);
				mp += 3*c;
				r -= c;
			}
		} else if (dtype == DTdouble) {
			double	dc[3];
			r = cm->ncols*cm->nrows;
			while (r--) {
				copycolor(dc, mp);
				if (fwrite(dc, sizeof(double), 3, fp) != 3)
					return(0);
				mp += 3;
			}
		} else /* dtype == DTfloat */ {
			float	fc[3];
			r = cm->ncols*cm->nrows;
			while (r--) {
				copycolor(fc, mp);
				if (fwrite(fc, sizeof(float), 3, fp) != 3)
					return(0);
				mp += 3;
			}
		}
		break;
	case DTrgbe:
	case DTxyze:
		fprtresolu(cm->ncols, cm->nrows, fp);
		for (r = 0; r < cm->nrows; r++, mp += 3*cm->ncols)
			if (fwritescan((COLOR *)mp, cm->ncols, fp) < 0)
				return(0);
		break;
	default:
		fputs("Unsupported data type in cm_write()!\n", stderr);
		return(0);
	}
	return(fflush(fp) == 0);
}
