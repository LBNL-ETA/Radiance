#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Color matrix routines.
 *
 *	G. Ward
 */

#include <ctype.h>
#include "platform.h"
#include "standard.h"
#include "cmatrix.h"
#include "platform.h"
#include "paths.h"
#include "resolu.h"
#ifdef INTEL_DCTOPT
#include <immintrin.h> // AVX2 Instruction
#include <omp.h>
#endif /* INTEL_DCTOPT */

/* Allocate a color coefficient matrix */
CMATRIX *
cm_alloc(int nrows, int ncols)
{
	CMATRIX	*cm;

	if ((nrows <= 0) | (ncols <= 0))
		error(USER, "attempt to create empty matrix");
	cm = (CMATRIX *)malloc(sizeof(CMATRIX) +
				sizeof(COLOR)*((size_t)nrows*ncols - 1));
	if (!cm)
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

	if (!cm)
		return(NULL);
	if (nrows == cm->nrows)
		return(cm);
	if (nrows <= 0) {
		cm_free(cm);
		return(NULL);
	}
	old_size = sizeof(CMATRIX) + sizeof(COLOR)*((size_t)cm->nrows*cm->ncols - 1);
	adjacent_ra_sizes(ra_bounds, old_size);
	new_size = sizeof(CMATRIX) + sizeof(COLOR)*((size_t)nrows*cm->ncols - 1);
	if (nrows < cm->nrows ? new_size <= ra_bounds[0] :
				new_size > ra_bounds[1]) {
		adjacent_ra_sizes(ra_bounds, new_size);
		cm = (CMATRIX *)realloc(cm, ra_bounds[1]);
		if (!cm)
			error(SYSTEM, "out of memory in cm_resize()");
	}
	cm->nrows = nrows;
	return(cm);
}

typedef struct {
	int	dtype;		/* data type */
	int	need2swap;	/* need byte swap? */
	int	nrows, ncols;	/* matrix size */
	COLOR	expos;		/* exposure value */
	char	*err;		/* error message */
} CMINFO;		/* header info record */

static int
get_cminfo(char *s, void *p)
{
	CMINFO	*ip = (CMINFO *)p;
	char	fmt[MAXFMTLEN];
	int	i;

	if (isncomp(s) && ncompval(s) != 3) {
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
	if ((i = isbigendian(s)) >= 0) {
		ip->need2swap = (nativebigendian() != i);
		return(0);
	}
	if (isexpos(s)) {
		double	d = exposval(s);
		scalecolor(ip->expos, d);
		return(0);
	}
	if (iscolcor(s)) {
		COLOR	ctmp;
		colcorval(ctmp, s);
		multcolor(ip->expos, ctmp);
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
cm_getheader(int *dt, int *nr, int *nc, int *swp, COLOR scale, FILE *fp)
{
	CMINFO	cmi;
						/* read header */
	cmi.dtype = DTfromHeader;
	cmi.need2swap = 0;
	cmi.expos[0] = cmi.expos[1] = cmi.expos[2] = 1.f;
	cmi.nrows = cmi.ncols = 0;
	cmi.err = "unexpected EOF in header";
	if (getheader(fp, get_cminfo, &cmi) < 0)
		return(cmi.err);
	if (dt) {				/* get/check data type? */
		if (cmi.dtype == DTfromHeader) {
			if (*dt == DTfromHeader)
				return("missing/unknown data format in header");
		} else if (*dt == DTfromHeader)
			*dt = cmi.dtype;
		else if (*dt != cmi.dtype)
			return("unexpected data format in header");
	}
	if (nr) {				/* get/check #rows? */
		if (*nr <= 0)
			*nr = cmi.nrows;
		else if ((cmi.nrows > 0) & (*nr != cmi.nrows))
			return("unexpected row count in header");
	}
	if (nc) {				/* get/check #columns? */
		if (*nc <= 0)
			*nc = cmi.ncols;
		else if ((cmi.ncols > 0) & (*nc != cmi.ncols))
			return("unexpected column count in header");
	}
	if (swp)				/* get/check swap? */
		*swp = cmi.need2swap;
	if (scale) {				/* transfer exposure comp. */
		scale[0] = 1.f/cmi.expos[0];
		scale[1] = 1.f/cmi.expos[1];
		scale[2] = 1.f/cmi.expos[2];
	}
	return(NULL);
}

/* Allocate and load image data into matrix */
static CMATRIX *
cm_load_rgbe(FILE *fp, int nrows, int ncols)
{
	CMATRIX	*cm;
	COLORV	*mp;
						/* header already loaded */
	cm = cm_alloc(nrows, ncols);
	if (!cm)
		return(NULL);
	mp = cm->cmem;
	while (nrows--) {
		if (freadscan((COLOR *)mp, ncols, fp) < 0) {
			error(USER, "error reading color picture as matrix");
			cm_free(cm);
			return(NULL);
		}
		mp += 3*ncols;
	}					/* caller closes stream */
	return(cm);
}

#ifdef INTEL_DCTOPT

CMATRIX *
cm_alloc_u(int rows, int cols) {
	size_t total_elements = (size_t)rows * cols;
	size_t data_size = sizeof(float) * 3 * total_elements;

	size_t total_size = sizeof(CMATRIX) + data_size + 32;

#ifdef  _WIN32
	CMATRIX* m = (CMATRIX*)_aligned_malloc((total_size + 31) & ~31, 32);
#else
	CMATRIX* m = (CMATRIX*)aligned_alloc(32, (total_size + 31) & ~31);
#endif //  _WIN32
	if (m) {
		m->nrows = rows;
		m->ncols = cols;

		memset(m->cmem, 0, data_size + 32);
	}
	return m;
}


CMATRIX *
cm_resize_aligned(CMATRIX* cm, int nrows)
{
	CMATRIX* new_cm;
	size_t copy_bytes;

	if (!cm) return NULL;
	if (nrows == cm->nrows) return cm;
	if (nrows <= 0) {
		cm_free_u(cm);
		return NULL;
	}

	new_cm = cm_alloc_u(nrows, cm->ncols);
	if (!new_cm) {
		error(SYSTEM, "out of memory in cm_resize_aligned()");
		return NULL;
	}
	int copy_rows = (nrows < cm->nrows) ? nrows : cm->nrows;
	size_t total_copy_elements = (size_t)copy_rows * cm->ncols;
	size_t copy_data_size = sizeof(float) * 3 * total_copy_elements;
	copy_bytes = sizeof(CMATRIX) + copy_data_size - sizeof(float) * 3;
	memcpy(new_cm, cm, copy_bytes);
	cm_free_u(cm);
	return new_cm;
}


static CMATRIX*
cm_load_rgbe_u(FILE* fp, int nrows, int ncols)
{
	CMATRIX* cm;
	COLORV* mp;
	/* header already loaded */
	cm = cm_alloc_u(nrows, ncols);
	if (!cm)
		return(NULL);
	mp = cm->cmem;
	while (nrows--) {
		if (freadscan((COLOR*)mp, ncols, fp) < 0) {
			error(USER, "error reading color picture as matrix");
			cm_free(cm);
			return(NULL);
		}
		mp += 3 * ncols;
	}					/* caller closes stream */
	return(cm);
}


CMATRIX*
cm_load_u(const char* inspec, int nrows, int ncols, int dtype)
{
	const int	ROWINC = 2048;
	int		dimsOK = (dtype == DTascii) | (nrows > 0) && ncols;
	int		swap = 0;
	FILE* fp;
	COLOR		scale;
	CMATRIX* cm;

	if (!inspec)
		inspec = stdin_name;
	else if (!*inspec)
		return(NULL);
	if (inspec == stdin_name) {		/* reading from stdin? */
		fp = stdin;
	}
	else if (inspec[0] == '!') {
		fp = popen(inspec + 1, "r");
		if (!fp) {
			sprintf(errmsg, "cannot start command '%s'", inspec);
			error(SYSTEM, errmsg);
		}
	}
	else if (!(fp = fopen(inspec, "r"))) {
		sprintf(errmsg, "cannot open file '%s'", inspec);
		error(SYSTEM, errmsg);
	}
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (dtype != DTascii)
		SET_FILE_BINARY(fp);		/* doesn't really work */
	if (!dtype | !dimsOK) {			/* expecting header? */
		char* err = cm_getheader(&dtype, &nrows, &ncols, &swap, scale, fp);
		if (err)
			error(USER, err);
		dimsOK = ncols > 0 && (nrows > 0 ||
			(dtype != DTrgbe) & (dtype != DTxyze));
	}
	if (!dimsOK && !fscnresolu(&ncols, &nrows, fp))
		error(USER, "unspecified matrix size");
	switch (dtype) {
	case DTascii:
	case DTfloat:
	case DTdouble:
		break;
	case DTrgbe:
	case DTxyze:
		cm = cm_load_rgbe_u(fp, nrows, ncols);
		goto cleanup;
	default:
		error(USER, "unexpected data type in cm_load()");
	}
	if (nrows <= 0) {			/* don't know length? */
		int	guessrows = 147;	/* usually big enough */
		if (cm_elem_size[dtype] && (fp != stdin) & (inspec[0] != '!')) {
			long	startpos = ftell(fp);
			if (fseek(fp, 0L, SEEK_END) == 0) {
				long	rowsiz = (long)ncols * cm_elem_size[dtype];
				long	endpos = ftell(fp);

				if ((endpos - startpos) % rowsiz) {
					sprintf(errmsg,
						"improper length for binary file '%s'",
						inspec);
					error(USER, errmsg);
				}
				guessrows = (endpos - startpos) / rowsiz;
				if (fseek(fp, startpos, SEEK_SET) < 0) {
					sprintf(errmsg,
						"fseek() error on file '%s'",
						inspec);
					error(SYSTEM, errmsg);
				}
				nrows = guessrows;	/* we're confident */
			}
		}
		cm = cm_alloc_u(guessrows, ncols);
	}
	else
		cm = cm_alloc_u(nrows, ncols);
	if (!cm)					/* XXX never happens */
		return(NULL);
	if (dtype == DTascii) {				/* read text file */
		int	maxrow = (nrows > 0 ? nrows : 32000);
		int	r, c;
		for (r = 0; r < maxrow; r++) {
			if (r >= cm->nrows)			/* need more space? */
				cm = cm_resize_aligned(cm, cm->nrows + ROWINC);
			for (c = 0; c < ncols; c++) {
				COLORV* cv = cm_lval(cm, r, c);
				if (fscanf(fp, COLSPEC, cv, cv + 1, cv + 2) != 3) {
					if ((nrows <= 0) & (r > 0) & !c) {
						cm = cm_resize_aligned(cm, maxrow = r);
						break;
					}
					else
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
	}
	else {					/* read binary file */
		if (sizeof(COLOR) == cm_elem_size[dtype]) {
			size_t	nread = 0;
			do {				/* read all we can */
				nread += getbinary(cm->cmem + 3 * nread,
					sizeof(COLOR),
					(size_t)cm->nrows * cm->ncols - nread,
					fp);
				if (nrows <= 0) {	/* unknown length */
					if (nread == (size_t)cm->nrows * cm->ncols)
						/* need more space? */
						cm = cm_resize_aligned(cm, cm->nrows + ROWINC);
					else if (nread && !(nread % cm->ncols))
						/* seem to be  done */
						cm = cm_resize_aligned(cm, nread / cm->ncols);
					else		/* ended mid-row */
						goto EOFerror;
				}
				else if (nread < (size_t)cm->nrows * cm->ncols)
					goto EOFerror;
			} while (nread < (size_t)cm->nrows * cm->ncols);

			if (swap) {
				if (sizeof(COLORV) == 4)
					swap32((char*)cm->cmem,
						3 * (size_t)cm->nrows * cm->ncols);
				else /* sizeof(COLORV) == 8 */
					swap64((char*)cm->cmem,
						3 * (size_t)cm->nrows * cm->ncols);
			}
		}
		else if (dtype == DTdouble) {
			double	dc[3];			/* load from double */
			COLORV* cvp = cm->cmem;
			size_t	n = (size_t)nrows * ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (getbinary(dc, sizeof(double), 3, fp) != 3)
					goto EOFerror;
				if (swap) swap64((char*)dc, 3);
				copycolor(cvp, dc);
				cvp += 3;
			}
		}
		else /* dtype == DTfloat */ {
			float	fc[3];			/* load from float */
			COLORV* cvp = cm->cmem;
			size_t	n = (size_t)nrows * ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (getbinary(fc, sizeof(float), 3, fp) != 3)
					goto EOFerror;
				if (swap) swap32((char*)fc, 3);
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
cleanup:
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
	if ((scale[0] < .99) | (scale[0] > 1.01) |
		(scale[1] < .99) | (scale[1] > 1.01) |
		(scale[2] < .99) | (scale[2] > 1.01)) {
		size_t	n = (size_t)ncols * nrows;
		COLORV* mp = cm->cmem;
		while (n--) {		/* apply exposure scaling */
			*mp++ *= scale[0];
			*mp++ *= scale[1];
			*mp++ *= scale[2];
		}
	}
	return(cm);
EOFerror:
	sprintf(errmsg, "unexpected EOF reading %s", inspec);
	error(USER, errmsg);
	return(NULL);
not_handled:
	error(INTERNAL, "unhandled data size or length in cm_load()");
	return(NULL);	/* gratis return */
}

CMATRIX *
cm_column_u(const CMATRIX* cm, int c)
{
	CMATRIX* cvr;
	int	dr;

	if (!cm)
		return(NULL);
	if ((c < 0) | (c >= cm->ncols))
		error(INTERNAL, "column requested outside matrix");
	cvr = cm_alloc_u(cm->nrows, 1);
	if (!cvr)
		return(NULL);
	for (dr = 0; dr < cm->nrows; dr++) {
		const COLORV* sp = cm_lval(cm, dr, c);
		COLORV* dp = cv_lval(cvr, dr);
		dp[0] = sp[0];
		dp[1] = sp[1];
		dp[2] = sp[2];
	}
	return(cvr);
}

CMATRIX *
cm_get_batch(const CMATRIX* src, int start_col, int num_cols) {
	if (start_col + num_cols > src->ncols)
		num_cols = src->ncols - start_col;

	CMATRIX* dst = cm_alloc_u(src->nrows, num_cols);
	if (!dst) return NULL;

	for (int i = 0; i < src->nrows; i++) {

		float* src_ptr = (float*)cm_lval(src, i, start_col);
		float* dst_ptr = (float*)cm_lval(dst, i, 0);

		memcpy(dst_ptr, src_ptr, sizeof(float) * 3 * num_cols);
	}
	return dst;
}

CMATRIX *
cm_multiply_avx2(const CMATRIX* cm1, const CMATRIX* cm2) {
	int R = cm1->nrows, K = cm1->ncols, C = cm2->ncols;
	CMATRIX* cmr = cm_alloc_u(R, C);
	if (!cmr) return NULL;
	int i = 0;
#pragma omp parallel for schedule(dynamic, 4)
	for (i = 0; i < R - 1; i += 2) {
		float* r_a0 = (float*)cm_lval(cm1, i, 0);
		float* r_a1 = (float*)cm_lval(cm1, i + 1, 0);
		float* r_c0 = (float*)cm_lval(cmr, i, 0);
		float* r_c1 = (float*)cm_lval(cmr, i + 1, 0);

		for (int k = 0; k < K; k++) {
			float ra0 = r_a0[k * 3 + 0], ga0 = r_a0[k * 3 + 1], ba0 = r_a0[k * 3 + 2];
			float ra1 = r_a1[k * 3 + 0], ga1 = r_a1[k * 3 + 1], ba1 = r_a1[k * 3 + 2];

			if (ra0 == 0 && ga0 == 0 && ba0 == 0 && ra1 == 0 && ga1 == 0 && ba1 == 0) continue;

			__m256 v0_1 = _mm256_setr_ps(ra0, ga0, ba0, ra0, ga0, ba0, ra0, ga0);
			__m256 v0_2 = _mm256_setr_ps(ba0, ra0, ga0, ba0, ra0, ga0, ba0, ra0);
			__m256 v0_3 = _mm256_setr_ps(ga0, ba0, ra0, ga0, ba0, ra0, ga0, ba0);

			__m256 v1_1 = _mm256_setr_ps(ra1, ga1, ba1, ra1, ga1, ba1, ra1, ga1);
			__m256 v1_2 = _mm256_setr_ps(ba1, ra1, ga1, ba1, ra1, ga1, ba1, ra1);
			__m256 v1_3 = _mm256_setr_ps(ga1, ba1, ra1, ga1, ba1, ra1, ga1, ba1);

			float* row_b = (float*)cm_lval(cm2, k, 0);
			int j = 0, limit = C * 3;

			for (; j <= limit - 24; j += 24) {
				__m256 vb1 = _mm256_loadu_ps(&row_b[j]);
				__m256 vb2 = _mm256_loadu_ps(&row_b[j + 8]);
				__m256 vb3 = _mm256_loadu_ps(&row_b[j + 16]);

				__m256 vc0_1 = _mm256_loadu_ps(&r_c0[j]);
				__m256 vc0_2 = _mm256_loadu_ps(&r_c0[j + 8]);
				__m256 vc0_3 = _mm256_loadu_ps(&r_c0[j + 16]);
				_mm256_storeu_ps(&r_c0[j], _mm256_fmadd_ps(vb1, v0_1, vc0_1));
				_mm256_storeu_ps(&r_c0[j + 8], _mm256_fmadd_ps(vb2, v0_2, vc0_2));
				_mm256_storeu_ps(&r_c0[j + 16], _mm256_fmadd_ps(vb3, v0_3, vc0_3));


				__m256 vc1_1 = _mm256_loadu_ps(&r_c1[j]);
				__m256 vc1_2 = _mm256_loadu_ps(&r_c1[j + 8]);
				__m256 vc1_3 = _mm256_loadu_ps(&r_c1[j + 16]);
				_mm256_storeu_ps(&r_c1[j], _mm256_fmadd_ps(vb1, v1_1, vc1_1));
				_mm256_storeu_ps(&r_c1[j + 8], _mm256_fmadd_ps(vb2, v1_2, vc1_2));
				_mm256_storeu_ps(&r_c1[j + 16], _mm256_fmadd_ps(vb3, v1_3, vc1_3));
			}

			for (; j < limit; j += 3) {
				r_c0[j + 0] += ra0 * row_b[j + 0]; r_c0[j + 1] += ga0 * row_b[j + 1]; r_c0[j + 2] += ba0 * row_b[j + 2];
				r_c1[j + 0] += ra1 * row_b[j + 0]; r_c1[j + 1] += ga1 * row_b[j + 1]; r_c1[j + 2] += ba1 * row_b[j + 2];
			}
		}
	}


	if (R % 2 != 0) {
		int last_i = R - 1;
		float* r_a = (float*)cm_lval(cm1, last_i, 0);
		float* r_c = (float*)cm_lval(cmr, last_i, 0);

		for (int k = 0; k < K; k++) {
			float ra = r_a[k * 3 + 0], ga = r_a[k * 3 + 1], ba = r_a[k * 3 + 2];
			if (ra == 0 && ga == 0 && ba == 0) continue;

			__m256 va1 = _mm256_setr_ps(ra, ga, ba, ra, ga, ba, ra, ga);
			__m256 va2 = _mm256_setr_ps(ba, ra, ga, ba, ra, ga, ba, ra);
			__m256 va3 = _mm256_setr_ps(ga, ba, ra, ga, ba, ra, ga, ba);

			float* row_b = (float*)cm_lval(cm2, k, 0);
			int j = 0, limit = C * 3;
			for (; j <= limit - 24; j += 24) {
				__m256 vb1 = _mm256_loadu_ps(&row_b[j]);
				__m256 vb2 = _mm256_loadu_ps(&row_b[j + 8]);
				__m256 vb3 = _mm256_loadu_ps(&row_b[j + 16]);

				_mm256_storeu_ps(&r_c[j], _mm256_fmadd_ps(vb1, va1, _mm256_loadu_ps(&r_c[j])));
				_mm256_storeu_ps(&r_c[j + 8], _mm256_fmadd_ps(vb2, va2, _mm256_loadu_ps(&r_c[j + 8])));
				_mm256_storeu_ps(&r_c[j + 16], _mm256_fmadd_ps(vb3, va3, _mm256_loadu_ps(&r_c[j + 16])));
			}

			for (; j < limit; j += 3) {
				r_c[j + 0] += ra * row_b[j + 0];
				r_c[j + 1] += ga * row_b[j + 1];
				r_c[j + 2] += ba * row_b[j + 2];
			}
		}
	}

	return cmr;
}

CMATRIX *
cm_multiply_smart(const CMATRIX* basis, const CMATRIX* cv_batch) {
	int R = basis->nrows;
	int K = basis->ncols;
	int C = cv_batch->ncols; // Batch Size

	int* is_zero_map = (int*)malloc(C * sizeof(int));
	int* nonzero_map = (int*)malloc(C * sizeof(int));
	int num_nonzero = 0;

	for (int c = 0; c < C; c++) {
		int zero_flag = 1;
		for (int k = 0; k < K; k++) {
			float* val = (float*)cm_lval(cv_batch, k, c);
			if (val[0] != 0.0f || val[1] != 0.0f || val[2] != 0.0f) {
				zero_flag = 0;
				break;
			}
		}
		is_zero_map[c] = zero_flag;
		if (!zero_flag) {
			nonzero_map[num_nonzero++] = c;
		}
	}

	if (num_nonzero == 0) {
		CMATRIX* final_res = cm_alloc_u(R, C);
		free(is_zero_map); free(nonzero_map);
		return final_res;
	}

	if (num_nonzero == C) {
		free(is_zero_map); free(nonzero_map);
		return cm_multiply_avx2(basis, cv_batch);
	}
	CMATRIX* compact_cv = cm_alloc_u(K, num_nonzero);
	for (int nz = 0; nz < num_nonzero; nz++) {
		int original_col = nonzero_map[nz];
		for (int k = 0; k < K; k++) {
			float* src = (float*)cm_lval(cv_batch, k, original_col);
			float* dst = (float*)cm_lval(compact_cv, k, nz);
			dst[0] = src[0]; dst[1] = src[1]; dst[2] = src[2];
		}
	}

	CMATRIX* compact_res = cm_multiply_avx2(basis, compact_cv);
	cm_free_u(compact_cv);
	if (!compact_res) {
		free(is_zero_map); free(nonzero_map);
		return NULL;
	}
	CMATRIX* final_res = cm_alloc_u(R, C);
	if (final_res) {
		for (int r = 0; r < R; r++) {
			float* src_row = (float*)cm_lval(compact_res, r, 0);
			float* dst_row = (float*)cm_lval(final_res, r, 0);

			int nz_idx = 0;
			for (int c = 0; c < C; c++) {
				if (!is_zero_map[c]) {

					dst_row[c * 3 + 0] = src_row[nz_idx * 3 + 0];
					dst_row[c * 3 + 1] = src_row[nz_idx * 3 + 1];
					dst_row[c * 3 + 2] = src_row[nz_idx * 3 + 2];
					nz_idx++;
				}
			}
		}
	}

	cm_free_u(compact_res);
	free(is_zero_map);
	free(nonzero_map);
	return final_res;
}

#endif /* INTEL_DCTOPT */

/* Allocate and load a matrix from the given input (or stdin if NULL) */
CMATRIX *
cm_load(const char *inspec, int nrows, int ncols, int dtype)
{
	const int	ROWINC = 2048;
	int		dimsOK = (dtype == DTascii) | (nrows > 0) && ncols;
	int		swap = 0;
	FILE		*fp;
	COLOR		scale;
	CMATRIX		*cm;

	if (!inspec)
		inspec = stdin_name;
	else if (!*inspec)
		return(NULL);
	if (inspec == stdin_name) {		/* reading from stdin? */
		fp = stdin;
	} else if (inspec[0] == '!') {
		fp = popen(inspec+1, "r");
		if (!fp) {
			sprintf(errmsg, "cannot start command '%s'", inspec);
			error(SYSTEM, errmsg);
		}
	} else if (!(fp = fopen(inspec, "r"))) {
		sprintf(errmsg, "cannot open file '%s'", inspec);
		error(SYSTEM, errmsg);
	}
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (dtype != DTascii)
		SET_FILE_BINARY(fp);		/* doesn't really work */
	if (!dtype | !dimsOK) {			/* expecting header? */
		char	*err = cm_getheader(&dtype, &nrows, &ncols, &swap, scale, fp);
		if (err)
			error(USER, err);
		dimsOK = ncols > 0 && ( nrows > 0 ||
				(dtype != DTrgbe) & (dtype != DTxyze) );
	}
	if (!dimsOK && !fscnresolu(&ncols, &nrows, fp))
		error(USER, "unspecified matrix size");
	switch (dtype) {
	case DTascii:
	case DTfloat:
	case DTdouble:
		break;
	case DTrgbe:
	case DTxyze:
		cm = cm_load_rgbe(fp, nrows, ncols);
		goto cleanup;
	default:
		error(USER, "unexpected data type in cm_load()");
	}
	if (nrows <= 0) {			/* don't know length? */
		int	guessrows = 147;	/* usually big enough */
		if (cm_elem_size[dtype] && (fp != stdin) & (inspec[0] != '!')) {
			long	startpos = ftell(fp);
			if (fseek(fp, 0L, SEEK_END) == 0) {
				long	rowsiz = (long)ncols*cm_elem_size[dtype];
				long	endpos = ftell(fp);

				if ((endpos - startpos) % rowsiz) {
					sprintf(errmsg,
					"improper length for binary file '%s'",
							inspec);
					error(USER, errmsg);
				}
				guessrows = (endpos - startpos)/rowsiz;
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
	if (!cm)					/* XXX never happens */
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
			size_t	nread = 0;
			do {				/* read all we can */
				nread += getbinary(cm->cmem + 3*nread,
						sizeof(COLOR),
						(size_t)cm->nrows*cm->ncols - nread,
						fp);
				if (nrows <= 0) {	/* unknown length */
					if (nread == (size_t)cm->nrows*cm->ncols)
							/* need more space? */
						cm = cm_resize(cm, cm->nrows+ROWINC);
					else if (nread && !(nread % cm->ncols))
							/* seem to be  done */
						cm = cm_resize(cm, nread/cm->ncols);
					else		/* ended mid-row */
						goto EOFerror;
				} else if (nread < (size_t)cm->nrows*cm->ncols)
					goto EOFerror;
			} while (nread < (size_t)cm->nrows*cm->ncols);

			if (swap) {
				if (sizeof(COLORV) == 4)
					swap32((char *)cm->cmem,
							3*(size_t)cm->nrows*cm->ncols);
				else /* sizeof(COLORV) == 8 */
					swap64((char *)cm->cmem,
							3*(size_t)cm->nrows*cm->ncols);
			}
		} else if (dtype == DTdouble) {
			double	dc[3];			/* load from double */
			COLORV	*cvp = cm->cmem;
			size_t	n = (size_t)nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (getbinary(dc, sizeof(double), 3, fp) != 3)
					goto EOFerror;
				if (swap) swap64((char *)dc, 3);
				copycolor(cvp, dc);
				cvp += 3;
			}
		} else /* dtype == DTfloat */ {
			float	fc[3];			/* load from float */
			COLORV	*cvp = cm->cmem;
			size_t	n = (size_t)nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (getbinary(fc, sizeof(float), 3, fp) != 3)
					goto EOFerror;
				if (swap) swap32((char *)fc, 3);
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
cleanup:
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
	if ((scale[0] < .99) | (scale[0] > 1.01) |
			(scale[1] < .99) | (scale[1] > 1.01) |
			(scale[2] < .99) | (scale[2] > 1.01)) {
		size_t	n = (size_t)ncols*nrows;
		COLORV	*mp = cm->cmem;
		while (n--) {		/* apply exposure scaling */
			*mp++ *= scale[0];
			*mp++ *= scale[1];
			*mp++ *= scale[2];
		}
	}
	return(cm);
EOFerror:
	sprintf(errmsg, "unexpected EOF reading %s", inspec);
	error(USER, errmsg);
	return(NULL);
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

	if (!cm)
		return(NULL);
	if ((c < 0) | (c >= cm->ncols))
		error(INTERNAL, "column requested outside matrix");
	cvr = cm_alloc(cm->nrows, 1);
	if (!cvr)
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


/* Multiply two matrices (or a matrix and a vector) and allocate the result */
CMATRIX *
cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2)
{
	char	*rowcheck=NULL, *colcheck=NULL;
	CMATRIX	*cmr;
	int	dr, dc, i;

	if (!cm1 | !cm2)
		return(NULL);
	if ((cm1->ncols <= 0) | (cm1->ncols != cm2->nrows))
		error(INTERNAL, "matrix dimension mismatch in cm_multiply()");
	cmr = cm_alloc(cm1->nrows, cm2->ncols);
	if (!cmr)
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
		if (rowcheck && !rowcheck[dr])
			continue;
		if (colcheck && !colcheck[dc])
			continue;
		res[0] = res[1] = res[2] = 0;
		for (i = 0; i < cm1->ncols; i++) {
		    const COLORV	*cp1 = cm_lval(cm1,dr,i);
		    const COLORV	*cp2 = cm_lval(cm2,i,dc);
		    res[0] += (double)cp1[0] * cp2[0];
		    res[1] += (double)cp1[1] * cp2[1];
		    res[2] += (double)cp1[2] * cp2[2];
		}
		copycolor(dp, res);
	    }
	if (rowcheck) free(rowcheck);
	if (colcheck) free(colcheck);
	return(cmr);
}


/* write out matrix to file (precede by resolution string if picture) */
int
cm_write(const CMATRIX *cm, int dtype, FILE *fp)
{
	static const char	tabEOL[2] = {'\t','\n'};
	const COLORV		*mp;
	int			r, c;
	size_t			n, rv;

	if (!cm)
		return(0);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	mp = cm->cmem;
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
			n = (size_t)cm->ncols*cm->nrows;
			while (n > 0) {
				rv = fwrite(mp, sizeof(COLOR), n, fp);
				if (rv <= 0)
					return(0);
				mp += 3*rv;
				n -= rv;
			}
		} else if (dtype == DTdouble) {
			double	dc[3];
			n = (size_t)cm->ncols*cm->nrows;
			while (n--) {
				copycolor(dc, mp);
				if (putbinary(dc, sizeof(double), 3, fp) != 3)
					return(0);
				mp += 3;
			}
		} else /* dtype == DTfloat */ {
			float	fc[3];
			n = (size_t)cm->ncols*cm->nrows;
			while (n--) {
				copycolor(fc, mp);
				if (putbinary(fc, sizeof(float), 3, fp) != 3)
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
#ifdef getc_unlocked
	funlockfile(fp);
#endif
	return(fflush(fp) == 0);
}
