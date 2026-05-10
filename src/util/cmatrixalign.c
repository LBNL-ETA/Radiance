#ifndef lint
static const char RCSid[] = "$Id$";
#endif

/*
* matrix routines (avx2).
* 
* Yongqing
*/


#include "cmatrix.h"
#include "cmmap.h"
#include "platform.h"
#include "standard.h"
#include "resolu.h"
#include "mmavx2.h"
#include <ctype.h>
#include <stdbool.h>
#include "cmatrixalign.h"
#ifdef _WIN32
    #include <intrin.h>
	static inline void get_cpuid(int leaf, int subleaf, unsigned int* eax, unsigned int* ebx, unsigned int* ecx, unsigned int* edx) {
        int cpuInfo[4];
        __cpuidex(cpuInfo, leaf, subleaf);
        *eax = cpuInfo[0];
        *ebx = cpuInfo[1];
        *ecx = cpuInfo[2];
        *edx = cpuInfo[3];
    }
    
    static inline unsigned long long get_xcr0(void) {
        return _xgetbv(0);
    }


#else
    #include <cpuid.h>
    
    static inline void get_cpuid(int leaf, int subleaf, unsigned int* eax, unsigned int* ebx, unsigned int* ecx, unsigned int* edx) {
        __cpuid_count(leaf, subleaf, *eax, *ebx, *ecx, *edx);
    }
    
    static inline unsigned long long get_xcr0(void) {
        unsigned int eax, edx;
        __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
        return ((unsigned long long)edx << 32) | eax;
    }

#endif

int cpu_supports_avx2_fma(void) {
    unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;

    // get Leaf 1
    get_cpuid(1, 0, &eax, &ebx, &ecx, &edx);
    
    // check FMA (ECX bit 12)
    if ((ecx & (1 << 12)) == 0) return 0;
    
    // check OSXSAVE (ECX bit 27) - XGETBV
    if ((ecx & (1 << 27)) == 0) return 0;


    // getLeaf 7
    get_cpuid(7, 0, &eax, &ebx, &ecx, &edx);
    
    // check AVX2 (EBX bit 5)
    if ((ebx & (1 << 5)) == 0) return 0;


    // OS (XCR0)
    // Make sure that the operating system has enabled the context save mechanism for SSE (bit 1) and AVX (bit 2) (1 | 2 = 3, that is, decimal 6).
    unsigned long long xcr0 = get_xcr0();
    if ((xcr0 & 6) != 6) return 0;

    return 1;
}

CMATRIX* cm_alloc_u(int rows, int cols) {
	size_t total_elements = (size_t)rows * cols;
	size_t data_size = sizeof(float) * 3 * total_elements;

	size_t total_size = sizeof(CMATRIX) + data_size + 32;

#ifdef  _WIN32
	CMATRIX* m = (CMATRIX*)_aligned_malloc((total_size + 31) & ~31, 32);
#else
	CMATRIX* m;
	if (posix_memalign((void **)&m, 32, (total_size + 31) & ~31) < 0)
		m = NULL;
#endif //  _WIN32
	if (m) {
		m->nrows = rows;
		m->ncols = cols;

		memset(m->cmem, 0, data_size + 32);
	}
	return m;
}

CMATRIX*
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

CMATRIX*
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


CMATRIX* cm_multiply_sd(const CMATRIX* cm1, const CMATRIX* cm2) {
	int R = cm1->nrows, K = cm1->ncols, C = cm2->ncols;
	CMATRIX* cmr = cm_alloc_u(R, C);
	if (!cmr) return NULL;
	int i = 0;
	//#pragma omp parallel for schedule(dynamic, 8)
	for (i = 0; i < R - 1; i += 2) {
		float* r_a0 = (float*)cm_lval(cm1, i, 0);
		float* r_a1 = (float*)cm_lval(cm1, i + 1, 0);
		float* r_c0 = (float*)cm_lval(cmr, i, 0);
		float* r_c1 = (float*)cm_lval(cmr, i + 1, 0);

		for (int k = 0; k < K; k++) {
			float ra0 = r_a0[k * 3 + 0], ga0 = r_a0[k * 3 + 1], ba0 = r_a0[k * 3 + 2];
			float ra1 = r_a1[k * 3 + 0], ga1 = r_a1[k * 3 + 1], ba1 = r_a1[k * 3 + 2];

			if (ra0 == 0 && ga0 == 0 && ba0 == 0 && ra1 == 0 && ga1 == 0 && ba1 == 0) continue;

			float* row_b = (float*)cm_lval(cm2, k, 0);
			int j = 0, limit = C * 3;
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
			float* row_b = (float*)cm_lval(cm2, k, 0);
			int j = 0, limit = C * 3;
			for (; j < limit; j += 3) {
				r_c[j + 0] += ra * row_b[j + 0];
				r_c[j + 1] += ga * row_b[j + 1];
				r_c[j + 2] += ba * row_b[j + 2];
			}
		}
	}

	return cmr;
}

int*
stat_numn_onzerocol(CMATRIX* cv, int* size) {
	int K = cv->nrows;
	int C = cv->ncols;

	int* pixel_count_map = (int*)calloc(C, sizeof(int));
	int num_nonzero = 0;

	for (int c = 0; c < C; c++) {
		int col_pixel_count = 0;
		for (int k = 0; k < K; k++) {
			float* val = (float*)cm_lval(cv, k, c);
			if (val[0] != 0.0f || val[1] != 0.0f || val[2] != 0.0f) {
				col_pixel_count++;
			}
		}

		if (col_pixel_count > 0) {
			pixel_count_map[num_nonzero] = col_pixel_count;
			num_nonzero++;
		}
	}

	*size = num_nonzero;
	return pixel_count_map;
}

CMATRIX* cm_multiply_smart(const CMATRIX* basis, const CMATRIX* cv_batch) {
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
		CMATRIX* res = NULL;
		if(cpu_supports_avx2_fma()) {
			res = cm_multiply_avx2(basis, cv_batch);
		}
		else {
			res = cm_multiply_sd(basis, cv_batch);
		}
		return res;
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
	CMATRIX* compact_res = NULL;
	if(cpu_supports_avx2_fma()) {
		compact_res = cm_multiply_avx2(basis, compact_cv);
	}
	else {
		compact_res = cm_multiply_sd(basis, compact_cv);
	}
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

CMATRIX*
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




