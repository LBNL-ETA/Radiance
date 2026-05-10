#ifndef lint
static const char RCSid[] = "$Id$";
#endif

/*
* matrix multiplication (avx2).
* 
* Yongqing
*/

#include "cmatrix.h"
#include "platform.h"
#include "standard.h"

#include <ctype.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#include "cmmap.h"
#include "cmatrixalign.h"
#include "mmavx2.h"

CMATRIX* cm_multiply_avx2(const CMATRIX* cm1, const CMATRIX* cm2) {
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

			// =========================================================================
			// Intel/AMD x86_64 256 AVX2
#if defined(__AVX2__)
			__m256 v0_1 = _mm256_setr_ps(ra0, ga0, ba0, ra0, ga0, ba0, ra0, ga0);
			__m256 v0_2 = _mm256_setr_ps(ba0, ra0, ga0, ba0, ra0, ga0, ba0, ra0);
			__m256 v0_3 = _mm256_setr_ps(ga0, ba0, ra0, ga0, ba0, ra0, ga0, ba0);

			__m256 v1_1 = _mm256_setr_ps(ra1, ga1, ba1, ra1, ga1, ba1, ra1, ga1);
			__m256 v1_2 = _mm256_setr_ps(ba1, ra1, ga1, ba1, ra1, ga1, ba1, ra1);
			__m256 v1_3 = _mm256_setr_ps(ga1, ba1, ra1, ga1, ba1, ra1, ga1, ba1);

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
#endif 
			// =========================================================================


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

			// =========================================================================
#if defined(__AVX2__)
			__m256 va1 = _mm256_setr_ps(ra, ga, ba, ra, ga, ba, ra, ga);
			__m256 va2 = _mm256_setr_ps(ba, ra, ga, ba, ra, ga, ba, ra);
			__m256 va3 = _mm256_setr_ps(ga, ba, ra, ga, ba, ra, ga, ba);

			for (; j <= limit - 24; j += 24) {

				__m256 vb1 = _mm256_loadu_ps(&row_b[j]);
				__m256 vb2 = _mm256_loadu_ps(&row_b[j + 8]);
				__m256 vb3 = _mm256_loadu_ps(&row_b[j + 16]);

				_mm256_storeu_ps(&r_c[j], _mm256_fmadd_ps(vb1, va1, _mm256_loadu_ps(&r_c[j])));
				_mm256_storeu_ps(&r_c[j + 8], _mm256_fmadd_ps(vb2, va2, _mm256_loadu_ps(&r_c[j + 8])));
				_mm256_storeu_ps(&r_c[j + 16], _mm256_fmadd_ps(vb3, va3, _mm256_loadu_ps(&r_c[j + 16])));
			}
#endif
			// =========================================================================

			for (; j < limit; j += 3) {
				r_c[j + 0] += ra * row_b[j + 0];
				r_c[j + 1] += ga * row_b[j + 1];
				r_c[j + 2] += ba * row_b[j + 2];
			}
		}
	}

	return cmr;
}
