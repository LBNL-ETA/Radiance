
/* RCSid $Id$ */
/*
* matrix routines declarations (avx2).
*
* Yongqing
*/

#include "cmatrix.h"

#ifdef _WIN32
#define cm_free_u(cm)	_aligned_free(cm)
#else
#define cm_free_u(cm)	free(cm)
#endif

extern int MAX_CHUNK_SIZE;
extern CMATRIX* cm_alloc_u(int rows, int cols);
extern CMATRIX* cm_column_u(const CMATRIX* cm, int c);

extern CMATRIX* cm_get_batch(const CMATRIX* src, int start_col, int num_cols);

static CMATRIX* cm_load_rgbe_u(FILE* fp, int nrows, int ncols);

extern CMATRIX* cm_load_u(const char* inspec, int nrows, int ncols, int dtype);
extern CMATRIX* cm_multiply_avx2(const CMATRIX* cm1, const CMATRIX* cm2);
extern CMATRIX* cm_multiply_smart(const CMATRIX* basis, const CMATRIX* cv_batch);

CMATRIX* cm_resize_aligned(CMATRIX* cm, int nrows);

extern int* stat_numn_onzerocol(CMATRIX* cv, int* size);
