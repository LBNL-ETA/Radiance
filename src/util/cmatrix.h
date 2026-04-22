/* RCSid $Id$ */
/*
 * Color matrix routine declarations.
 *
 *	G. Ward
 */

#ifndef _RAD_CMATRIX_H_
#define _RAD_CMATRIX_H_

#include  <sys/types.h>
#include "color.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Data types for file loading (used to be an enum) */
#define DTfromHeader	0
#define DTrgbe		1
#define DTxyze		2
#define DTspec		3
#define DTfloat		4
#define DTascii		5
#define	DTdouble	6
#define DTend		7

#ifdef _WIN32
#define cm_free_u(cm)	_aligned_free(cm)
#else
#define cm_free_u(cm)	free(cm)
#endif

/* Defined in cmconst.c */
extern const char	stdin_name[];
extern const char	*cm_fmt_id[];
extern const int	cm_elem_size[];

/* A color coefficient matrix -- vectors have ncols==1 */
typedef struct {
	int	nrows, ncols;
	COLORV	cmem[3];		/* extends struct */
} CMATRIX;

#define COLSPEC	(sizeof(COLORV)==sizeof(float) ? "%f %f %f" : "%lf %lf %lf")

#define cm_lval(cm,r,c)	((cm)->cmem + 3*((size_t)(r)*(cm)->ncols + (c)))

#define cv_lval(cm,i)	((cm)->cmem + 3*(i))

/* Allocate a color coefficient matrix */
extern CMATRIX	*cm_alloc(int nrows, int ncols);
#ifdef INTEL_DCTOPT
extern CMATRIX* cm_alloc_u(int rows, int cols);
extern CMATRIX* cm_resize_aligned(CMATRIX* cm, int nrows);
extern CMATRIX* cm_load_u(const char* inspec, int nrows, int ncols, int dtype);
extern CMATRIX* cm_column_u(const CMATRIX* cm, int c);
extern CMATRIX* cm_get_batch(const CMATRIX* src, int start_col, int num_cols);
extern CMATRIX* cm_multiply_avx2(const CMATRIX* cm1, const CMATRIX* cm2);
extern CMATRIX* cm_multiply_smart(const CMATRIX* basis, const CMATRIX* cv_batch);
#endif /* INTEL_DCTOPT */

/* Resize color coefficient matrix */
extern CMATRIX	*cm_resize(CMATRIX *cm, int nrows);

#define cm_free(cm)	free(cm)

/* Load header to obtain/check data type and matrix dimensions */
extern char	*cm_getheader(int *dt, int *nr, int *nc,
				int *swp, COLOR scale, FILE *fp);

/* Allocate and load a matrix from the given input (or stdin if NULL) */
extern CMATRIX	*cm_load(const char *inspec, int nrows, int ncols, int dtype);

/* Extract a column vector from a matrix */
extern CMATRIX	*cm_column(const CMATRIX *cm, int c);

/* Multiply two matrices (or a matrix and a vector) and allocate the result */
extern CMATRIX	*cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2);

/* write out matrix to file (precede by resolution string if picture) */
extern int	cm_write(const CMATRIX *cm, int dtype, FILE *fp);

/* Load and convert a matrix BTDF from the given XML file */
extern CMATRIX	*cm_loadBTDF(const char *fname);

/* Load and convert a matrix BRDF from the given XML file */
extern CMATRIX	*cm_loadBRDF(const char *fname, int backside);

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_CMATRIX_H_ */
