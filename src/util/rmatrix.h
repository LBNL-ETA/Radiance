/* RCSid $Id: rmatrix.h,v 2.1 2014/05/31 05:02:37 greg Exp $ */
/*
 * Header file for general matrix routines.
 */

#ifndef _RAD_RMATRIX_H_
#define _RAD_RMATRIX_H_

#include "cmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* General plane-ordered component matrix */
typedef struct {
	int	nrows, ncols, ncomp;
	double	mtx[1];			/* extends struct */
} RMATRIX;

#define rmx_indx(rm,r,c,i)	(((i)*(rm)->nrows+(r))*(rm)->ncols+(c))

#define rmx_lval(rm,r,c,i)	(rm)->mtx[rmx_indx(rm,r,c,i)]

#define rmx_free(rm)		free(rm)

/* Allocate a nr x nc matrix with n components */
extern RMATRIX	*rmx_alloc(int nr, int nc, int n);

/* Load matrix from supported file type */
extern RMATRIX	*rmx_load(const char *fname);

/* Write matrix to file type indicated by dt */
extern long	rmx_write(const RMATRIX *rm, int dtype, FILE *fp);

/* Allocate and assign square identity matrix with n components */
extern RMATRIX	*rmx_identity(int dim, int n);

/* Duplicate the given matrix */
extern RMATRIX	*rmx_copy(const RMATRIX *rm);

/* Swap rows and columns in the given matrix in situ */
extern int	rmx_transpose(RMATRIX *rm);

/* Multiply (concatenate) two matrices and allocate the result */
extern RMATRIX	*rmx_multiply(const RMATRIX *m1, const RMATRIX *m2);

/* Sum second matrix into first, applying scale factor beforehand */
extern int	rmx_sum(RMATRIX *msum, const RMATRIX *madd, const double sf[]);

/* Scale the given matrix by the indicated scalar component vector */
extern int	rmx_scale(RMATRIX *rm, const double sf[]);

/* Allocate new matrix and apply component transformation */
extern RMATRIX	*rmx_transform(const RMATRIX *msrc, int n, const double cmat[]);

/* Convert a color matrix to newly allocated RMATRIX buffer */
extern RMATRIX	*rmx_from_cmatrix(const CMATRIX *cm);

/* Convert general matrix to newly allocated CMATRIX buffer */
extern CMATRIX	*cm_from_rmatrix(const RMATRIX *rm);

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_RMATRIX_H_ */
