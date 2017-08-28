/* RCSid $Id: rmatrix.h,v 2.8 2017/08/28 15:59:46 greg Exp $ */
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
	int	dtype;
	char	*info;
	double	mtx[1];			/* extends struct */
} RMATRIX;

#define rmx_lval(rm,r,c,i)	(rm)->mtx[(i)+(rm)->ncomp*((c)+(rm)->ncols*(r))]

/* Allocate a nr x nc matrix with n components */
extern RMATRIX	*rmx_alloc(int nr, int nc, int n);

/* Free a RMATRIX array */
extern void	rmx_free(RMATRIX *rm);

/* Resolve data type based on two input types (returns 0 for mismatch) */
extern int	rmx_newtype(int dtyp1, int dtyp2);

/* Load matrix from supported file type (NULL for stdin, '!' with command) */
extern RMATRIX	*rmx_load(const char *inspec);

/* Append header information associated with matrix data */
extern int	rmx_addinfo(RMATRIX *rm, const char *info);

/* Write matrix to file type indicated by dtype */
extern int	rmx_write(const RMATRIX *rm, int dtype, FILE *fp);

/* Allocate and assign square identity matrix with n components */
extern RMATRIX	*rmx_identity(int dim, int n);

/* Duplicate the given matrix */
extern RMATRIX	*rmx_copy(const RMATRIX *rm);

/* Allocate and assign transposed matrix */
extern RMATRIX	*rmx_transpose(const RMATRIX *rm);

/* Multiply (concatenate) two matrices and allocate the result */
extern RMATRIX	*rmx_multiply(const RMATRIX *m1, const RMATRIX *m2);

/* Element-wise multiplication (or division) of m2 into m1 */
extern int	rmx_elemult(RMATRIX *m1, const RMATRIX *m2, int divide);

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
