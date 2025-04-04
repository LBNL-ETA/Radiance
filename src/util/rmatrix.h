/* RCSid $Id: rmatrix.h,v 2.26 2025/04/04 18:06:48 greg Exp $ */
/*
 * Header file for general matrix routines.
 */

#ifndef _RAD_RMATRIX_H_
#define _RAD_RMATRIX_H_

#include "cmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Preferred BSDF component:
	none, transmission, reflection front (normal side), reflection back */
typedef enum {RMPnone=-1, RMPtrans=0, RMPreflF, RMPreflB} RMPref;

/* RMATRIX flags (usually private):
	need to swap input, we should free memory */
#define	RMF_SWAPIN	1
#define RMF_FREEMEM	2

#ifndef MAXCOMP
#define MAXCOMP		MAXCSAMP	/* #components we support */
#endif
					/* Set in-core data type */
#if !defined(DTrmx_native) || DTrmx_native==DTfloat
#define	DTrmx_native	DTfloat
#define rmx_dtype	float
#define rmx_scanfmt	"%f"
#elif DTrmx_native==DTdouble
#define	rmx_dtype	double
#define rmx_scanfmt	"%lf"
#endif

/* General [row][col][cmp] component matrix */
typedef struct {
	char		*info;
	void		*mapped;
	rmx_dtype	*mtx;
	COLOR		cexp;
	float		wlpart[4];
	int		nrows, ncols;
	short		ncomp;
	uby8		dtype;
	uby8		pflags;
} RMATRIX;

#define rmx_lval(rm,r,c)	((rm)->mtx + (rm)->ncomp*((c)+(size_t)(rm)->ncols*(r)))
#define rmx_val			rmx_lval

#define rmx_array_size(rm)	(sizeof(rmx_dtype)*(rm)->nrows*(rm)->ncols*(rm)->ncomp)
#define rmx_mapped_size(rm)	((char *)(rm)->mtx + rmx_array_size(rm) - (char *)(rm)->mapped)

/* Initialize a RMATRIX struct but don't allocate array space */
extern RMATRIX	*rmx_new(int nr, int nc, int n);

/* Prepare a RMATRIX for writing (allocate array if needed) */
extern int	rmx_prepare(RMATRIX *rm);

/* Call rmx_new() and rmx_prepare() */
extern RMATRIX	*rmx_alloc(int nr, int nc, int n);

/* Clear state by freeing info and matrix data */
extern void	rmx_reset(RMATRIX *rm);

/* Free an RMATRIX struct and data */
extern void	rmx_free(RMATRIX *rm);

/* Resolve data type based on two input types (returns 0 for mismatch) */
extern int	rmx_newtype(int dtyp1, int dtyp2);

/* Read matrix header from input stream (cannot be XML) */
extern int	rmx_load_header(RMATRIX *rm, FILE *fp);

/* Load next row as rmx_dtype (cannot be XML) */
extern int	rmx_load_row(rmx_dtype *drp, const RMATRIX *rm, FILE *fp);

/* Allocate & load post-header data from stream given type set in rm->dtype */
extern int	rmx_load_data(RMATRIX *rm, FILE *fp);

/* Load matrix from supported file type (NULL for stdin, '!' with command) */
extern RMATRIX	*rmx_load(const char *inspec, RMPref rmp);

/* Append header information associated with matrix data */
extern int	rmx_addinfo(RMATRIX *rm, const char *info);

/* Finish writing header data with resolution and format, returning type used */
extern int	rmx_write_header(const RMATRIX *rm, int dtype, FILE *fp);

/* Write out matrix data (usually by row) */
extern int	rmx_write_data(const rmx_dtype *dp, int nc, int len,
				int dtype, FILE *fp);

/* Write matrix using file format indicated by dtype */
extern int	rmx_write(const RMATRIX *rm, int dtype, FILE *fp);

/* Allocate and assign square identity matrix with n components */
extern RMATRIX	*rmx_identity(int dim, int n);

/* Duplicate the given matrix */
extern RMATRIX	*rmx_copy(const RMATRIX *rm);

/* Replace data in first matrix with data from second */
extern int	rmx_transfer_data(RMATRIX *rdst, RMATRIX *rsrc, int dometa);

/* Transpose the given matrix */
extern int	rmx_transpose(RMATRIX *rm);

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
