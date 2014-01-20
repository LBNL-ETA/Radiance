/* RCSid $Id: cmatrix.h,v 2.2 2014/01/20 21:30:34 greg Exp $ */
/*
 * Color matrix routine declarations.
 *
 *	G. Ward
 */

#ifndef _RAD_CMATRIX_H_
#define _RAD_CMATRIX_H_

#include "color.h"

#ifdef __cplusplus
extern "C" {
#endif

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
extern CMATRIX *cm_alloc(int nrows, int ncols);

/* Resize color coefficient matrix */
extern CMATRIX *cm_resize(CMATRIX *cm, int nrows);

#define cm_free(cm)	free(cm)

/* Load header to obtain data type */
extern int getDTfromHeader(FILE *fp);

/* Allocate and load a matrix from the given file (or stdin if NULL) */
extern CMATRIX *cm_load(const char *fname, int nrows, int ncols, int dtype);

/* Extract a column vector from a matrix */
extern CMATRIX *cm_column(const CMATRIX *cm, int c);

/* Scale a matrix by a single value */
extern CMATRIX *cm_scale(const CMATRIX *cm1, const COLOR sca);

/* Multiply two matrices (or a matrix and a vector) and allocate the result */
extern CMATRIX *cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2);

/* print out matrix as ASCII text -- no header */
extern void cm_print(const CMATRIX *cm, FILE *fp);

/* Load and convert a matrix BSDF from the given XML file */
extern CMATRIX *cm_loadBSDF(char *fname, COLOR cLamb);

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_CMATRIX_H_ */
