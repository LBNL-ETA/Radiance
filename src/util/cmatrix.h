/* RCSid $Id: cmatrix.h,v 2.5 2014/05/30 00:00:54 greg Exp $ */
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
enum {DTfromHeader=0, DTascii, DTfloat, DTdouble, DTrgbe, DTxyze, DTend};

extern const char	*cm_fmt_id[];
extern const int	cm_elem_size[];

/* A color coefficient matrix -- vectors have ncols==1 */
typedef struct {
	int	nrows, ncols;
	COLORV	cmem[3];		/* extends struct */
} CMATRIX;

#define COLSPEC	(sizeof(COLORV)==sizeof(float) ? "%f %f %f" : "%lf %lf %lf")

#define cm_lval(cm,r,c)	((cm)->cmem + 3*((r)*(cm)->ncols + (c)))

#define cv_lval(cm,i)	((cm)->cmem + 3*(i))

/* Allocate a color coefficient matrix */
extern CMATRIX	*cm_alloc(int nrows, int ncols);

/* Resize color coefficient matrix */
extern CMATRIX	*cm_resize(CMATRIX *cm, int nrows);

#define cm_free(cm)	free(cm)

/* Load header to obtain/check data type and matrix dimensions */
extern char	*cm_getheader(int *dt, int *nr, int *nc, FILE *fp);

/* Allocate and load a matrix from the given file (or stdin if NULL) */
extern CMATRIX	*cm_load(const char *fname, int nrows, int ncols, int dtype);

/* Extract a column vector from a matrix */
extern CMATRIX	*cm_column(const CMATRIX *cm, int c);

/* Scale a matrix by a single value */
extern CMATRIX	*cm_scale(const CMATRIX *cm1, const COLOR sca);

/* Multiply two matrices (or a matrix and a vector) and allocate the result */
extern CMATRIX	*cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2);

/* write out matrix to file (precede by resolution string if picture) */
extern int	cm_write(const CMATRIX *cm, int dtype, FILE *fp);

/* Load and convert a matrix BTDF from the given XML file */
extern CMATRIX	*cm_loadBTDF(char *fname);

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_CMATRIX_H_ */
