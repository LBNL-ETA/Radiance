/* RCSid $Id: rtmath.h,v 3.5 2013/02/07 20:19:52 greg Exp $ */
/*
 * Header for Radiance vector and math routines
 */

#ifndef _RAD_RTMATH_H_
#define _RAD_RTMATH_H_

#include  <math.h>

#include  "tiff.h"		/* needed for int32, etc. */
#include  "mat4.h"

#ifdef __cplusplus
extern "C" {
#endif

				/* regular transformation */
typedef struct {
	MAT4  xfm;				/* transform matrix */
	RREAL  sca;				/* scalefactor */
}  XF;
				/* complemetary tranformation */
typedef struct {
	XF  f;					/* forward */
	XF  b;					/* backward */
}  FULLXF;

#ifndef  PI
#ifdef	M_PI
#define	 PI		((double)M_PI)
#else
#define	 PI		3.14159265358979323846
#endif
#endif

#ifdef  FASTMATH
#define  tcos			cos
#define  tsin			sin
#define  ttan			tan
#else
					/* table-based cosine approximation */
#define  tsin(x)		tcos((x)-(PI/2.))
#define  ttan(x)		(tsin(x)/tcos(x))
#endif
					/* defined in tcos.c */
extern double	tcos(double x);
extern double	aatan2(double y, double x);
					/* defined in xf.c */
extern int	xf(XF *ret, int ac, char *av[]);
extern int	invxf(XF *ret, int ac, char *av[]);
extern int	fullxf(FULLXF *fx, int ac, char *av[]);
					/* defined in zeroes.c */
extern int	quadratic(double *r, double a, double b, double c);
					/* defined in dircode.c */
extern int32	encodedir(FVECT dv);
extern void	decodedir(FVECT dv, int32 dc);
extern double	dir2diff(int32 dc1, int32 dc2);
extern double	fdir2diff(int32 dc1, FVECT v2);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMATH_H_ */

