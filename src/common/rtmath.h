/* RCSid $Id: rtmath.h,v 3.2 2003/07/12 09:56:24 schorsch Exp $ */
/*
 * Header for Radiance vector and math routines
 */

#ifndef _RAD_RTMATH_H_
#define _RAD_RTMATH_H_
#ifdef __cplusplus
extern "C" {
#endif

#include  <math.h>

#include  "tifftypes.h"

#include  "mat4.h"

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

