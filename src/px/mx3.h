/* RCSid: $Id: mx3.h,v 2.6 2004/03/28 20:33:14 schorsch Exp $ */
#ifndef _RAD_MX3_H_
#define _RAD_MX3_H_

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define mx3_copy(a, b) memcpy(b, a, 3*3*sizeof(double))

/* FLOAT VERSIONS */
extern float  mx3f_adjoint(/* a, b */);
extern void  mx3f_mul(/* a, b, c */);
extern void  mx3f_print(/* str, a */);
extern void  mx3f_translate_mat(/* tx, ty, a */);
extern void  mx3f_translate(/* a, tx, ty */);
extern void  mx3f_scale_mat(/* sx, sy, a */);
extern void  mx3f_scale(/* a, sx, sy */);
extern void  mx3f_transform(/* p, a, q */);
extern float  mx3f_transform_div(/* p, a, q */);
#define mx3f_copy(a, b) memcpy(b, a, 3*3*sizeof(float))

/* DOUBLE VERSIONS */
extern double  mx3d_adjoint(double a[3][3], double b[3][3]);
extern void  mx3d_mul(double a[3][3], double b[3][3], double c[3][3]);
extern void  mx3d_print(/* str, a */);  /* XXX doesn't exist */
extern void  mx3d_translate_mat(double tx, double ty, double m[3][3]);
extern void  mx3d_translate(double m[3][3], double tx, double ty);
extern void  mx3d_scale_mat(double sx, double sy, double m[3][3]);
extern void  mx3d_scale(double m[3][3], double sx, double sy);
extern void  mx3d_transform(double p[3], double a[3][3], double q[3]);
extern double  mx3d_transform_div(double p[3], double a[3][3], double q[2]);
#define mx3d_copy(a, b) memcpy(b, a, 3*3*sizeof(double))

#ifdef __cplusplus
}
#endif
#endif /* _RAD_MX3_H_ */
