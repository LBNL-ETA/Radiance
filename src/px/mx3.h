/* RCSid: $Id: mx3.h,v 2.5 2003/07/14 22:24:00 schorsch Exp $ */
#ifndef _RAD_MX3_H_
#define _RAD_MX3_H_

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define mx3_copy(a, b) memcpy(b, a, 3*3*sizeof(double))

/* FLOAT VERSIONS */
float  mx3f_adjoint(/* a, b */);
void  mx3f_mul(/* a, b, c */);
void  mx3f_print(/* str, a */);
void  mx3f_translate_mat(/* tx, ty, a */);
void  mx3f_translate(/* a, tx, ty */);
void  mx3f_scale_mat(/* sx, sy, a */);
void  mx3f_scale(/* a, sx, sy */);
void  mx3f_transform(/* p, a, q */);
float  mx3f_transform_div(/* p, a, q */);
#define mx3f_copy(a, b) memcpy(b, a, 3*3*sizeof(float))

/* DOUBLE VERSIONS */
double  mx3d_adjoint(/* a, b */);
void  mx3d_mul(/* a, b, c */);
void  mx3d_print(/* str, a */);
void  mx3d_translate_mat(/* tx, ty, a */);
void  mx3d_translate(/* a, tx, ty */);
void  mx3d_scale_mat(/* sx, sy, a */);
void  mx3d_scale(/* a, sx, sy */);
void  mx3d_transform(/* p, a, q */);
double  mx3d_transform_div(/* p, a, q */);
#define mx3d_copy(a, b) memcpy(b, a, 3*3*sizeof(double))

#ifdef __cplusplus
}
#endif
#endif /* _RAD_MX3_H_ */
