#ifndef lint
static const char	RCSid[] = "$Id: mx3.c,v 2.3 2004/03/28 20:33:14 schorsch Exp $";
#endif

#include "mx3.h"

extern double
mx3d_adjoint(
	register double a[3][3],
	register double b[3][3]
)
{
    b[0][0] =  ((a[1][1])*( a[2][2]) - ( a[1][2])*( a[2][1]));
    b[1][0] =  ((a[1][2])*( a[2][0]) - ( a[1][0])*( a[2][2]));
    b[2][0] =  ((a[1][0])*( a[2][1]) - ( a[1][1])*( a[2][0]));
    b[0][1] =  ((a[2][1])*( a[0][2]) - ( a[2][2])*( a[0][1]));
    b[1][1] =  ((a[2][2])*( a[0][0]) - ( a[2][0])*( a[0][2]));
    b[2][1] =  ((a[2][0])*( a[0][1]) - ( a[2][1])*( a[0][0]));
    b[0][2] =  ((a[0][1])*( a[1][2]) - ( a[0][2])*( a[1][1]));
    b[1][2] =  ((a[0][2])*( a[1][0]) - ( a[0][0])*( a[1][2]));
    b[2][2] =  ((a[0][0])*( a[1][1]) - ( a[0][1])*( a[1][0]));
    return a[0][0]*b[0][0] + a[1][0]*b[0][1] + a[2][0]*b[0][2];
}


extern void
mx3d_mul(
	register double a[3][3],
	register double b[3][3],
	register double c[3][3]
)
{
    c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
    c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
    c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
    c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
    c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
    c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
    c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
    c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
    c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
}


extern void
mx3d_transform(
	register double p[3],
	register double a[3][3],
	register double q[3]
)
{
    q[0] = p[0]*a[0][0] + p[1]*a[1][0] + p[2]*a[2][0];
    q[1] = p[0]*a[0][1] + p[1]*a[1][1] + p[2]*a[2][1];
    q[2] = p[0]*a[0][2] + p[1]*a[1][2] + p[2]*a[2][2];
}


extern double
mx3d_transform_div(
	register double p[3],
	register double a[3][3],
	register double q[2]
)
{
    double q2;
    q2 = p[0]*a[0][2] + p[1]*a[1][2] + p[2]*a[2][2];
    if (q2!=0.) {
	q[0] = (p[0]*a[0][0] + p[1]*a[1][0] + p[2]*a[2][0]) / q2;
	q[1] = (p[0]*a[0][1] + p[1]*a[1][1] + p[2]*a[2][1]) / q2;
    }
    return q2;
}


extern void
mx3d_translate_mat(
	double tx,
	double ty,
	register double m[3][3]
)
{
    m[0][0] = 1.; m[0][1] = 0.; m[0][2] = 0.;
    m[1][0] = 0.; m[1][1] = 1.; m[1][2] = 0.;
    m[2][0] = tx; m[2][1] = ty; m[2][2] = 1.;
}


extern void
mx3d_translate(
	register double m[3][3],
	double tx,
	double ty
)
{
    m[3][0] += tx*m[0][0]+ty*m[1][0];
    m[3][1] += tx*m[0][1]+ty*m[1][1];
    m[3][2] += tx*m[0][2]+ty*m[1][2];
}


extern void
mx3d_scale_mat(
	double sx,
	double sy,
	register double m[3][3]
)
{
    m[0][0] = sx; m[0][1] = 0.; m[0][2] = 0.;
    m[1][0] = 0.; m[1][1] = sy; m[1][2] = 0.;
    m[2][0] = 0.; m[2][1] = 0.; m[2][2] = 1.;
}


extern void
mx3d_scale(
	register double m[3][3],
	double sx,
	double sy
)
{
    m[0][0] *= sx; m[0][1] *= sx; m[0][2] *= sx;
    m[1][0] *= sy; m[1][1] *= sy; m[1][2] *= sy;
}
