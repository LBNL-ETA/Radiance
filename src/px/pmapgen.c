#ifndef lint
static const char	RCSid[] = "$Id: pmapgen.c,v 2.4 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 * pmapgen.c: general routines for 2-D perspective mappings.
 * These routines are independent of the poly structure,
 * so we do not think in terms of texture and screen space.
 *
 * Paul Heckbert	5 Nov 85, 12 Dec 85
 */

#include <stdio.h>
#include "pmap.h"
#include "mx3.h"

#define TOLERANCE 1e-13
#define ZERO(x) ((x)<TOLERANCE && (x)>-TOLERANCE)

#define X(i) qdrl[i][0]		/* quadrilateral x and y */
#define Y(i) qdrl[i][1]

/*
 * pmap_quad_rect: find mapping between quadrilateral and rectangle.
 * The correspondence is:
 *
 *	qdrl[0] --> (u0,v0)
 *	qdrl[1] --> (u1,v0)
 *	qdrl[2] --> (u1,v1)
 *	qdrl[3] --> (u0,v1)
 *
 * This method of computing the adjoint numerically is cheaper than
 * computing it symbolically.
 */
	
extern int
pmap_quad_rect(
	double u0,		/* bounds of rectangle */
	double v0,
	double u1,
	double v1,
	double qdrl[4][2],		/* vertices of quadrilateral */
	double QR[3][3]		/* qdrl->rect transform (returned) */
)
{
    int ret;
    double du, dv;
    double RQ[3][3];		/* rect->qdrl transform */

    du = u1-u0;
    dv = v1-v0;
    if (du==0. || dv==0.) {
	fprintf(stderr, "pmap_quad_rect: null rectangle\n");
	return PMAP_BAD;
    }

    /* first find mapping from unit uv square to xy quadrilateral */
    ret = pmap_square_quad(qdrl, RQ);
    if (ret==PMAP_BAD) return PMAP_BAD;

    /* concatenate transform from uv rectangle (u0,v0,u1,v1) to unit square */
    RQ[0][0] /= du;
    RQ[1][0] /= dv;
    RQ[2][0] -= RQ[0][0]*u0 + RQ[1][0]*v0;
    RQ[0][1] /= du;
    RQ[1][1] /= dv;
    RQ[2][1] -= RQ[0][1]*u0 + RQ[1][1]*v0;
    RQ[0][2] /= du;
    RQ[1][2] /= dv;
    RQ[2][2] -= RQ[0][2]*u0 + RQ[1][2]*v0;

    /* now RQ is transform from uv rectangle to xy quadrilateral */
    /* QR = inverse transform, which maps xy to uv */
    if (mx3d_adjoint(RQ, QR)==0.)
	fprintf(stderr, "pmap_quad_rect: warning: determinant=0\n");
    return ret;
}

/*
 * pmap_square_quad: find mapping between unit square and quadrilateral.
 * The correspondence is:
 *
 *	(0,0) --> qdrl[0]
 *	(1,0) --> qdrl[1]
 *	(1,1) --> qdrl[2]
 *	(0,1) --> qdrl[3]
 */

extern int
pmap_square_quad(
	register double qdrl[4][2],	/* vertices of quadrilateral */
	register double SQ[3][3]	/* square->qdrl transform */
)
{
    double px, py;

    px = X(0)-X(1)+X(2)-X(3);
    py = Y(0)-Y(1)+Y(2)-Y(3);

    if (ZERO(px) && ZERO(py)) {		/* affine */
	SQ[0][0] = X(1)-X(0);
	SQ[1][0] = X(2)-X(1);
	SQ[2][0] = X(0);
	SQ[0][1] = Y(1)-Y(0);
	SQ[1][1] = Y(2)-Y(1);
	SQ[2][1] = Y(0);
	SQ[0][2] = 0.;
	SQ[1][2] = 0.;
	SQ[2][2] = 1.;
	return PMAP_LINEAR;
    }
    else {				/* perspective */
	double dx1, dx2, dy1, dy2, del;

	dx1 = X(1)-X(2);
	dx2 = X(3)-X(2);
	dy1 = Y(1)-Y(2);
	dy2 = Y(3)-Y(2);
	del = DET2(dx1,dx2, dy1,dy2);
	if (del==0.) {
	    fprintf(stderr, "pmap_square_quad: bad mapping\n");
	    return PMAP_BAD;
	}
	SQ[0][2] = DET2(px,dx2, py,dy2)/del;
	SQ[1][2] = DET2(dx1,px, dy1,py)/del;
	SQ[2][2] = 1.;
	SQ[0][0] = X(1)-X(0)+SQ[0][2]*X(1);
	SQ[1][0] = X(3)-X(0)+SQ[1][2]*X(3);
	SQ[2][0] = X(0);
	SQ[0][1] = Y(1)-Y(0)+SQ[0][2]*Y(1);
	SQ[1][1] = Y(3)-Y(0)+SQ[1][2]*Y(3);
	SQ[2][1] = Y(0);
	return PMAP_PERSP;
    }
}
