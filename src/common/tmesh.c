#ifndef lint
static const char RCSid[] = "$Id: tmesh.c,v 2.4 2003/11/14 17:22:06 schorsch Exp $";
#endif
/*
 * Compute and print barycentric coordinates for triangle meshes
 */

#include <stdio.h>

#include "fvect.h"

#include "tmesh.h"

#define ABS(x)		((x) >= 0 ? (x) : -(x))


int
flat_tri(v1, v2, v3, n1, n2, n3)	/* determine if triangle is flat */
FVECT	v1, v2, v3, n1, n2, n3;
{
	double	d1, d2, d3;
	FVECT	vt1, vt2, vn;
					/* compute default normal */
	VSUB(vt1, v2, v1);
	VSUB(vt2, v3, v2);
	VCROSS(vn, vt1, vt2);
	if (normalize(vn) == 0.0)
		return(DEGEN);
					/* compare to supplied normals */
	d1 = DOT(vn, n1); d2 = DOT(vn, n2); d3 = DOT(vn, n3);
	if (d1 < 0 && d2 < 0 && d3 < 0) {
		if (d1 > -COSTOL || d2 > -COSTOL || d3 > -COSTOL)
			return(RVBENT);
		return(RVFLAT);
	}
	if (d1 < COSTOL || d2 < COSTOL || d3 < COSTOL)
		return(ISBENT);
	return(ISFLAT);
}


int
comp_baryc(bcm, v1, v2, v3)		/* compute barycentric vectors */
register BARYCCM	*bcm;
RREAL			*v1, *v2, *v3;
{
	RREAL	*vt;
	FVECT	va, vab, vcb;
	double	d;
	int	ax0, ax1;
	register int	i;
					/* compute major axis */
	VSUB(vab, v1, v2);
	VSUB(vcb, v3, v2);
	VCROSS(va, vab, vcb);
	bcm->ax = ABS(va[0]) > ABS(va[1]) ? 0 : 1;
	bcm->ax = ABS(va[bcm->ax]) > ABS(va[2]) ? bcm->ax : 2;
	ax0 = (bcm->ax + 1) % 3;
	ax1 = (bcm->ax + 2) % 3;
	for (i = 0; i < 2; i++) {
		vab[0] = v1[ax0] - v2[ax0];
		vcb[0] = v3[ax0] - v2[ax0];
		vab[1] = v1[ax1] - v2[ax1];
		vcb[1] = v3[ax1] - v2[ax1];
		d = vcb[0]*vcb[0] + vcb[1]*vcb[1];
		if (d <= FTINY*FTINY)
			return(-1);
		d = (vcb[0]*vab[0]+vcb[1]*vab[1])/d;
		va[0] = vab[0] - vcb[0]*d;
		va[1] = vab[1] - vcb[1]*d;
		d = va[0]*va[0] + va[1]*va[1];
		if (d <= FTINY*FTINY)
			return(-1);
		d = 1.0/d;
		bcm->tm[i][0] = va[0] *= d;
		bcm->tm[i][1] = va[1] *= d;
		bcm->tm[i][2] = -(v2[ax0]*va[0]+v2[ax1]*va[1]);
					/* rotate vertices */
		vt = v1;
		v1 = v2;
		v2 = v3;
		v3 = vt;
	}
	return(0);
}


void
eval_baryc(wt, p, bcm)		/* evaluate barycentric weights at p */
RREAL	wt[3];
FVECT	p;
register BARYCCM	*bcm;
{
	double	u, v;
	
	u = p[(bcm->ax + 1) % 3];
	v = p[(bcm->ax + 2) % 3];
	wt[0] = u*bcm->tm[0][0] + v*bcm->tm[0][1] + bcm->tm[0][2];
	wt[1] = u*bcm->tm[1][0] + v*bcm->tm[1][1] + bcm->tm[1][2];
	wt[2] = 1. - wt[1] - wt[0];
}


int
get_baryc(wt, p, v1, v2, v3)	/* compute barycentric weights at p */
RREAL	wt[3];
FVECT	p;
FVECT	v1, v2, v3;
{
	BARYCCM	bcm;
	
	if (comp_baryc(&bcm, v1, v2, v3) < 0)
		return(-1);
	eval_baryc(wt, p, &bcm);
	return(0);
}


#if 0
int
get_baryc(wt, p, v1, v2, v3)	/* compute barycentric weights at p */
RREAL	wt[3];
FVECT	p;
FVECT	v1, v2, v3;
{
	FVECT	ac, bc, pc, cros;
	double	normf;
				/* area formula w/o 2-D optimization */
	VSUB(ac, v1, v3);
	VSUB(bc, v2, v3);
	VSUB(pc, p, v3);
	VCROSS(cros, ac, bc);
	normf = DOT(cros,cros)
	if (normf <= 0.0)
		return(-1);
	normf = 1./sqrt(normf);
	VCROSS(cros, bc, pc);
	wt[0] = VLEN(cros) * normf;
	VCROSS(cros, ac, pc);
	wt[1] = VLEN(cros) * normf;
	wt[2] = 1. - wt[1] - wt[0];
	return(0);
}
#endif


void
put_baryc(bcm, com, n)		/* put barycentric coord. vectors */
register BARYCCM	*bcm;
register RREAL		com[][3];
int			n;
{
	double	a, b;
	register int	i;

	printf("%d\t%d\n", 1+3*n, bcm->ax);
	for (i = 0; i < n; i++) {
		a = com[i][0] - com[i][2];
		b = com[i][1] - com[i][2];
		printf("%14.8f %14.8f %14.8f\n",
			bcm->tm[0][0]*a + bcm->tm[1][0]*b,
			bcm->tm[0][1]*a + bcm->tm[1][1]*b,
			bcm->tm[0][2]*a + bcm->tm[1][2]*b + com[i][2]);
	}
}
