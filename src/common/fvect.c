#ifndef lint
static const char	RCSid[] = "$Id: fvect.c,v 2.10 2008/06/24 02:16:14 greg Exp $";
#endif
/*
 *  fvect.c - routines for floating-point vector calculations
 */

#include "copyright.h"

#include  <math.h>
#include  "fvect.h"


double
fdot(				/* return the dot product of two vectors */
register FVECT v1,
register FVECT v2
)
{
	return(DOT(v1,v2));
}


double
dist2(				/* return square of distance between points */
register FVECT p1,
register FVECT p2
)
{
	FVECT  delta;

	delta[0] = p2[0] - p1[0];
	delta[1] = p2[1] - p1[1];
	delta[2] = p2[2] - p1[2];

	return(DOT(delta, delta));
}


double
dist2line(			/* return square of distance to line */
FVECT p,		/* the point */
FVECT ep1,
FVECT ep2		/* points on the line */
)
{
	register double  d, d1, d2;

	d = dist2(ep1, ep2);
	d1 = dist2(ep1, p);
	d2 = d + d1 - dist2(ep2, p);

	return(d1 - 0.25*d2*d2/d);
}


double
dist2lseg(			/* return square of distance to line segment */
FVECT p,		/* the point */
FVECT ep1,
FVECT ep2		/* the end points */
)
{
	register double  d, d1, d2;

	d = dist2(ep1, ep2);
	d1 = dist2(ep1, p);
	d2 = dist2(ep2, p);

	if (d2 > d1) {			/* check if past endpoints */
		if (d2 - d1 > d)
			return(d1);
	} else {
		if (d1 - d2 > d)
			return(d2);
	}
	d2 = d + d1 - d2;

	return(d1 - 0.25*d2*d2/d);	/* distance to line */
}


void
fcross(				/* vres = v1 X v2 */
register FVECT vres,
register FVECT v1,
register FVECT v2
)
{
	vres[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vres[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vres[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


void
fvsum(				/* vres = v0 + f*v1 */
register FVECT vres,
register FVECT v0,
register FVECT v1,
register double f
)
{
	vres[0] = v0[0] + f*v1[0];
	vres[1] = v0[1] + f*v1[1];
	vres[2] = v0[2] + f*v1[2];
}


double
normalize(			/* normalize a vector, return old magnitude */
register FVECT  v
)
{
	register double  len, d;
	
	d = DOT(v, v);
	
	if (d == 0.0)
		return(0.0);
	
	if (d <= 1.0+FTINY && d >= 1.0-FTINY)
		len = 0.5 + 0.5*d;	/* first order approximation */
	else
		len = sqrt(d);

	v[0] *= d = 1.0/len;
	v[1] *= d;
	v[2] *= d;

	return(len);
}


int
closestapproach(			/* closest approach of two rays */
RREAL t[2],		/* returned distances along each ray */
FVECT rorg0,		/* first origin */
FVECT rdir0,		/* first direction (normalized) */
FVECT rorg1,		/* second origin */
FVECT rdir1		/* second direction (normalized) */
)
{
	double	dotprod = DOT(rdir0, rdir1);
	double	denom = 1. - dotprod*dotprod;
	double	o1o2_d1;
	FVECT	o0o1;

	if (denom <= FTINY) {		/* check if lines are parallel */
		t[0] = t[1] = 0.0;
		return(0);
	}
	VSUB(o0o1, rorg0, rorg1);
	o1o2_d1 = DOT(o0o1, rdir1);
	t[0] = (o1o2_d1*dotprod - DOT(o0o1,rdir0)) / denom;
	t[1] = o1o2_d1 + t[0]*dotprod;
	return(1);
}


#if 0
int
closestapproach(			/* closest approach of two rays */
RREAL t[2],		/* returned distances along each ray */
FVECT rorg0,		/* first origin */
FVECT rdir0,		/* first direction (unnormalized) */
FVECT rorg1,		/* second origin */
FVECT rdir1		/* second direction (unnormalized) */
)
{
	double	dotprod = DOT(rdir0, rdir1);
	double	d0n2 = DOT(rdir0, rdir0);
	double  d1n2 = DOT(rdir1, rdir1);
	double	denom = d0n2*d1n2 - dotprod*dotprod;
	double	o1o2_d1;
	FVECT	o0o1;

	if (denom <= FTINY) {		/* check if lines are parallel */
		t[0] = t[1] = 0.0;
		return(0);
	}
	VSUB(o0o1, rorg0, rorg1);
	o1o2_d1 = DOT(o0o1, rdir1);
	t[0] = (o1o2_d1*dotprod - DOT(o0o1,rdir0)*d1n2) / denom;
	t[1] = (o1o2_d1 + t[0]*dotprod) / d1n2;
	return(1);
}
#endif


void
spinvector(				/* rotate vector around normal */
FVECT vres,		/* returned vector */
FVECT vorig,		/* original vector */
FVECT vnorm,		/* normalized vector for rotation */
double theta		/* left-hand radians */
)
{
	double  sint, cost, normprod;
	FVECT  vperp;
	register int  i;
	
	if (theta == 0.0) {
		if (vres != vorig)
			VCOPY(vres, vorig);
		return;
	}
	cost = cos(theta);
	sint = sin(theta);
	normprod = DOT(vorig, vnorm)*(1.-cost);
	fcross(vperp, vnorm, vorig);
	for (i = 0; i < 3; i++)
		vres[i] = vorig[i]*cost + vnorm[i]*normprod + vperp[i]*sint;
}
