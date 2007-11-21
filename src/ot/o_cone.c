#ifndef lint
static const char	RCSid[] = "$Id: o_cone.c,v 2.6 2007/11/21 18:51:04 greg Exp $";
#endif
/*
 *  o_cone.c - routines for intersecting cubes with cones.
 *
 *     2/3/86
 */

#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "cone.h"
#include  "plocate.h"

#define  ROOT3		1.732050808

/*
 *     The algorithm used to detect cube intersection with cones is
 *  recursive.  First, we approximate the cube to be a sphere.  Then
 *  we test for cone intersection with the sphere by testing the
 *  segment of the cone which is nearest the sphere's center.
 *     If the cone has points within the cube's bounding sphere,
 *  we must check for intersection with the cube.  This is done with
 *  the 3D line clipper.  The same cone segment is used in this test.
 *  If the clip fails, we still cannot be sure there is no intersection,
 *  so we subdivide the cube and recurse.
 *     If none of the sub-cubes intersect, then our cube does not intersect.
 */

extern double  mincusize;		/* minimum cube size */

static int findcseg(FVECT ep0, FVECT ep1, CONE *co, FVECT p);



extern int
o_cone(		/* determine if cone intersects cube */
	OBJREC  *o,
	register CUBE  *cu
)
{
	CONE  *co;
	FVECT  ep0, ep1;
#ifdef STRICT
	FVECT  cumin, cumax;
	CUBE  cukid;
	register int  j;
#endif
	double  r;
	FVECT  p;
	register int  i;
					/* get cone arguments */
	co = getcone(o, 0);
					/* get cube center */
	r = cu->cusize * 0.5;
	for (i = 0; i < 3; i++)
		p[i] = cu->cuorg[i] + r;
	r *= ROOT3;			/* bounding radius for cube */

	if (findcseg(ep0, ep1, co, p)) {
					/* check min. distance to cone */
		if (dist2lseg(p, ep0, ep1) > (r+FTINY)*(r+FTINY))
			return(O_MISS);
#ifdef  STRICT
					/* get cube boundaries */
		for (i = 0; i < 3; i++)
			cumax[i] = (cumin[i] = cu->cuorg[i]) + cu->cusize;
					/* closest segment intersects? */
		if (clip(ep0, ep1, cumin, cumax))
			return(O_HIT);
	}
					/* check sub-cubes */
	cukid.cusize = cu->cusize * 0.5;
	if (cukid.cusize < mincusize)
		return(O_HIT);		/* cube too small */
	cukid.cutree = EMPTY;

	for (j = 0; j < 8; j++) {
		for (i = 0; i < 3; i++) {
			cukid.cuorg[i] = cu->cuorg[i];
			if (1<<i & j)
				cukid.cuorg[i] += cukid.cusize;
		}
		if (o_cone(o, &cukid))
			return(O_HIT);	/* sub-cube intersects */
	}
	return(O_MISS);			/* no intersection */
#else
	}
	return(O_HIT);			/* assume intersection */
#endif
}


static int
findcseg(	/* find line segment from cone closest to p */
	FVECT  ep0,
	FVECT  ep1,
	register CONE  *co,
	FVECT  p
)
{
	double  d;
	FVECT  v;
	register int  i;
					/* find direction from axis to point */
	VSUB(v, p, CO_P0(co));
	d = DOT(v, co->ad);
	for (i = 0; i < 3; i++)
		v[i] -= d*co->ad[i];
	if (normalize(v) == 0.0)
		return(0);
					/* find endpoints of segment */
	for (i = 0; i < 3; i++) {
		ep0[i] = CO_R0(co)*v[i] + CO_P0(co)[i];
		ep1[i] = CO_R1(co)*v[i] + CO_P1(co)[i];
	}
	return(1);			/* return distance from axis */
}
