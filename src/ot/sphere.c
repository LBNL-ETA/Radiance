#ifndef lint
static const char	RCSid[] = "$Id: sphere.c,v 2.2 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  sphere.c - routines for creating octrees for spheres.
 *
 *     7/28/85
 */

#include  "standard.h"

#include  "octree.h"

#include  "object.h"

#include  "otypes.h"

#define  ROOT3		1.732050808

/*
 *	Regrettably, the algorithm for determining a cube's location
 *  with respect to a sphere is not simple.  First, a quick test is 
 *  made to determine if the sphere and the bounding sphere of the cube 
 *  are disjoint.  This of course means no intersection.  Failing this,
 *  we determine if the cube lies inside the sphere.  The cube is
 *  entirely inside if the bounding sphere on the cube is
 *  contained within our sphere.  This means no intersection.  Otherwise,
 *  if the cube radius is smaller than the sphere's and the cube center is
 *  inside the sphere, we assume intersection.  If these tests fail,
 *  we proceed as follows.
 *	The sphere center is located in relation to the 6 cube faces,
 *  and one of four things is done depending on the number of
 *  planes the center lies between:
 *
 *	0:  The sphere is closest to a cube corner, find the
 *		distance to that corner.
 *
 *	1:  The sphere is closest to a cube edge, find this
 *		distance.
 *
 *	2:  The sphere is closest to a cube face, find the distance.
 *
 *	3:  The sphere has its center inside the cube.
 *
 *	In cases 0-2, if the closest part of the cube is within
 *  the radius distance from the sphere center, we have intersection.
 *  If it is not, the cube must be outside the sphere.
 *	In case 3, there must be intersection, and no further
 *  tests are necessary.
 */


o_sphere(o, cu)			/* determine if sphere intersects cube */
OBJREC  *o;
register CUBE  *cu;
{
	FVECT  v1;
	double  d1, d2;
	register FLOAT  *fa;
	register int  i;
#define  cent		fa
#define  rad		fa[3]
					/* get arguments */
	if (o->oargs.nfargs != 4)
		objerror(o, USER, "bad # arguments");
	fa = o->oargs.farg;
	if (rad < -FTINY) {
		objerror(o, WARNING, "negative radius");
		o->otype = o->otype == OBJ_SPHERE ?
				OBJ_BUBBLE : OBJ_SPHERE;
		rad = -rad;
	} else if (rad <= FTINY)
		objerror(o, USER, "zero radius");

	d1 = ROOT3/2.0 * cu->cusize;	/* bounding radius for cube */

	d2 = cu->cusize * 0.5;		/* get distance between centers */
	for (i = 0; i < 3; i++)
		v1[i] = cu->cuorg[i] + d2 - cent[i];
	d2 = DOT(v1,v1);

	if (d2 > (rad+d1+FTINY)*(rad+d1+FTINY))	/* quick test */
		return(O_MISS);			/* cube outside */
	
					/* check sphere interior */
	if (d1 < rad) {
		if (d2 < (rad-d1-FTINY)*(rad-d1-FTINY))
			return(O_MISS);		/* cube inside sphere */
		if (d2 < (rad+FTINY)*(rad+FTINY))
			return(O_HIT);		/* cube center inside */
	}
					/* find closest distance */
	for (i = 0; i < 3; i++)
		if (cent[i] < cu->cuorg[i])
			v1[i] = cu->cuorg[i] - cent[i];
		else if (cent[i] > cu->cuorg[i] + cu->cusize)
			v1[i] = cent[i] - (cu->cuorg[i] + cu->cusize);
		else
			v1[i] = 0;
					/* final intersection check */
	if (DOT(v1,v1) <= (rad+FTINY)*(rad+FTINY))
		return(O_HIT);
	else
		return(O_MISS);
}
