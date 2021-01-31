#ifndef lint
static const char RCSid[] = "$Id: o_cone.c,v 2.10 2021/01/31 18:08:04 greg Exp $";
#endif
/*
 *  o_cone.c - routine to determine ray intersection with cones.
 */

#include "copyright.h"

#include  "ray.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "cone.h"


int
o_cone(			/* intersect ray with cone */
	OBJREC  *o,
	RAY  *r
)
{
	FVECT  rox, rdx;
	double  a, b, c;
	double  root[2];
	int  nroots, rn;
	CONE  *co;
	int  i;

						/* get cone structure */
	co = getcone(o, 1);
	if (co == NULL)
		objerror(o, INTERNAL, "unexpected illegal");

	/*
	 *     To intersect a ray with a cone, we transform the
	 *  ray into the cone's normalized space.  This greatly
	 *  simplifies the computation.
	 *     For a cone or cup, normalization results in the
	 *  equation:
	 *
	 *		x*x + y*y - z*z == 0
	 *
	 *     For a cylinder or tube, the normalized equation is:
	 *
	 *		x*x + y*y - r*r == 0
	 *
	 *     A normalized ring obeys the following set of equations:
	 *
	 *		z == 0			&&
	 *		x*x + y*y >= r0*r0	&&
	 *		x*x + y*y <= r1*r1
	 */

					/* transform ray */
	multp3(rox, r->rorg, co->tm);
	multv3(rdx, r->rdir, co->tm);
					/* compute intersection */

	if (o->otype == OBJ_CONE || o->otype == OBJ_CUP) {

		a = rdx[0]*rdx[0] + rdx[1]*rdx[1] - rdx[2]*rdx[2];
		b = 2.0*(rdx[0]*rox[0] + rdx[1]*rox[1] - rdx[2]*rox[2]);
		c = rox[0]*rox[0] + rox[1]*rox[1] - rox[2]*rox[2];

	} else if (o->otype == OBJ_CYLINDER || o->otype == OBJ_TUBE) {

		a = rdx[0]*rdx[0] + rdx[1]*rdx[1];
		b = 2.0*(rdx[0]*rox[0] + rdx[1]*rox[1]);
		c = rox[0]*rox[0] + rox[1]*rox[1] - CO_R0(co)*CO_R0(co);

	} else { /* OBJ_RING */

		if ((rdx[2] <= FTINY) & (rdx[2] >= -FTINY))
			return(0);			/* parallel */
		root[0] = -rox[2]/rdx[2];
		if (rayreject(o, r, root[0]))
			return(0);			/* have better */
		b = root[0]*rdx[0] + rox[0];
		c = root[0]*rdx[1] + rox[1];
		a = b*b + c*c;
		if (a < CO_R0(co)*CO_R0(co) || a > CO_R1(co)*CO_R1(co))
			return(0);			/* outside radii */
		r->ro = o;
		r->rot = root[0];
		VSUM(r->rop, r->rorg, r->rdir, r->rot);
		VCOPY(r->ron, co->ad);
		r->rod = -rdx[2];
		r->rox = NULL;
		return(1);				/* good */
	}
					/* roots for cone, cup, cyl., tube */
	nroots = quadratic(root, a, b, c);

	for (rn = 0; rn < nroots; rn++) {	/* check real roots */
		if (root[rn] <= FTINY)
			continue;		/* too small */
		if (root[rn] > r->rot + FTINY)
			break;			/* too big */
						/* check endpoints */
		VSUM(rox, r->rorg, r->rdir, root[rn]);
		VSUB(rdx, rox, CO_P0(co));
		b = DOT(rdx, co->ad); 
		if (b < 0.0)
			continue;		/* before p0 */
		if (b > co->al)
			continue;		/* after p1 */
		if (rayreject(o, r, root[rn]))
			break;			/* previous hit better */
		r->ro = o;
		r->rot = root[rn];
		VCOPY(r->rop, rox);
						/* get normal */
		if (o->otype == OBJ_CYLINDER)
			a = CO_R0(co);
		else if (o->otype == OBJ_TUBE)
			a = -CO_R0(co);
		else { /* OBJ_CONE || OBJ_CUP */
			c = CO_R1(co) - CO_R0(co);
			a = CO_R0(co) + b*c/co->al;
			if (o->otype == OBJ_CUP) {
				c = -c;
				a = -a;
			}
		}
		for (i = 0; i < 3; i++)
			r->ron[i] = (rdx[i] - b*co->ad[i])/a;
		if (o->otype == OBJ_CONE || o->otype == OBJ_CUP)
			for (i = 0; i < 3; i++)
				r->ron[i] = (co->al*r->ron[i] - c*co->ad[i])
						/ co->sl;
		a = DOT(r->ron, r->ron);
		if (a > 1.+FTINY || a < 1.-FTINY) {
			c = 1./(.5 + .5*a);     /* avoid numerical error */
			r->ron[0] *= c; r->ron[1] *= c; r->ron[2] *= c;
		}
		r->rod = -DOT(r->rdir, r->ron);
		r->pert[0] = r->pert[1] = r->pert[2] = 0.0;
		r->uv[0] = r->uv[1] = 0.0;
		r->rox = NULL;
		return(1);			/* good */
	}
	return(0);
}
