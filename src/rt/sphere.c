#ifndef lint
static const char RCSid[] = "$Id: sphere.c,v 2.5 2003/03/12 04:59:05 greg Exp $";
#endif
/*
 *  sphere.c - compute ray intersection with spheres.
 */

#include "copyright.h"

#include  "ray.h"

#include  "otypes.h"


o_sphere(so, r)			/* compute intersection with sphere */
OBJREC  *so;
register RAY  *r;
{
	double  a, b, c;	/* coefficients for quadratic equation */
	double  root[2];	/* quadratic roots */
	int  nroots;
	double  t;
	register FLOAT  *ap;
	register int  i;

	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad # arguments");
	ap = so->oargs.farg;
	if (ap[3] < -FTINY) {
		objerror(so, WARNING, "negative radius");
		so->otype = so->otype == OBJ_SPHERE ?
				OBJ_BUBBLE : OBJ_SPHERE;
		ap[3] = -ap[3];
	} else if (ap[3] <= FTINY)
		objerror(so, USER, "zero radius");

	/*
	 *	We compute the intersection by substituting into
	 *  the surface equation for the sphere.  The resulting
	 *  quadratic equation in t is then solved for the
	 *  smallest positive root, which is our point of
	 *  intersection.
	 *	Since the ray is normalized, a should always be
	 *  one.  We compute it here to prevent instability in the
	 *  intersection calculation.
	 */
				/* compute quadratic coefficients */
	a = b = c = 0.0;
	for (i = 0; i < 3; i++) {
		a += r->rdir[i]*r->rdir[i];
		t = r->rorg[i] - ap[i];
		b += 2.0*r->rdir[i]*t;
		c += t*t;
	}
	c -= ap[3] * ap[3];

	nroots = quadratic(root, a, b, c);	/* solve quadratic */
	
	for (i = 0; i < nroots; i++)		/* get smallest positive */
		if ((t = root[i]) > FTINY)
			break;
	if (i >= nroots)
		return(0);			/* no positive root */

	if (t >= r->rot)
		return(0);			/* other is closer */

	r->ro = so;
	r->rot = t;
					/* compute normal */
	a = ap[3];
	if (so->otype == OBJ_BUBBLE)
		a = -a;			/* reverse */
	for (i = 0; i < 3; i++) {
		r->rop[i] = r->rorg[i] + r->rdir[i]*t;
		r->ron[i] = (r->rop[i] - ap[i]) / a;
	}
	r->rod = -DOT(r->rdir, r->ron);
	r->rox = NULL;
	r->pert[0] = r->pert[1] = r->pert[2] = 0.0;
	r->uv[0] = r->uv[1] = 0.0;

	return(1);			/* hit */
}
