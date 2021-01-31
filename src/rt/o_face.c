#ifndef lint
static const char RCSid[] = "$Id: o_face.c,v 2.8 2021/01/31 18:08:04 greg Exp $";
#endif
/*
 *  o_face.c - compute ray intersection with faces.
 */

#include "copyright.h"

#include  "ray.h"
#include  "face.h"
#include  "rtotypes.h"


int
o_face(		/* compute intersection with polygonal face */
	OBJREC  *o,
	RAY  *r
)
{
	double  rdot;		/* direction . normal */
	double  t;		/* distance to intersection */
	FVECT  pisect;		/* intersection point */
	FACE  *f;	/* face record */

	f = getface(o);
		
	/*
	 *  First, we find the distance to the plane containing the
	 *  face.  If the plane is parallel to our ray, or the
	 *  previous hit was better, we return.  Otherwise,
	 *  we determine whether in fact the ray intersects the
	 *  face.  The ray intersects the face if the
	 *  point of intersection with the plane of the face
	 *  is inside the face.
	 */
						/* compute dist. to plane */
	rdot = -DOT(r->rdir, f->norm);
	if ((rdot <= FTINY) & (rdot >= -FTINY))	/* ray parallels plane */
		return(0);

	t = (DOT(r->rorg, f->norm) - f->offset) / rdot;
	
	if (rayreject(o, r, t))			/* previous hit is better? */
		return(0);
						/* compute intersection */
	VSUM(pisect, r->rorg, r->rdir, t);

	if (!inface(pisect, f))			/* ray intersects face? */
		return(0);

	r->ro = o;
	r->rot = t;
	VCOPY(r->rop, pisect);
	VCOPY(r->ron, f->norm);
	r->rod = rdot;
	r->pert[0] = r->pert[1] = r->pert[2] = 0.0;
	r->uv[0] = r->uv[1] = 0.0;
	r->rox = NULL;

	return(1);				/* hit */
}
