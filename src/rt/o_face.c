#ifndef lint
static const char RCSid[] = "$Id: o_face.c,v 2.4 2003/03/11 17:08:55 greg Exp $";
#endif
/*
 *  o_face.c - compute ray intersection with faces.
 */

#include "copyright.h"

#include  "ray.h"

#include  "face.h"


o_face(o, r)		/* compute intersection with polygonal face */
OBJREC  *o;
register RAY  *r;
{
	double  rdot;		/* direction . normal */
	double  t;		/* distance to intersection */
	FVECT  pisect;		/* intersection point */
	register FACE  *f;	/* face record */
	register int  i;

	f = getface(o);
		
	/*
	 *  First, we find the distance to the plane containing the
	 *  face.  If this distance is less than zero or greater
	 *  than a previous intersection, we return.  Otherwise,
	 *  we determine whether in fact the ray intersects the
	 *  face.  The ray intersects the face if the
	 *  point of intersection with the plane of the face
	 *  is inside the face.
	 */
						/* compute dist. to plane */
	rdot = -DOT(r->rdir, f->norm);
	if (rdot <= FTINY && rdot >= -FTINY)	/* ray parallels plane */
		t = FHUGE;
	else
		t = (DOT(r->rorg, f->norm) - f->offset) / rdot;
	
	if (t <= FTINY || t >= r->rot)		/* not good enough */
		return(0);
						/* compute intersection */
	for (i = 0; i < 3; i++)
		pisect[i] = r->rorg[i] + r->rdir[i]*t;

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
