#ifndef lint
static const char	RCSid[] = "$Id: o_face.c,v 2.3 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
 *  o_face.c - routines for creating octrees for polygonal faces.
 *
 *     8/27/85
 */

#include  "standard.h"

#include  "octree.h"

#include  "object.h"

#include  "face.h"

#include  "plocate.h"

/*
 *	The algorithm for determining a face's intersection
 *  with a cube is relatively straightforward:
 *
 *	1) Check to see if any vertices are inside the cube
 *	   (intersection).
 *
 *	2) Check to see if all vertices are to one side of
 *	   cube (no intersection).
 *
 *	3) Check to see if any portion of any edge is inside
 *	   cube (intersection).
 *
 *	4) Check to see if the cube cuts the plane of the
 *	   face and one of its edges passes through
 *	   the face (intersection).
 *
 *	5) If test 4 fails, we have no intersection.
 */


o_face(o, cu)			/* determine if face intersects cube */
OBJREC  *o;
CUBE  *cu;
{
	FVECT  cumin, cumax;
	FVECT  v1, v2;
	double  d1, d2;
	int  vloc;
	register FACE  *f;
	register int  i, j;
				/* get face arguments */
	f = getface(o);
	if (f->area == 0.0)		/* empty face */
		return(O_MISS);
					/* compute cube boundaries */
	for (j = 0; j < 3; j++)
		cumax[j] = (cumin[j] = cu->cuorg[j]-FTINY)
				+ cu->cusize + 2.0*FTINY;

	vloc = ABOVE | BELOW;		/* check vertices */
	for (i = 0; i < f->nv; i++)
		if ( (j = plocate(VERTEX(f,i), cumin, cumax)) )
			vloc &= j;
		else
			return(O_HIT);	/* vertex inside */

	if (vloc)			/* all to one side */
		return(O_MISS);
	
	for (i = 0; i < f->nv; i++) {	/* check edges */
		if ((j = i + 1) >= f->nv)
			j = 0;			/* wrap around */
		VCOPY(v1, VERTEX(f,i));		/* clip modifies */
		VCOPY(v2, VERTEX(f,j));		/* the vertices! */
		if (clip(v1, v2, cumin, cumax))
			return(O_HIT);		/* edge inside */
	}
					/* see if cube cuts plane */
	for (j = 0; j < 3; j++)
		if (f->norm[j] > 0.0) {
			v1[j] = cumin[j];
			v2[j] = cumax[j];
		} else {
			v1[j] = cumax[j];
			v2[j] = cumin[j];
		}
	if ((d1 = DOT(v1, f->norm) - f->offset) > FTINY)
		return(O_MISS);
	if ((d2 = DOT(v2, f->norm) - f->offset) < -FTINY)
		return(O_MISS);
					/* intersect face */
	for (j = 0; j < 3; j++)
		v1[j] = (v1[j]*d2 - v2[j]*d1)/(d2 - d1);
	if (inface(v1, f))
		return(O_HIT);
	
	return(O_MISS);		/* no intersection */
}
