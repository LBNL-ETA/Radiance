#ifndef lint
static const char	RCSid[] = "$Id: plocate.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 *  plocate.c - routine to locate 3D vector w.r.t. box.
 */

#include "copyright.h"

#include  "fvect.h"

#include  "plocate.h"


int
plocate(p, min, max)		/* return location of p w.r.t. min & max */
register FVECT  p;
FVECT  min, max;
{
	register int  loc = 0;

	if (p[0] < min[0] - EPSILON)
		loc |= XPOS & BELOW;
	else if (p[0] > max[0] + EPSILON)
		loc |= XPOS & ABOVE;
	if (p[1] < min[1] - EPSILON)
		loc |= YPOS & BELOW;
	else if (p[1] > max[1] + EPSILON)
		loc |= YPOS & ABOVE;
	if (p[2] < min[2] - EPSILON)
		loc |= ZPOS & BELOW;
	else if (p[2] > max[2] + EPSILON)
		loc |= ZPOS & ABOVE;
	
	return(loc);
}
