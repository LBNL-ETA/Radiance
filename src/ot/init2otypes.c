#ifndef lint
static const char	RCSid[] = "$Id: init2otypes.c,v 2.2 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 * Initialize ofun[] list for bounding box checker
 */

#include  "standard.h"

#include  "octree.h"

#include  "otypes.h"

FUN  ofun[NUMOTYPE] = INIT_OTYPE;


o_default()			/* default action is no intersection */
{
	return(O_MISS);
}
