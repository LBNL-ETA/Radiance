#ifndef lint
static const char	RCSid[] = "$Id: init2otypes.c,v 2.3 2004/03/27 12:41:45 schorsch Exp $";
#endif
/*
 * Initialize ofun[] list for bounding box checker
 */

#include  "standard.h"

#include  "octree.h"

#include  "otypes.h"

FUN  ofun[NUMOTYPE] = INIT_OTYPE;


int
o_default(void)			/* default action is no intersection */
{
	return(O_MISS);
}
