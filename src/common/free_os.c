#ifndef lint
static const char	RCSid[] = "$Id: free_os.c,v 3.3 2003/07/10 03:47:00 greg Exp $";
#endif
/*
 * Free memory associated with object(s)
 *
 *   External symbols declared in object.h
 */

#include "copyright.h"

#include "standard.h"
#include "octree.h"
#include "object.h"
#include "otypes.h"
#include "face.h"
#include "cone.h"
#include "instance.h"


int
free_os(op)			/* free unneeded memory for object */
register OBJREC	*op;
{
	if (op->os == NULL)
		return(0);
	switch (op->otype) {
	case OBJ_FACE:		/* polygon */
		freeface(op);
		return(1);
	case OBJ_CONE:		/* cone */
	case OBJ_RING:		/* disk */
	case OBJ_CYLINDER:	/* cylinder */
	case OBJ_CUP:		/* inverted cone */
	case OBJ_TUBE:		/* inverted cylinder */
		freecone(op);
		return(1);
	case OBJ_INSTANCE:	/* octree instance */
		freeinstance(op);
		return(1);
	case OBJ_MESH:		/* mesh instance */
		freemeshinst(op);
		return(1);
	}
				/* don't really know */
	free((void *)op->os);
	op->os = NULL;
	return(1);
}
