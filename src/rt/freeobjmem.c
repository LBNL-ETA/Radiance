#ifndef lint
static const char	RCSid[] = "$Id: freeobjmem.c,v 2.6 2003/07/10 03:47:00 greg Exp $";
#endif
/*
 * Free memory associated with object(s)
 *
 * External symbols declared in ray.h
 */

#include "copyright.h"

#include "ray.h"
#include "otypes.h"
#include "face.h"
#include "cone.h"
#include "instance.h"
#include "data.h"
#include "font.h"


int
free_os(op)			/* free unneeded memory for object */
register OBJREC	*op;
{
	if (op->os == NULL)
		return(0);
	if (hasfunc(op->otype)) {
		freefunc(op);
		return(1);
	}
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
	case PAT_BTEXT:		/* monochromatic text */
	case PAT_CTEXT:		/* colored text */
	case MIX_TEXT:		/* mixing text */
		freetext(op);
		return(1);
	case MAT_CLIP:		/* clipping surface */
	case MAT_SPOT:		/* spot light source */
		free((void *)op->os);
		op->os = NULL;
		return(1);
	}
#ifdef DEBUG
	objerror(op, WARNING, "cannot free structure");
#endif
	return(0);
}


int
free_objs(on, no)		/* free some object structures */
register OBJECT	on;
OBJECT	no;
{
	int	nfreed;
	register OBJREC	*op;

	for (nfreed = 0; no-- > 0; on++) {
		op = objptr(on);
		if (op->os != NULL)
			nfreed += free_os(op);
	}
	return(nfreed);
}


void
free_objmem()			/* free all object cache memory */
{
	free_objs(0, nobjects);
	freedata(NULL);
	freefont(NULL);
}
