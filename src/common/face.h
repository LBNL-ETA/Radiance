/* RCSid $Id: face.h,v 2.4 2003/06/06 16:38:47 schorsch Exp $ */
/*
 *  face.h - header for routines using polygonal faces.
 */
#ifndef _RAD_FACE_H_
#define _RAD_FACE_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#define  VERTEX(f,n)	((f)->va + 3*(n))

typedef struct {	/* a polygonal face */
	FVECT  norm;		/* the plane's unit normal */
	FLOAT  offset;		/* plane equation:  DOT(norm, v) == offset */
	FLOAT  area;		/* area of face */
	FLOAT  *va;		/* vertex array (o->oargs.farg) */
	short  nv;		/* # of vertices */
	short  ax;		/* axis closest to normal */
} FACE;


extern FACE  *getface(OBJREC *o);
extern void  freeface(OBJREC *o);
extern int  inface(FVECT p, FACE *f);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_FACE_H_ */

