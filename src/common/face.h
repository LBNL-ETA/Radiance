/* RCSid $Id: face.h,v 2.3 2003/02/25 02:47:21 greg Exp $ */
/*
 *  face.h - header for routines using polygonal faces.
 */

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

#ifdef NOPROTO

extern FACE  *getface();
extern void  freeface();
extern int  inface();

#else

extern FACE  *getface(OBJREC *o);
extern void  freeface(OBJREC *o);
extern int  inface(FVECT p, FACE *f);

#endif
