/* RCSid $Id: face.h,v 2.6 2003/06/27 06:53:21 greg Exp $ */
/*
 *  face.h - header for routines using polygonal faces.
 */
#ifndef _RAD_FACE_H_
#define _RAD_FACE_H_
#ifdef __cplusplus
extern "C" {
#endif

#define  VERTEX(f,n)	((f)->va + 3*(n))

typedef struct {	/* a polygonal face */
	FVECT  norm;		/* the plane's unit normal */
	RREAL  offset;		/* plane equation:  DOT(norm, v) == offset */
	RREAL  area;		/* area of face */
	RREAL  *va;		/* vertex array (o->oargs.farg) */
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

