/* RCSid $Id: cone.h,v 2.3 2003/02/25 02:47:21 greg Exp $ */
/*
 *  cone.h - header file for cones (cones, cylinders, rings, cups, tubes).
 *
 *	Storage of arguments in the cone structure is a little strange.
 *  To save space, we use an index into the real arguments of the
 *  object structure through ca.  The indices are for the axis
 *  endpoints and radii:  p0, p1, r0 and r1.
 *
 *     2/12/86
 */

#include "copyright.h"

typedef struct cone {
	FLOAT  *ca;		/* cone arguments (o->oargs.farg) */
	char  p0, p1;		/* indices for endpoints */
	char  r0, r1;		/* indices for radii */
	FVECT  ad;		/* axis direction vector */
	FLOAT  al;		/* axis length */
	FLOAT  sl;		/* side length */
	FLOAT  (*tm)[4];	/* pointer to transformation matrix */
}  CONE;

#define  CO_R0(co)	((co)->ca[(co)->r0])
#define  CO_R1(co)	((co)->ca[(co)->r1])
#define  CO_P0(co)	((co)->ca+(co)->p0)
#define  CO_P1(co)	((co)->ca+(co)->p1)

#ifdef NOPROTO

extern CONE  *getcone();
extern void  freecone();
extern void  conexform();

#else

extern CONE  *getcone(OBJREC *o, int getxf);
extern void  freecone(OBJREC *o);
extern void  conexform(CONE *co);

#endif
