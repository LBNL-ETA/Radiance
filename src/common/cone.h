/* RCSid $Id: cone.h,v 2.8 2005/09/23 19:04:52 greg Exp $ */
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
#ifndef _RAD_CONE_H_
#define _RAD_CONE_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct cone {
	FVECT  ad;		/* axis direction vector */
	RREAL  al;		/* axis length */
	RREAL  sl;		/* side length */
	RREAL  *ca;		/* cone arguments (o->oargs.farg) */
	RREAL  (*tm)[4];	/* pointer to transformation matrix */
	char  p0, p1;		/* indices for endpoints */
	char  r0, r1;		/* indices for radii */
}  CONE;

#define  CO_R0(co)	((co)->ca[(int)((co)->r0)])
#define  CO_R1(co)	((co)->ca[(int)((co)->r1)])
#define  CO_P0(co)	((co)->ca+(co)->p0)
#define  CO_P1(co)	((co)->ca+(co)->p1)


extern CONE  *getcone(OBJREC *o, int getxf);
extern void  freecone(OBJREC *o);
extern void  conexform(CONE *co);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_CONE_H_ */

