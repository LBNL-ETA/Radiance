/* RCSid $Id: instance.h,v 2.6 2003/06/27 06:53:21 greg Exp $ */
/*
 *  instance.h - header file for routines using octree objects.
 *
 *  Include after object.h and octree.h
 */
#ifndef _RAD_INSTANCE_H_
#define _RAD_INSTANCE_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct scene {
	char  *name;			/* octree name */
	int  nref;			/* number of references */
	int  ldflags;			/* what was loaded */
	CUBE  scube;			/* scene cube */
	OBJECT  firstobj, nobjs;	/* first object and count */
	struct scene  *next;		/* next in list */
}  SCENE;			/* loaded octree */

typedef struct {
	FULLXF  x;			/* forward and backward transforms */
	SCENE  *obj;			/* loaded object */
}  INSTANCE;			/* instance of octree */


extern SCENE  *getscene(char *sname, int flags);
extern INSTANCE  *getinstance(OBJREC *o, int flags);
extern void  freescene(SCENE *sc);
extern void  freeinstance(OBJREC *o);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_INSTANCE_H_ */

