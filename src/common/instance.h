/* RCSid $Id: instance.h,v 2.4 2003/02/25 02:47:21 greg Exp $ */
/*
 *  instance.h - header file for routines using octree objects.
 *
 *  Include after object.h and octree.h
 */

#include "copyright.h"

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

#ifdef NOPROTO

extern SCENE  *getscene();
extern INSTANCE  *getinstance();
extern void  freescene();
extern void  freeinstance();

#else

extern SCENE  *getscene(char *sname, int flags);
extern INSTANCE  *getinstance(OBJREC *o, int flags);
extern void  freescene(SCENE *sc);
extern void  freeinstance(OBJREC *o);

#endif
