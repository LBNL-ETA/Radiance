/* RCSid $Id: object.h,v 2.15 2003/06/26 00:58:09 schorsch Exp $ */
/*
 *  object.h - header file for routines using objects and object sets.
 *
 *  Include after "standard.h"
 */
#ifndef _RAD_OBJECT_H_
#define _RAD_OBJECT_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#ifndef  OCTREE
#define  OCTREE		int
#endif

/*
 *	Object definitions require general specifications
 *	which may include a number of different argument types.
 *	The following structure aids in the storage of such
 *	argument lists.
 */

typedef struct {
	short  nsargs;			/* # of string arguments */
	short  nfargs;			/* # of real arguments */
	char  **sarg;			/* string arguments */
	RREAL  *farg;			/* real arguments */
#ifdef  IARGS
	short  niargs;			/* # of integer arguments */
	long  *iarg;			/* integer arguments */
#endif
}  FUNARGS;

#define  MAXSTR		128		/* maximum string length */

/*
 *	An object is defined as an index into an array of
 *	structures containing the object type and specification
 *	and the modifier index.
 */

#ifndef  OBJECT
#ifdef  SMLMEM
#define  OBJECT		int16		/* index to object array */
#else
#define  OBJECT		int32		/* index to object array */
#endif
#endif

typedef struct {
	OBJECT  omod;			/* modifier number */
	short  otype;			/* object type number */
	char  *oname;			/* object name */
	FUNARGS  oargs;			/* object specification */
	char  *os;			/* object structure */
}  OBJREC;

#ifndef  MAXOBJBLK
#ifdef  SMLMEM
#define  MAXOBJBLK	63		/* maximum number of object blocks */
#else
#define  MAXOBJBLK	65535		/* maximum number of object blocks */
#endif
#endif

extern OBJREC  *objblock[MAXOBJBLK];	/* the object blocks */
extern OBJECT  nobjects;		/* # of objects */

#define  OBJBLKSHFT	9
#define  OBJBLKSIZ	(1<<OBJBLKSHFT)	/* object block size */
#define  objptr(obj)	(objblock[(obj)>>OBJBLKSHFT]+((obj)&(OBJBLKSIZ-1)))

#define  OVOID		(-1)		/* void object */
#define  VOIDID		"void"		/* void identifier */

/*
 *     Object sets begin with the number of objects and proceed with
 *  the objects in ascending order.
 */

#define  MAXSET		511		/* maximum object set size */

#define setfree(os)	free((void *)(os))

extern void  (*addobjnotify[])();        /* people to notify of new objects */

					/* defined in modobject.c */
extern OBJECT	objndx(OBJREC *op);
extern OBJECT	lastmod(OBJECT obj, char *mname);
extern OBJECT	modifier(char *name);
extern OBJECT	object(char *oname);
extern void	insertobject(OBJECT obj);
extern void	clearobjndx(void);
					/* defined in objset.c */
extern void	insertelem(OBJECT *os, OBJECT obj);
extern void	deletelem(OBJECT *os, OBJECT obj);
extern int	inset(OBJECT *os, OBJECT obj);
extern int	setequal(OBJECT *os1, OBJECT *os2);
extern void	setcopy(OBJECT *os1, OBJECT *os2);
extern OBJECT *	setsave(OBJECT *os);
extern void	setunion(OBJECT *osr, OBJECT *os1, OBJECT *os2);
extern void	setintersect(OBJECT *osr, OBJECT *os1, OBJECT *os2);
extern OCTREE	fullnode(OBJECT *oset);
extern void	objset(OBJECT *oset, OCTREE ot);
extern int	dosets(int (*f)());
extern void	donesets(void);

					/* defined in otypes.c */
extern int	otype(char *ofname);
extern void	objerror(OBJREC *o, int etyp, char *msg);
					/* defined in readfargs.c */
extern int	readfargs(FUNARGS *fa, FILE *fp);
extern void	freefargs(FUNARGS *fa);
					/* defined in readobj.c */
extern void	readobj(char *inpspec);
extern void	getobject(char *name, FILE *fp);
extern OBJECT	newobject();
extern void	freeobjects(int firstobj, int nobjs);
					/* defined in free_os.c */
extern int	free_os(OBJREC *op);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_OBJECT_H_ */

