/* RCSid: $Id: object.h,v 2.8 2003/02/22 02:07:22 greg Exp $ */
/*
 *  object.h - header file for routines using objects and object sets.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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
	FLOAT  *farg;			/* real arguments */
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
#ifdef  BIGMEM
#define  OBJECT		int4		/* index to object array */
#else
#define  OBJECT		int2		/* index to object array */
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
#ifdef  BIGMEM
#define  MAXOBJBLK	65535		/* maximum number of object blocks */
#else
#define  MAXOBJBLK	63		/* maximum number of object blocks */
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

#define  MAXSET		127		/* maximum object set size */

#define setfree(os)	free((void *)(os))

extern void  (*addobjnotify[])();        /* people to notify of new objects */

#ifdef NOPROTO

extern int	objndx();
extern int	lastmod();
extern int	modifier();
extern int	object();
extern void	insertobject();
extern void	clearobjndx();
extern void	insertelem();
extern void	deletelem();
extern int	inset();
extern int	setequal();
extern void	setcopy();
extern OBJECT *	setsave();
extern void	setunion();
extern void	setintersect();
extern OCTREE	fullnode();
extern void	objset();
extern int	dosets();
extern void	donesets();
extern int	otype();
extern void	objerror();
extern int	readfargs();
extern void	freefargs();
extern void	readobj();
extern void	getobject();
extern int	newobject();
extern void	freeobjects();
extern int	free_os();

#else
					/* defined in modobject.c */
extern int	objndx(OBJREC *op);
extern int	lastmod(OBJECT obj, char *mname);
extern int	modifier(char *name);
extern int	object(char *oname);
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
extern int	newobject();
extern void	freeobjects(int firstobj, int nobjs);
					/* defined in free_os.c */
extern int	free_os(OBJREC *op);

#endif
