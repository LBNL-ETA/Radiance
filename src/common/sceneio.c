#ifndef lint
static const char RCSid[] = "$Id: sceneio.c,v 2.1 2003/03/14 21:27:46 greg Exp $";
#endif
/*
 *  Portable, binary Radiance i/o routines.
 *
 *  Called from octree and mesh i/o routines.
 */

#include "copyright.h"
#include "standard.h"
#include "octree.h"
#include "object.h"
#include "otypes.h"

static OBJECT  object0;			/* zeroeth object */
static short  otypmap[NUMOTYPE+32];	/* object type map */


static int
getobj(fp, objsiz)			/* get next object */
FILE	*fp;
int	objsiz;
{
	char  sbuf[MAXSTR];
	int  obj;
	register int  i;
	register long  m;
	register OBJREC	 *objp;

	i = getint(1, fp);
	if (i == -1)
		return(OVOID);		/* terminator */
	if ((obj = newobject()) == OVOID)
		error(SYSTEM, "out of object space");
	objp = objptr(obj);
	if ((objp->otype = otypmap[i]) < 0)
		error(USER, "reference to unknown type");
	if ((m = getint(objsiz, fp)) != OVOID) {
		m += object0;
		if ((OBJECT)m != m)
			error(INTERNAL, "too many objects in getobj");
	}
	objp->omod = m;
	objp->oname = savqstr(getstr(sbuf, fp));
	if ((objp->oargs.nsargs = getint(2, fp)) > 0) {
		objp->oargs.sarg = (char **)malloc
				(objp->oargs.nsargs*sizeof(char *));
		if (objp->oargs.sarg == NULL)
			goto memerr;
		for (i = 0; i < objp->oargs.nsargs; i++)
			objp->oargs.sarg[i] = savestr(getstr(sbuf, fp));
	} else
		objp->oargs.sarg = NULL;
#ifdef	IARGS
	if ((objp->oargs.niargs = getint(2, fp)) > 0) {
		objp->oargs.iarg = (long *)malloc
				(objp->oargs.niargs*sizeof(long));
		if (objp->oargs.iarg == NULL)
			goto memerr;
		for (i = 0; i < objp->oargs.niargs; i++)
			objp->oargs.iarg[i] = getint(4, fp);
	} else
		objp->oargs.iarg = NULL;
#endif
	if ((objp->oargs.nfargs = getint(2, fp)) > 0) {
		objp->oargs.farg = (FLOAT *)malloc
				(objp->oargs.nfargs*sizeof(FLOAT));
		if (objp->oargs.farg == NULL)
			goto memerr;
		for (i = 0; i < objp->oargs.nfargs; i++)
			objp->oargs.farg[i] = getflt(fp);
	} else
		objp->oargs.farg = NULL;
	if (feof(fp))
		error(SYSTEM, "unexpected EOF in getobj");
						/* initialize */
	objp->os = NULL;
						/* insert */
	insertobject(obj);
	return(obj);
memerr:
	error(SYSTEM, "out of memory in getobj");
}


void
readscene(fp, objsiz)			/* read binary scene description */
FILE	*fp;
int	objsiz;
{
	char  sbuf[32];
	int  i;
					/* record starting object */
	object0 = nobjects;
					/* read type map */
	for (i = 0; getstr(sbuf, fp) != NULL && sbuf[0]; i++)
		if ((otypmap[i] = otype(sbuf)) < 0) {
			sprintf(errmsg, "unknown object type \"%s\"",
					sbuf);
			error(WARNING, errmsg);
		}
					/* read objects */
	while (getobj(fp, objsiz) != OVOID)
		;
}


static void
putobj(o, fp)			/* write out object */
FILE	*fp;
register OBJREC  *o;
{
	register int  i;

	if (o == NULL) {		/* terminator */
		putint(-1L, 1, fp);
		return;
	}
	putint((long)o->otype, 1, fp);
	putint((long)o->omod, sizeof(OBJECT), fp);
	putstr(o->oname, fp);
	putint((long)o->oargs.nsargs, 2, fp);
	for (i = 0; i < o->oargs.nsargs; i++)
		putstr(o->oargs.sarg[i], fp);
#ifdef  IARGS
	putint((long)o->oargs.niargs, 2, fp);
	for (i = 0; i < o->oargs.niargs; i++)
		putint((long)o->oargs.iarg[i], 4, fp);
#endif
	putint((long)o->oargs.nfargs, 2, fp);
	for (i = 0; i < o->oargs.nfargs; i++)
		putflt(o->oargs.farg[i], fp);
}


void
writescene(firstobj, nobjs, fp)		/* write binary scene description */
OBJECT	firstobj, nobjs;
FILE	*fp;
{
	int	i;
					/* write out type list */
	for (i = 0; i < NUMOTYPE; i++)
		putstr(ofun[i].funame, fp);
	putstr("", fp);
					/* write objects */
	for (i = firstobj; i < firstobj+nobjs; i++)
		putobj(objptr(i), fp);
	putobj(NULL, fp);		/* terminator */
	if (ferror(fp))
		error(SYSTEM, "write error in writescene");
}
