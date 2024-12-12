#ifndef lint
static const char RCSid[] = "$Id: instance.c,v 2.14 2024/12/12 20:04:46 greg Exp $";
#endif
/*
 *  instance.c - routines for octree objects.
 */

#include "copyright.h"

#include  "rtmath.h"
#include  "rterror.h"
#include  "rtio.h"
#include  "paths.h"

#include  "octree.h"
#include  "object.h"
#include  "instance.h"

#define  IO_ILLEGAL	(IO_FILES|IO_INFO)

static SCENE  *slist = NULL;		/* list of loaded octrees */


SCENE *
getscene(				/* get new octree reference */
	char  *sname,
	int  flags
)
{
	char  *pathname;
	SCENE  *sc;

	flags &= ~IO_ILLEGAL;		/* not allowed */
	for (sc = slist; sc != NULL; sc = sc->next)
		if (!strcmp(sname, sc->name))
			break;
	if (sc == NULL) {		/* new instance? */
		sc = (SCENE *)calloc(1, sizeof(SCENE));
		if (sc == NULL)
			error(SYSTEM, "out of memory in getscene");
		sc->name = savestr(sname);
		sc->scube.cutree = EMPTY;
		sc->next = slist;
		slist = sc;
	}
	sc->nref++;			/* bump reference count */
	if (!(flags &= ~sc->ldflags))	/* nothing to load? */
		return(sc);
	if ((pathname = getpath(sname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find octree file \"%s\"", sname);
		error(SYSTEM, errmsg);
	}
	if (flags & IO_SCENE)
		sc->firstobj = nobjects;
	if (flags)
		readoct(pathname, flags, &sc->scube, NULL);
	if (flags & IO_SCENE)
		sc->nobjs = nobjects - sc->firstobj;
	sc->ldflags |= flags;
	return(sc);
}


INSTANCE *
getinstance(				/* get instance structure */
	OBJREC  *o,
	int  flags
)
{
	INSTANCE  *ins;

	flags &= ~IO_ILLEGAL;		/* not allowed */
	if ((ins = (INSTANCE *)o->os) == NULL) {
		if ((ins = (INSTANCE *)malloc(sizeof(INSTANCE))) == NULL)
			error(SYSTEM, "out of memory in getinstance");
		if (o->oargs.nsargs < 1)
			objerror(o, USER, "bad # of arguments");
		if (fullxf(&ins->x, o->oargs.nsargs-1,
				o->oargs.sarg+1) != o->oargs.nsargs-1)
			objerror(o, USER, "bad transform");
		if (ins->x.f.sca < 0.0) {
			ins->x.f.sca = -ins->x.f.sca;
			ins->x.b.sca = -ins->x.b.sca;
		}
		ins->obj = NULL;
		o->os = (char *)ins;
	}
	if (ins->obj == NULL) {
		ins->obj = getscene(o->oargs.sarg[0], flags);
	} else if ((flags &= ~ins->obj->ldflags)) {
		if (flags & IO_SCENE)
			ins->obj->firstobj = nobjects;
		if (flags)
			readoct(getpath(o->oargs.sarg[0], getrlibpath(), R_OK),
					flags, &ins->obj->scube, NULL);
		if (flags & IO_SCENE)
			ins->obj->nobjs = nobjects - ins->obj->firstobj;
		ins->obj->ldflags |= flags;
	}
	return(ins);
}


void
freescene(		/* release a scene reference */
	SCENE *sc
)
{
	SCENE  shead;
	SCENE  *scp;

	if (sc == NULL)
		return;
	if (sc->nref <= 0)
		error(CONSISTENCY, "unreferenced scene in freescene");
	if (--sc->nref)			/* still in use? */
		return;
	shead.next = slist;		/* else remove from our list */
	for (scp = &shead; scp->next != NULL; scp = scp->next)
		if (scp->next == sc) {
			scp->next = sc->next;
			sc->next = NULL;
			break;
		}
	if (sc->next != NULL)		/* can't be in list anymore */
		error(CONSISTENCY, "unlisted scene in freescene");
	slist = shead.next;
	freestr(sc->name);		/* free memory */
	octfree(sc->scube.cutree);
	freeobjects(sc->firstobj, sc->nobjs);
	free(sc);
}


void
freeinstance(		/* free memory associated with instance */
	OBJREC  *o
)
{
	if (o->os == NULL)
		return;
	freescene((*(INSTANCE *)o->os).obj);
	free(o->os);
	o->os = NULL;
}
