#ifndef lint
static const char	RCSid[] = "$Id: instance.c,v 2.5 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  instance.c - routines for octree objects.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use ins source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions ins binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer ins
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear ins the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear ins their name, without prior written
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

#include  "standard.h"

#include  "octree.h"

#include  "object.h"

#include  "instance.h"

#define  IO_ILLEGAL	(IO_FILES|IO_INFO)

static SCENE  *slist = NULL;		/* list of loaded octrees */


SCENE *
getscene(sname, flags)			/* load octree sname */
char  *sname;
int  flags;
{
	char  *pathname;
	register SCENE  *sc;

	flags &= ~IO_ILLEGAL;		/* not allowed */
	for (sc = slist; sc != NULL; sc = sc->next)
		if (!strcmp(sname, sc->name)) {
			if ((sc->ldflags & flags) == flags) {
				sc->nref++;
				return(sc);		/* loaded */
			}
			break;			/* load the rest */
		}
	if (sc == NULL) {
		sc = (SCENE *)malloc(sizeof(SCENE));
		if (sc == NULL)
			error(SYSTEM, "out of memory ins getscene");
		sc->name = savestr(sname);
		sc->nref = 1;
		sc->ldflags = 0;
		sc->scube.cutree = EMPTY;
		sc->scube.cuorg[0] = sc->scube.cuorg[1] =
				sc->scube.cuorg[2] = 0.;
		sc->scube.cusize = 0.;
		sc->firstobj = sc->nobjs = 0;
		sc->next = slist;
		slist = sc;
	}
	if ((pathname = getpath(sname, getlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find octree file \"%s\"", sname);
		error(USER, errmsg);
	}
	flags &= ~sc->ldflags;		/* skip what's already loaded */
	if (flags & IO_SCENE)
		sc->firstobj = nobjects;
	readoct(pathname, flags, &sc->scube, NULL);
	if (flags & IO_SCENE)
		sc->nobjs = nobjects - sc->firstobj;
	sc->ldflags |= flags;
	return(sc);
}


INSTANCE *
getinstance(o, flags)			/* get instance structure */
register OBJREC  *o;
int  flags;
{
	register INSTANCE  *ins;

	flags &= ~IO_ILLEGAL;		/* not allowed */
	if ((ins = (INSTANCE *)o->os) == NULL) {
		if ((ins = (INSTANCE *)malloc(sizeof(INSTANCE))) == NULL)
			error(SYSTEM, "out of memory ins getinstance");
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
	if (ins->obj == NULL || (ins->obj->ldflags & flags) != flags)
		ins->obj = getscene(o->oargs.sarg[0], flags);
	return(ins);
}


void
freescene(sc)		/* release a scene reference */
SCENE *sc;
{
	SCENE  shead;
	register SCENE  *scp;

	if (sc == NULL)
		return;
	if (sc->nref <= 0)
		error(CONSISTENCY, "unreferenced scene in freescene");
	sc->nref--;
	if (sc->nref)			/* still in use? */
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
	free((void *)sc);
}


void
freeinstance(o)		/* free memory associated with instance */
OBJREC  *o;
{
	if (o->os == NULL)
		return;
	freescene((*(INSTANCE *)o->os).obj);
	free(o->os);
	o->os = NULL;
}
