#ifndef lint
static const char	RCSid[] = "$Id: readobj.c,v 2.8 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  readobj.c - routines for reading in object descriptions.
 *
 *  External symbols declared in object.h
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

#include  "standard.h"

#include  "object.h"

#include  "otypes.h"

#include  <ctype.h>

OBJREC  *objblock[MAXOBJBLK];		/* our objects */
OBJECT  nobjects = 0;			/* # of objects */


void
readobj(inpspec)		/* read in an object file or stream */
char  *inpspec;
{
	OBJECT  lastobj;
	FILE  *infp;
	char  buf[1024];
	register int  c;

	lastobj = nobjects;
	if (inpspec == NULL) {
		infp = stdin;
		inpspec = "standard input";
	} else if (inpspec[0] == '!') {
		if ((infp = popen(inpspec+1, "r")) == NULL) {
			sprintf(errmsg, "cannot execute \"%s\"", inpspec);
			error(SYSTEM, errmsg);
		}
	} else if ((infp = fopen(inpspec, "r")) == NULL) {
		sprintf(errmsg, "cannot open scene file \"%s\"", inpspec);
		error(SYSTEM, errmsg);
	}
	while ((c = getc(infp)) != EOF) {
		if (isspace(c))
			continue;
		if (c == '#') {				/* comment */
			fgets(buf, sizeof(buf), infp);
		} else if (c == '!') {			/* command */
			ungetc(c, infp);
			fgetline(buf, sizeof(buf), infp);
			readobj(buf);
		} else {				/* object */
			ungetc(c, infp);
			getobject(inpspec, infp);
		}
	}
	if (inpspec[0] == '!')
		pclose(infp);
	else
		fclose(infp);
	if (nobjects == lastobj) {
		sprintf(errmsg, "(%s): empty file", inpspec);
		error(WARNING, errmsg);
	}
}


void
getobject(name, fp)			/* read the next object */
char  *name;
FILE  *fp;
{
	OBJECT  obj;
	char  sbuf[MAXSTR];
	int  rval;
	register OBJREC  *objp;

	if ((obj = newobject()) == OVOID)
		error(SYSTEM, "out of object space");
	objp = objptr(obj);
					/* get modifier */
	strcpy(sbuf, "EOF");
	fgetword(sbuf, MAXSTR, fp);
	if (!strcmp(sbuf, VOIDID))
		objp->omod = OVOID;
	else if ((objp->omod = modifier(sbuf)) == OVOID) {
		sprintf(errmsg, "(%s): undefined modifier \"%s\"", name, sbuf);
		error(USER, errmsg);
	}
					/* get type */
	strcpy(sbuf, "EOF");
	fgetword(sbuf, MAXSTR, fp);
	if (!strcmp(sbuf, ALIASID))
		objp->otype = -1;
	else if ((objp->otype = otype(sbuf)) < 0) {
		sprintf(errmsg, "(%s): unknown type \"%s\"", name, sbuf);
		error(USER, errmsg);
	}
					/* get identifier */
	sbuf[0] = '\0';
	fgetword(sbuf, MAXSTR, fp);
	objp->oname = savqstr(sbuf);
					/* get arguments */
	if (objp->otype == -1) {
		register OBJECT  alias;
		strcpy(sbuf, "EOF");
		fgetword(sbuf, MAXSTR, fp);
		if ((alias = modifier(sbuf)) == OVOID) {
			sprintf(errmsg,
			"(%s): bad reference \"%s\" for %s \"%s\"",
					name, sbuf, ALIASID, objp->oname);
			error(USER, errmsg);
		}
		objp->otype = objptr(alias)->otype;
		copystruct(&objp->oargs, &objptr(alias)->oargs);
	} else if ((rval = readfargs(&objp->oargs, fp)) == 0) {
		sprintf(errmsg, "(%s): bad arguments", name);
		objerror(objp, USER, errmsg);
	} else if (rval < 0) {
		sprintf(errmsg, "(%s): error reading scene", name);
		error(SYSTEM, errmsg);
	}
					/* initialize */
	objp->os = NULL;

	insertobject(obj);		/* add to global structure */
}


int
newobject()				/* get a new object */
{
	register int  i;

	if ((nobjects & (OBJBLKSIZ-1)) == 0) {	/* new block */
		errno = 0;
		i = nobjects >> OBJBLKSHFT;
		if (i >= MAXOBJBLK)
			return(OVOID);
		objblock[i] = (OBJREC *)calloc(OBJBLKSIZ, sizeof(OBJREC));
		if (objblock[i] == NULL)
			return(OVOID);
	}
	return(nobjects++);
}

void
freeobjects(firstobj, nobjs)		/* free a range of objects */
OBJECT  firstobj, nobjs;
{
	register int  obj;
					/* check bounds */
	if (firstobj < 0)
		return;
	if (nobjs <= 0)
		return;
	if (firstobj + nobjs > nobjects)
		return;
					/* clear objects */
	for (obj = firstobj+nobjs; obj-- > firstobj; ) {
		register OBJREC  *o = objptr(obj);
		free_os(o);		/* free client memory */
		freeqstr(o->oname);
		freefargs(&o->oargs);
		bzero(o, sizeof(OBJREC));
	}
	clearobjndx();
					/* free objects off end */
	for (obj = nobjects; obj-- > 0; )
		if (objptr(obj)->oname != NULL)
			break;
	++obj;
	while (nobjects > obj)		/* free empty end blocks */
		if ((--nobjects & (OBJBLKSIZ-1)) == 0) {
			int	i = nobjects >> OBJBLKSHFT;
			free((void *)objblock[i]);
			objblock[i] = NULL;
		}
}
