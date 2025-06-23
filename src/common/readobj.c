#ifndef lint
static const char RCSid[] = "$Id: readobj.c,v 2.30 2025/06/23 19:56:47 greg Exp $";
#endif
/*
 *  readobj.c - routines for reading in object descriptions.
 *
 *  External symbols declared in object.h
 */

#include "copyright.h"

#include  <ctype.h>
#include  <string.h>
#include  <stdio.h>

#include  "platform.h"
#include  "standard.h"
#include  "object.h"
#include  "otypes.h"

#ifndef OBJMEMOPT
#define OBJMEMOPT	1		/* optimize object block memory? */
#endif

OBJREC  *objblock[MAXOBJBLK];		/* our objects */
OBJECT  nobjects = 0;			/* # of objects */


void
readobj(				/* read in an object file or stream */
	char  *inpspec
)
{
	OBJECT  lastobj;
	FILE  *infp;
	char  buf[2048];
	int  c;

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
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(infp);
#endif
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
	if (inpspec[0] == '!') {
		if (pclose(infp) != 0) {
			sprintf(errmsg, "bad status from \"%s\"", inpspec);
			error(WARNING, errmsg);
		}
	} else if (infp != stdin)
		fclose(infp);
#ifdef getc_unlocked
	else
		funlockfile(infp);
#endif
	if (nobjects == lastobj) {
		sprintf(errmsg, "(%s): empty file", inpspec);
		error(WARNING, errmsg);
	}
}


void
getobject(				/* read the next object */
	char  *name,
	FILE  *fp
)
{
#define	OALIAS	-2
	OBJECT  obj;
	char  sbuf[MAXSTR];
	int  rval;
	OBJREC  *objp;

	if ((obj = newobject()) == OVOID)
		error(SYSTEM, "out of object space");
	objp = objptr(obj);
					/* get modifier */
	strcpy(sbuf, "EOF");
	fgetword(sbuf, MAXSTR, fp);
	if (strchr(sbuf, '\t')) {
		sprintf(errmsg, "(%s): illegal tab in modifier \"%s\"",
					name, sbuf);
		error(USER, errmsg);
	}
	if (!strcmp(sbuf, VOIDID))
		objp->omod = OVOID;
	else if (!strcmp(sbuf, ALIASMOD))
		objp->omod = OALIAS;
	else if ((objp->omod = modifier(sbuf)) == OVOID) {
		sprintf(errmsg, "(%s): undefined modifier \"%s\"", name, sbuf);
		error(USER, errmsg);
	}
					/* get type */
	strcpy(sbuf, "EOF");
	fgetword(sbuf, MAXSTR, fp);
	if ((objp->otype = otype(sbuf)) < 0) {
		sprintf(errmsg, "(%s): unknown type \"%s\"", name, sbuf);
		error(USER, errmsg);
	}
					/* get identifier */
	sbuf[0] = '\0';
	fgetword(sbuf, MAXSTR, fp);
	if (strchr(sbuf, '\t')) {
		sprintf(errmsg, "(%s): illegal tab in identifier \"%s\"",
					name, sbuf);
		error(USER, errmsg);
	}
	objp->oname = savqstr(sbuf);
					/* get arguments */
	if (objp->otype == MOD_ALIAS) {
		OBJECT  ref;
		OBJREC  *rfp;
		strcpy(sbuf, "EOF");
		fgetword(sbuf, MAXSTR, fp);
		if ((ref = modifier(sbuf)) == OVOID) {
			sprintf(errmsg, "(%s): bad reference \"%s\"",
					name, sbuf);
			objerror(objp, USER, errmsg);
		}			/* skip pass-thru aliases */
		while ((rfp=objptr(ref))->otype == MOD_ALIAS &&
				!rfp->oargs.nsargs & (rfp->omod != OVOID))
			ref = rfp->omod;

		if ((objp->omod == OALIAS) | (objp->omod == rfp->omod)) {
			objp->omod = ref;
		} else {
			objp->oargs.sarg = (char **)malloc(sizeof(char *));
			if (objp->oargs.sarg == NULL)
				error(SYSTEM, "out of memory in getobject");
			objp->oargs.nsargs = 1;
			objp->oargs.sarg[0] = savestr(sbuf);
		}
	} else if ((rval = readfargs(&objp->oargs, fp)) == 0) {
		sprintf(errmsg, "(%s): bad arguments", name);
		objerror(objp, USER, errmsg);
	} else if (rval < 0) {
		sprintf(errmsg, "(%s): error reading scene", name);
		error(SYSTEM, errmsg);
	}
	if (objp->omod == OALIAS) {
		sprintf(errmsg, "(%s): inappropriate use of '%s' modifier",
				name, ALIASMOD);
		objerror(objp, USER, errmsg);
	}
					/* initialize */
	objp->os = NULL;

	insertobject(obj);		/* add to global structure */
#undef OALIAS
}


static void
optimize_objblock(int i)		/* consolidate memory in object block */
{
#if OBJMEMOPT
	OBJREC		*o, *co;
	int		n = 0;
	unsigned long	sargcnt = 0, iargcnt = 0, fargcnt = 0, namecnt = 0;

	if (i < 0 || objblock[i] == NULL || objblock[i][OBJBLKSIZ].otype < 0)
		return;			/* invalid or already flagged */

	for (o = objblock[i]+OBJBLKSIZ; o-- > objblock[i]; ) {
		if (o->oname == NULL)	/* too early to optimize? */
			return;
		if (o->os != NULL)	/* too late to optimize? */
			return;
		n += (o->oargs.nsargs > 0) | (o->oargs.nfargs > 0);
		sargcnt += o->oargs.nsargs;
		fargcnt += o->oargs.nfargs;
#ifdef  IARGS
		iargcnt += o->oargs.niargs;
#endif
		namecnt += strlen(o->oname)+1;
	}
	if (n < OBJBLKSIZ/10)	/* never happens? */
		return;
					/* prep consolidation object */
	co = objblock[i]+OBJBLKSIZ;
	co->oargs.nsargs = sargcnt;
	co->oargs.nfargs = fargcnt;
	if ((co->oargs.nsargs != sargcnt) | (co->oargs.nfargs != fargcnt))
		return;			/* overrun condition */

	co->oname = (char *)malloc(sizeof(char)*namecnt);
	co->oargs.sarg = (char **)malloc(sizeof(char *)*sargcnt);
	co->oargs.farg = (RREAL *)malloc(sizeof(RREAL)*fargcnt);
	if ((co->oname == NULL) | (co->oargs.sarg == NULL) |
			(co->oargs.farg == NULL)) {
		free(co->oname);
		free(co->oargs.sarg); free(co->oargs.farg);
		return;			/* insufficient memory */
	}
#ifdef  IARGS
	co->oargs.niargs = iargcnt;
	co->oargs.iarg = (long *)malloc(sizeof(long)*iargcnt);
	if (co->oargs.iarg == NULL) {
		free(co->oname);
		free(co->oargs.sarg); free(co->oargs.farg);
		return;			/* insufficient memory */
	}
	iargcnt = 0;
#endif
	namecnt = sargcnt = fargcnt = 0;
	for (o = objblock[i]+OBJBLKSIZ; o-- > objblock[i]; ) {
		n = strlen(o->oname)+1;
		memcpy(co->oname + namecnt, o->oname, n);
		freeqstr(o->oname);
		o->oname = co->oname + namecnt;
		namecnt += n;
		if (o->oargs.nsargs) {
			memcpy(co->oargs.sarg+sargcnt, o->oargs.sarg,
					sizeof(char *)*o->oargs.nsargs);
			free(o->oargs.sarg);
			o->oargs.sarg = co->oargs.sarg + sargcnt;
			sargcnt += o->oargs.nsargs;
		}
		if (o->oargs.nfargs) {
			memcpy(co->oargs.farg+fargcnt, o->oargs.farg,
					sizeof(RREAL)*o->oargs.nfargs);
			free(o->oargs.farg);
			o->oargs.farg = co->oargs.farg + fargcnt;
			fargcnt += o->oargs.nfargs;
		}
#ifdef  IARGS
		if (o->oargs.niargs) {
			memcpy(co->oargs.iarg+iargcnt, o->oargs.iarg,
					sizeof(long)*o->oargs.niargs);
			free(o->oargs.iarg);
			o->oargs.iarg = co->oargs.iarg + iargcnt;
			iargcnt += o->oargs.niargs;
		}
#endif
	}
	co->otype = -1;		/* flag for optimized block */
#endif
}


OBJECT
newobject(void)				/* get a new object */
{
	int  i;

	if ((nobjects & (OBJBLKSIZ-1)) == 0) {	/* new block */
		i = nobjects >> OBJBLKSHFT;
		optimize_objblock(i-1);		/* optimize previous block */
		errno = 0;
		if (i >= MAXOBJBLK)
			return(OVOID);
		objblock[i] = (OBJREC *)calloc(OBJBLKSIZ+OBJMEMOPT,
						sizeof(OBJREC));
		if (objblock[i] == NULL)
			return(OVOID);
	}
	return(nobjects++);
}

void
freeobjects(				/* free a range of objects */
	int firstobj,
	int nobjs
)
{
	int  obj;
					/* check bounds */
	if (firstobj < 0)
		return;
	if (nobjs <= 0)
		return;
	if (firstobj + nobjs > nobjects)
		return;
					/* clear objects */
	for (obj = firstobj+nobjs; obj-- > firstobj; ) {
		OBJREC  *o = objptr(obj);
		free_os(o);		/* free client memory */
		if (!OBJMEMOPT || !objblock[obj>>OBJBLKSHFT][OBJBLKSIZ].otype) {
			freeqstr(o->oname);
			freefargs(&o->oargs);
		}
		memset(o, 0, sizeof(OBJREC));
	}
					/* free objects off end */
	for (obj = nobjects; obj-- > 0; )
		if (objptr(obj)->oname != NULL)
			break;
	if (++obj >= nobjects)
		return;
	while (nobjects > obj)		/* free empty end blocks */
		if ((--nobjects & (OBJBLKSIZ-1)) == 0) {
			int	i = nobjects >> OBJBLKSHFT;
					/* consolidated block? */
			if (OBJMEMOPT && objblock[i][OBJBLKSIZ].otype < 0) {
				free(objblock[i][OBJBLKSIZ].oname);
				freefargs(&objblock[i][OBJBLKSIZ].oargs);
			}
			free(objblock[i]);
			objblock[i] = NULL;
		}
	truncobjndx();			/* truncate modifier look-up */
}
