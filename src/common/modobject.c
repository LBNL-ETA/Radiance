#ifndef lint
static const char	RCSid[] = "$Id: modobject.c,v 2.8 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 *  Routines for tracking object modifiers
 *
 *  External symbols declared in object.h
 */

#include "copyright.h"

#include  "standard.h"

#include  "object.h"

#include  "otypes.h"


static struct ohtab {
	int  hsiz;			/* current table size */
	OBJECT  *htab;			/* table, if allocated */
}  modtab = {100, NULL}, objtab = {1000, NULL};	/* modifiers and objects */

static int  otndx();


int
objndx(op)			/* get object number from pointer */
register OBJREC  *op;
{
	register int  i, j;

	for (i = nobjects>>OBJBLKSHFT; i >= 0; i--) {
		j = op - objblock[i];
		if (j >= 0 && j < OBJBLKSIZ)
			return((i<<OBJBLKSHFT) + j);
	}
	return(OVOID);
}


int
lastmod(obj, mname)		/* find modifier definition before obj */
OBJECT  obj;
char  *mname;
{
	register OBJREC  *op;
	register int  i;

	i = modifier(mname);		/* try hash table first */
	if (i < obj)
		return(i);
	for (i = obj; i-- > 0; ) {	/* need to search */
		op = objptr(i);
		if (ismodifier(op->otype) && !strcmp(op->oname, mname))
			return(i);
	}
	return(OVOID);
}


int
modifier(mname)			/* get a modifier number from its name */
char  *mname;
{
	register int  ndx;

	ndx = otndx(mname, &modtab);
	return(modtab.htab[ndx]);
}


#ifdef  GETOBJ
int
object(oname)			/* get an object number from its name */
char  *oname;
{
	register int  ndx;

	ndx = otndx(oname, &objtab);
	return(objtab.htab[ndx]);
}
#endif


void
insertobject(obj)		/* insert new object into our list */
register OBJECT  obj;
{
	register int  i;

	if (ismodifier(objptr(obj)->otype)) {
		i = otndx(objptr(obj)->oname, &modtab);
		modtab.htab[i] = obj;
	}
#ifdef  GETOBJ
	else {
		i = otndx(objptr(obj)->oname, &objtab);
		objtab.htab[i] = obj;
	}
#endif
	for (i = 0; addobjnotify[i] != NULL; i++)
		(*addobjnotify[i])(obj);
}


void
clearobjndx()			/* clear object hash tables */
{
	if (modtab.htab != NULL) {
		free((void *)modtab.htab);
		modtab.htab = NULL;
		modtab.hsiz = 100;
	}
	if (objtab.htab != NULL) {
		free((void *)objtab.htab);
		objtab.htab = NULL;
		objtab.hsiz = 100;
	}
}


static int
nexthsiz(oldsiz)		/* return next hash table size */
int  oldsiz;
{
	static int  hsiztab[] = {
		251, 509, 1021, 2039, 4093, 8191, 16381, 0
	};
	register int  *hsp;

	for (hsp = hsiztab; *hsp; hsp++)
		if (*hsp > oldsiz)
			return(*hsp);
	return(oldsiz*2 + 1);		/* not always prime */
}


static int
otndx(name, tab)		/* get object table index for name */
char  *name;
register struct ohtab  *tab;
{
	OBJECT  *oldhtab;
	int  hval, i;
	register int  ndx;

	if (tab->htab == NULL) {		/* new table */
		tab->hsiz = nexthsiz(tab->hsiz);
		tab->htab = (OBJECT *)malloc(tab->hsiz*sizeof(OBJECT));
		if (tab->htab == NULL)
			error(SYSTEM, "out of memory in otndx");
		ndx = tab->hsiz;
		while (ndx--)			/* empty it */
			tab->htab[ndx] = OVOID;
	}
					/* look up object */
	hval = shash(name);
tryagain:
	for (i = 0; i < tab->hsiz; i++) {
		ndx = (hval + i*i) % tab->hsiz;
		if (tab->htab[ndx] == OVOID ||
				!strcmp(objptr(tab->htab[ndx])->oname, name))
			return(ndx);
	}
					/* table is full, reallocate */
	oldhtab = tab->htab;
	ndx = tab->hsiz;
	tab->htab = NULL;
	while (ndx--)
		if (oldhtab[ndx] != OVOID) {
			i = otndx(objptr(oldhtab[ndx])->oname, tab);
			tab->htab[i] = oldhtab[ndx];
		}
	free((void *)oldhtab);
	goto tryagain;			/* should happen only once! */
}
