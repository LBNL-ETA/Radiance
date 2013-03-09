#ifndef lint
static const char RCSid[] = "$Id: modobject.c,v 2.13 2013/03/09 18:53:05 greg Exp $";
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

static int  otndx(char *, struct ohtab *);


OBJECT
objndx(				/* get object number from pointer */
	OBJREC  *op
)
{
	int  i, j;

	for (i = nobjects>>OBJBLKSHFT; i >= 0; i--) {
		j = op - objblock[i];
		if ((j >= 0) & (j < OBJBLKSIZ))
			return((i<<OBJBLKSHFT) + j);
	}
	return(OVOID);
}


OBJECT
lastmod(			/* find modifier definition before obj */
	OBJECT  obj,
	char  *mname
)
{
	OBJREC  *op;
	int  i;

	i = modifier(mname);		/* try hash table first */
	if ((obj == OVOID) | (i < obj))
		return(i);
	for (i = obj; i-- > 0; ) {	/* need to search */
		op = objptr(i);
		if (ismodifier(op->otype) && op->oname[0] == mname[0] &&
					!strcmp(op->oname, mname))
			return(i);
	}
	return(OVOID);
}


OBJECT
modifier(			/* get a modifier number from its name */
	char  *mname
)
{
	int  ndx;

	ndx = otndx(mname, &modtab);
	return(modtab.htab[ndx]);
}


#ifdef  GETOBJ
OBJECT
object(				/* get an object number from its name */
	char  *oname
)
{
	int  ndx;

	ndx = otndx(oname, &objtab);
	return(objtab.htab[ndx]);
}
#endif


void
insertobject(			/* insert new object into our list */
	OBJECT  obj
)
{
	int  i;

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
clearobjndx(void)		/* clear object hash tables */
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
nexthsiz(			/* return next hash table size */
	int  oldsiz
)
{
	static int  hsiztab[] = {
		251, 509, 1021, 2039, 4093, 8191, 16381, 0
	};
	int  *hsp;

	for (hsp = hsiztab; *hsp; hsp++)
		if (*hsp > oldsiz)
			return(*hsp);
	return(oldsiz*2 + 1);		/* not always prime */
}


static int
otndx(				/* get object table index for name */
	char  *name,
	struct ohtab  *tab
)
{
	OBJECT  *oldhtab;
	int  hval, i;
	int  ndx;

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
		ndx = (hval + (unsigned long)i*i) % tab->hsiz;
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
