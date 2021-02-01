#ifndef lint
static const char RCSid[] = "$Id: modobject.c,v 2.20 2021/02/01 17:36:45 greg Exp $";
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
	int	i;
	long	j;

	for (i = (nobjects-1)>>OBJBLKSHFT; i >= 0; i--) {
		j = op - objblock[i];
		if ((j >= 0) & (j < OBJBLKSIZ))
			return(((OBJECT)i<<OBJBLKSHFT) + (OBJECT)j);
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


int
eqreal(				/* are two real values close enough to equal? */
	double	d1,
	double	d2
)
{
	if (d2 != 0.0)
		d1 = d1/d2 - 1.0;
	return((-FTINY <= d1) & (d1 <= FTINY));
}


int
eqobjects(			/* check if two objects are equal */
	OBJECT	obj1,
	OBJECT	obj2
)
{
	OBJREC	*op1, *op2;
	int	i, n;

	while (obj1 != obj2) {
		if (obj1 == OVOID)
			return(0);
		if (obj2 == OVOID)
			return(0);
		op1 = objptr(obj1);
		op2 = objptr(obj2);
		if (op1->otype != op2->otype)
			return(0);
		if (op1->oargs.nsargs != op2->oargs.nsargs)
			return(0);
		if (op1->oargs.nfargs != op2->oargs.nfargs)
			return(0);
#ifdef IARGS
		if (op1->oargs.niargs != op2->oargs.niargs)
			return(0);
		for (i = op1->oargs.niargs; i-- > 0; )
			if (op1->oargs.iarg[i] != op2->oargs.iarg[i])
				return(0);
#endif
		for (i = op1->oargs.nfargs; i-- > 0; )
			if (!eqreal(op1->oargs.farg[i], op2->oargs.farg[i]))
				return(0);
		n = 0;
		switch (op1->otype) {	/* special cases (KEEP CONSISTENT!) */
		case MOD_ALIAS:
		case MAT_ILLUM:
		case MAT_MIRROR:
			n = (op1->oargs.nsargs > 0);
			break;
		case MIX_FUNC:
		case MIX_DATA:
		case MIX_TEXT:
		case MIX_PICT:
			n = 2*(op1->oargs.nsargs >= 2);
			break;
		case MAT_CLIP:
			n = op1->oargs.nsargs;
			break;
		}
					/* check other string arguments */
		for (i = op1->oargs.nsargs; i-- > n; )
			if (strcmp(op1->oargs.sarg[i], op2->oargs.sarg[i]))
				return(0);
		while (n-- > 0)		/* check modifier references */
			if (!eqobjects( lastmod(obj1, op1->oargs.sarg[n]),
					lastmod(obj2, op2->oargs.sarg[n]) ))
				return(0);
		obj1 = op1->omod;
		obj2 = op2->omod;
	}
	return(1);
}


void
insertobject(			/* insert new object into our list */
	OBJECT  obj
)
{
	int  i;

	if (ismodifier(objptr(obj)->otype)) {
		i = otndx(objptr(obj)->oname, &modtab);
		if (eqobjects(obj, modtab.htab[i]))
			return;	/* don't index if same as earlier def. */
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
truncobjndx(void)		/* remove bogus table entries past end */
{
	int	ndx;

	if (nobjects <= 0) {
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
		return;
	}
	for (ndx = modtab.hsiz*(modtab.htab != NULL); ndx--; )
		if (modtab.htab[ndx] >= nobjects)
			modtab.htab[ndx] = OVOID;

	for (ndx = objtab.hsiz*(objtab.htab != NULL); ndx--; )
		if (objtab.htab[ndx] >= nobjects)
			objtab.htab[ndx] = OVOID;
}


static int
nexthsiz(			/* return next hash table size */
	int  oldsiz
)
{
	static int  hsiztab[] = {
		251, 509, 1021, 2039, 4093, 8191, 16381,
		32749, 65521, 131071, 262139, 0
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
	char	*onm;
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
		if (tab->htab[ndx] == OVOID)
			return(ndx);
		onm = objptr(tab->htab[ndx])->oname;
		if (onm != NULL && !strcmp(onm, name))
			return(ndx);
	}
					/* table is full, reallocate */
	oldhtab = tab->htab;
	ndx = tab->hsiz;
	tab->htab = NULL;
	while (ndx--)
		if (oldhtab[ndx] != OVOID) {
			onm = objptr(oldhtab[ndx])->oname;
			if (onm == NULL)
				continue;
			i = otndx(onm, tab);
			tab->htab[i] = oldhtab[ndx];
		}
	free((void *)oldhtab);
	goto tryagain;			/* should happen only once! */
}
