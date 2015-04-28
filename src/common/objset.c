#ifndef lint
static const char	RCSid[] = "$Id: objset.c,v 2.17 2015/04/28 16:56:10 greg Exp $";
#endif
/*
 *  objset.c - routines for maintaining object sets.
 *
 *  External symbols declared in object.h
 */

#include "copyright.h"

#include  "standard.h"

#include  "octree.h"

#include  "object.h"

#ifndef  OSTSIZ
#ifdef  SMLMEM
#define  OSTSIZ		32749		/* object table size (a prime!) */
#else
#define  OSTSIZ		262139		/* object table size (a prime!) */
#endif
#endif

#undef xtra_long
#ifdef _WIN64
typedef unsigned long long int	xtra_long;
#else
typedef unsigned long int	xtra_long;
#endif

static OBJECT  *ostable[OSTSIZ];	/* the object set table */


void
insertelem(			/* insert obj into os, no questions */
	OBJECT  *os,
	OBJECT  obj
)
{
	int  i;
	
	for (i = os[0]++; i > 0; i--)
		if (os[i] > obj)
			os[i+1] = os[i];
		else
			break;
	os[i+1] = obj;
}


void
deletelem(			/* delete obj from os, no questions */
	OBJECT  *os,
	OBJECT  obj
)
{
	int  i;

	i = (*os)--;
	os++;
	while (i > 0 && *os < obj) {
		i--;
		os++;
	}
	while (--i > 0) {
		os[0] = os[1];
		os++;
	}
}


int
inset(				/* determine if object is in set */
	OBJECT  *os,
	OBJECT  obj
)
{
	int  upper, lower;
	int  cm, i;

	if ((i = os[0]) <= 12) {	/* linear search algorithm */
		cm = obj;
		while (i-- > 0)
			if (*++os == cm)
				return(1);
		return(0);
	}
	lower = 1;
	upper = cm = i + 1;
					/* binary search algorithm */
	while ((i = (lower + upper) >> 1) != cm) {
		cm = obj - os[i];
		if (cm > 0)
			lower = i;
		else if (cm < 0)
			upper = i;
		else
			return(1);
		cm = i;
	}
	return(0);
}


int
setequal(			/* determine if two sets are equal */
	OBJECT  *os1,
	OBJECT	*os2
)
{
	int  i;

	for (i = *os1; i-- >= 0; )
		if (*os1++ != *os2++)
			return(0);
	return(1);
}


void
setcopy(			/* copy object set os2 into os1 */
	OBJECT  *os1,
	OBJECT	*os2
)
{
	int  i;

	for (i = *os2; i-- >= 0; )
		*os1++ = *os2++;
}


OBJECT *
setsave(			/* allocate space and save set */
	OBJECT  *os
)
{
	OBJECT  *osnew;
	OBJECT  *oset;
	int  i;

	if ((osnew = oset = (OBJECT *)malloc((*os+1)*sizeof(OBJECT))) == NULL)
		error(SYSTEM, "out of memory in setsave\n");
	for (i = *os; i-- >= 0; )	/* inline setcopy */
		*oset++ = *os++;
	return(osnew);
}


void
setunion(			/* osr = os1 Union os2 */
	OBJECT  *osr,
	OBJECT  *os1,
	OBJECT	*os2
)
{
	int	i1, i2;

	osr[0] = 0;
	for (i1 = i2 = 1; i1 <= os1[0] || i2 <= os2[0]; ) {
		while (i1 <= os1[0] && (i2 > os2[0] || os1[i1] <= os2[i2])) {
			osr[++osr[0]] = os1[i1];
			if (i2 <= os2[0] && os2[i2] == os1[i1])
				i2++;
			i1++;
		}
		while (i2 <= os2[0] && (i1 > os1[0] || os2[i2] < os1[i1]))
			osr[++osr[0]] = os2[i2++];
	}
}


void
setintersect(			/* osr = os1 Intersect os2 */
	OBJECT  *osr,
	OBJECT  *os1,
	OBJECT	*os2
)
{
	int	i1, i2;

	osr[0] = 0;
	if (os1[0] <= 0)
		return;
	for (i1 = i2 = 1; i2 <= os2[0]; ) {
		while (os1[i1] < os2[i2])
			if (++i1 > os1[0])
				return;
		while (os2[i2] < os1[i1])
			if (++i2 > os2[0])
				return;
		if (os1[i1] == os2[i2])
			osr[++osr[0]] = os2[i2++];
	}
}


OCTREE
fullnode(			/* return octree for object set */
	OBJECT  *oset
)
{
	unsigned int  ntries;
	xtra_long  hval;
	OCTREE  ot;
	int  osentry, i;
	OBJECT  *os;
					/* hash on set */
	hval = 0;
	os = oset;
	i = *os++;
	while (i-- > 0)
		hval += *os++;
	ntries = 0;
tryagain:
	osentry = (hval + (xtra_long)ntries*ntries) % OSTSIZ;
	os = ostable[osentry];
	if (os == NULL) {
		os = ostable[osentry] = (OBJECT *)malloc(
				(unsigned)(oset[0]+2)*sizeof(OBJECT));
		if (os == NULL)
			goto memerr;
		ot = oseti(osentry);
	} else {
						/* look for set */
		for (i = 0; *os > 0; i++, os += *os + 1)
			if (setequal(os, oset))
				break;
		ot = oseti(i*OSTSIZ + osentry);
		if (*os > 0)			/* found it */
			return(ot);
		if (!isfull(ot)) {		/* entry overflow */
			if (++ntries < OSTSIZ)
				goto tryagain;
			else
				error(INTERNAL, "hash table overflow in fullnode");
		}
						/* remember position */
		i = os - ostable[osentry];
		os = ostable[osentry] = (OBJECT *)realloc(
				(void *)ostable[osentry],
				(unsigned)(i+oset[0]+2)*sizeof(OBJECT));
		if (os == NULL)
			goto memerr;
		os += i;			/* last entry */
	}
	setcopy(os, oset);		/* add new set */
	os += *os + 1;
	*os = 0;			/* terminator */
	return(ot);
memerr:
	error(SYSTEM, "out of memory in fullnode");
	return 0; /* pro forma return */
}


void
objset(			/* get object set for full node */
	OBJECT  *oset,
	OCTREE  ot
)
{
	OBJECT  *os;
	int  i;

	if (!isfull(ot))
		goto noderr;
	ot = oseti(ot);
	if ((os = ostable[ot%OSTSIZ]) == NULL)
		goto noderr;
	for (i = ot/OSTSIZ; i--; os += *os + 1)
		if (*os <= 0)
			goto noderr;
	for (i = *os; i-- >= 0; )		/* copy set here */
		*oset++ = *os++;
	return;
noderr:
	error(CONSISTENCY, "bad node in objset");
}


void
donesets(void)			/* free ALL SETS in our table */
{
	int  n;

	for (n = 0; n < OSTSIZ; n++) 
		if (ostable[n] != NULL) {
			free((void *)ostable[n]);
			ostable[n] = NULL;
		}
}
