#ifndef lint
static const char	RCSid[] = "$Id: octree.c,v 2.8 2009/11/01 04:41:55 greg Exp $";
#endif
/*
 *  octree.c - routines dealing with octrees and cubes.
 */

#include "copyright.h"

#include  "standard.h"

#include  "octree.h"

OCTREE  *octblock[MAXOBLK];		/* our octree */
static OCTREE  ofreelist = EMPTY;	/* freed octree nodes */
static OCTREE  treetop = 0;		/* next free node */


OCTREE
octalloc()			/* allocate an octree */
{
	OCTREE  freet;

	if ((freet = ofreelist) != EMPTY) {
		ofreelist = octkid(freet, 0);
		return(freet);
	}
	freet = treetop;
	if (octti(freet) == 0) {
		errno = 0;
		if (octbi(freet) >= MAXOBLK)
			return(EMPTY);
		if ((octblock[octbi(freet)] = (OCTREE *)malloc(
				(unsigned)OCTBLKSIZ*8*sizeof(OCTREE))) == NULL)
			return(EMPTY);
	}
	treetop++;
	return(freet);
}


void
octfree(ot)			/* free an octree */
OCTREE  ot;
{
	int  i;

	if (!istree(ot))
		return;
	for (i = 0; i < 8; i++)
		octfree(octkid(ot, i));
	octkid(ot, 0) = ofreelist;
	ofreelist = ot;
}


void
octdone()			/* free EVERYTHING */
{
	int	i;

	for (i = 0; i < MAXOBLK; i++) {
		if (octblock[i] == NULL)
			break;
		free((void *)octblock[i]);
		octblock[i] = NULL;
	}
	ofreelist = EMPTY;
	treetop = 0;
}


OCTREE
combine(ot)			/* recursively combine nodes */
OCTREE  ot;
{
	int  i;
	OCTREE  ores;

	if (!istree(ot))	/* not a tree */
		return(ot);
	ores = octkid(ot, 0) = combine(octkid(ot, 0));
	for (i = 1; i < 8; i++)
		if ((octkid(ot, i) = combine(octkid(ot, i))) != ores)
			ores = ot;
	if (!istree(ores)) {	/* all were identical leaves */
		octkid(ot, 0) = ofreelist;
		ofreelist = ot;
	}
	return(ores);
}


void
culocate(cu, pt)		/* locate point within cube */
CUBE  *cu;
FVECT  pt;
{
	int  i;
	int  branch;

	while (istree(cu->cutree)) {
		cu->cusize *= 0.5;
		branch = 0;
		for (i = 0; i < 3; i++)
			if (cu->cuorg[i] + cu->cusize <= pt[i]) {
				cu->cuorg[i] += cu->cusize;
				branch |= 1 << i;
			}
		cu->cutree = octkid(cu->cutree, branch);
	}
}


int
incube(cu, pt)			/* determine if a point is inside a cube */
CUBE  *cu;
FVECT  pt;
{
	if (cu->cuorg[0] > pt[0] || pt[0] >= cu->cuorg[0] + cu->cusize)
		return(0);
	if (cu->cuorg[1] > pt[1] || pt[1] >= cu->cuorg[1] + cu->cusize)
		return(0);
	if (cu->cuorg[2] > pt[2] || pt[2] >= cu->cuorg[2] + cu->cusize)
		return(0);
	return(1);
}
