#ifndef lint
static const char	RCSid[] = "$Id: octree.c,v 2.5 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  octree.c - routines dealing with octrees and cubes.
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

#include  "octree.h"

OCTREE  *octblock[MAXOBLK];		/* our octree */
static OCTREE  ofreelist = EMPTY;	/* freed octree nodes */
static OCTREE  treetop = 0;		/* next free node */


OCTREE
octalloc()			/* allocate an octree */
{
	register OCTREE  freet;

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
	treetop += 8;
	return(freet);
}


void
octfree(ot)			/* free an octree */
register OCTREE  ot;
{
	register int  i;

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
	register int	i;

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
register OCTREE  ot;
{
	register int  i;
	register OCTREE  ores;

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
register CUBE  *cu;
register FVECT  pt;
{
	register int  i;
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


void
cucopy(cu1, cu2)		/* copy cu2 into cu1 */
register CUBE  *cu1, *cu2;
{
	cu1->cutree = cu2->cutree;
	cu1->cusize = cu2->cusize;
	VCOPY(cu1->cuorg, cu2->cuorg);
}


int
incube(cu, pt)			/* determine if a point is inside a cube */
register CUBE  *cu;
register FVECT  pt;
{
	if (cu->cuorg[0] > pt[0] || pt[0] >= cu->cuorg[0] + cu->cusize)
		return(0);
	if (cu->cuorg[1] > pt[1] || pt[1] >= cu->cuorg[1] + cu->cusize)
		return(0);
	if (cu->cuorg[2] > pt[2] || pt[2] >= cu->cuorg[2] + cu->cusize)
		return(0);
	return(1);
}
