#ifndef lint
static const char	RCSid[] = "$Id: objset.c,v 2.8 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  objset.c - routines for maintaining object sets.
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

#include  "octree.h"

#include  "object.h"

#ifndef  OSTSIZ
#ifdef  BIGMEM
#define  OSTSIZ		262139		/* object table size (a prime!) */
#else
#define  OSTSIZ		32749		/* object table size (a prime!) */
#endif
#endif

static OBJECT  *ostable[OSTSIZ];	/* the object set table */


void
insertelem(os, obj)		/* insert obj into os, no questions */
register OBJECT  *os;
OBJECT  obj;
{
	register int  i;
	
	for (i = os[0]++; i > 0; i--)
		if (os[i] > obj)
			os[i+1] = os[i];
		else
			break;
	os[i+1] = obj;
}


void
deletelem(os, obj)		/* delete obj from os, no questions */
register OBJECT  *os;
OBJECT  obj;
{
	register int  i;

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
inset(os, obj)			/* determine if object is in set */
register OBJECT  *os;
OBJECT  obj;
{
	int  upper, lower;
	register int  cm, i;

	if ((i = os[0]) <= 6) {		/* linear search algorithm */
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
setequal(os1, os2)		/* determine if two sets are equal */
register OBJECT  *os1, *os2;
{
	register int  i;

	for (i = *os1; i-- >= 0; )
		if (*os1++ != *os2++)
			return(0);
	return(1);
}


void
setcopy(os1, os2)		/* copy object set os2 into os1 */
register OBJECT  *os1, *os2;
{
	register int  i;

	for (i = *os2; i-- >= 0; )
		*os1++ = *os2++;
}


OBJECT *
setsave(os)			/* allocate space and save set */
register OBJECT  *os;
{
	OBJECT  *osnew;
	register OBJECT  *oset;
	register int  i;

	if ((osnew = oset = (OBJECT *)malloc((*os+1)*sizeof(OBJECT))) == NULL)
		error(SYSTEM, "out of memory in setsave\n");
	for (i = *os; i-- >= 0; )	/* inline setcopy */
		*oset++ = *os++;
	return(osnew);
}


void
setunion(osr, os1, os2)		/* osr = os1 Union os2 */
register OBJECT  *osr, *os1, *os2;
{
	register int	i1, i2;

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
setintersect(osr, os1, os2)	/* osr = os1 Intersect os2 */
register OBJECT  *osr, *os1, *os2;
{
	register int	i1, i2;

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
fullnode(oset)			/* return octree for object set */
OBJECT  *oset;
{
	int  osentry, ntries;
	long  hval;
	OCTREE  ot;
	register int  i;
	register OBJECT  *os;
					/* hash on set */
	hval = 0;
	os = oset;
	i = *os++;
	while (i-- > 0)
		hval += *os++;
	ntries = 0;
tryagain:
	osentry = (hval + (long)ntries*ntries) % OSTSIZ;
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
		if (!isfull(ot))		/* entry overflow */
			if (++ntries < OSTSIZ)
				goto tryagain;
			else
				error(INTERNAL, "hash table overflow in fullnode");
						/* remember position */
		i = os - ostable[osentry];
		os = ostable[osentry] = (OBJECT *)realloc(
				(char *)ostable[osentry],
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
}


void
objset(oset, ot)		/* get object set for full node */
register OBJECT  *oset;
OCTREE  ot;
{
	register OBJECT  *os;
	register int  i;

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


int
dosets(f)				/* loop through all sets */
int	(*f)();
{
	int  res = 0;
	int  n;
	register OBJECT  *os;

	for (n = 0; n < OSTSIZ; n++) {
		if ((os = ostable[n]) == NULL)
			continue;
		while (*os > 0) {
			res += (*f)(os);
			os += *os + 1;
		}
	}
	return(res);
}


void
donesets()			/* free ALL SETS in our table */
{
	register int  n;

	for (n = 0; n < OSTSIZ; n++) 
		if (ostable[n] != NULL) {
			free((void *)ostable[n]);
			ostable[n] = NULL;
		}
}
