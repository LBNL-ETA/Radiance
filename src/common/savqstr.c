#ifndef lint
static const char	RCSid[] = "$Id: savqstr.c,v 2.3 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  Save unshared strings.
 *
 *  External symbols declared in standard.h
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

#include <stdlib.h>

extern void	eputs();
extern void	quit();

#if 1

char *
savqstr(s)			/* save a private string */
register char  *s;
{
	register char  *cp;
	char  *newp;

	for (cp = s; *cp++; )			/* compute strlen()+1 */
		;
	newp = (char *)malloc(cp-s);
	if (newp == NULL) {
		eputs("out of memory in savqstr");
		quit(1);
	}
	for (cp = newp; *cp++ = *s++; )		/* inline strcpy() */
		;
	return(newp);				/* return new location */
}


void
freeqstr(s)			/* free a private string */
char  *s;
{
	free((void *)s);
}

#else

/*
 *  Save unshared strings, packing them together into
 *  large blocks to optimize paging in VM environments.
 */

#ifdef  BIGMEM
#ifndef  MINBLOCK
#define  MINBLOCK	(1<<12)		/* minimum allocation block size */
#endif
#ifndef  MAXBLOCK
#define  MAXBLOCK	(1<<16)		/* maximum allocation block size */
#endif
#else
#ifndef  MINBLOCK
#define  MINBLOCK	(1<<10)		/* minimum allocation block size */
#endif
#ifndef  MAXBLOCK
#define  MAXBLOCK	(1<<14)		/* maximum allocation block size */
#endif
#endif

extern char  *bmalloc();


char *
savqstr(s)			/* save a private string */
register char  *s;
{
	static char  *curp = NULL;		/* allocated memory pointer */
	static unsigned  nrem = 0;		/* bytes remaining in block */
	static unsigned  nextalloc = MINBLOCK;	/* next block size */
	register char  *cp;
	register unsigned  n;

	for (cp = s; *cp++; )			/* compute strlen()+1 */
		;
	if ((n = cp-s) > nrem) {		/* do we need more core? */
		bfree(curp, nrem);			/* free remnant */
		while (n > nextalloc)
			nextalloc <<= 1;
		if ((curp = bmalloc(nrem=nextalloc)) == NULL) {
			eputs("out of memory in savqstr");
			quit(1);
		}
		if ((nextalloc <<= 1) > MAXBLOCK)	/* double block size */
			nextalloc = MAXBLOCK;
	}
	for (cp = curp; *cp++ = *s++; )		/* inline strcpy() */
		;
	s = curp;				/* update allocation info. */
	curp = cp;
	nrem -= n;
	return(s);				/* return new location */
}


void
freeqstr(s)			/* free a private string (not recommended) */
char  *s;
{
	bfree(s, strlen(s)+1);
}

#endif
