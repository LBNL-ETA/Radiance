#ifndef lint
static const char	RCSid[] = "$Id: savqstr.c,v 2.4 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 *  Save unshared strings.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

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
