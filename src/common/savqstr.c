#ifndef lint
static const char RCSid[] = "$Id: savqstr.c,v 2.12 2024/10/29 00:35:06 greg Exp $";
#endif
/*
 *  Save unshared strings.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdlib.h>

#include "rtio.h"
#include "rterror.h"


#if 1

char *
savqstr(const char *s)			/* save a private string */
{
	char  *cp;
	char  *newp;

	if (s == NULL)
		return(NULL);
	if (!*s)
		return("");
	for (cp = s; *cp++; )			/* compute strlen()+1 */
		;
	newp = (char *)malloc(cp-s);
	if (newp == NULL) {
		eputs("out of memory in savqstr");
		quit(1);
	}
	for (cp = newp; (*cp++ = *s++); )	/* inline strcpy() */
		;
	return(newp);				/* return new location */
}


void
freeqstr(char *s)			/* free a private string */
{
	if (s != NULL && *s)
		free(s);
}

#else

/*
 *  Save unshared strings, packing them together into
 *  large blocks to optimize paging in VM environments.
 */

#include "rtmisc.h"

#ifdef  SMLMEM
#ifndef  MINBLOCK
#define  MINBLOCK	(1<<10)		/* minimum allocation block size */
#endif
#ifndef  MAXBLOCK
#define  MAXBLOCK	(1<<14)		/* maximum allocation block size */
#endif
#else
#ifndef  MINBLOCK
#define  MINBLOCK	(1<<12)		/* minimum allocation block size */
#endif
#ifndef  MAXBLOCK
#define  MAXBLOCK	(1<<16)		/* maximum allocation block size */
#endif
#endif


char *
savqstr(const char *s)			/* save a private string */
{
	static char  *curp = NULL;		/* allocated memory pointer */
	static unsigned  nrem = 0;		/* bytes remaining in block */
	static unsigned  nextalloc = MINBLOCK;	/* next block size */
	char  *cp;
	unsigned  n;

	if (s == NULL)
		return(NULL);
	if (!*s)
		return("");
	for (cp = s; *cp++; )			/* compute strlen()+1 */
		;
	if ((n = cp-s) > nrem) {		/* do we need more core? */
		bfree(curp, nrem);			/* free remnant */
		while (n > nextalloc)
			nextalloc <<= 1;
		if ((curp = (char *)bmalloc(nrem=nextalloc)) == NULL) {
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
freeqstr(char *s)		/* free a private string (not recommended) */
{
	if (s != NULL && *s)
		bfree(s, strlen(s)+1);
}

#endif
