#ifndef lint
static const char	RCSid[] = "$Id: ealloc.c,v 2.4 2003/04/23 00:52:33 greg Exp $";
#endif
/*
 *  ealloc.c - memory routines which call quit on error.
 */

#include "copyright.h"


#include  <stdio.h>
#include  <stdlib.h>


char *
emalloc(n)			/* return pointer to n uninitialized bytes */
unsigned int  n;
{
	register char  *cp;
	
	if (n == 0)
		return(NULL);

	if ((cp = (char *)malloc(n)) != NULL)
		return(cp);

	eputs("Out of memory in emalloc\n");
	quit(1);
}


char *
ecalloc(ne, es)			/* return pointer to initialized memory */
register unsigned int  ne;
unsigned int  es;
{
	register char  *cp;
	
	ne *= es;
	if (ne == 0)
		return(NULL);

	if ((cp = (char *)malloc(ne)) == NULL) {
		eputs("Out of memory in ecalloc\n");
		quit(1);
	}
	cp += ne;
	while (ne--)
		*--cp = 0;
	return(cp);
}


char *
erealloc(cp, n)			/* reallocate cp to size n */
register char  *cp;
unsigned int  n;
{
	if (n == 0) {
		if (cp != NULL)
			free((void *)cp);
		return(NULL);
	}

	if (cp == NULL)
		cp = (char *)malloc(n);
	else 
		cp = (char *)realloc((void *)cp, n);

	if (cp != NULL)
		return(cp);

	eputs("Out of memory in erealloc\n");
	quit(1);
}


void
efree(cp)			/* free memory allocated by above */
char  *cp;
{
	free((void *)cp);
}
