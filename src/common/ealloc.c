#ifndef lint
static const char	RCSid[] = "$Id: ealloc.c,v 2.7 2003/10/27 10:19:31 schorsch Exp $";
#endif
/*
 *  ealloc.c - memory routines which call quit on error.
 */

#include "copyright.h"


#include  <stdio.h>
#include  <stdlib.h>

#include "rterror.h"

extern char *	/* return pointer to n uninitialized bytes */
emalloc(unsigned int  n)
{
	register char  *cp;
	
	if (n == 0)
		return(NULL);

	if ((cp = (char *)malloc(n)) != NULL)
		return(cp);

	eputs("Out of memory in emalloc\n");
	quit(1);
}


extern char *			/* return pointer to initialized memory */
ecalloc(register unsigned int ne, unsigned int es)
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


extern char *			/* reallocate cp to size n */
erealloc(register char  *cp, unsigned int  n)
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


extern void			/* free memory allocated by above */
efree(char  *cp)
{
	free((void *)cp);
}
