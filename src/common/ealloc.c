#ifndef lint
static const char	RCSid[] = "$Id: ealloc.c,v 2.9 2004/03/28 20:33:12 schorsch Exp $";
#endif
/*
 *  ealloc.c - memory routines which call quit on error.
 */

#include "copyright.h"


#include  <stdio.h>
#include  <stdlib.h>

#include "rterror.h"
#include "rtmisc.h"

extern void *	/* return pointer to n uninitialized bytes */
emalloc(size_t  n)
{
	register void  *cp;
	
	if (n == 0)
		return(NULL);

	if ((cp = malloc(n)) != NULL)
		return(cp);

	eputs("Out of memory in emalloc\n");
	quit(1);
	return NULL; /* pro forma return */
}


extern void *			/* return pointer to initialized memory */
ecalloc(register size_t ne, size_t es)
{
	register char  *cp;
	
	ne *= es;
	if (ne == 0)
		return(NULL);

	if ((cp = malloc(ne)) == NULL) {
		eputs("Out of memory in ecalloc\n");
		quit(1);
	}
	cp += ne;
	while (ne--)
		*--cp = 0;
	return(cp);
}


extern void *			/* reallocate cp to size n */
erealloc(register void  *cp, size_t  n)
{
	if (n == 0) {
		if (cp != NULL)
			free(cp);
		return(NULL);
	}

	if (cp == NULL)
		cp = malloc(n);
	else 
		cp = realloc(cp, n);

	if (cp != NULL)
		return(cp);

	eputs("Out of memory in erealloc\n");
	quit(1);
	return NULL; /* pro forma return */
}


extern void			/* free memory allocated by above */
efree(void  *cp)
{
	free(cp);
}
