#ifndef lint
static const char	RCSid[] = "$Id: ealloc.c,v 2.11 2022/01/19 00:59:33 greg Exp $";
#endif
/*
 *  ealloc.c - memory routines which call quit on error.
 */

#include "copyright.h"


#include  <stdio.h>
#include  <stdlib.h>

#include "rterror.h"
#include "rtmisc.h"

void *			/* return pointer to n uninitialized bytes */
emalloc(size_t  n)
{
	void  *cp;
	
	if (n == 0)
		return(NULL);

	if ((cp = malloc(n)) != NULL)
		return(cp);

	eputs("Out of memory in emalloc\n");
	quit(1);
	return NULL; /* pro forma return */
}


void *			/* return pointer to initialized memory */
ecalloc(size_t ne, size_t es)
{
	void  *cp;
	
	if (!ne | !es)
		return(NULL);

	if ((cp = calloc(ne, es)) != NULL)
		return(cp);

	eputs("Out of memory in ecalloc\n");
	quit(1);
	return(NULL); /* pro forma return */
}


void *			/* reallocate cp to size n */
erealloc(void  *cp, size_t  n)
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


void			/* free memory allocated by above */
efree(void *cp)
{
	if (cp)
		free(cp);
}
