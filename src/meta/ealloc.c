#ifndef lint
static const char	RCSid[] = "$Id: ealloc.c,v 1.2 2003/04/23 00:52:34 greg Exp $";
#endif
/*
 *  ealloc.c - memory routines which call quit on error.
 */


#include  <stdio.h>
#include  <stdlib.h>



char *
emalloc(n)			/* return pointer to n uninitialized bytes */
unsigned  n;
{
	register char  *cp;
	
	if (n == 0)
		return(NULL);

	if ((cp = malloc(n)) != NULL)
		return(cp);

	eputs("Out of memory in emalloc\n");
	quit(1);
}


char *
ecalloc(ne, es)			/* return pointer to initialized memory */
register unsigned  ne;
unsigned  es;
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


char *
erealloc(cp, n)			/* reallocate cp to size n */
register char  *cp;
unsigned  n;
{
	if (n == 0) {
		if (cp != NULL)
			free(cp);
		return(NULL);
	}

	if (cp == NULL)
		cp = malloc(n);
	else 
		cp = realloc((void *)cp, n);

	if (cp != NULL)
		return(cp);

	eputs("Out of memory in erealloc\n");
	quit(1);
}


efree(cp)			/* free memory allocated by above */
char  *cp;
{
	free((void *)cp);
}
