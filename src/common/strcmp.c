#ifndef lint
static const char	RCSid[] = "$Id: strcmp.c,v 2.3 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * String comparison routine optimized for use with savestr.c
 */

#include "copyright.h"


int
strcmp(s1, s2)				/* check for s1==s2 */
register char  *s1, *s2;
{
	if (s1 == s2)
		return(0);

	while (*s1 == *s2++)
		if (!*s1++)
			return(0);

	return(*s1 - *--s2);
}
