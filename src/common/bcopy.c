#ifndef lint
static const char	RCSid[] = "$Id: bcopy.c,v 1.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * bcopy.c - substitutes for library routines.
 */

#include "copyright.h"


bcopy(src, dest, nbytes)
register char	*src, *dest;
register int	nbytes;
{
	while (nbytes-- > 0)
		*dest++ = *src++;
}


bzero(b, nbytes)
register char	*b;
register int	nbytes;
{
	while (nbytes-- > 0)
		*b++ = 0;
}


int
bcmp(b1, b2, nbytes)
register unsigned char	*b1, *b2;
register int	nbytes;
{
	while (nbytes-- > 0)
		if (*b1++ - *b2++)
			return(*--b1 - *--b2);
	return(0);
}
