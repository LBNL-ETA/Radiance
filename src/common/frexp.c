#ifndef lint
static const char	RCSid[] = "$Id: frexp.c,v 1.3 2003/02/25 02:47:21 greg Exp $";
#endif
#include "copyright.h"

frexp(x, ip)		/* call it paranoia, I've seen the lib version */
register double  x;
int  *ip;
{
	int  neg;
	register int  i;

	if (neg = (x < 0.0))
		x = -x;
	else if (x == 0.0) {
		*ip = 0;
		return(0.0);
	}
	if (x < 0.5)
		for (i = 0; x < 0.5; i--)
			x *= 2.0;
	else
		for (i = 0; x >= 1.0; i++)
			x /= 2.0;
	*ip = i;
	if (neg)
		return(-x);
	else
		return(x);
}
