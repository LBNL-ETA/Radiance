#ifndef lint
static const char	RCSid[] = "$Id: brandom.c,v 1.2 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  brandom.c - blue noise function.
 *
 *	11/8/87
 */

#include  "random.h"


double
brandom(l)			/* blue noise function */
int  l;				/* length between 1 and 8 */
{
	static float  his[8] = {.5,.5,.5,.5,.5,.5,.5,.5};
	static double  avg = .5;
	static int  x = 0;
	double  y = frandom();
	
	if (avg < .5)
		y = 1 - 2*y*avg;
	else
		y = 2*y*(1 - avg);
					/* update */
	avg += (y - his[x])/l;
	his[x] = y;
	if (++x >= l) x = 0;
	
	return(y);
}
