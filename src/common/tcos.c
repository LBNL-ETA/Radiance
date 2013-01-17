#ifndef lint
static const char	RCSid[] = "$Id: tcos.c,v 3.6 2013/01/17 22:48:21 greg Exp $";
#endif
/*
 * Table-based cosine approximation.
 *
 * Use doubles in table even though we're not nearly that accurate just
 * to avoid conversion and guarantee that tsin(x)^2 + tcos(x)^2 == 1.
 *
 * No interpolation in this version.
 *
 * External symbols declared in rtmath.h
 */

#include "copyright.h"

#include <math.h>

#include "rtmath.h"

#ifndef NCOSENTRY
#define NCOSENTRY	1024
#endif


double
tcos(double x)				/* approximate cosine */
{
	static double	costab[NCOSENTRY+1];
	register int	i;

	if (costab[0] < 0.5)		/* initialize table */
		for (i = 0; i <= NCOSENTRY; i++)
			costab[i] = cos((PI/2./NCOSENTRY)*i);
					/* normalize angle */
	if (x < 0.)
		x = -x;
	i = (NCOSENTRY*2./PI) * x  +  0.5;
	while (i >= 4*NCOSENTRY)
		i -= 4*NCOSENTRY;
	switch (i / NCOSENTRY) {
	case 0:
		return(costab[i]);
	case 1:
		return(-costab[(2*NCOSENTRY)-i]);
	case 2:
		return(-costab[i-(2*NCOSENTRY)]);
	case 3:
		return(costab[(4*NCOSENTRY)-i]);
	}
	return(0.);		/* should never be reached */
}
