#ifndef lint
static const char	RCSid[] = "$Id: multisamp.c,v 2.6 2024/04/17 15:07:29 greg Exp $";
#endif
/*
 * Binary space partitioning curve for multidimensional sampling.
 *
 *	Written by Christophe Schlick
 */

#include "copyright.h"

#include <stdlib.h>

#include "random.h"


/* Convert 1-dimensional sample to N dimensions */
void
multisamp(double t[], int n, double r)
{
	int	j;
	int	i, k;
	int	ti[8];
	double	s;

	i = n;
	while (i-- > 0)
		ti[i] = 0;
	j = 8;
	while (j--) {
		k = s = r*(1<<n);
		r = s - k;
		i = n;
		while (i-- > 0)
			ti[i] += ti[i] + ((k>>i) & 1);
	}
	i = n;
	while (i-- > 0)
		t[i] = (1./256.) * (ti[i] + frandom());
}
