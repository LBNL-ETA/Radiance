#ifndef lint
static const char	RCSid[] = "$Id: urand.c,v 2.6 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Anticorrelated random function due to Christophe Schlick
 */

#include "copyright.h"

#include  <stdlib.h>
#include  "random.h"

#undef initurand

#define  MAXORDER	(8*sizeof(unsigned short))

unsigned short  *urperm = NULL;	/* urand() permutation */
int  urmask;		/* bits used in permutation */

int
initurand(size)		/* initialize urand() for size entries */
int  size;
{
	int  order, n;
	register int  i, offset;

	if (urperm != NULL)
		free((void *)urperm);
	if (--size <= 0) {
		urperm = NULL;
		urmask = 0;
		return(0);
	}
	for (i = 1; size >>= 1; i++)
		;
	order = i>MAXORDER ? MAXORDER : i;
	urmask = (1<<i) - 1;
	urperm = (unsigned short *)malloc((urmask+1)*sizeof(unsigned short));
	if (urperm == NULL) {
		eputs("out of memory in initurand\n");
		quit(1);
	}
	urperm[0] = 0;
	for (n = 1, offset = 1; n <= order; n++, offset <<= 1)
		for (i = offset; i--; ) {
			urperm[i] =
			urperm[i+offset] = 2*urperm[i];
			if (random() & 0x4000)
				urperm[i]++;
			else
				urperm[i+offset]++;
		}
	return(1<<order);
}


int
ilhash(d, n)			/* hash a set of integer values */
register int  *d;
register int  n;
{
	static int  tab[8] = {13623,353,1637,5831,2314,3887,5832,8737};
	register int  hval;

	hval = 0;
	while (n-- > 0)
		hval += *d++ * tab[n&7];
	return(hval);
}
