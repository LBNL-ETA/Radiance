#ifndef lint
static const char	RCSid[] = "$Id: urand.c,v 2.8 2004/05/14 20:49:13 greg Exp $";
#endif
/*
 * Anticorrelated random function due to Christophe Schlick
 */

#include "copyright.h"

#include  <stdlib.h>

#include  "standard.h"
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
	for (i = 1; (size >>= 1); i++)
		if (i == MAXORDER)
			break;
	order = i;
	urmask = (1<<order) - 1;
	urperm = (unsigned short *)malloc((urmask+1)*sizeof(unsigned short));
	if (urperm == NULL) {
		eputs("out of memory in initurand\n");
		quit(1);
	}
	urperm[0] = (random() & 0x4000) != 0;
	for (n = 1, offset = 1; n <= order; n++, offset <<= 1)
		for (i = offset; i--; ) {
			urperm[i+offset] = urperm[i] <<= 1;
			if (random() & 0x4000)
				urperm[i]++;
			else
				urperm[i+offset]++;
		}
	return(urmask+1);
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
