#ifndef lint
static const char	RCSid[] = "$Id: urind.c,v 2.4 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Compute pseudo-asyncronous entry point for urand(3)
 */

#include "copyright.h"
#include "random.h"

#define NBITS		32		/* number of bits in an integer */

static char  bctab[256] = {
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
	};

#if NBITS==32
#define  bitcount(i)	(bctab[(i)>>24&0xff]+bctab[(i)>>16&0xff]+ \
				bctab[(i)>>8&0xff]+bctab[(i)&0xff])
#endif
#if  NBITS==16
#define  bitcount(i)	(bctab[(i)>>8&0xff]+bctab[(i)&0xff])
#endif


int
urind(s, i)			/* compute i'th index from seed s */
int	s, i;
{
	register int  ss, k;
	int  left;

	ss = s*1103515245 + 12345;
	left = 0;
	for (k = i/NBITS; k--; ) {
		left += bitcount(ss);
		ss = ss*1103515245 + 12345;
	}
	for (k = i&(NBITS-1); k--; ss >>= 1)
		left += ss & 1;
	if (ss & 1)
		return(s-left-1);
	return(s-left+i);
}
