#ifndef lint
static const char	RCSid[] = "$Id: bmalloc.c,v 2.10 2025/04/12 01:43:53 greg Exp $";
#endif
/*
 * Bmalloc provides basic memory allocation without overhead (no free lists).
 * Use only to take the load off of malloc for all those
 * piddling little requests that you never expect to free.
 * Bmalloc defers to malloc for big requests.
 * Bfree should hand memory to bmalloc, but it usually fails here.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdlib.h>

#include "rtmisc.h"

#ifndef  ALIGNT
#define  ALIGNT		double		/* type for alignment */
#endif
#define  BYTES_WORD	sizeof(ALIGNT)

#ifndef  MBLKSIZ			/* size of memory allocation block */
#define  MBLKSIZ	((1<<15)-BYTES_WORD)
#endif
#define  WASTEFRAC	12		/* don't waste more than a fraction */

static char  *bposition = NULL;
static size_t  nremain = 0;

void *
bmalloc(		/* quickly allocate a block of n bytes */
size_t  n
)
{
	if ((n > nremain) & ((n > MBLKSIZ) | (nremain > MBLKSIZ/WASTEFRAC)))
		return(malloc(n));			/* too big */

	n = (n+(BYTES_WORD-1)) & ~(BYTES_WORD-1);	/* word align */

	if (n > nremain && (bposition = malloc(nremain = MBLKSIZ)) == NULL) {
		nremain = 0;
		return(NULL);
	}
	bposition += n;
	nremain -= n;
	return(bposition - n);
}

void
bfree(			/* free some memory; anything in process space */
void	*pp,
size_t	n
)
{
	char *p = pp;
	size_t	bsiz;
					/* check alignment */
	bsiz = BYTES_WORD - ((size_t)p&(BYTES_WORD-1));
	if (bsiz < BYTES_WORD) {
		p += bsiz;
		n -= bsiz;
	}
	if (p + n == bposition) {	/* just allocated? */
		bposition = p;
		nremain += n;
		return;
	}
	if (n > nremain) {		/* better than what we've got? */
		bposition = p;
		nremain = n;
		return;
	}
				/* just throw it away, then */
}
