#ifndef lint
static const char	RCSid[] = "$Id: bmalloc.c,v 2.5 2003/07/17 09:21:29 schorsch Exp $";
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

#ifndef  MBLKSIZ
#define  MBLKSIZ	16376		/* size of memory allocation block */
#endif
#define  WASTEFRAC	12		/* don't waste more than a fraction */
#ifndef  ALIGNT
#define  ALIGNT		double		/* type for alignment */
#endif
#define  BYTES_WORD	sizeof(ALIGNT)

static char  *bposition = NULL;
static unsigned int  nremain = 0;


char *
bmalloc(n)		/* allocate a block of n bytes */
register unsigned int  n;
{
	if (n > nremain && (n > MBLKSIZ || nremain > MBLKSIZ/WASTEFRAC))
		return(malloc(n));			/* too big */

	n = (n+(BYTES_WORD-1))&~(BYTES_WORD-1);		/* word align */

	if (n > nremain && (bposition = (char *)malloc(nremain = MBLKSIZ)) == NULL) {
		nremain = 0;
		return(NULL);
	}
	bposition += n;
	nremain -= n;
	return(bposition - n);
}


void
bfree(p, n)			/* free random memory */
register char	*p;
register unsigned int	n;
{
	register unsigned int	bsiz;
					/* check alignment */
	bsiz = BYTES_WORD - ((unsigned int)p&(BYTES_WORD-1));
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
