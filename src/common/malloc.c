#ifndef lint
static const char	RCSid[] = "$Id: malloc.c,v 2.16 2003/06/30 14:59:11 schorsch Exp $";
#endif
/*
 * Fast malloc for memory hogs and VM environments.
 * Performs a minimum of searching through free lists.
 * malloc(), realloc() and free() use buckets
 *	containing blocks in powers of two, similar to CIT routines.
 * bfree() puts memory from any source into malloc free lists.
 * bmalloc() doesn't keep track of free lists -- it's usually
 *	just a buffered call to sbrk(2).  However, bmalloc will
 *	call mscrounge() if sbrk fails.
 * mscrounge() returns whatever free memory it can find to satisfy the
 *	request along with the number of bytes (modified).
 * mcompact() takes the same arguments as mscrounge() (and is in fact
 *	called by mscrounge() under dire circumstances), but it uses
 *	memory compaction to return a larger block.  (MCOMP)
 *
 *	Greg Ward	Lawrence Berkeley Laboratory
 */

#include "copyright.h"

#include  <errno.h>
#include  <string.h>


#ifdef MSTATS
#include  <stdio.h>
static unsigned	b_nsbrked = 0;
static unsigned	b_nalloced = 0;
static unsigned	b_nfreed = 0;
static unsigned	b_nscrounged = 0;
static unsigned b_ncompacted = 0;
static unsigned	m_nalloced = 0;
static unsigned	m_nfreed = 0;
static unsigned	m_nwasted = 0;
#else
#define  NULL		0
#endif

#ifndef ALIGNT
#define  ALIGNT		int			/* align type */
#endif
#define  BYTES_WORD	sizeof(ALIGNT)

#ifndef  MAXINCR
#define  MAXINCR	(1<<16)			/* largest sbrk(2) increment */
#endif
					/* malloc free lists */
typedef union m_head {
	union m_head	*next;
	struct {
		short		magic;
		short		bucket;
	}	a;
	ALIGNT		dummy;
} M_HEAD;

#define MAGIC		0x1a2		/* magic number for allocated memory */

#define FIRSTBUCKET	3
#define NBUCKETS	30

static M_HEAD	*free_list[NBUCKETS];

static ALIGNT	dummy_mem;

static char	*memlim[2];

#define DUMMYLOC	((char *)&dummy_mem)

#define BADPTR(p)	((p) < memlim[0] | (p) >= memlim[1])

#ifdef  MCOMP			/* memory compaction routines */
static char	seedtab[1024];		/* seed for compaction table */

static struct mblk {
	char		*ptr;	/* pointer to memory block */
	unsigned	siz;	/* size of memory block */
}	cptab = {seedtab, sizeof(seedtab)};

#define mtab(mb)	((struct mblk *)(mb)->ptr)
#define mtablen(mb)	((mb)->siz/sizeof(struct mblk))


static int
mbcmp(m1, m2)		/* compare memory block locations */
register struct mblk	*m1, *m2;
{
				/* put empty blocks at end */
	if (m1->siz == 0)
		return(m2->siz == 0 ? 0 : 1);
	if (m2->siz == 0)
		return(-1);
				/* otherwise, ascending order */
	return(m1->ptr - m2->ptr);
}


int
compactfree()		/* compact free lists (moving to table) */
{
	register struct mblk	*tp;
	register int	bucket;
	unsigned	n, bsiz, ncomp;
				/* fill table from free lists */
	tp = mtab(&cptab); n = mtablen(&cptab);
	bucket = NBUCKETS-1; bsiz = 1<<(NBUCKETS-1);
	ncomp = 0;
	for ( ; ; ) {
		while (tp->siz) {	/* find empty slot */
			if (--n == 0)
				goto loopexit;
			tp++;
		}
		while (free_list[bucket] == NULL) {	/* find block */
			if (--bucket < FIRSTBUCKET)
				goto loopexit;
			bsiz >>= 1;
		}
		tp->ptr = (char *)free_list[bucket];
		ncomp += (tp->siz = bsiz + sizeof(M_HEAD));
		free_list[bucket] = free_list[bucket]->next;
	}
loopexit:
	if (ncomp == 0)
		return(0);		/* nothing new to compact */
#ifdef MSTATS
	b_ncompacted += ncomp;
#endif
	tp = mtab(&cptab); n = mtablen(&cptab);		/* sort table */
	qsort((char *)tp, n, sizeof(struct mblk), mbcmp);
	ncomp = 0;					/* collect blocks */
	while (--n && tp[1].siz) {
		if (tp[0].ptr + tp[0].siz == tp[1].ptr) {
			tp[1].ptr = tp[0].ptr;
			tp[1].siz += tp[0].siz;
			ncomp += tp[0].siz;
			tp[0].siz = 0;
		}
		tp++;
	}
	return(ncomp);
}


char *
mcompact(np)		/* compact free lists to satisfy request */
unsigned	*np;
{
	struct mblk	*tab;
	unsigned	tablen;
	register struct mblk	*bp, *big;

	for ( ; ; ) {
					/* compact free lists */
		while (compactfree())
			;
					/* find largest block */
		tab = mtab(&cptab); tablen = mtablen(&cptab);
		big = tab;
		bp = tab + tablen;
		while (--bp > tab)
			if (bp->siz > big->siz)
				big = bp;
		if (big->siz >= *np) {		/* big enough? */
			*np = big->siz;
			big->siz = 0;		/* remove from table */
			return(big->ptr);	/* return it */
		}
		if (mtablen(big) <= tablen) {
			*np = 0;		/* cannot grow table */
			return(NULL);		/* report failure */
		}
				/* enlarge table, free old one */
		mtab(big)[0].ptr = cptab.ptr;
		mtab(big)[0].siz = cptab.siz;
		cptab.ptr = big->ptr;
		cptab.siz = big->siz;
		big->siz = 0;			/* clear and copy */
		memcpy((char *)(mtab(&cptab)+1), (char *)tab,
				tablen*sizeof(struct mblk));
		memset((char *)(mtab(&cptab)+tablen+1), '\0',
			(mtablen(&cptab)-tablen-1)*sizeof(struct mblk));
	}			/* next round */
}
#endif		/* MCOMP */


char *
mscrounge(np)		/* search free lists to satisfy request */
register unsigned	*np;
{
	char	*p;
	register int	bucket;
	register unsigned	bsiz;

	for (bucket = FIRSTBUCKET, bsiz = 1<<FIRSTBUCKET;
			bucket < NBUCKETS; bucket++, bsiz <<= 1)
		if (free_list[bucket] != NULL && bsiz+sizeof(M_HEAD) >= *np) {
			*np = bsiz+sizeof(M_HEAD);
			p = (char *)free_list[bucket];
			free_list[bucket] = free_list[bucket]->next;
#ifdef MSTATS
			b_nscrounged += *np;
#endif
			return(p);
		}
#ifdef MCOMP
	return(mcompact(np));		/* try compaction */
#else
	*np = 0;
	return(NULL);
#endif
}


char *
bmalloc(n)		/* allocate a block of n bytes, sans header */
register unsigned  n;
{
	extern char  *sbrk();
	static unsigned  pagesz = 0;
	static unsigned  amnt;
	static char  *bpos;
	static unsigned  nrem;
	unsigned  thisamnt;
	register char	*p;

	if (pagesz == 0) {				/* initialize */
		pagesz = amnt = getpagesize();
		nrem = (int)sbrk(0);			/* page align break */
		nrem = pagesz - (nrem&(pagesz-1));
		bpos = sbrk(nrem);
		if ((int)bpos == -1)
			return(NULL);
#ifdef MSTATS
		b_nsbrked += nrem;
#endif
		bpos += nrem & (BYTES_WORD-1);		/* align pointer */
		nrem &= ~(BYTES_WORD-1);
		memlim[0] = bpos;
		memlim[1] = bpos + nrem;
	}

	n = (n+(BYTES_WORD-1))&~(BYTES_WORD-1);		/* word align rqst. */

	if (n > nrem) {					/* need more core */
	tryagain:
		if (n > amnt) {				/* big chunk */
			thisamnt = (n+(pagesz-1))&~(pagesz-1);
			if (thisamnt <= MAXINCR)	/* increase amnt */
				amnt = thisamnt;
		} else
			thisamnt = amnt;
		p = sbrk(thisamnt);
		if ((int)p == -1) {			/* uh-oh, ENOMEM */
			errno = 0;			/* call cavalry */
			if (thisamnt >= n+pagesz) {
				amnt = pagesz;		/* minimize request */
				goto tryagain;
			}
			thisamnt = n;
			p = mscrounge(&thisamnt);	/* search free lists */
			if (p == NULL) {		/* we're really out */
				errno = ENOMEM;
				return(NULL);
			}
		}
#ifdef MSTATS
		else b_nsbrked += thisamnt;
#endif
		if (p != bpos+nrem) {			/* not an increment */
			bfree(bpos, nrem);		/* free what we have */
			bpos = p;
			nrem = thisamnt;
		} else					/* otherwise tack on */
			nrem += thisamnt;
		if (bpos < memlim[0])
			memlim[0] = bpos;
		if (bpos + nrem > memlim[1])
			memlim[1] = bpos + nrem;
	}
	p = bpos;
	bpos += n;					/* advance */
	nrem -= n;
#ifdef MSTATS
	b_nalloced += n;
#endif
	return(p);
}


bfree(p, n)			/* free n bytes of random memory */
register char  *p;
register unsigned  n;
{
	register int	bucket;
	register unsigned	bsiz;

	if (n < 1<<FIRSTBUCKET || p == NULL)
		return;
#ifdef MSTATS
	b_nfreed += n;
#endif
					/* align pointer */
	bsiz = BYTES_WORD - ((unsigned)p&(BYTES_WORD-1));
	if (bsiz < BYTES_WORD) {
		p += bsiz;
		n -= bsiz;
	}
	if (p < memlim[0])
		memlim[0] = p;
	if (p + n > memlim[1])
		memlim[1] = p + n;
					/* fill big buckets first */
	for (bucket = NBUCKETS-1, bsiz = 1<<(NBUCKETS-1);
			bucket >= FIRSTBUCKET; bucket--, bsiz >>= 1)
		if (n >= bsiz+sizeof(M_HEAD)) {
			((M_HEAD *)p)->next = free_list[bucket];
			free_list[bucket] = (M_HEAD *)p;
			p += bsiz+sizeof(M_HEAD);
			n -= bsiz+sizeof(M_HEAD);
		}
}


char *
malloc(n)			/* allocate n bytes of memory */
unsigned	n;
{
	register M_HEAD	*mp;
	register int	bucket;
	register unsigned	bsiz;
					/* don't return NULL on 0 request */
	if (n == 0)
		return(DUMMYLOC);
					/* find first bucket that fits */
	for (bucket = FIRSTBUCKET, bsiz = 1<<FIRSTBUCKET;
			bucket < NBUCKETS; bucket++, bsiz <<= 1)
		if (bsiz >= n)
			break;
	if (bucket >= NBUCKETS) {
		errno = EINVAL;
		return(NULL);
	}
	if (free_list[bucket] == NULL) {	/* need more core */
		mp = (M_HEAD *)bmalloc(bsiz+sizeof(M_HEAD));
		if (mp == NULL)
			return(NULL);
	} else {				/* got it */
		mp = free_list[bucket];
		free_list[bucket] = mp->next;
	}
#ifdef MSTATS
	m_nalloced += bsiz + sizeof(M_HEAD);
	m_nwasted += bsiz + sizeof(M_HEAD) - n;
#endif
	mp->a.magic = MAGIC;			/* tag block */
	mp->a.bucket = bucket;
	return((char *)(mp+1));
}


char *
realloc(op, n)			/* reallocate memory using malloc() */
register char	*op;
unsigned	n;
{
	char	*p;
	register unsigned	on;
					/* get old size */
	if (op != DUMMYLOC && !BADPTR(op) &&
			((M_HEAD *)op-1)->a.magic == MAGIC)
		on = 1 << ((M_HEAD *)op-1)->a.bucket;
	else
		on = 0;
	if (n <= on && (n > on>>1 || on == 1<<FIRSTBUCKET))
		return(op);		/* same bucket */
	if ((p = malloc(n)) == NULL)
		return(n<=on ? op : NULL);
	if (on) {
		memcpy(p, op, n>on ? on : n);
		free(op);
	}
	return(p);
}


free(p)			/* free memory allocated by malloc */
char	*p;
{
	register M_HEAD	*mp;
	register int	bucket;

	if (p == DUMMYLOC)
		return(1);
	if (BADPTR(p))
		goto invalid;
	mp = (M_HEAD *)p - 1;
	if (mp->a.magic != MAGIC)		/* sanity check */
		goto invalid;
	bucket = mp->a.bucket;
	if (bucket < FIRSTBUCKET | bucket >= NBUCKETS)
		goto invalid;
	mp->next = free_list[bucket];
	free_list[bucket] = mp;
#ifdef MSTATS
	m_nfreed += (1 << bucket) + sizeof(M_HEAD);
#endif
	return(1);
invalid:
	errno = EINVAL;
	return(0);
}


#ifdef MSTATS
printmemstats(fp)		/* print memory statistics to stream */
FILE	*fp;
{
	register int	i;
	int	n;
	register M_HEAD	*mp;
	unsigned int	total = 0;

	fprintf(fp, "Memory statistics:\n");
	fprintf(fp, "\tbmalloc: %d bytes from sbrk\n", b_nsbrked);
	fprintf(fp, "\tbmalloc: %d bytes allocated\n", b_nalloced);
	fprintf(fp, "\tbmalloc: %d bytes freed\n", b_nfreed);
	fprintf(fp, "\tbmalloc: %d bytes scrounged\n", b_nscrounged);
	fprintf(fp, "\tbmalloc: %d bytes compacted\n", b_ncompacted);
	fprintf(fp, "\tmalloc: %d bytes allocated\n", m_nalloced);
	fprintf(fp, "\tmalloc: %d bytes wasted (%.1f%%)\n", m_nwasted,
			100.0*m_nwasted/m_nalloced);
	fprintf(fp, "\tmalloc: %d bytes freed\n", m_nfreed);
	for (i = NBUCKETS-1; i >= FIRSTBUCKET; i--) {
		n = 0;
		for (mp = free_list[i]; mp != NULL; mp = mp->next)
			n++;
		if (n) {
			fprintf(fp, "\t%d * %u\n", n, 1<<i);
			total += n * ((1<<i) + sizeof(M_HEAD));
		}
	}
	fprintf(fp, "\t %u total bytes in free list\n", total);
}
#endif
