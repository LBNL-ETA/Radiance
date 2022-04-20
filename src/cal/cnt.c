#ifndef lint
static const char	RCSid[] = "$Id: cnt.c,v 1.5 2022/04/20 20:56:01 greg Exp $";
#endif
/*
 *  cnt.c - simple counting program.
 *
 *	2/1/88
 *
 *  Added -s (shuffle) option April 2022
 */

#include  <stdio.h>
#include  <time.h>
#include  "random.h"

#ifndef uint16
#define uint16	unsigned short		/* 16-bit unsigned integer */
#endif
#undef uby8
#define uby8  unsigned char		/* 8-bit unsigned integer */

#define	MAXDIM	50

#define NLEVELS	9			/* number of tree levels */
#define BRORDER	6			/* branches/level */

/* Tree branch structure for quick occupancy search */
/* with 9 levels & 6 branches per level, we can store 1.94 Gbits in 259 MBytes (4.5% overhead) */
const struct {
	long	capacity;		/* slots/branch this level */
	long	skip_bytes;		/* bytes until next branch */
	int	cntr_siz;		/* occupancy counter size */
} tree_br[NLEVELS] = {
	{248L, 32L, 1},
	{248L*6, 32L*6+2, 2},
	{248L*6*6, (32L*6+2)*6+2, 2},
	{248L*6*6*6, ((32L*6+2)*6+2)*6+2, 2},
	{248L*6*6*6*6, (((32L*6+2)*6+2)*6+2)*6+3, 3},
	{248L*6*6*6*6*6, ((((32L*6+2)*6+2)*6+2)*6+3)*6+3, 3},
	{248L*6*6*6*6*6*6, (((((32L*6+2)*6+2)*6+2)*6+3)*6+3)*6+3, 3},
	{248L*6*6*6*6*6*6*6, ((((((32L*6+2)*6+2)*6+2)*6+3)*6+3)*6+3)*6+4, 4},
	{248L*6*6*6*6*6*6*6*6, (((((((32L*6+2)*6+2)*6+2)*6+3)*6+3)*6+3)*6+4)*6+4, 4},
};

char  buf[256];				/* buffer for ordered array output */


/* Encode integer in string and return pointer to end */
static char *
tack(char *b, long i)
{
	char  *cp;
	char  *res;

	*b++ = '\t';
	cp = b;
	if (i == 0)
		*cp++ = '0';
	else
		do {
			*cp++ = i%10L + '0';
			i /= 10L;
		} while (i);
	res = cp--;
#define c i
	while (cp > b) {		/* reverse string */
		c = *cp;
		*cp-- = *b;
		*b++ = c;
	}
#undef c
	return(res);
}


/* Loop over dimensions, spitting out buffer after each increment */
static void
loop(long *n, char *b)
{
	long  i;

	if (n[0] == 0) {
		*b = '\0';
		puts(buf);
		return;
	}
	for (i = 0; i < n[0]; i++)
		loop(n+1, tack(b, i));
}


/* Print out shuffled value */
static void
print_shuf(long *n, long aval)
{
	int	i;

	for (i = 0; n[i+1]; i++) {
		printf("\t%ld", aval % n[i]);
		aval /= n[i];
	}
	printf("\t%ld\n", aval);
}


/* Allocate and prepare occupancy tree */
static uby8 *
tree_alloc(long alen)
{
	uby8  *troot;
	double  bytes_per_bit;
	int  i;
	int  ht = 0;
					/* how tall does our tree need to be? */
	while (tree_br[ht].capacity*BRORDER < alen)
		if (++ht >= NLEVELS) {
			fputs("Array too large to shuffle\n", stderr);
			exit(1);
		}
	bytes_per_bit = 1.;		/* figure out tree size (with overhead) */
	for (i = ht; i >= 0; i--)
		bytes_per_bit += (double)tree_br[i].cntr_siz;
	bytes_per_bit += (double)tree_br[ht].skip_bytes;
	bytes_per_bit /= (double)tree_br[ht].capacity;
	troot = (uby8 *)calloc((long)(alen*bytes_per_bit)+2, 1);
	if (troot == NULL) {
		fputs("Not enough memory for shuffle\n", stderr);
		exit(1);
	}
	*troot = ht;			/* first byte is tree height */
	for (i = 256; i--; ) {		/* assign 0-bit count table */
		int	b;
		buf[i] = 8;
		for (b = i; b; b >>= 1)
			buf[i] -= b&1;
	}
	return(troot);
}


/* Get number of slots available at this branch location */
static long
get_avail(const uby8 *ctrp, int lvl)
{
	long	cnt = 0;
	int	n = tree_br[lvl].cntr_siz;

	while (--n > 0) {		/* LSB first */
		cnt |= ctrp[n];
		cnt <<= 8;
	}
	cnt |= ctrp[0];

	return(tree_br[lvl].capacity - cnt);
}


/* Increment branch occupancy counter */
static void
incr_counter(uby8 *ctrp, int n)
{
	n = tree_br[n].cntr_siz;

	while (n-- > 0)			/* LSB first */
		if (++(*ctrp++))
			break;
}


/* Skip to and allocate a leaf from tree */
static long
eat_nth_leaf(uby8 *brp, long ski)
{
	int	lvl = *brp++;			/* tree height in first byte */
	long	pos = 0;
	int	b;

	while (lvl >= 0) {			/* descend to leaves */
		long  navail;
		b = 0;				/* select each branch */
		while (ski >= (navail = get_avail(brp, lvl))) {
			if (++b >= BRORDER) {
				fputs("Shuffle tree error!\n", stderr);
				exit(1);
			}
			pos += tree_br[lvl].capacity;
			ski -= navail;
			brp += tree_br[lvl].skip_bytes;
		}
		incr_counter(brp, lvl);		/* we intend to eat one */
		brp += tree_br[lvl--].cntr_siz;	/* drop a level */
	}
	while (ski >= buf[*brp]) {		/* browse the leaves */
		pos += 8;
		ski -= buf[*brp++];		/* buf contains 0-bit counts */
	}
	b = 0;					/* find target bit in byte */
	while ((ski -= !(*brp & 1<<b)) >= 0) {
		pos++;
		b++;
	}
	*brp |= 1<<b;				/* eat it */
	return(pos);				/* & return leaf's slot# */
}


/* Shuffle all possible output strings and spit out randomly (tree version) */
static void
big_shuffle(long *n, long alen)
{
	uby8  *tree_root;
					/* size and allocate holder tree */
	tree_root = tree_alloc(alen);

	while (alen > 0)		/* allocate and print random array entries */
		print_shuf(n, eat_nth_leaf(tree_root, random() % alen--));

	free(tree_root);		/* all done */
}


/* Shuffle all possible output strings and spit out randomly */
static void
shuffle(long *n)
{
	long  alen;
	uint16  *myshuf;
	int  i;

	alen = 1;		/* compute shuffle size */
	for (i = 0; n[i]; i++) {
		if (alen*n[i] <= alen) {
			fputs("Array too large to count!\n", stderr);
			exit(1);
		}
		alen *= n[i];
	}
				/* get unique starting point */
	srandom((long)time(0));

	if (alen > 1L<<16) {	/* use large shuffle method? */
		big_shuffle(n, alen);
		return;
	}
	myshuf = (uint16 *)malloc(alen*sizeof(uint16));
	if (myshuf == NULL) {
		fputs("Insufficient memory for shuffle\n", stderr);
		exit(1);
	}
	for (i = alen; i--; )	/* initialize in any order */
		myshuf[i] = i;
				/* perform Fisher-Yates shuffle */
	for (i = 0; i < alen-1; i++) {
		int	ix = random()%(alen-i) + i;
		int	ndx = myshuf[i];
		myshuf[i] = myshuf[ix];
		myshuf[ix] = ndx;
	}
				/* put randomly indexed output */
	for (i = alen; i--; )
		print_shuf(n, (long)myshuf[i]);

	free(myshuf);		/* all done */
}


int
main(int argc, char *argv[])
{
	char  *prog = argv[0];
	int  doshuffle = 0;
	long  n[MAXDIM];
	int  a;

	argv++; argc--;
	if (argc <= 0)
		goto userr;
	if (argv[0][0] == '-' && argv[0][1] == 's') {
		doshuffle = 1;
		argv++; argc--;
	}
	for (a = 0; a < argc; a++)
		if ((n[a] = atol(argv[a])) <= 1)
			goto userr;
	n[a] = 0;
	if (!a)
		goto userr;

	if (doshuffle)
		shuffle(n);
	else
		loop(n, buf);

	return(0);
userr:
	fputs("Usage: ", stderr);
	fputs(prog, stderr);
	fputs(" [-s] N0 [N1 ..]\n", stderr);
	return(1);
}
