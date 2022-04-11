#ifndef lint
static const char	RCSid[] = "$Id: cnt.c,v 1.3 2022/04/11 18:08:19 greg Exp $";
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

#define	MAXDIM	50

char  buf[256];


/* Encode integer in string and return pointer to end */
static char *
tack(
char  *b,
int  i
)
{
	char  *cp;
	char  *res;

	*b++ = '\t';
	cp = b;
	if (i == 0)
		*cp++ = '0';
	else
		do {
			*cp++ = i%10 + '0';
			i /= 10;
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
loop(
int  *n,
char  *b
)
{
	int  i;

	if (n[0] == 0) {
		*b = '\0';
		puts(buf);
		return;
	}
	for (i = 0; i < n[0]; i++)
		loop(n+1, tack(b, i));
}


/* Shuffle all possible output strings and spit out randomly */
static void
shuffle(
int *n
)
{
	int  sub[MAXDIM];
	int  ndim;
	int  alen;
	int  *myshuf;
	int  i, j;

	alen = 1;		/* allocate shuffle index array */
	for (j = 0; n[j]; j++)
		if ((alen *= n[j]) < 0)
			exit(1);

	myshuf = (int *)malloc(alen*sizeof(int));
	if (myshuf == NULL) {
		fputs("Insufficient memory for shuffle\n", stderr);
		exit(1);
	}
	for (i = alen; i--; )	/* initialize in any order */
		myshuf[i] = i;
				/* get unique starting point */
	srandom((long)time(0));
				/* perform Fisher-Yates shuffle */
	for (i = 0; i < alen-1; i++) {
		int	ix = random()%(alen-i) + i;
		int	ndx = myshuf[i];
		myshuf[i] = myshuf[ix];
		myshuf[ix] = ndx;
	}
				/* put randomly indexed output */
	for (i = alen; i--; ) {
		int	aval = myshuf[i];
		for (j = 0; n[j+1]; j++) {
			printf("\t%d", aval % n[j]);
			aval /= n[j];
		}
		printf("\t%d\n", aval);
	}
	free(myshuf);
}


int
main(
int  argc,
char  *argv[]
)
{
	char  *prog = argv[0];
	int  doshuffle = 0;
	int  n[MAXDIM];
	int  a;

	argv++; argc--;
	if (argc <= 0)
		goto userr;
	if (argv[0][0] == '-' && argv[0][1] == 's') {
		doshuffle = 1;
		argv++; argc--;
	}
	for (a = 0; a < argc; a++)
		if ((n[a] = atoi(argv[a])) <= 1)
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
