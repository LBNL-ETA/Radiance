#ifndef lint
static const char	RCSid[] = "$Id: cnt.c,v 1.2 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 *  cnt.c - simple counting program.
 *
 *	2/1/88
 */

#include  <stdlib.h>
#include  <stdio.h>


int  n[50];
char  buf[256];

static void loop(int *n, char *b);

int
main(
int  argc,
char  *argv[]
)
{
	int  a;

	argv++; argc--;
	for (a = 0; a < argc; a++)
		n[a] = atoi(argv[a]);
	n[a] = 0;
	loop(n, buf);

	exit(0);
}


char *
tack(
register char  *b,
register int  i
)
{
	register char  *cp;
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
