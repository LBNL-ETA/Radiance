#ifndef lint
static const char	RCSid[] = "$Id: neat.c,v 1.3 2003/07/27 22:12:01 schorsch Exp $";
#endif
/*
 *  neat.c - program to tidy up columns.
 *
 *	10/24/86
 */

#include  <stdlib.h>
#include  <stdio.h>

#include  <ctype.h>


char  *format = "8.8";			/* default format */


main(argc, argv)
int  argc;
char  *argv[];
{
	int  left, anchor, right;
	char  buf[512];
	register char  *cp, *word;
	register int  i;

	if (argc == 2)
		format = argv[1];
	else if (argc > 2)
		goto userror;
	left = 0;
	for (cp = format; isdigit(*cp); cp++)
		left = left*10 + *cp - '0';
	right = 0;
	if ( (anchor = *cp) )
		for (cp++; isdigit(*cp); cp++)
			right = right*10 + *cp - '0';
	if (*cp)
		goto userror;
	
	while ((cp = fgets(buf, sizeof(buf), stdin)) != NULL)
		for ( ; ; ) {
			while (isspace(*cp))
				cp++;
			if (!*cp) {
				putchar('\n');
				break;
			}
			word = cp;
			while (*cp && *cp != anchor && !isspace(*cp))
				cp++;
			i = left-(cp-word);
			while (i-- > 0)
				putchar(' ');
			while (word < cp)
				putchar(*word++);
			i = right+1;
			if (*cp == anchor)
				do {
					putchar(*cp++);
					i--;
				} while (*cp && !isspace(*cp));
			while (i-- > 0)
				putchar(' ');
		}
	exit(0);
userror:
	fputs("Usage: ", stderr);
	fputs(argv[0], stderr);
	fputs(" [format]\n", stderr);
	exit(1);
}
