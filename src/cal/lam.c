#ifndef lint
static const char	RCSid[] = "$Id: lam.c,v 1.6 2005/06/13 22:40:47 greg Exp $";
#endif
/*
 *  lam.c - simple program to laminate files.
 *
 *	7/14/88		Greg Ward
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "platform.h"
#include "rtprocess.h"

#define MAXFILE		32		/* maximum number of files */

#define MAXLINE		512		/* maximum input line */

FILE	*input[MAXFILE];
char	*tabc[MAXFILE];
int	nfiles;

char	buf[MAXLINE];

int
main(argc, argv)
int	argc;
char	*argv[];
{
	register int	i;
	char	*curtab;
	int	running;

	curtab = "\t";
	nfiles = 0;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 't':
				curtab = argv[i]+2;
				break;
			case '\0':
				tabc[nfiles] = curtab;
				input[nfiles++] = stdin;
				break;
			default:
				fputs(argv[0], stderr);
				fputs(": bad option\n", stderr);
				exit(1);
			}
		} else if (argv[i][0] == '!') {
			tabc[nfiles] = curtab;
			if ((input[nfiles++] = popen(argv[i]+1, "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot start command\n", stderr);
				exit(1);
			}
		} else {
			tabc[nfiles] = curtab;
			if ((input[nfiles++] = fopen(argv[i], "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot open file\n", stderr);
				exit(1);
			}
		}
		if (nfiles >= MAXFILE) {
			fputs(argv[0], stderr);
			fputs(": too many input streams\n", stderr);
			exit(1);
		}
	}

	do {
		running = 0;
		for (i = 0; i < nfiles; i++) {
			if (fgets(buf, MAXLINE, input[i]) != NULL) {
				if (i)
					fputs(tabc[i], stdout);
				buf[strlen(buf)-1] = '\0';
				fputs(buf, stdout);
				running++;
			}
		}
		if (running)
			putchar('\n');
	} while (running);

	exit(0);
}
