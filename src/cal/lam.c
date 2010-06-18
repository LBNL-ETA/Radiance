#ifndef lint
static const char	RCSid[] = "$Id: lam.c,v 1.11 2010/06/18 21:22:49 greg Exp $";
#endif
/*
 *  lam.c - simple program to laminate files.
 *
 *	7/14/88		Greg Ward
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "platform.h"
#include "rtprocess.h"

#define MAXFILE		512		/* maximum number of files */

#define MAXLINE		4096		/* maximum input line */

FILE	*input[MAXFILE];
int	bytsiz[MAXFILE];
char	*tabc[MAXFILE];
int	nfiles;

char	buf[MAXLINE];

int
main(argc, argv)
int	argc;
char	*argv[];
{
	int	unbuff = 0;
	int	i;
	char	*curtab;
	int	curbytes;
	int	puteol;

	curtab = "\t";
	curbytes = 0;
	nfiles = 0;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 't':
				curtab = argv[i]+2;
				break;
			case 'u':
				unbuff = !unbuff;
				break;
			case 'i':
				switch (argv[i][2]) {
				case 'f':
					curbytes = sizeof(float);
					break;
				case 'd':
					curbytes = sizeof(double);
					break;
				case 'i':
					curbytes = sizeof(int);
					break;
				case 'w':
					curbytes = 2;
					break;
				case 'b':
					curbytes = 1;
					break;
				case 'a':
					curbytes = argv[i][3] ? 1 : 0;
					break;
				default:
					goto badopt;
				}
				if (isdigit(argv[i][3]))
					curbytes *= atoi(argv[i]+3);
				if (curbytes < 0 || curbytes > MAXLINE) {
					fputs(argv[0], stderr);
					fputs(": illegal input size\n", stderr);
					exit(1);
				}
				if (curbytes)
					curtab = "";
				break;
			case '\0':
				tabc[nfiles] = curtab;
				bytsiz[nfiles] = curbytes;
				input[nfiles++] = stdin;
				break;
			badopt:;
			default:
				fputs(argv[0], stderr);
				fputs(": bad option\n", stderr);
				exit(1);
			}
		} else if (argv[i][0] == '!') {
			tabc[nfiles] = curtab;
			bytsiz[nfiles] = curbytes;
			if ((input[nfiles] = popen(argv[i]+1, "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot start command\n", stderr);
				exit(1);
			}
			if (bytsiz[nfiles])
				SET_FILE_BINARY(input[nfiles]);
			++nfiles;
		} else {
			tabc[nfiles] = curtab;
			bytsiz[nfiles] = curbytes;
			if ((input[nfiles] = fopen(argv[i], "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot open file\n", stderr);
				exit(1);
			}
			if (bytsiz[nfiles])
				SET_FILE_BINARY(input[nfiles]);
			++nfiles;
		}
		if (nfiles >= MAXFILE) {
			fputs(argv[0], stderr);
			fputs(": too many input streams\n", stderr);
			exit(1);
		}
	}
	puteol = 0;				/* check for ASCII output */
	for (i = nfiles; i--; )
		if (!bytsiz[i] || isprint(tabc[i][0]) || tabc[i][0] == '\t') {
			puteol++;
			break;
		}
	for ( ; ; ) {				/* main loop */
		for (i = 0; i < nfiles; i++) {
			if (bytsiz[i]) {		/* binary file */
				if (fread(buf, bytsiz[i], 1, input[i]) < 1)
					break;
				if (i)
					fputs(tabc[i], stdout);
				fwrite(buf, bytsiz[i], 1, stdout);
			} else {
				if (fgets(buf, MAXLINE, input[i]) == NULL)
					break;
				if (i)
					fputs(tabc[i], stdout);
				buf[strlen(buf)-1] = '\0';
				fputs(buf, stdout);
			}
		}
		if (i < nfiles)
			break;
		if (puteol)
			putchar('\n');
		if (unbuff)
			fflush(stdout);
	}
	return(0);
}
