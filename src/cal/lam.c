#ifndef lint
static const char	RCSid[] = "$Id: lam.c,v 1.26 2024/06/07 18:19:22 greg Exp $";
#endif
/*
 *  lam.c - simple program to laminate files.
 *
 *	7/14/88		Greg Ward
 */

#include <ctype.h>

#include "rtio.h"
#include "platform.h"
#include "paths.h"

#define MAXLINE		262144		/* maximum input line */

struct instream {		/* structure to hold input stream info */
	FILE		*input;		/* input stream */
	int		bytsiz;		/* bytes/component if binary */
	char		*tabc;		/* data separation string */
}	*rifile = NULL;

int	nfiles = 0;

char	buf[MAXLINE];

int
main(int argc, char *argv[])
{
	long	incnt = 0;
	int	unbuff = 0;
	int	binout = 0;
	char	*curtab = "\t";
	int	curbytes = 0;
	int	puteol;
	int	i;

	rifile = (struct instream *)calloc(argc-1, sizeof(struct instream));
	if (!rifile) {
		fputs(argv[0], stderr);
		fputs(": not enough memory\n", stderr);
		return(1);
	}
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 't':
				curtab = argv[i]+2;
				if (!*curtab) curtab = "\n";
				break;
			case 'u':
				unbuff = !unbuff;
				break;
			case 'i':
				switch (argv[i][2]) {
				case 'n':
					incnt = atol(argv[++i]);
					break;
				case 'f':
				case 'F':
					curbytes = sizeof(float);
					break;
				case 'd':
				case 'D':
					curbytes = sizeof(double);
					break;
				case 'i':
				case 'I':
					curbytes = sizeof(int);
					break;
				case 'w':
				case 'W':
					curbytes = 2;
					break;
				case 'b':
					curbytes = 1;
					break;
				case 'a':
					curbytes = argv[i][3] ? -1 : 0;
					break;
				default:
					goto badopt;
				}
				if (isdigit(argv[i][3]))
					curbytes *= atoi(argv[i]+3);
				curbytes += (curbytes == -1);
				if (curbytes > MAXLINE) {
					fputs(argv[0], stderr);
					fputs(": input size too big\n", stderr);
					return(1);
				}
				if (curbytes > 0) {
					curtab = "";
					++binout;
				}
				break;
			case '\0':
				rifile[nfiles].tabc = curtab;
				rifile[nfiles].input = stdin;
				if (curbytes > 0)
					SET_FILE_BINARY(rifile[nfiles].input);
				rifile[nfiles++].bytsiz = curbytes;
				break;
			badopt:;
			default:
				fputs(argv[0], stderr);
				fputs(": bad option\n", stderr);
				return(1);
			}
		} else if (argv[i][0] == '!') {
			rifile[nfiles].tabc = curtab;
			if ((rifile[nfiles].input = popen(argv[i]+1, "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot start command\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(rifile[nfiles].input);
			rifile[nfiles++].bytsiz = curbytes;
		} else {
			rifile[nfiles].tabc = curtab;
			if ((rifile[nfiles].input = fopen(argv[i], "r")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot open file\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(rifile[nfiles].input);
			rifile[nfiles++].bytsiz = curbytes;
		}
	}
	if (!nfiles) {
		fputs(argv[0], stderr);
		fputs(": no input streams\n", stderr);
		return(1);
	}					/* reduce array to size we need */
	rifile = (struct instream *)realloc(rifile, nfiles*sizeof(struct instream));
	if (!rifile) {
		fputs(argv[0], stderr);
		fputs(": realloc() failed!\n", stderr);
		return(1);
	}
	if (binout)				/* binary output? */
		SET_FILE_BINARY(stdout);
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
	for (i = nfiles; i--; )
		flockfile(rifile[i].input);
	flockfile(stdout);
#endif
	puteol = 0;				/* any ASCII output at all? */
	for (i = nfiles; i--; )
		puteol += (rifile[i].bytsiz <= 0);
	do {					/* main loop */
		for (i = 0; i < nfiles; i++) {
			if (rifile[i].bytsiz > 0) {		/* binary input */
				if (getbinary(buf, rifile[i].bytsiz, 1, rifile[i].input) < 1)
					break;
				if (putbinary(buf, rifile[i].bytsiz, 1, stdout) != 1)
					break;
			} else if (rifile[i].bytsiz < 0) {	/* multi-line input */
				int	n = -rifile[i].bytsiz;
				while (n--) {
					if (fgets(buf, MAXLINE, rifile[i].input) == NULL)
						break;
					if ((i > 0) | (n < -rifile[i].bytsiz-1))
						fputs(rifile[i].tabc, stdout);
					buf[strlen(buf)-1] = '\0';
					if (fputs(buf, stdout) == EOF)
						break;
				}
				if (n >= 0)		/* fell short? */
					break;
			} else {			/* single-line input */
				if (fgets(buf, MAXLINE, rifile[i].input) == NULL)
					break;
				if (i)
					fputs(rifile[i].tabc, stdout);
				buf[strlen(buf)-1] = '\0';
				if (fputs(buf, stdout) == EOF)
					break;
			}
		}
		if (i < nfiles)
			break;
		if (puteol)
			putchar('\n');
		if (unbuff)
			fflush(stdout);
	} while (--incnt);
							/* check ending */
	if (fflush(stdout) == EOF) {
		fputs(argv[0], stderr);
		fputs(": write error on standard output\n", stderr);
		return(1);
	}
	if (incnt > 0) {
		fputs(argv[0], stderr);
		fputs(": warning: premature EOD\n", stderr);
	}
	return(0);
}
