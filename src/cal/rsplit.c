#ifndef lint
static const char	RCSid[] = "$Id: rsplit.c,v 1.2 2019/07/05 00:46:23 greg Exp $";
#endif
/*
 *  rsplit.c - split input into multiple output streams
 *
 *	7/4/19		Greg Ward
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "platform.h"
#include "resolu.h"
#include "rtio.h"

#define	DOHEADER	1
#define DORESOLU	2

#define MAXFILE		512		/* maximum number of files */

FILE		*output[MAXFILE];
int		bytsiz[MAXFILE];
short		hdrflags[MAXFILE];
const char	*format[MAXFILE];
int		termc[MAXFILE];
int		nfiles = 0;

int		outheader = 0;		/* output header to each stream? */

RESOLU		ourres = {PIXSTANDARD, 0, 0};

char		buf[16384];


/* process header line */
static int
headline(char *s, void *p)
{
	extern const char	FMTSTR[];
	int			i = nfiles;

	if (strstr(s, FMTSTR) == s)
		return(0);		/* don't copy format */
	if (!strncmp(s, "NROWS=", 6))
		return(0);
	if (!strncmp(s, "NCOLS=", 6))
		return(0);
	if (!strncmp(s, "NCOMP=", 6))
		return(0);
	while (i--)			/* else copy line to output streams */
		if (hdrflags[i] & DOHEADER)
			fputs(s, output[i]);
	return(1);
}


/* scan field into buffer up to and including terminating byte */
static int
scanOK(int termc)
{
	char	*cp = buf;
	int	c;

	while ((c = getchar()) != EOF) {
		*cp++ = c;
		if (cp-buf >= sizeof(buf))
			break;
		if (c == termc) {
			*cp = '\0';
			return(cp-buf);
		}
	}
	return(0);
}


int
main(int argc, char *argv[])
{
	int		inpheader = 0;
	int		inpres = 0;
	int		outres = 0;
	int		force = 0;
	int		append = 0;
	long		outcnt = 0;
	int		bininp = 0;
	int		curterm = '\n';
	int		curbytes = 0;
	int		curflags = 0;
	const char	*curfmt = "ascii";
	int	i;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 't':
				curterm = argv[i][2];
				break;
			case 'i':
				switch (argv[i][2]) {
				case 'h':
					inpheader = !inpheader;
					break;
				case 'H':
					inpres = !inpres;
					break;
				default:
					goto badopt;
				}
				break;
			case 'f':
				++force;
				append = 0;
				break;
			case 'a':
				if (outheader | outres) {
					fputs(argv[0], stderr);
					fputs(": -a option incompatible with -oh and -oH\n",
							stderr);
					return(1);
				}
				append = 1;
				break;
			case 'x':
				ourres.xr = atoi(argv[++i]);
				outres = 1;
				break;
			case 'y':
				ourres.yr = atoi(argv[++i]);
				outres = 1;
				break;
			case 'o':
				switch (argv[i][2]) {
				case 'n':
					outcnt = atol(argv[++i]);
					continue;
				case 'h':
					curflags ^= DOHEADER;
					continue;
				case 'H':
					curflags ^= DORESOLU;
					continue;
				case 'f':
					curfmt = "float";
					curbytes = sizeof(float);
					break;
				case 'd':
					curfmt = "double";
					curbytes = sizeof(double);
					break;
				case 'i':
					curfmt = "int";
					curbytes = sizeof(int);
					break;
				case 'w':
					curfmt = "16-bit";
					curbytes = 2;
					break;
				case 'b':
					curfmt = "byte";
					curbytes = 1;
					break;
				case 'a':
					curfmt = "ascii";
					curbytes = argv[i][3] ? -1 : 0;
					break;
				default:
					goto badopt;
				}
				if (isdigit(argv[i][3]))
					curbytes *= atoi(argv[i]+3);
				curbytes += (curbytes == -1);
				if (curbytes > (int)sizeof(buf)) {
					fputs(argv[0], stderr);
					fputs(": output size too big\n", stderr);
					return(1);
				}
				if (curbytes > 0) {
					curterm = '\0';
					++bininp;
				}
				break;
			case '\0':
				outres |= (curflags & DORESOLU);
				termc[nfiles] = curterm;
				hdrflags[nfiles] = curflags;
				output[nfiles] = stdout;
				if (curbytes > 0)
					SET_FILE_BINARY(output[nfiles]);
				format[nfiles] = curfmt;
				bytsiz[nfiles++] = curbytes;
				break;
			badopt:;
			default:
				fputs(argv[0], stderr);
				fputs(": bad option\n", stderr);
				return(1);
			}
		} else if (argv[i][0] == '!') {
			outres |= (curflags & DORESOLU);
			termc[nfiles] = curterm;
			hdrflags[nfiles] = curflags;
			if ((output[nfiles] = popen(argv[i]+1, "w")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot start command\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(output[nfiles]);
			format[nfiles] = curfmt;
			bytsiz[nfiles++] = curbytes;
		} else {
			outres |= (curflags & DORESOLU);
			termc[nfiles] = curterm;
			hdrflags[nfiles] = curflags;
			if (!append & !force && access(argv[i], F_OK) == 0) {
				fputs(argv[i], stderr);
				fputs(": file exists -- use -f to overwrite\n",
						stderr);
				return(1);
			}
			output[nfiles] = fopen(argv[i], append ? "a" : "w");
			if (output[nfiles] == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot open for output\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(output[nfiles]);
			format[nfiles] = curfmt;
			bytsiz[nfiles++] = curbytes;
		}
		if (nfiles >= MAXFILE) {
			fputs(argv[0], stderr);
			fputs(": too many output streams\n", stderr);
			return(1);
		}
	}
	if (!nfiles) {
		fputs(argv[0], stderr);
		fputs(": no output streams\n", stderr);
		return(1);
	}
	if (bininp)				/* binary input? */
		SET_FILE_BINARY(stdin);
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
	flockfile(stdin);
	for (i = nfiles; i--; )
		flockfile(output[i]);
#endif
						/* load/copy header */
	if (inpheader && getheader(stdin, headline, NULL) < 0) {
		fputs(argv[0], stderr);
		fputs(": cannot get header from standard input\n",
				stderr);
		return(1);
	}
						/* handle resolution string */
	if (inpres && !fgetsresolu(&ourres, stdin)) {
		fputs(argv[0], stderr);
		fputs(": cannot get resolution string from standard input\n",
				stderr);
		return(1);
	}
	if (outres && (ourres.xr <= 0) | (ourres.yr <= 0)) {
		fputs(argv[0], stderr);
		fputs(": -oH option requires -iH or -x and -y options\n", stderr);
		return(1);
	}
	if ((ourres.xr > 0) & (ourres.yr > 0)) {
		if (outcnt <= 0) {
			outcnt = ourres.xr * ourres.yr;
		} else if (outcnt != ourres.xr*ourres.yr) {
			fputs(argv[0], stderr);
			fputs(": warning: -on option does not agree with resolution\n",
					stderr);
		}
	}
	for (i = 0; i < nfiles; i++) {		/* complete headers */
		if (hdrflags[i] & DOHEADER) {
			if (!inpheader)
				newheader("RADIANCE", output[i]);
			printargs(argc, argv, output[i]);
			if (format[i] != NULL)
				fputformat(format[i], output[i]);
			fputc('\n', output[i]);
		}
		if (hdrflags[i] & DORESOLU)
			fputsresolu(&ourres, output[i]);
	}
	do {					/* main loop */
		for (i = 0; i < nfiles; i++) {
			if (bytsiz[i] > 0) {		/* binary output */
				if (getbinary(buf, bytsiz[i], 1, stdin) < 1)
					break;
				putbinary(buf, bytsiz[i], 1, output[i]);
			} else if (bytsiz[i] < 0) {	/* N-field output */
				int	n = -bytsiz[i];
				while (n--) {
					if (!scanOK(termc[i]))
						break;
					fputs(buf, output[i]);
				}
				if (n >= 0)		/* fell short? */
					break;
			} else {			/* 1-field output */
				if (!scanOK(termc[i]))
					break;
				fputs(buf, output[i]);
			}
		}
		if (i < nfiles)
			break;
	} while (--outcnt);
							/* check ending */
	if (outcnt > 0) {
		fputs(argv[0], stderr);
		fputs(": warning: premature EOD\n", stderr);
	}
	return(0);
}
