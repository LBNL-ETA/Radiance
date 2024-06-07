#ifndef lint
static const char	RCSid[] = "$Id: rsplit.c,v 1.17 2024/06/07 20:26:53 greg Exp $";
#endif
/*
 *  rsplit.c - split input into multiple output streams
 *
 *	7/4/19		Greg Ward
 */

#include <ctype.h>

#include "rtio.h"
#include "platform.h"
#include "paths.h"
#include "resolu.h"

#define	DOHEADER	1
#define DORESOLU	2

int			swapped = 0;	/* input is byte-swapped */

struct outstream {		/* structure to hold output stream info */
	const char	*outspec;	/* output specification */
	FILE		*output;	/* output stream */
	int		ncomp;		/* component count */
	int		bytsiz;		/* bytes/component if binary */
	int		hdrflags;	/* header output flags */
	int		termc;		/* data separation character */
	const char	*format;	/* data format */
}	*rofile = NULL;

int		nfiles = 0;	/* output file count */

RESOLU		ourres = {PIXSTANDARD, 0, 0};

char		buf[16384];		/* input buffer used in scanOK() */


/* process header line */
int
headline(char *s, void *p)
{
	extern const char	FMTSTR[];
	int			i;

	if (strstr(s, FMTSTR) == s)
		return(0);		/* don't copy format */
	if (!strncmp(s, "NROWS=", 6))
		return(0);
	if (!strncmp(s, "NCOLS=", 6))
		return(0);
	if (!strncmp(s, "NCOMP=", 6))
		return(0);
	if ((i = isbigendian(s)) >= 0) {
		swapped = (nativebigendian() != i);
		return(0);
	}
	i = nfiles;
	while (i--)			/* else copy line to output streams */
		if (rofile[i].hdrflags & DOHEADER)
			fputs(s, rofile[i].output);
	return(1);
}


/* scan field into buffer up to and including terminating byte */
int
scanOK(int termc)
{
	int	skip_white = (termc == ' ');
	char	*cp = buf;
	int	c;

	while ((c = getchar()) != EOF) {
		if (skip_white && isspace(c))
			continue;
		skip_white = 0;
		if (c == '\n' && isspace(termc))
			c = termc;	/* forgiving assumption */
		*cp++ = c;
		if (cp-buf >= sizeof(buf))
			break;
		if ((termc == ' ') ? isspace(c) : (c == termc)) {
			*cp = '\0';
			return(cp-buf);
		}
	}
	return(0);
}


int
main(int argc, char *argv[])
{
	int		inpflags = 0;
	int		needres = 0;
	int		force = 0;
	int		append = 0;
	long		outcnt = 0;
	int		bininp = 0;
	int		nstdoutcomp = 0;
	int		curterm = '\n';
	int		curncomp = 1;
	int		curbytes = 0;
	int		curflags = 0;
	const char	*curfmt = "ascii";
	int		i;

	rofile = (struct outstream *)calloc(argc-1, sizeof(struct outstream));
	if (!rofile) {
		fputs(argv[0], stderr);
		fputs(": not enough memory\n", stderr);
		return(1);
	}
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 't':
				curterm = argv[i][2];
				if (!curterm) curterm = '\n';
				break;
			case 'i':
				switch (argv[i][2]) {
				case 'h':
					inpflags ^= DOHEADER;
					break;
				case 'H':
					inpflags ^= DORESOLU;
					break;
				default:
					goto badopt;
				}
				break;
			case 'f':
				force = !force;
				break;
			case 'a':
				append = !append;
				break;
			case 'x':
				ourres.xr = atoi(argv[++i]);
				break;
			case 'y':
				ourres.yr = atoi(argv[++i]);
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
					curbytes = 0;
					break;
				default:
					goto badopt;
				}
				curncomp = isdigit(argv[i][3]) ?
							atoi(argv[i]+3) : 1 ;
				if (curbytes*curncomp > (int)sizeof(buf)) {
					fputs(argv[0], stderr);
					fputs(": output size too big\n", stderr);
					return(1);
				}
				bininp += (curbytes > 0);
				break;
			case '\0':			/* stdout */
				if (!nstdoutcomp) {	/* first use? */
					needres |= (curflags & DORESOLU);
					rofile[nfiles].hdrflags = curflags;
				}
				rofile[nfiles].output = stdout;
				if (curbytes > 0)
					SET_FILE_BINARY(rofile[nfiles].output);
				rofile[nfiles].termc = curterm;
				rofile[nfiles].format = curfmt;
				nstdoutcomp +=
					rofile[nfiles].ncomp = curncomp;
				rofile[nfiles].bytsiz = curbytes;
				rofile[nfiles++].outspec = argv[i];
				break;
			badopt:;
			default:
				fputs(argv[0], stderr);
				fputs(": bad option\n", stderr);
				return(1);
			}
		} else if (argv[i][0] == '.' && !argv[i][1]) {
			rofile[nfiles].output = NULL;		/* discard data */
			rofile[nfiles].termc = curterm;
			rofile[nfiles].format = curfmt;
			rofile[nfiles].ncomp = curncomp;
			rofile[nfiles].bytsiz = curbytes;
			rofile[nfiles++].outspec = argv[i];
		} else if (argv[i][0] == '!') {
			needres |= (curflags & DORESOLU);
			rofile[nfiles].hdrflags = curflags;
			rofile[nfiles].termc = curterm;
			if ((rofile[nfiles].output = popen(argv[i]+1, "w")) == NULL) {
				fputs(argv[i], stderr);
				fputs(": cannot start command\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(rofile[nfiles].output);
			rofile[nfiles].format = curfmt;
			rofile[nfiles].ncomp = curncomp;
			rofile[nfiles].bytsiz = curbytes;
			rofile[nfiles++].outspec = argv[i];
		} else {
			int	j = nfiles;
			while (j--)			/* check duplicates */
				if (!strcmp(argv[i], rofile[j].outspec)) {
					fputs(argv[0], stderr);
					fputs(": duplicate output: ", stderr);
					fputs(argv[i], stderr);
					fputc('\n', stderr);
					return(1);
				}
			if (append & (curflags != 0)) {
				fputs(argv[0], stderr);
				fputs(": -a option incompatible with -oh and -oH\n",
						stderr);
				return(1);
			}
			needres |= (curflags & DORESOLU);
			rofile[nfiles].hdrflags = curflags;
			rofile[nfiles].termc = curterm;
			if (!append & !force && access(argv[i], F_OK) == 0) {
				fputs(argv[i], stderr);
				fputs(": file exists -- use -f to overwrite\n",
						stderr);
				return(1);
			}
			rofile[nfiles].output = fopen(argv[i], append ? "a" : "w");
			if (!rofile[nfiles].output) {
				fputs(argv[i], stderr);
				fputs(": cannot open for output\n", stderr);
				return(1);
			}
			if (curbytes > 0)
				SET_FILE_BINARY(rofile[nfiles].output);
			rofile[nfiles].format = curfmt;
			rofile[nfiles].ncomp = curncomp;
			rofile[nfiles].bytsiz = curbytes;
			rofile[nfiles++].outspec = argv[i];
		}
	}
	if (!nfiles) {
		fputs(argv[0], stderr);
		fputs(": no output streams\n", stderr);
		return(1);
	}					/* reduce array to size we need */
	rofile = (struct outstream *)realloc(rofile, nfiles*sizeof(struct outstream));
	if (!rofile) {
		fputs(argv[0], stderr);
		fputs(": realloc() failed!\n", stderr);
		return(1);
	}
	if (bininp)				/* binary input? */
		SET_FILE_BINARY(stdin);
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
	flockfile(stdin);
	for (i = nfiles; i--; )
		if (rofile[i].output != NULL)
			ftrylockfile(rofile[i].output);
#endif
						/* load/copy header */
	if (inpflags & DOHEADER && getheader(stdin, headline, NULL) < 0) {
		fputs(argv[0], stderr);
		fputs(": cannot read header from standard input\n",
				stderr);
		return(1);
	}
						/* handle resolution string */
	if (inpflags & DORESOLU && !fgetsresolu(&ourres, stdin)) {
		fputs(argv[0], stderr);
		fputs(": bad resolution string on standard input\n", stderr);
		return(1);
	}
	if (needres && (ourres.xr <= 0) | (ourres.yr <= 0)) {
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
		if (rofile[i].hdrflags & DOHEADER) {
			if (!(inpflags & DOHEADER))
				newheader("RADIANCE", rofile[i].output);
			printargs(argc, argv, rofile[i].output);
			fprintf(rofile[i].output, "NCOMP=%d\n", rofile[i].output==stdout ?
						nstdoutcomp : rofile[i].ncomp);
			if (rofile[i].format != NULL) {
				extern const char  BIGEND[];
				if (rofile[i].bytsiz > 1) {
					fputs(BIGEND, rofile[i].output);
					fputs(nativebigendian() ^ swapped ?
						"1\n" : "0\n", rofile[i].output);
				}
				fputformat(rofile[i].format, rofile[i].output);
			}
			fputc('\n', rofile[i].output);
		}
		if (rofile[i].hdrflags & DORESOLU)
			fputsresolu(&ourres, rofile[i].output);
	}
	do {					/* main loop */
		for (i = 0; i < nfiles; i++) {
			if (rofile[i].bytsiz > 0) {		/* binary output */
				if (getbinary(buf, rofile[i].bytsiz, rofile[i].ncomp,
							stdin) != rofile[i].ncomp)
					break;
				if (rofile[i].output != NULL &&
					    putbinary(buf, rofile[i].bytsiz, rofile[i].ncomp,
							rofile[i].output) != rofile[i].ncomp)
					break;
			} else if (rofile[i].ncomp > 1) {	/* N-field output */
				int	n = rofile[i].ncomp;
				while (n--) {
					if (!scanOK(rofile[i].termc))
						break;
					if (rofile[i].output != NULL &&
						    fputs(buf, rofile[i].output) == EOF)
						break;
				}
				if (n >= 0)		/* fell short? */
					break;
				if ((rofile[i].output != NULL) &  /* add EOL if none */
						(rofile[i].termc != '\n'))
					fputc('\n', rofile[i].output);
			} else {			/* 1-field output */
				if (!scanOK(rofile[i].termc))
					break;
				if (rofile[i].output != NULL) {
					if (fputs(buf, rofile[i].output) == EOF)
						break;
					if (rofile[i].termc != '\n')	/* add EOL? */
						fputc('\n', rofile[i].output);
				}
			}
							/* skip input EOL? */
			if (!bininp && rofile[nfiles-1].termc != '\n') {
				int	c = getchar();
				if ((c != '\n') & (c != EOF))
					ungetc(c, stdin);
			}
		}
		if (i < nfiles)
			break;
	} while (--outcnt);
							/* check ending */
	if (fflush(NULL) == EOF) {
		fputs(argv[0], stderr);
		fputs(": write error on one or more outputs\n", stderr);
		return(1);
	}
	if (outcnt > 0) {
		fputs(argv[0], stderr);
		fputs(": warning: premature EOD\n", stderr);
	}
	return(0);
}
