#ifndef lint
static const char	RCSid[] = "$Id: ra_bn.c,v 2.9 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  ra_bn.c - program to convert between RADIANCE and barneyscan picture format.
 *
 *	10/12/88
 */

#include  <stdio.h>
#include  <string.h>
#include  <time.h>
#include  <math.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

double	gamcor = 2.0;			/* gamma correction */
int  bradj = 0;				/* brightness adjustment */
char  *progname;
char  errmsg[128];
FILE	*rafp, *bnfp[3];
int  xmax, ymax;

static void quiterr(char *err);
static int openbarney(char *fname, char *mode);
static int getint(FILE *fp);
static void putint(int val, FILE *fp);
static void ra2bn(void);
static void bn2ra(void);


int
main(
	int  argc,
	char  *argv[]
)
{
	int  reverse = 0;
	int  i;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			case 'r':
				reverse = !reverse;
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			default:
				goto userr;
			}
		else
			break;
						/* set gamma correction */
	setcolrgam(gamcor);

	if (reverse) {
		if (i > argc-1 || i < argc-2)
			goto userr;
		if (openbarney(argv[i], "r") < 0) {
			sprintf(errmsg, "cannot open input \"%s\"", argv[i]);
			quiterr(errmsg);
		}
		if (i == argc-1 || !strcmp(argv[i+1], "-"))
			rafp = stdout;
		else if ((rafp = fopen(argv[i+1], "w")) == NULL) {
			sprintf(errmsg, "cannot open output \"%s\"",
					argv[i+1]);
			quiterr(errmsg);
		}
					/* put header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, rafp);
		fputformat(COLRFMT, rafp);
		putc('\n', rafp);
		fprtresolu(xmax, ymax, rafp);
					/* convert file */
		bn2ra();
	} else {
		if (i != argc-2)
			goto userr;
		if (!strcmp(argv[i], "-"))
			rafp = stdin;
		else if ((rafp = fopen(argv[i], "r")) == NULL) {
			sprintf(errmsg, "cannot open input \"%s\"",
					argv[i]);
			quiterr(errmsg);
		}
					/* get header */
		if (checkheader(rafp, COLRFMT, NULL) < 0 ||
				fgetresolu(&xmax, &ymax, rafp) < 0)
			quiterr("bad RADIANCE format");
		if (openbarney(argv[i+1], "w") < 0) {
			sprintf(errmsg, "cannot open output \"%s\"", argv[i+1]);
			quiterr(errmsg);
		}
					/* convert file */
		ra2bn();
	}
	quiterr(NULL);
userr:
	fprintf(stderr, "Usage: %s [-g gamma][-e +/-stops] {input|-} output\n",
			progname);
	fprintf(stderr, "   or: %s -r [-g gamma][-e +/-stops] input [output|-]\n",
			progname);
	exit(1);
}


static void
quiterr(		/* print message and exit */
	char  *err
)
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


static int
openbarney(			/* open barneyscan files */
	char	*fname,
	char	*mode
)
{
	static char	suffix[3][4] = {"red", "grn", "blu"};
	int	i;
	char	file[128];

	for (i = 0; i < 3; i++) {
		sprintf(file, "%s.%s", fname, suffix[i]);
		if ((bnfp[i] = fopen(file, mode)) == NULL)
			return(-1);
	}
	if (mode[0] == 'r') {		/* get header */
		xmax = getint(bnfp[0]);
		ymax = getint(bnfp[0]);
		for (i = 1; i < 3; i++)
			if (getint(bnfp[i]) != xmax || getint(bnfp[i]) != ymax)
				return(-1);
	} else {			/* put header */
		for (i = 0; i < 3; i++) {
			putint(xmax, bnfp[i]);
			putint(ymax, bnfp[i]);
		}
	}
	return(0);
}


static int
getint(				/* get short int from barneyscan file */
	register FILE	*fp
)
{
	register short	val;

	val = getc(fp);
	val |= getc(fp) << 8;

	return(val);
}


static void
putint(				/* put short int to barneyscan file */
	register int	val,
	register FILE	*fp
)
{
	putc(val&0xff, fp);
	putc((val >> 8)&0xff, fp);
}


static void
ra2bn(void)					/* convert radiance to barneyscan */
{
	register int	i;
	register COLR	*inl;
	int	j;

	if ((inl = (COLR *)malloc(xmax*sizeof(COLR))) == NULL)
		quiterr("out of memory");
	for (j = 0; j < ymax; j++) {
		if (freadcolrs(inl, xmax, rafp) < 0)
			quiterr("error reading RADIANCE file");
		if (bradj)
			shiftcolrs(inl, xmax, bradj);
		colrs_gambs(inl, xmax);
		for (i = 0; i < xmax; i++) {
			putc(inl[i][RED], bnfp[0]);
			putc(inl[i][GRN], bnfp[1]);
			putc(inl[i][BLU], bnfp[2]);
		}
		if (ferror(bnfp[0]) || ferror(bnfp[1]) || ferror(bnfp[2]))
			quiterr("error writing Barney files");
	}
	free((void *)inl);
}


static void
bn2ra(void)					/* convert barneyscan to radiance */
{
	register int	i;
	register COLR	*outline;
	int	j;

	if ((outline = (COLR *)malloc(xmax*sizeof(COLR))) == NULL)
		quiterr("out of memory");
	for (j = 0; j < ymax; j++) {
		for (i = 0; i < xmax; i++) {
			outline[i][RED] = getc(bnfp[0]);
			outline[i][GRN] = getc(bnfp[1]);
			outline[i][BLU] = getc(bnfp[2]);
		}
		if (feof(bnfp[0]) || feof(bnfp[1]) || feof(bnfp[2]))
			quiterr("error reading barney file");
		gambs_colrs(outline, xmax);
		if (bradj)
			shiftcolrs(outline, xmax, bradj);
		if (fwritecolrs(outline, xmax, rafp) < 0)
			quiterr("error writing RADIANCE file");
	}
	free((void *)outline);
}
