#ifndef lint
static const char	RCSid[] = "$Id: ra_avs.c,v 2.9 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  Convert Radiance file to/from AVS file.
 */

#include  <stdio.h>
#include  <math.h>
#ifdef MSDOS
#include  <fcntl.h>
#endif
#include  <time.h>
#include  "color.h"
#include  "resolu.h"

double	gamcor = 2.2;			/* gamma correction */

int  bradj = 0;				/* brightness adjustment */

char  *progname;

int  xmax, ymax;


main(argc, argv)
int  argc;
char  *argv[];
{
	extern long  getint();
	int  reverse = 0;
	int  i;
	
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':		/* gamma correction */
				gamcor = atof(argv[++i]);
				break;
			case 'e':		/* exposure adjustment */
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'r':		/* reverse conversion */
				reverse = 1;
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				progname, argv[i]);
		exit(1);
	}
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		fprintf(stderr, "%s: can't open output \"%s\"\n",
				progname, argv[i+1]);
		exit(1);
	}
#ifdef MSDOS
	setmode(fileno(stdin), O_BINARY);
	setmode(fileno(stdout), O_BINARY);
#endif
	setcolrgam(gamcor);		/* set up gamma correction */
	if (reverse) {
					/* get their image resolution */
		xmax = getint(4, stdin);
		ymax = getint(4, stdin);
		if (feof(stdin))
			quiterr("empty input file");
					/* put our header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fprtresolu(xmax, ymax, stdout);
					/* convert file */
		avs2ra();
	} else {
					/* get our header */
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				fgetresolu(&xmax, &ymax, stdin) < 0)
			quiterr("bad picture format");
					/* write their header */
		putint((long)xmax, 4, stdout);
		putint((long)ymax, 4, stdout);
					/* convert file */
		ra2avs();
	}
	exit(0);
userr:
	fprintf(stderr,
		"Usage: %s [-r][-g gamma][-e +/-stops] [input [output]]\n",
			progname);
	exit(1);
}


quiterr(err)		/* print message and exit */
char  *err;
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


avs2ra()		/* convert 24-bit scanlines to Radiance picture */
{
	COLR	*scanout;
	register int	x;
	int	y;
						/* allocate scanline */
	scanout = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanout == NULL)
		quiterr("out of memory in avs2ra");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		for (x = 0; x < xmax; x++) {
			(void)getc(stdin);			/* toss alpha */
			scanout[x][RED] = getc(stdin);
			scanout[x][GRN] = getc(stdin);
			scanout[x][BLU] = getc(stdin);
		}
		if (feof(stdin) | ferror(stdin))
			quiterr("error reading AVS image");
						/* undo gamma */
		gambs_colrs(scanout, xmax);
		if (bradj)			/* adjust exposure */
			shiftcolrs(scanout, xmax, bradj);
		if (fwritecolrs(scanout, xmax, stdout) < 0)
			quiterr("error writing Radiance picture");
	}
						/* free scanline */
	free((void *)scanout);
}


ra2avs()		/* convert Radiance scanlines to 24-bit */
{
	COLR	*scanin;
	register int	x;
	int	y;
						/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in ra2avs");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		if (bradj)			/* adjust exposure */
			shiftcolrs(scanin, xmax, bradj);
		colrs_gambs(scanin, xmax);	/* gamma correction */
		for (x = 0; x < xmax; x++) {
			putc(0, stdout);		/* no alpha */
			putc(scanin[x][RED], stdout);
			putc(scanin[x][GRN], stdout);
			putc(scanin[x][BLU], stdout);
		}
		if (ferror(stdout))
			quiterr("error writing AVS file");
	}
						/* free scanline */
	free((void *)scanin);
}
