#ifndef lint
static const char	RCSid[] = "$Id: ra_hexbit.c,v 3.2 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 * Create a 4x1 hex bitmap from a Radiance picture.
 */

#include  <stdio.h>
#include  <time.h>
#include  "color.h"
#include  "resolu.h"

char  *progname;

int  xmax, ymax;

double	thresh = 0.5;		/* threshold value */

COLR	threshclr;

#define abovethresh(c)	((c)[EXP]>threshclr[EXP] || \
			((c)[EXP]==threshclr[EXP] && (c)[GRN]>threshclr[GRN]))


main(argc, argv)
int  argc;
char  *argv[];
{
	int  i;
	
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 't':		/* threshold value */
				thresh = atof(argv[++i]);
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
				/* assign threshold color */
	setcolr(threshclr, thresh, thresh, thresh);
				/* get our header */
	if (checkheader(stdin, COLRFMT, NULL) < 0 ||
			fgetresolu(&xmax, &ymax, stdin) < 0)
		quiterr("bad picture format");
				/* convert file */
	ra2hex();
	exit(0);
userr:
	fprintf(stderr,
		"Usage: %s [-t thresh] [input [output]]\n", progname);
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


ra2hex()		/* convert Radiance scanlines to 4x1 bit hex */
{
	static char	cmap[] = "0123456789ABCDEF";
	COLR	*scanin;
	register int	x, c, t;
	int	y;
						/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in ra2skel");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		c = 0;
		for (x = 0; x < xmax; x++)
			if ((t = 03 - (x&03)))
				c |= abovethresh(scanin[x]) << t;
			else {
				c |= abovethresh(scanin[x]);
				putchar(cmap[c]);
				c = 0;
			}
		if (t)
			fputc(cmap[c], stdout);
		fputc('\n', stdout);
		if (ferror(stdout))
			quiterr("error writing hex bit file");
	}
						/* free scanline */
	free((void *)scanin);
}
