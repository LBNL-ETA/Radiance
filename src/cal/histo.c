#ifndef lint
static const char	RCSid[] = "$Id: histo.c,v 1.2 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 * Compute a histogram from input data
 *
 *	9/5/96	Greg Ward
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAXCOL		64		/* maximum number of input columns */
#define MAXDIV		1024

#define isint(x)	(floor((x)+1e-6) != floor((x)-1e-6))

char	*progname;

int	cumulative = 0;

long	histo[MAXCOL][MAXDIV];
double	minv, maxv;
int	ndiv;

int	ncols;


static void
readinp(void)			/* gather statistics on input */
{
	char	buf[16*MAXCOL];
	double	d;
	register int	c;
	register char	*cp;

	while ((cp = fgets(buf, sizeof(buf), stdin)) != NULL) {
		for (c = 0; c < MAXCOL; c++) {
			while (isspace(*cp))
				cp++;
			if (!*cp)
				break;
			d = atof(cp);
			while (*cp && !isspace(*cp))
				cp++;
			if (d >= minv && d < maxv)
				histo[c][(int)(ndiv*(d-minv)/(maxv-minv))]++;
		}
		if (c > ncols)
			ncols = c;
	}
}


static void
printcumul(void)			/* print cumulative histogram results */
{
	long		ctot[MAXCOL];
	register int	i, c;

	for (c = ncols; c--; )
		ctot[c] = 0L;

	for (i = 0; i < ndiv; i++) {
		printf("%g", minv + (maxv-minv)*(i+1)/ndiv);
		for (c = 0; c < ncols; c++) {
			ctot[c] += histo[c][i];
			printf("\t%ld", ctot[c]);
		}
		putchar('\n');
	}
}


static void
printhisto(void)			/* print histogram results */
{
	register int	i, c;

	for (i = 0; i < ndiv; i++) {
		printf("%g", minv + (maxv-minv)*(i+.5)/ndiv);
		for (c = 0; c < ncols; c++)
			printf("\t%ld", histo[c][i]);
		putchar('\n');
	}
}


int
main(
int	argc,
char	*argv[]
)
{
	progname = argv[0];
	if (argc > 1 && !strcmp(argv[1], "-c")) {
		cumulative++;
		argc--; argv++;
	}
	if (argc < 3)
		goto userr;
	minv = atof(argv[1]);
	maxv = atof(argv[2]);
	if (argc == 4)
		ndiv = atoi(argv[3]);
	else {
		if (argc > 4 || !isint(minv) || !isint(maxv))
			goto userr;
		maxv += 0.5;
		minv -= 0.5;
		ndiv = maxv - minv + 0.5;
	}
	if (minv >= maxv | ndiv <= 0)
		goto userr;
	if (ndiv > MAXDIV) {
		fprintf(stderr, "%s: maximum number of divisions: %d\n",
				progname, MAXDIV);
		goto userr;
	}
	readinp();
	if (cumulative)
		printcumul();
	else
		printhisto();
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-c] min max n\n", progname);
	fprintf(stderr, "   Or: %s [-c] imin imax\n", progname);
	exit(1);
}


