#ifndef lint
static const char	RCSid[] = "$Id: histo.c,v 1.4 2003/12/14 16:33:37 greg Exp $";
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
int	percentile = 0;

long	outrange[MAXCOL][2];
long	histo[MAXCOL][MAXDIV];
long	ctotal[MAXCOL];
double	minv, maxv;
int	ndiv;

int	ncols;


static void
readinp(void)			/* gather statistics on input */
{
	char	buf[16*MAXCOL];
	double	d;
	int	i;
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
			if (d <= minv)
				outrange[c][0]++;
			else if (d >= maxv)
				outrange[c][1]++;
			else
				histo[c][(int)(ndiv*(d-minv)/(maxv-minv))]++;
		}
		if (c > ncols)
			ncols = c;
	}
	for (c = 0; c < ncols; c++) {
		ctotal[c] += outrange[c][0] + outrange[c][1];
		for (i = 0; i < ndiv; i++)
			ctotal[c] += histo[c][i];
	}
}


static void
printcumul(			/* print cumulative histogram results */
int	pctl
)
{
	long		csum[MAXCOL];
	register int	i, c;

	for (c = ncols; c--; )
		csum[c] = outrange[c][0];

	for (i = 0; i <= ndiv; i++) {
		printf("%g", minv + (maxv-minv)*i/ndiv);
		for (c = 0; c < ncols; c++) {
			if (pctl)
				printf("\t%f", 100.*csum[c]/ctotal[c]);
			else
				printf("\t%ld", csum[c]);
			if (i < ndiv)
				csum[c] += histo[c][i];
		}
		putchar('\n');
	}
}


static void
printhisto(			/* print histogram results */
int	pctl
)
{
	register int	i, c;

	for (i = 0; i < ndiv; i++) {
		printf("%g", minv + (maxv-minv)*(i+.5)/ndiv);
		for (c = 0; c < ncols; c++)
			if (pctl)
				printf("\t%f", 100.*histo[c][i]/ctotal[c]);
			else
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
	while (argc > 1 && argv[1][0] == '-') {
		if (argv[1][1] == 'c')
			cumulative++;
		else if (argv[1][1] == 'p')
			percentile++;
		else
			break;
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
	if ((minv >= maxv) | (ndiv <= 0))
		goto userr;
	if (ndiv > MAXDIV) {
		fprintf(stderr, "%s: maximum number of divisions: %d\n",
				progname, MAXDIV);
		goto userr;
	}
	readinp();
	if (cumulative)
		printcumul(percentile);
	else
		printhisto(percentile);
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-c][-p] min max n\n", progname);
	fprintf(stderr, "   Or: %s [-c][-p] imin imax\n", progname);
	exit(1);
}


