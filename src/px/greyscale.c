#ifndef lint
static const char	RCSid[] = "$Id: greyscale.c,v 2.6 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  greyscale.c - program to produce grey test levels.
 *
 *     4/21/86
 */

#include  <stdio.h>

#include  <math.h>

#include  "color.h"


double  minlog = 0.0;		/* minimum for log scale (0 == linear) */


main(argc, argv)
int  argc;
char  *argv[];
{
	COLOR  col;
	double  d1,d2,d3;
	int  i;

	printargs(argc, argv, stdout);

	setcolor(col, 1.0, 1.0, 1.0);

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'c':
			d1 = atof(argv[++i]);
			d2 = atof(argv[++i]);
			d3 = atof(argv[++i]);
			setcolor(col, d1, d2, d3);
			break;
		case 'l':
			d1 = atof(argv[++i]);
			if (d1 <= 0.0)
				minlog = 0.0;
			else
				minlog = log(d1);
			break;
		default:
			fprintf(stderr, "%s: unknown option \"%s\"\n",
					argv[0], argv[i]);
			exit(1);
		}

	printf("\n");
	printf("-Y 512 +X 512\n");
	greyscale(col);
}


greyscale(col0)			/* output our grey scale */
COLOR  col0;
{
	COLOR  col1, col2, scanline[512];
	double  x;
	int  j;
	register int  i, k;

	for (j = 0; j < 512; j += 32) {
		for (k = 0; k < 512; k++)
			setcolor(scanline[k], 0.0, 0.0, 0.0);
		for (k = 0; k < 4; k++)
			if (fwritescan(scanline, 512, stdout) < 0)
				goto writerr;
		x = j/32 / 16.0;
		if (minlog != 0.0)
			x = exp((1.0-x)*minlog);
		setcolor(col1, x, x, x);

		for (i = 0; i < 512; i += 32) {
			for (k = 0; k < 4; k++)
				setcolor(scanline[i+k], 0.0, 0.0, 0.0);
			x = i/32 / 255.0;
			if (minlog != 0.0) {
				x = exp(-x*minlog);
				setcolor(col2, x, x, x);
				multcolor(col2, col1);
			} else {
				setcolor(col2, x, x, x);
				addcolor(col2, col1);
			}
			multcolor(col2, col0);
			for (k = 4; k < 32; k++)
				copycolor(scanline[i+k], col2);
		}
		for (i = 0; i < 28; i++)
			if (fwritescan(scanline, 512, stdout) < 0)
				goto writerr;
	}
	return;
writerr:
	fprintf(stderr, "write error in greyscale\n");
	exit(1);
}


printargs(ac, av, fp)		/* print arguments to a file */
int  ac;
char  **av;
FILE  *fp;
{
	while (ac-- > 0) {
		fputs(*av++, fp);
		putc(' ', fp);
	}
	putc('\n', fp);
}
