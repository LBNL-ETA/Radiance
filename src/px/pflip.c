#ifndef lint
static const char	RCSid[] = "$Id: pflip.c,v 2.6 2003/06/05 19:29:34 schorsch Exp $";
#endif
/*
 * flip picture file horizontally and/or vertically
 */

#include <stdio.h>
#include  <time.h>

#include "platform.h"
#include "color.h"
#include "resolu.h"

int	order;				/* input orientation */
int	xres, yres;			/* resolution (scanlen, nscans) */

long	*scanpos;			/* scanline positions */

int	fhoriz=0, fvert=0;		/* flip flags */

int	correctorder = 0;		/* correcting orientation? */

FILE	*fin;				/* input file */

char	*progname;


int
neworder()			/* figure out new order from old */
{
	register int  no;

	if (correctorder)
		return(order);		/* just leave it */
	if ((no = order) & YMAJOR) {
		if (fhoriz) no ^= XDECR;
		if (fvert) no ^= YDECR;
	} else {
		if (fhoriz) no ^= YDECR;
		if (fvert) no ^= XDECR;
	}
	return(no);
}


main(argc, argv)
int	argc;
char	*argv[];
{
	static char	picfmt[LPICFMT+1] = PICFMT;
	int	i, rval;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdout);
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (!strcmp(argv[i], "-h"))
			fhoriz++;
		else if (!strcmp(argv[i], "-v"))
			fvert++;
		else if (!strcmp(argv[i], "-c"))
			correctorder++;
		else
			break;
	if (i >= argc || argv[i][0] == '-') {
		fprintf(stderr, "Usage: %s [-h][-v][-c] infile [outfile]\n",
				progname);
		exit(1);
	}
	if ((fin = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", argv[1]);
		exit(1);
	}
	if (i < argc-1 && freopen(argv[i+1], "w", stdout) == NULL) {
		fprintf(stderr, "%s: cannot open\n", argv[i+1]);
		exit(1);
	}
					/* transfer header */
	if ((rval = checkheader(fin, picfmt, stdout)) < 0) {
		fprintf(stderr, "%s: input not a Radiance picture\n",
				progname);
		exit(1);
	}
	if (rval)
		fputformat(picfmt, stdout);
					/* add new header info. */
	printargs(i, argv, stdout);
	putchar('\n');
					/* get picture size */
	if ((order = fgetresolu(&xres, &yres, fin)) < 0) {
		fprintf(stderr, "%s: bad picture size\n", progname);
		exit(1);
	}
					/* write new picture size */
	fputresolu(neworder(), xres, yres, stdout);
					/* goto end if vertical flip */
	if (fvert)
		scanfile();
	flip();				/* flip the image */
	exit(0);
}


memerr()
{
	fprintf(stderr, "%s: out of memory\n", progname);
	exit(1);
}


scanfile()				/* scan to the end of file */
{
	extern long	ftell();
	COLR	*scanin;
	int	y;

	if ((scanpos = (long *)malloc(yres*sizeof(long))) == NULL)
		memerr();
	if ((scanin = (COLR *)malloc(xres*sizeof(COLR))) == NULL)
		memerr();
	for (y = yres-1; y > 0; y--) {
		scanpos[y] = ftell(fin);
		if (freadcolrs(scanin, xres, fin) < 0) {
			fprintf(stderr, "%s: read error\n", progname);
			exit(1);
		}
	}
	scanpos[0] = ftell(fin);
	free((void *)scanin);
}


flip()					/* flip the picture */
{
	COLR	*scanin, *scanout;
	int	y;
	register int	x;

	if ((scanin = (COLR *)malloc(xres*sizeof(COLR))) == NULL)
		memerr();
	if (fhoriz) {
		if ((scanout = (COLR *)malloc(xres*sizeof(COLR))) == NULL)
			memerr();
	} else
		scanout = scanin;
	for (y = yres-1; y >= 0; y--) {
		if (fvert && fseek(fin, scanpos[yres-1-y], 0) == EOF) {
			fprintf(stderr, "%s: seek error\n", progname);
			exit(1);
		}
		if (freadcolrs(scanin, xres, fin) < 0) {
			fprintf(stderr, "%s: read error\n", progname);
			exit(1);
		}
		if (fhoriz)
			for (x = 0; x < xres; x++)
				copycolr(scanout[x], scanin[xres-1-x]);
		if (fwritecolrs(scanout, xres, stdout) < 0) {
			fprintf(stderr, "%s: write error\n", progname);
			exit(1);
		}
	}
	free((void *)scanin);
	if (fhoriz)
		free((void *)scanout);
}
