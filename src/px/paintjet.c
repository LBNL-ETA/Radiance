#ifndef lint
static const char	RCSid[] = "$Id: paintjet.c,v 2.6 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  paintjet.c - program to dump pixel file to HP PaintJet color printer.
 *
 *     10/6/89
 */

#include  <stdio.h>
#include  <time.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

#define	 NCOLS		1440		/* 8" at 180 dpi */

static int printp(char  *fname);
static void plotscan(COLR  scan[], int  len, int  y);
static int colbit(COLR  col, int  x, int  a);


int
main(
	int  argc,
	char  *argv[]
)
{
	int  i, status = 0;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	if (argc < 2)
		status = printp(NULL) == -1;
	else
		for (i = 1; i < argc; i++)
			status += printp(argv[i]) == -1;
	exit(status);
}


static int
printp(				/* print a picture */
	char  *fname
)
{
	FILE  *input;
	int  xres, yres;
	COLR  scanline[NCOLS];
	int  i;

	if (fname == NULL) {
		input = stdin;
		fname = "<stdin>";
	} else if ((input = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fname);
		return(-1);
	}
				/* discard header */
	if (checkheader(input, COLRFMT, NULL) < 0) {
		fprintf(stderr, "%s: not a Radiance picture\n", fname);
		return(-1);
	}
				/* get picture dimensions */
	if (fgetresolu(&xres, &yres, input) < 0) {
		fprintf(stderr, "%s: bad picture size\n", fname);
		return(-1);
	}
	if (xres > NCOLS) {
		fprintf(stderr, "%s: resolution mismatch\n", fname);
		return(-1);
	}
	/*
	 * Set 180 dpi, 3 planes, left margin and width.
	 */
	printf("\033*t180R\033*r3U\033&a%dH\033*r%dS\033*r1A",
			((NCOLS-xres)>>4)<<5, xres);
	
	for (i = yres-1; i >= 0; i--) {
		if (freadcolrs(scanline, xres, input) < 0)
			return(-1);
		normcolrs(scanline, xres, 0);
		plotscan(scanline, xres, i);
	}

	printf("\033*r0B\f");
	
	fclose(input);

	return(0);
}


static void
plotscan(			/* plot a scanline */
	COLR  scan[],
	int  len,
	int  y
)
{
	int  c;
	register int  x, b;

	for (c = 0; c < 3; c++) {
		printf("\033*b%d%c", (len+7)>>3, c<2 ? 'V' : 'W');
		b = 0;
		for (x = 0; x < len; x++) {
			b = b<<1 | colbit(scan[x], x, c);
			if ((x&7) == 7) {
				putc(b, stdout);
				b = 0;
			}
		}
		if (x & 7)
			putc(b, stdout);
	}
}


static int
colbit(		/* determine bit value for primary at x */
	COLR  col,
	register int  x,
	register int  a
)
{
	static int  cerr[NCOLS][3];
	static int  err[3];
	int  b, errp;
	register int  ison;

	b = col[a];
	errp = err[a];
	err[a] += b + cerr[x][a];
	ison = err[a] > 128;
	if (ison) err[a] -= 256;
	err[a] /= 3;
	cerr[x][a] = err[a] + errp;
	return(ison);
}
