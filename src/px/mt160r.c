#ifndef lint
static const char	RCSid[] = "$Id: mt160r.c,v 2.8 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  mt160r.c - program to dump pixel file to Mannesman-Tally 160.
 *
 *     8/16/85
 */

#include  <stdio.h>
#include  <time.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

#define	 NCOLS		880		/* for wide carriage */

static int printp(char  *fname);
static void plotscan(COLR  scan[], int  len, int  y);
static int bit(COLR  clr, int  x);


int
main(
	int  argc,
	char  *argv[]
)
{
	int  i;
	int  status = 0;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	if (argc < 2)
		status += printp(NULL) == -1;
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

	fputs("\033[6~\033[7z", stdout);
	
	for (i = yres-1; i >= 0; i--) {
		if (freadcolrs(scanline, xres, input) < 0) {
			fprintf(stderr, "%s: read error (y=%d)\n", fname, i);
			return(-1);
		}
		normcolrs(scanline, xres, 0);
		plotscan(scanline, xres, i);
	}

	fputs("\f\033[6~", stdout);

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
	static BYTE  pat[NCOLS];
	int  bpos;
	register int  i;

	if ((bpos = y & 7)) {

		for (i = 0; i < len; i++)
			pat[i] |= bit(scan[i], i) << bpos;

	} else {

		fputs("\033%5", stdout);
		putchar(len & 255);
		putchar(len >> 8);

		for (i = 0; i < len; i++) {
			pat[i] |= bit(scan[i], i);
			putchar(pat[i]);
			pat[i] = 0;
		}
		putchar('\r');
		putchar('\n');
		fflush(stdout);
	}
}


static int
bit(				/* return bit for color at x */
	COLR  clr,
	register int  x
)
{
	static int  cerr[NCOLS];
	static int  err;
	int  b, errp;
	register int  isblack;

	b = normbright(clr);
	errp = err;
	err += b + cerr[x];
	isblack = err < 128;
	if (!isblack) err -= 256;
	err /= 3;
	cerr[x] = err + errp;
	return(isblack);
}
