#ifndef lint
static const char	RCSid[] = "$Id: oki20.c,v 2.11 2003/10/27 10:24:51 schorsch Exp $";
#endif
/*
 *  oki20.c - program to dump pixel file to OkiMate 20 printer.
 */

#include  <stdio.h>
#include  <time.h>

#include  "platform.h"
#include  "rtprocess.h"
#include  "color.h"
#include  "resolu.h"

#define	 NROWS		1440		/* 10" at 144 dpi */
#define	 NCOLS		960		/* 8" at 120 dpi */

#define	 ASPECT		(120./144.)	/* pixel aspect ratio */

#define	 FILTER		"pfilt -1 -x %d -y %d -p %f %s",NCOLS,NROWS,ASPECT

long  lpat[NCOLS];

int  dofilter = 0;		/* filter through pfilt first? */


main(argc, argv)
int  argc;
char  *argv[];
{
	int  i, status = 0;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	if (argc > 1 && !strcmp(argv[1], "-p")) {
		dofilter++;
		argv++; argc--;
	}
	if (argc < 2)
		status = printp(NULL) == -1;
	else
		for (i = 1; i < argc; i++)
			status += printp(argv[i]) == -1;
	exit(status);
}


printp(fname)				/* print a picture */
char  *fname;
{
	char  buf[64];
	FILE  *input;
	int  xres, yres;
	COLR  scanline[NCOLS];
	int  i;

	if (dofilter) {
		if (fname == NULL) {
			sprintf(buf, FILTER, "");
			fname = "<stdin>";
		} else
			sprintf(buf, FILTER, fname);
		if ((input = popen(buf, "r")) == NULL) {
			fprintf(stderr, "Cannot execute: %s\n", buf);
			return(-1);
		}
	} else if (fname == NULL) {
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
				/* set line spacing (overlap for knitting) */
	fputs("\0333\042\022", stdout);
				/* put out scanlines */
	for (i = yres-1; i >= 0; i--) {
		if (freadcolrs(scanline, xres, input) < 0) {
			fprintf(stderr, "%s: read error (y=%d)\n", fname, i);
			return(-1);
		}
		normcolrs(scanline, xres, 0);
		plotscan(scanline, xres, i);
	}
				/* advance page */
	putchar('\f');
	
	if (dofilter)
		pclose(input);
	else
		fclose(input);

	return(0);
}


plotscan(scan, len, y)			/* plot a scanline */
COLR  scan[];
int  len;
int  y;
{
	int  bpos, start, end;
	register long  c;
	register int  i;

	bpos = y % 23;
	for (i = 0; i < len; i++)
		lpat[i] |= (long)bit(scan[i],i) << bpos;

	if (bpos)
		return;
				/* find limits of non-zero print buffer */
	for (i = 0; lpat[i] == 0; i++)
		if (i == len-1) {
			putchar('\n');
			return;
		}
	start = i - i%12;
	i = len;
	while (lpat[--i] == 0)
		;
	end = i;
				/* skip to start position */
	for (i = start/12; i-- > 0; )
		putchar(' ');
				/* print non-zero portion of buffer */
	fputs("\033%O", stdout);
	i = end+1-start;
	putchar(i & 255);
	putchar(i >> 8);
	for (i = start; i <= end; i++) {
		c = lpat[i];
		putchar((int)(c>>16));
		putchar((int)(c>>8 & 255));
		putchar((int)(c & 255));
		if (y)			/* repeat this row next time */
			lpat[i] = (c & 1) << 23;
		else			/* or clear for next image */
			lpat[i] = 0L;
	}
	putchar('\r');
	putchar('\n');
	fflush(stdout);
}


bit(col, x)		/* determine bit value for pixel at x */
COLR  col;
register int  x;
{
	static int  cerr[NCOLS];
	static int  err;
	int  b, errp;
	register int  ison;

	b = normbright(col);
	errp = err;
	err += b + cerr[x];
	ison = err < 128;
	if (!ison) err -= 256;
	err /= 3;
	cerr[x] = err + errp;
	return(ison);
}
