#ifndef lint
static const char	RCSid[] = "$Id: oki20c.c,v 2.11 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  oki20c.c - program to dump pixel file to OkiMate 20 color printer.
 */

#include  <stdio.h>
#ifdef MSDOS
#include  <fcntl.h>
#endif
#include  <time.h>

#include  "color.h"
#include  "resolu.h"

#define	 NROWS		1440		/* 10" at 144 dpi */
#define	 NCOLS		960		/* 8" at 120 dpi */

#define	 ASPECT		(120./144.)	/* pixel aspect ratio */

#define	 FILTER		"pfilt -1 -x %d -y %d -p %f %s",NCOLS,NROWS,ASPECT

/*
 *  Subtractive primaries are ordered:	Yellow, Magenta, Cyan.
 */

#define	 sub_add(sub)	(2-(sub))	/* map subtractive to additive pri. */

long  lpat[NCOLS][3];

int  dofilter = 0;		/* filter through pfilt first? */

extern FILE  *popen();


main(argc, argv)
int  argc;
char  *argv[];
{
	int  i, status = 0;
#ifdef MSDOS
	extern int  _fmode;
	_fmode = O_BINARY;
	setmode(fileno(stdin), O_BINARY);
	setmode(fileno(stdout), O_BINARY);
#endif
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
	fputs("\0333\042", stdout);
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
	int  bpos;
	register long  c;
	register int  i, j;

	if (bpos = y % 23) {

		for (j = 0; j < 3; j++)
			for (i = 0; i < len; i++)
				lpat[i][j] |= (long)colbit(scan[i],i,j) << bpos;
		return;
	}
	fputs("\033\031", stdout);

	for (j = 0; j < 3; j++) {
		i = (NCOLS + len)/2;		/* center image */
		fputs("\033%O", stdout);
		putchar(i & 255);
		putchar(i >> 8);
		while (i-- > len) {
			putchar(0);
			putchar(0);
			putchar(0);
		}
		for (i = 0; i < len; i++) {
			c = lpat[i][j] | colbit(scan[i],i,j);
			putchar((int)(c>>16));
			putchar((int)(c>>8 & 255));
			putchar((int)(c & 255));
			if (y)			/* repeat this row */
				lpat[i][j] = (c & 1) << 23;
			else			/* or clear for next image */
				lpat[i][j] = 0L;
		}
		putchar('\r');
	}
	putchar('\n');
	fflush(stdout);
}


colbit(col, x, s)		/* determine bit value for primary at x */
COLR  col;
register int  x;
int  s;
{
	static int  cerr[NCOLS][3];
	static int  err[3];
	int  b, errp;
	register int  a, ison;

	a = sub_add(s);			/* use additive primary */
	b = col[a];
	errp = err[a];
	err[a] += b + cerr[x][a];
	ison = err[a] < 128;
	if (!ison) err[a] -= 256;
	err[a] /= 3;
	cerr[x][a] = err[a] + errp;
	return(ison);
}
