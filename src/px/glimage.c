#ifndef lint
static const char	RCSid[] = "$Id: glimage.c,v 2.2 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  glimage.c - program to dump pixel file NeWS window
 *
 *     1/24/89
 */

#include  <stdio.h>

#include  <gl.h>

#include  <time.h>

#include  "color.h"

#include  "resolu.h"


main(argc, argv)				/* print a picture */
int  argc;
char  *argv[];
{
	char  *fname;
	FILE  *input;
	int  xres, yres;
	int  winorgx, winorgy;
	COLR  *scanline;
	long  *winout;
	int  i;
	register int  j;

	fname = argv[1];
	if (fname == NULL) {
		input = stdin;
		fname = "<stdin>";
	} else if ((input = fopen(fname, "r")) == NULL) {
		perror(fname);
		exit(1);
	}
				/* discard header */
	if (checkheader(input, COLRFMT, NULL) < 0) {
		fprintf(stderr, "%s: not a Radiance picture\n", fname);
		exit(1);
	}
				/* get picture dimensions */
	if (fgetresolu(&xres, &yres, input) < 0) {
		fprintf(stderr, "%s: bad picture size\n", fname);
		exit(1);
	}
	prefsize(xres, yres);
	winopen("newsimage");
	RGBmode();
	gconfig();
	RGBcolor(0,0,0);
	clear();
	getorigin(&winorgx,&winorgy);
	winout = (long *)malloc(xres*sizeof(long));
	scanline = (COLR *)malloc(xres*sizeof(COLR));
	if (winout == NULL || scanline == NULL) {
		perror(argv[0]);
		exit(1);
	}
	for (i = yres-1; i >= 0; i--) {
		if (freadcolrs(scanline, xres, input) < 0) {
			fprintf(stderr, "%s: read error (y=%d)\n", fname, i);
			exit(1);
		}
		normcolrs(scanline, xres, 0);
		for (j = 0; j < xres; j++) {
			winout[j] = 255L<<24;
			winout[j] |= (long)scanline[j][BLU]<<16;
			winout[j] |= (long)scanline[j][GRN]<<8;
			winout[j] |= (long)scanline[j][RED];
		}
		lrectwrite(0,i,xres-1,i,winout);
	}
	pause();
	exit(0);
}
