#ifndef lint
static const char	RCSid[] = "$Id: ttyimage.c,v 2.4 2003/06/05 19:29:34 schorsch Exp $";
#endif
/*
 *  ttyimage.c - program to dump pixel file to dumb terminal.
 *
 *     8/15/85
 */

#include  <stdio.h>
#include  <time.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"


#define	 NCOLS		133


main(argc, argv)
int  argc;
char  **argv;
{
	FILE  *input;
	int  xres, yres;
	COLR  scanline[NCOLS];
	register int  i, j;
	
	if (argc < 2)
		input = stdin;
	else if ((input = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "%s: can't open file \"%s\"\n", argv[0], argv[1]);
		exit(1);
	}
	SET_FILE_BINARY(input);
				/* get picture dimensions */
	if (checkheader(input, COLRFMT, NULL) < 0 ||
			fgetresolu(&xres, &yres, input) < 0) {
		fprintf(stderr, "%s: bad picture format\n", argv[0]);
		exit(1);
	}
	if (xres > NCOLS) {
		fprintf(stderr, "%s: resolution mismatch\n", argv[0]);
		exit(1);
	}
	
	for (i = 0; i < yres; i++) {
		if (freadcolrs(scanline, xres, input) < 0) {
			fprintf(stderr, "%s: read error\n", argv[0]);
			exit(1);
		}
		normcolrs(scanline, xres, 0);
		for (j = 0; j < xres; j++)
			putchar(shade(scanline[j]));
		putchar('\n');
	}

	exit(0);
}


int
shade(clr)			/* return character for color */
COLR  clr;
{
#define NSHADES	 13

	static char  shadech[NSHADES+1] = " .,:;+?%&*$@#";

	return(shadech[normbright(clr)*NSHADES/256]);

#undef NSHADES
}
