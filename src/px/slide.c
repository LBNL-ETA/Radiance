#ifndef lint
static const char	RCSid[] = "$Id: slide.c,v 1.3 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  slide.c - driver for dicomed film recorder using GRAFPAC under VMS.
 *
 *     2/19/86
 *
 *	Version for color slide ouput.
 *	Link with grafpac driver "DS":
 *		link slide.obj,color.obj,sys_grafpac:ds/lib
 */

#include  <stdio.h>

#include  "color.h"

#include  "random.h"


#define  XSIZE		4096		/* device width */
#define  YSIZE		4096		/* device length */

#define  NCOLS		4018		/* maximum # columns for output */
#define  NROWS		2712		/* maximum # rows for output */

#define  XCENTER	2048		/* center for x */
#define  YCENTER	2048		/* center for y */

int  nrows;				/* number of rows for output */
int  ncols;				/* number of columns for output */

int  xoffset;				/* offset for x */
int  yoffset;				/* offset for y */

#define  RSIZ		17		/* size of response map */

float  rmap[2*RSIZ] = {				/* response map */
	.0,		.0,
	.0004,		.059,
	.0008,		.122,
	.0018,		.184,
	.0033,		.247,
	.0061,		.310,
	.0119,		.373,
	.0202,		.435,
	.0373,		.498,
	.065,		.561,
	.113,		.624,
	.199,		.686,
	.300,		.749,
	.432,		.812,
	.596,		.875,
	.767,		.937,
	1.000, 		1.000,
};


main(argc, argv)
int  argc;
char  **argv;
{
	FILE  *fin;
	COLOR  scanline[NCOLS];
	int  mult;
	char  sbuf[128];
	int  i, j;

/*
	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			default:
				fprintf(stderr, "unknown option: %s\n", argv[i]);
				exit(1);
				break;
			}
		else
			break;
			
*/
	argc = 2;
	argv[1] = "pixfile";
	i = 1;

	if (argc - i < 1)
		fin = 0;
	else if ((fin = fopen(argv[i], "r")) == NULL) {
			fprintf(stderr, "%s: can't open file \"%s\"\n",
						argv[0], argv[i]);
			exit(1);
		}
	
				/* discard header */
	while (fgets(sbuf, sizeof(sbuf), fin) != NULL && sbuf[0] != '\n')
		;
				/* get picture size */
	if (fgets(sbuf, sizeof(sbuf), fin) == NULL ||
			sscanf(sbuf, "-Y %d +X %d\n", &nrows, &ncols) != 2) {
		fprintf(stderr, "%s: bad picture size\n", argv[0]);
		exit(1);
	}

	if (ncols > NCOLS || nrows > NROWS) {
		fprintf(stderr, "%s: resolution mismatch\n", argv[0]);
		exit(1);
	}

	mult = NCOLS/ncols;
	if (NROWS/nrows < mult)
		mult = NROWS/nrows;

	xoffset = XCENTER/mult - ncols/2;
	yoffset = YCENTER/mult - nrows/2;

	tvinit();
	i = XSIZE/mult; j = YSIZE/mult; tvvrs(&i, &j);
	i = 12; tvset("ifile", &i);
	tvset("nrco", "rgb");

	for (i = nrows-1; i >= 0; i--) {
		if (freadscan(scanline, ncols, fin) < 0) {
			fprintf(stderr, "%s: read error in row %d\n",
					argv[0], i);
			exit(1);
		}
		scanout(scanline, i);
	}

	tvsend();
	tvend();

	return(0);
}


scanout(scan, y)			/* output scan line */
register COLOR  *scan;
int  y;
{
	static int  zero;
	static float  rgb[3][NCOLS];
	int  i;
	register int  j;
	
	zero = xoffset;
	y += yoffset;

	for (j = 0; j < ncols; j++) {
		mapcolor(scan[j], scan[j]);
		for (i = 0; i < 3; i++)
			rgb[i][j] = scan[j][i] + (.5-frandom())/255.0;
	}
	tvrstc(&zero, &y, rgb[0], rgb[1], rgb[2], &ncols);
}


mapcolor(c1, c2)		/* map c1 into c2 */
register COLOR  c1;
COLOR  c2;
{
	static float  *mp[3] = {rmap, rmap, rmap};
	double  x;
	register int  i;

	for (i = 0; i < 3; i++) {
		if (colval(c1,i) >= 1.0) {
			colval(c2,i) = 1.0;
			continue;
		}
		while (colval(c1,i) >= mp[i][2])
			mp[i] += 2;
		while (colval(c1,i) < mp[i][0])
			mp[i] -= 2;
		x = (colval(c1,i) - mp[i][0]) / (mp[i][2] - mp[i][0]);
		colval(c2,i) = (1.0-x)*mp[i][1] + x*mp[i][3];
	}
}
