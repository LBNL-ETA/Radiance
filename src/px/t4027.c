#ifndef lint
static const char	RCSid[] = "$Id: t4027.c,v 2.4 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  t4027.c - program to dump pixel file to Tektronix 4027.
 *
 *     8/15/85
 */

#include  <stdio.h>
#ifndef _WIN32
 #include  <unistd.h>
#endif

#include  "color.h"


#define  NROWS		462
#define  NCOLS		640

#define  COM		31		/* command character */

char  *progname;

static void plotscan(COLOR  scan[], int  len, int  y);
static int colormap(COLOR  col);
static void vector(int  x0, int  y0, int  x1, int  y1);
static void setdcolr(int  col);
static void printpat(int  p);
static int rorb(int  b);


int
main(
	int  argc,
	char  **argv
)
{
	int  xres, yres;
	COLOR  scanline[NCOLS];
	char  sbuf[128];
	FILE  *fp;
	register int  i;
	
	progname = argv[0];

	if (argc < 2)
		fp = stdin;
	else if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "%s: can't open file \"%s\"\n",
				progname, argv[1]);
		exit(1);
	}
				/* discard header */
	while (fgets(sbuf, sizeof(sbuf), fp) != NULL && sbuf[0] != '\n')
		;
				/* get picture dimensions */
	if (fgets(sbuf, sizeof(sbuf), fp) == NULL ||
			sscanf(sbuf, "-Y %d +X %d\n", &yres, &xres) != 2) {
		fprintf(stderr, "%s: bad picture size\n", progname);
		exit(1);
	}
	if (xres > NCOLS || yres > NROWS) {
		fprintf(stderr, "%s: resolution mismatch\n", progname);
		exit(1);
	}

	printf("%cMON1H", COM);
	printf("%cGRA1 33", COM);
	
	for (i = 0; i < yres; i++) {
		if (freadscan(scanline, xres, fp) < 0) {
			fprintf(stderr, "%s: read error\n", progname);
			exit(1);
		}
		plotscan(scanline, xres, NROWS-1-i);
	}
	fclose(fp);

#ifndef _WIN32  /* XXX extract console functions */
	if (isatty(fileno(stdout))) {
		printf("Hit return when done: ");
		fflush(stdout);
		fp = fopen("/dev/tty", "r");
		getc(fp);
		fclose(fp);
		printf("%cMON34", COM);
	}
	putchar('\r');
#endif

	exit(0);
	return 0; /* pro forma return */
}


static void
plotscan(			/* plot a scanline */
	COLOR  scan[],
	int  len,
	int  y
)
{
	int  color, lastcolor = -1;
	int  lastpos = -1;
	register int  i;

	for (i = 0; i < len; i++) {

		color = colormap(scan[i]);
		if (color != lastcolor) {
			if (lastpos >= 0) {
				setdcolr(lastcolor);
				vector(lastpos, y, i, y);
			}
			lastcolor = color;
			lastpos = i;
		}
	}
	setdcolr(lastcolor);
	vector(lastpos, y, len-1, y);
}


static int
colormap(			/* compute 4027 color value for col */
	COLOR  col
)
{
	register int  red, grn, blu;

	red = col[RED] >= 1.0 ? 3 : col[RED]*4;
	grn = col[GRN] >= 1.0 ? 3 : col[GRN]*4;
	blu = col[BLU] >= 1.0 ? 3 : col[BLU]*4;

	return(red<<4 | grn<<2 | blu);
}


static void
vector(		/* output a vector */
	int  x0,
	int  y0,
	int  x1,
	int  y1
)
{
	printf("%cVEC%d %d %d %d", COM, x0, y0, x1, y1);
}


static void
setdcolr(			/* set dithered pattern for terminal */
	int  col
)
{
	static int  cmask[3][4] = {
		{ 0, 0x88, 0xca, 0xff },
		{ 0, 0x44, 0xb2, 0xff },
		{ 0, 0x22, 0x2d, 0xff }
	};
	/* isset is a common bitmap related macro in <sys/param.h> */
	static short  cisset[64];
	int  pat;
	register int  r, g, b;
	
	if (cisset[col]) {
		printf("%cCOL P%d", COM, col);
		return;
	}

	printf("%cPAT P%d C7", COM, col);
	
	r = cmask[0][col>>4];
	g = cmask[1][(col>>2) & 03];
	b = cmask[2][col & 03];

	pat = r & ~(g|b);
	if (pat) {
		printf(" C1");		/* red */
		printpat(pat);
	}
	pat = g & ~(r|b);
	if (pat) {
		printf(" C2");		/* green */
		printpat(pat);
	}
	pat = b & ~(r|g);
	if (pat) {
		printf(" C3");		/* blue */
		printpat(pat);
	}
	pat = r & g & ~b;
	if (pat) {
		printf(" C4");		/* yellow */
		printpat(pat);
	}
	pat = r & b & ~g;
	if (pat) {
		printf(" C6");		/* magenta */
		printpat(pat);
	}
	pat = g & b & ~r;
	if (pat) {
		printf(" C5");		/* cyan */
		printpat(pat);
	}
	pat = r & g & b;
	if (pat) {
		printf(" C0");		/* white */
		printpat(pat);
	}
	cisset[col] = 1;
	printf("%cCOL P%d", COM, col);
}


static void
printpat(			/* print out a pattern */
	register int  p
)
{
	register int  i;
	
	for (i = 0; i < 14; i++) {		/* 14 rows */
	
		printf(",%d", p);
		p = rorb(rorb(rorb(p)));	/* rotate 3 bits */
	}
}


static int
rorb(				/* rotate a byte to the right */
	register int  b
)
{
	b |= (b&01) << 8;
	b >>= 1;
	return(b);
}
