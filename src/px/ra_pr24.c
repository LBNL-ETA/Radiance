#ifndef lint
static const char	RCSid[] = "$Id: ra_pr24.c,v 2.12 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  program to convert between RADIANCE and 24-bit rasterfiles.
 */

#include  <stdio.h>
#include  <time.h>
#include  <math.h>
#include  <string.h>

#include  "platform.h"
#include  "rasterfile.h"
#include  "color.h"
#include  "resolu.h"

double	gamcor = 2.2;			/* gamma correction */
int  bradj = 0;				/* brightness adjustment */
char  *progname;
int  xmax, ymax;

static void quiterr(char  *err);
static void pr2ra(int	rf, int	pad);
static void ra2pr(int	rf, int	pad);


int
main(int  argc, char  *argv[])
{
	struct rasterfile  head;
	int  reverse = 0;
	int  i;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	progname = argv[0];

	head.ras_type = RT_STANDARD;
	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'r':
				if (!strcmp(argv[i], "-rgb"))
					head.ras_type = RT_FORMAT_RGB;
				else
					reverse = 1;
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				progname, argv[i]);
		exit(1);
	}
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		fprintf(stderr, "%s: can't open output \"%s\"\n",
				progname, argv[i+1]);
		exit(1);
	}
	setcolrgam(gamcor);
	if (reverse) {
					/* get header */
		if (fread((char *)&head, sizeof(head), 1, stdin) != 1)
			quiterr("missing header");
		if (head.ras_magic != RAS_MAGIC)
			quiterr("bad raster format");
		xmax = head.ras_width;
		ymax = head.ras_height;
		if ((head.ras_type != RT_STANDARD
					&& head.ras_type != RT_FORMAT_RGB)
				|| head.ras_maptype != RMT_NONE
				|| head.ras_depth != 24)
			quiterr("incompatible format");
					/* put header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fprtresolu(xmax, ymax, stdout);
					/* convert file */
		pr2ra(head.ras_type, head.ras_length/ymax - xmax*3);
	} else {
					/* get header info. */
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				fgetresolu(&xmax, &ymax, stdin) < 0)
			quiterr("bad picture format");
					/* write rasterfile header */
		head.ras_magic = RAS_MAGIC;
		head.ras_width = xmax + (xmax&1);
		head.ras_height = ymax;
		head.ras_depth = 24;
		head.ras_length = head.ras_width*head.ras_height*3;
		head.ras_maptype = RMT_NONE;
		head.ras_maplength = 0;
		fwrite((char *)&head, sizeof(head), 1, stdout);
					/* convert file */
		ra2pr(head.ras_type, head.ras_length/ymax - xmax*3);
	}
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-r][-g gamma][-e +/-stops] [input [output]]\n",
			progname);
	exit(1);
}


static void
quiterr(		/* print message and exit */
	char  *err
)
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


static void
pr2ra(		/* convert 24-bit scanlines to Radiance picture */
	int	rf,
	int	pad
)
{
	COLR	*scanout;
	register int	x;
	int	y;
						/* allocate scanline */
	scanout = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanout == NULL)
		quiterr("out of memory in pr2ra");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (rf == RT_FORMAT_RGB)
			for (x = 0; x < xmax; x++) {
				scanout[x][RED] = getc(stdin);
				scanout[x][GRN] = getc(stdin);
				scanout[x][BLU] = getc(stdin);
			}
		else
			for (x = 0; x < xmax; x++) {
				scanout[x][BLU] = getc(stdin);
				scanout[x][GRN] = getc(stdin);
				scanout[x][RED] = getc(stdin);
			}
		for (x = pad; x--; getc(stdin));
		if (feof(stdin) || ferror(stdin))
			quiterr("error reading rasterfile");
		gambs_colrs(scanout, xmax);
		if (bradj)
			shiftcolrs(scanout, xmax, bradj);
		if (fwritecolrs(scanout, xmax, stdout) < 0)
			quiterr("error writing Radiance picture");
	}
						/* free scanline */
	free((void *)scanout);
}


static void
ra2pr(		/* convert Radiance scanlines to 24-bit rasterfile */
	int	rf,
	int	pad
)
{
	int	ord[3];
	COLR	*scanin;
	register int	x;
	int	y;
						/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in ra2pr");
	if (rf == RT_FORMAT_RGB) {
		ord[0] = RED; ord[1] = GRN; ord[2] = BLU;
	} else {
		ord[0] = BLU; ord[1] = GRN; ord[2] = RED;
	}
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		if (bradj)
			shiftcolrs(scanin, xmax, bradj);
		colrs_gambs(scanin, xmax);
		if (rf == RT_FORMAT_RGB)
			for (x = 0; x < xmax; x++) {
				putc(scanin[x][RED], stdout);
				putc(scanin[x][GRN], stdout);
				putc(scanin[x][BLU], stdout);
			}
		else
			for (x = 0; x < xmax; x++) {
				putc(scanin[x][BLU], stdout);
				putc(scanin[x][GRN], stdout);
				putc(scanin[x][RED], stdout);
			}
		for (x = 0; x < pad; x++)
			putc(scanin[xmax-1][ord[x%3]], stdout);
		if (ferror(stdout))
			quiterr("error writing rasterfile");
	}
						/* free scanline */
	free((void *)scanin);
}
