#ifndef lint
static const char	RCSid[] = "$Id: pflip.c,v 2.14 2025/06/07 05:09:46 greg Exp $";
#endif
/*
 * flip picture file horizontally and/or vertically
 */

#include "rtio.h"
#include "platform.h"
#include "color.h"
#include "resolu.h"

int	order;				/* input orientation */
int	xres, yres;			/* resolution (scanlen, nscans) */

long	*scanpos;			/* scanline positions */

int	fhoriz=0, fvert=0;		/* flip flags */

int	correctorder = 0;		/* correcting orientation? */

FILE	*fin;				/* input file */


static void memerr(void);
static void scanfile(void);
static void flip(void);


static int
neworder(void)			/* figure out new order from old */
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

int
main(
	int	argc,
	char	*argv[]
)
{
	static char	picfmt[MAXFMTLEN] = PICFMT;
	int	i, rval;
	
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdout);
	fixargv0(argv[0]);

	for (i = 1; i < argc; i++)
		if (!strcmp(argv[i], "-h"))
			fhoriz++;
		else if (!strcmp(argv[i], "-v"))
			fvert++;
		else if (!strcmp(argv[i], "-c"))
			correctorder++;
		else
			break;
	if (i < argc-2)
		goto userr;
	if (!fhoriz && !fvert)
		fprintf(stderr, "%s: warning - no operation\n", argv[0]);
	if (i >= argc || argv[i][0] == '-') {
		if (fvert)
			goto userr;
		SET_FILE_BINARY(stdin);
		fin = stdin;
	} else if ((fin = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", argv[i]);
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
userr:
	fprintf(stderr, "Usage: %s [-h][-v][-c] infile [outfile]\n",
			progname);
	exit(1);
}


static void
memerr(void)
{
	fprintf(stderr, "%s: out of memory\n", progname);
	exit(1);
}


static void
scanfile(void)				/* scan to the end of file */
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


static void
flip(void)					/* flip the picture */
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
