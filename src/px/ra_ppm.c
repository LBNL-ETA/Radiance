#ifndef lint
static const char	RCSid[] = "$Id: ra_ppm.c,v 2.11 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
 *  program to convert between RADIANCE and Poskanzer Pixmaps
 */

#include  <stdio.h>
#include  <math.h>
#include  <ctype.h>
#include  <time.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"


int  agryscan(), bgryscan(), aclrscan(), bclrscan();
int  agryscan2(), bgryscan2(), aclrscan2(), bclrscan2();
int  normval();
unsigned int	scanint(), intv(), getby2();

#define	fltv(i)		(((i)+.5)/(maxval+1.))

int  bradj = 0;				/* brightness adjustment */

double	gamcor = 2.2;			/* gamma correction value */

int  maxval = 255;			/* maximum primary value */

char  *progname;

int  xmax, ymax;


main(argc, argv)
int  argc;
char  *argv[];
{
	char  inpbuf[2];
	int  gryflag = 0;
	int  binflag = 1;
	int  reverse = 0;
	int  ptype;
	int  i;
	
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 's':
				maxval = atoi(argv[++i]) & 0xffff;
				break;
			case 'b':
				gryflag = 1;
				break;
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'a':
				binflag = 0;
				break;
			case 'r':
				reverse = !reverse;
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
	if (maxval < 256)
		setcolrgam(gamcor);
	if (reverse) {
					/* get header */
		if (read(fileno(stdin), inpbuf, 2) != 2 || inpbuf[0] != 'P')
			quiterr("input not a Poskanzer Pixmap");
		ptype = inpbuf[1];
#ifdef _WIN32
		if (ptype > 4)
			SET_FILE_BINARY(stdin);
		SET_FILE_BINARY(stdout);
#endif
		xmax = scanint(stdin);
		ymax = scanint(stdin);
		maxval = scanint(stdin);
					/* put header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fprtresolu(xmax, ymax, stdout);
					/* convert file */
		if (maxval >= 256)
			switch (ptype) {
			case '2':
				ppm2ra2(agryscan2);
				break;
			case '5':
				ppm2ra2(bgryscan2);
				break;
			case '3':
				ppm2ra2(aclrscan2);
				break;
			case '6':
				ppm2ra2(bclrscan2);
				break;
			default:
				quiterr("unsupported Pixmap type");
			}
		else
			switch (ptype) {
			case '2':
				ppm2ra(agryscan);
				break;
			case '5':
				ppm2ra(bgryscan);
				break;
			case '3':
				ppm2ra(aclrscan);
				break;
			case '6':
				ppm2ra(bclrscan);
				break;
			default:
				quiterr("unsupported Pixmap type");
			}
	} else {
#ifdef _WIN32
		SET_FILE_BINARY(stdin);
		if (binflag)
			SET_FILE_BINARY(stdout);
#endif
					/* get header info. */
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				fgetresolu(&xmax, &ymax, stdin) < 0)
			quiterr("bad picture format");
					/* write PPM header */
		printf("P%1d\n%d %d\n%u\n", (gryflag?2:3)+(binflag?3:0),
				xmax, ymax, maxval);
					/* convert file */
		if (maxval >= 256)
			ra2ppm2(binflag, gryflag);
		else
			ra2ppm(binflag, gryflag);
	}
	exit(0);
userr:
	fprintf(stderr,
		"Usage: %s [-r][-a][-b][-s maxv][-g gamma][-e +/-stops] [input [output]]\n",
			progname);
	exit(1);
}


quiterr(err)		/* print message and exit */
char  *err;
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


ppm2ra(getscan)		/* convert 1-byte Pixmap to Radiance picture */
int  (*getscan)();
{
	COLR	*scanout;
	int	y;
						/* allocate scanline */
	scanout = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanout == NULL)
		quiterr("out of memory in ppm2ra");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if ((*getscan)(scanout, xmax, stdin) < 0)
			quiterr("error reading Pixmap");
		gambs_colrs(scanout, xmax);
		if (bradj)
			shiftcolrs(scanout, xmax, bradj);
		if (fwritecolrs(scanout, xmax, stdout) < 0)
			quiterr("error writing Radiance picture");
	}
						/* free scanline */
	free((void *)scanout);
}


ra2ppm(binary, grey)	/* convert Radiance picture to Pixmap */
int  binary, grey;
{
	COLR	*scanin;
	register int	x;
	int	y;
						/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in ra2ppm");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		if (bradj)
			shiftcolrs(scanin, xmax, bradj);
		for (x = grey?xmax:0; x--; )
			scanin[x][GRN] = normbright(scanin[x]);
		colrs_gambs(scanin, xmax);
		if (grey)
			if (binary)
				for (x = 0; x < xmax; x++)
					putc(scanin[x][GRN], stdout);
			else
				for (x = 0; x < xmax; x++)
					printf("%d\n", scanin[x][GRN]);
		else
			if (binary)
				for (x = 0; x < xmax; x++) {
					putc(scanin[x][RED], stdout);
					putc(scanin[x][GRN], stdout);
					putc(scanin[x][BLU], stdout);
				}
			else
				for (x = 0; x < xmax; x++)
					printf("%d %d %d\n", scanin[x][RED],
							scanin[x][GRN],
							scanin[x][BLU]);
		if (ferror(stdout))
			quiterr("error writing Pixmap");
	}
						/* free scanline */
	free((void *)scanin);
}


ppm2ra2(getscan)	/* convert 2-byte Pixmap to Radiance picture */
int  (*getscan)();
{
	COLOR	*scanout;
	double	mult;
	int	y;
	register int	x;
						/* allocate scanline */
	scanout = (COLOR *)malloc(xmax*sizeof(COLOR));
	if (scanout == NULL)
		quiterr("out of memory in ppm2ra2");
	if (bradj)
		mult = pow(2., (double)bradj);
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if ((*getscan)(scanout, xmax, stdin) < 0)
			quiterr("error reading Pixmap");
		for (x = (gamcor>1.01)|(gamcor<0.99)?xmax:0; x--; ) {
			colval(scanout[x],RED) =
					pow(colval(scanout[x],RED), gamcor);
			colval(scanout[x],GRN) =
					pow(colval(scanout[x],GRN), gamcor);
			colval(scanout[x],BLU) =
					pow(colval(scanout[x],BLU), gamcor);
		}
		for (x = bradj?xmax:0; x--; )
			scalecolor(scanout[x], mult);
		if (fwritescan(scanout, xmax, stdout) < 0)
			quiterr("error writing Radiance picture");
	}
						/* free scanline */
	free((void *)scanout);
}


ra2ppm2(binary, grey)	/* convert Radiance picture to Pixmap (2-byte) */
int  binary, grey;
{
	COLOR	*scanin;
	double	mult, d;
	register int	x;
	int	y;
						/* allocate scanline */
	scanin = (COLOR *)malloc(xmax*sizeof(COLOR));
	if (scanin == NULL)
		quiterr("out of memory in ra2ppm2");
	if (bradj)
		mult = pow(2., (double)bradj);
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadscan(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		for (x = bradj?xmax:0; x--; )
			scalecolor(scanin[x], mult);
		for (x = grey?xmax:0; x--; )
			colval(scanin[x],GRN) = bright(scanin[x]);
		d = 1./gamcor;
		for (x = (d>1.01)|(d<0.99)?xmax:0; x--; ) {
			colval(scanin[x],GRN) = pow(colval(scanin[x],GRN), d);
			if (!grey) {
				colval(scanin[x],RED) =
						pow(colval(scanin[x],RED), d);
				colval(scanin[x],BLU) =
						pow(colval(scanin[x],BLU), d);
			}
		}
		if (grey)
			if (binary)
				for (x = 0; x < xmax; x++)
					putby2(intv(colval(scanin[x],GRN)),
							stdout);
			else
				for (x = 0; x < xmax; x++)
					printf("%u\n",
						intv(colval(scanin[x],GRN)));
		else
			if (binary)
				for (x = 0; x < xmax; x++) {
					putby2(intv(colval(scanin[x],RED)),
							stdout);
					putby2(intv(colval(scanin[x],GRN)),
							stdout);
					putby2(intv(colval(scanin[x],BLU)),
							stdout);
				}
			else
				for (x = 0; x < xmax; x++)
					printf("%u %u %u\n",
						intv(colval(scanin[x],RED)),
						intv(colval(scanin[x],GRN)),
						intv(colval(scanin[x],BLU)));
		if (ferror(stdout))
			quiterr("error writing Pixmap");
	}
						/* free scanline */
	free((void *)scanin);
}


agryscan(scan, len, fp)			/* get an ASCII greyscale scanline */
register COLR  *scan;
register int  len;
FILE  *fp;
{
	while (len-- > 0) {
		scan[0][RED] =
		scan[0][GRN] =
		scan[0][BLU] = normval(scanint(fp));
		scan++;
	}
	return(0);
}


bgryscan(scan, len, fp)			/* get a binary greyscale scanline */
register COLR  *scan;
int  len;
register FILE  *fp;
{
	register int  c;

	while (len-- > 0) {
		if ((c = getc(fp)) == EOF)
			return(-1);
		if (maxval != 255)
			c = normval(c);
		scan[0][RED] =
		scan[0][GRN] =
		scan[0][BLU] = c;
		scan++;
	}
	return(0);
}


aclrscan(scan, len, fp)			/* get an ASCII color scanline */
register COLR  *scan;
register int  len;
FILE  *fp;
{
	while (len-- > 0) {
		scan[0][RED] = normval(scanint(fp));
		scan[0][GRN] = normval(scanint(fp));
		scan[0][BLU] = normval(scanint(fp));
		scan++;
	}
	return(0);
}


bclrscan(scan, len, fp)			/* get a binary color scanline */
register COLR  *scan;
int  len;
register FILE  *fp;
{
	int  r, g, b;

	while (len-- > 0) {
		r = getc(fp);
		g = getc(fp);
		if ((b = getc(fp)) == EOF)
			return(-1);
		if (maxval == 255) {
			scan[0][RED] = r;
			scan[0][GRN] = g;
			scan[0][BLU] = b;
		} else {
			scan[0][RED] = normval(r);
			scan[0][GRN] = normval(g);
			scan[0][BLU] = normval(b);
		}
		scan++;
	}
	return(0);
}


agryscan2(scan, len, fp)		/* get an ASCII greyscale scanline */
register COLOR  *scan;
register int  len;
FILE  *fp;
{
	while (len-- > 0) {
		colval(scan[0],RED) =
		colval(scan[0],GRN) =
		colval(scan[0],BLU) = fltv(scanint(fp));
		scan++;
	}
	return(0);
}


bgryscan2(scan, len, fp)		/* get a binary greyscale scanline */
register COLOR  *scan;
int  len;
register FILE  *fp;
{
	register int  c;

	while (len-- > 0) {
		if ((c = getby2(fp)) == EOF)
			return(-1);
		colval(scan[0],RED) =
		colval(scan[0],GRN) =
		colval(scan[0],BLU) = fltv(c);
		scan++;
	}
	return(0);
}


aclrscan2(scan, len, fp)		/* get an ASCII color scanline */
register COLOR  *scan;
register int  len;
FILE  *fp;
{
	while (len-- > 0) {
		colval(scan[0],RED) = fltv(scanint(fp));
		colval(scan[0],GRN) = fltv(scanint(fp));
		colval(scan[0],BLU) = fltv(scanint(fp));
		scan++;
	}
	return(0);
}


bclrscan2(scan, len, fp)		/* get a binary color scanline */
register COLOR  *scan;
int  len;
register FILE  *fp;
{
	int  r, g, b;

	while (len-- > 0) {
		r = getby2(fp);
		g = getby2(fp);
		if ((b = getby2(fp)) == EOF)
			return(-1);
		scan[0][RED] = fltv(r);
		scan[0][GRN] = fltv(g);
		scan[0][BLU] = fltv(b);
		scan++;
	}
	return(0);
}


unsigned int
scanint(fp)			/* scan the next positive integer value */
register FILE  *fp;
{
	register int  c;
	register unsigned int  i;
tryagain:
	while (isspace(c = getc(fp)))
		;
	if (c == EOF)
		quiterr("unexpected end of file");
	if (c == '#') {		/* comment */
		while ((c = getc(fp)) != EOF && c != '\n')
			;
		goto tryagain;
	}
				/* should be integer */
	i = 0;
	do {
		if (!isdigit(c))
			quiterr("error reading integer");
		i = 10*i + c - '0';
		c = getc(fp);
	} while (c != EOF && !isspace(c));
	return(i);
}


int
normval(v)			/* normalize a value to [0,255] */
register unsigned int  v;
{
	if (v >= maxval)
		return(255);
	if (maxval == 255)
		return(v);
	return(v*255L/maxval);
}


unsigned int
getby2(fp)			/* return 2-byte quantity from fp */
register FILE	*fp;
{
	register int	lowerb, upperb;

	lowerb = getc(fp);
	upperb = getc(fp);
	if (upperb == EOF)
		return(EOF);
	return(upperb<<8 | lowerb);
}


putby2(w, fp)			/* put 2-byte quantity to fp */
register unsigned int	w;
register FILE	*fp;
{
	putc(w & 0xff, fp);
	putc(w>>8 & 0xff, fp);
	if (ferror(fp)) {
		fprintf(stderr, "%s: write error on PPM output\n", progname);
		exit(1);
	}
}


unsigned int
intv(v)				/* return integer quantity for v */
register double	v;
{
	if (v >= 0.99999)
		return(maxval);
	if (v <= 0.)
		return(0);
	return((int)(v*(maxval+1)));
}
