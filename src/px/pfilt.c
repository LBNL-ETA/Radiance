#ifndef lint
static const char RCSid[] = "$Id: pfilt.c,v 2.27 2003/10/22 02:06:35 greg Exp $";
#endif
/*
 *  pfilt.c - program to post-process picture file.
 *
 *     9/26/85
 *     6/23/93	Added additional buffers for value spreading
 */

#include  "copyright.h"

#include  <signal.h>
#include  <string.h>

#include  "platform.h"
#include  "standard.h"
#include  "color.h"
#include  "view.h"
#include  "paths.h"

extern float  *matchlamp();

#define	 FEQ(a,b)	((a) >= .98*(b) && (a) <= 1.02*(b))

double	 CHECKRAD = 2.0;	/* radius to check for filtering */

#define  THRESHRAD	5.0	/* maximum sample spread in output */

COLOR  exposure = WHTCOLOR;	/* exposure for the frame */

double	rad = 0.0;		/* output pixel radius for filtering */

double  thresh = 0.0;		/* maximum contribution for subpixel */

int  nrows = 0;			/* number of rows for output */
int  ncols = 0;			/* number of columns for output */

double	x_c = 1.0;		/* ratio of output x size to input */
double	y_r = 1.0;		/* ratio of output y size to input */

int  singlepass = 0;		/* true means skip first pass */

int  avghot = 0;		/* true means average in bright spots */

double	hotlvl = 100.0;		/* level considered "hot" */

int  npts = 0;			/* (half) number of points for stars */

double	spread = 1e-4;		/* spread for star points */

char  *tfname = NULL;

char  template[] = TEMPLATE;

char  *lampdat = "lamp.tab";	/* lamp data file */

int  order;			/* scanline ordering of input */
int  xres, yres;		/* resolution of input */
double	inpaspect = 1.0;	/* pixel aspect ratio of input */
int  correctaspect = 0;		/* aspect ratio correction? */

int  wrongformat = 0;

VIEW  ourview = STDVIEW;
int  gotview = 0;
int  wrapfilt = 0;		/* wrap filter horizontally? */

int  estatus = 0;		/* exit status (for non-fatal errors) */

int  xrad;			/* x search radius */
int  yrad;			/* y search radius */
int  xbrad;			/* x box size */
int  ybrad;			/* y box size */

int  barsize;			/* size of input scan bar */
COLOR  **scanin;		/* input scan bar */
COLOR  *scanout;		/* output scan line */
COLOR  **scoutbar;		/* output scan bar (if thresh > 0) */
float  **greybar;		/* grey-averaged input values */
int  obarsize = 0;		/* size of output scan bar */
int  orad = 0;			/* output window radius */

char  *progname;


main(argc, argv)
int  argc;
char  **argv;
{
	extern int  headline();
	FILE  *fin;
	float  *lampcolor;
	char  *lamptype = NULL;
	long  fpos;
	double	outaspect = 0.0;
	double	d;
	int  i, j;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	if (signal(SIGINT, quit) == SIG_IGN)
		signal(SIGINT, SIG_IGN);
#ifdef SIGHUP
	if (signal(SIGHUP, quit) == SIG_IGN)
		signal(SIGHUP, SIG_IGN);
#endif
	signal(SIGTERM, quit);
#ifdef SIGPIPE
	signal(SIGPIPE, quit);
#endif
#ifdef	SIGXCPU
	signal(SIGXCPU, quit);
	signal(SIGXFSZ, quit);
#endif

	progname = argv[0] = fixargv0(argv[0]);

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'x':
				i++;
				if (argv[i][0] == '/') {
					x_c = 1.0/atof(argv[i]+1);
					ncols = 0;
				} else
					ncols = atoi(argv[i]);
				break;
			case 'y':
				i++;
				if (argv[i][0] == '/') {
					y_r = 1.0/atof(argv[i]+1);
					nrows = 0;
				} else
					nrows = atoi(argv[i]);
				break;
			case 'c':
				correctaspect = !correctaspect;
				break;
			case 'p':
				i++;
				outaspect = atof(argv[i]);
				break;
			case 'e':
				if (argv[i+1][0] == '+' || argv[i+1][0] == '-')
					d = pow(2.0, atof(argv[i+1]));
				else
					d = atof(argv[i+1]);
				if (d < 1e-20 || d > 1e20) {
					fprintf(stderr,
						"%s: exposure out of range\n",
							argv[0]);
					quit(1);
				}
				switch (argv[i][2]) {
				case '\0':
					scalecolor(exposure, d);
					break;
				case 'r':
					colval(exposure,RED) *= d;
					break;
				case 'g':
					colval(exposure,GRN) *= d;
					break;
				case 'b':
					colval(exposure,BLU) *= d;
					break;
				default:
					goto badopt;
				}
				i++;
				break;
			case 'f':
				lampdat = argv[++i];
				break;
			case 't':
				lamptype = argv[++i];
				break;
			case '1':
				singlepass = 1;
				break;
			case '2':
				singlepass = 0;
				break;
			case 'n':
				npts = atoi(argv[++i]) / 2;
				break;
			case 's':
				spread = atof(argv[++i]);
				break;
			case 'a':
				avghot = !avghot;
				break;
			case 'h':
				hotlvl = atof(argv[++i]);
				break;
			case 'r':
				rad = atof(argv[++i]);
				break;
			case 'm':
				thresh = atof(argv[++i]);
				if (rad <= FTINY)
					rad = 0.6;
				break;
			case 'b':
				rad = thresh = 0.0;
				break;
			default:;
			badopt:
				fprintf(stderr, "%s: unknown option: %s\n",
						progname, argv[i]);
				quit(1);
				break;
			}
		else
			break;
					/* get lamp data (if necessary) */
	if (lamptype != NULL) {
		if (loadlamps(lampdat) < 0)
			quit(1);
		if ((lampcolor = matchlamp(lamptype)) == NULL) {
			fprintf(stderr, "%s: unknown lamp type\n", lamptype);
			quit(1);
		}
		for (j = 0; j < 3; j++)
			if (lampcolor[j] > 1e-4)
				colval(exposure,j) /= lampcolor[j];
		freelamps();
	}
					/* open input file */
	if (i == argc) {
		if (singlepass)
			fin = stdin;
		else {
			tfname = mktemp(template);
			if ((fin = fopen(tfname, "w+")) == NULL) {
				fprintf(stderr, "%s: can't create ", progname);
				fprintf(stderr, "temp file \"%s\"\n", tfname);
				quit(1);
			}
			copyfile(stdin, fin);
			if (fseek(fin, 0L, 0) == -1) {
				fprintf(stderr, "%s: seek fail\n", progname);
				quit(1);
			}
		}
	} else if (i == argc-1) {
		if ((fin = fopen(argv[i], "r")) == NULL) {
			fprintf(stderr, "%s: can't open file \"%s\"\n",
						progname, argv[i]);
			quit(1);
		}
	} else {
		fprintf(stderr, "%s: bad # file arguments\n", progname);
		quit(1);
	}
					/* get header */
	getheader(fin, headline, NULL);
	if (wrongformat) {
		fprintf(stderr, "%s: input must be a Radiance picture\n",
				progname);
		quit(1);
	}
					/* add new header info. */
	printargs(i, argv, stdout);
					/* get picture size */
	if ((order = fgetresolu(&xres, &yres, fin)) < 0) {
		fprintf(stderr, "%s: bad picture size\n", progname);
		quit(1);
	}
	if (!(order & YMAJOR))
		inpaspect = 1.0/inpaspect;
					/* wrap around for cylindrical view? */
	wrapfilt = gotview && ourview.type == VT_CYL &&
			ourview.horiz >= 360.-FTINY && order & YMAJOR;
					/* compute output resolution */
	if (ncols <= 0)
		ncols = x_c*xres + .5;
	if (nrows <= 0)
		nrows = y_r*yres + .5;
	if (outaspect > .01) {
		d = inpaspect * yres/xres / outaspect;
		if (d * ncols > nrows)
			ncols = nrows / d;
		else
			nrows = ncols * d;
	}
	x_c = (double)ncols/xres;
	y_r = (double)nrows/yres;

	if (singlepass) {		/* skip exposure, etc. */
		pass1default();
		pass2(fin);
		quit(0);
	}

	fpos = ftell(fin);		/* save input file position */
	
	pass1(fin);

	if (fseek(fin, fpos, 0) == -1) {
		fprintf(stderr, "%s: seek fail\n", progname);
		quit(1);
	}
	pass2(fin);

	quit(estatus);
}


double
rgb_bright(clr)
COLOR  clr;
{
	return(bright(clr));
}


double
xyz_bright(clr)
COLOR  clr;
{
	return(clr[CIEY]);
}


double	(*ourbright)() = rgb_bright;


int
headline(s)				/* process line from header */
char  *s;
{
	char  fmt[32];

	fputs(s, stdout);		/* copy to output */
	if (isaspect(s))		/* get aspect ratio */
		inpaspect *= aspectval(s);
	else if (isexpos(s))		/* get exposure */
		hotlvl *= exposval(s);
	else if (formatval(fmt, s)) {	/* get format */
		wrongformat = 0;
		if (!strcmp(COLRFMT, fmt))
			ourbright = rgb_bright;
		else if (!strcmp(CIEFMT, fmt))
			ourbright = xyz_bright;
		else
			wrongformat = !globmatch(PICFMT, fmt);
	} else if (isview(s) && sscanview(&ourview, s) > 0)
		gotview++;
	return(0);
}


copyfile(in, out)			/* copy a file */
register FILE  *in, *out;
{
	register int  c;

	while ((c = getc(in)) != EOF)
		putc(c, out);

	if (ferror(out)) {
		fprintf(stderr, "%s: write error in copyfile\n", progname);
		quit(1);
	}
}


pass1(in)				/* first pass of picture file */
FILE  *in;
{
	int  i;
	COLOR  *scan;

	pass1init();

	scan = (COLOR *)malloc(xres*sizeof(COLOR));
	if (scan == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		quit(1);
	}
	for (i = 0; i < yres; i++) {
		if (freadscan(scan, xres, in) < 0) {
			nrows = (long)nrows * i / yres;	/* adjust frame */
			if (nrows <= 0) {
				fprintf(stderr, "%s: empty frame\n", progname);
				quit(1);
			}
			fprintf(stderr, "%s: warning - partial frame (%d%%)\n",
					progname, (int)(100L*i/yres));
			yres = i;
			y_r = (double)nrows/yres;
			estatus++;
			break;
		}
		pass1scan(scan, i);
	}
	free((void *)scan);
}


pass2(in)			/* last pass on file, write to stdout */
FILE  *in;
{
	int  yread;
	int  ycent, xcent;
	int  r, c;
	
	pass2init();
	scan2init();
	yread = 0;
	for (r = 0; r < nrows; r++) {
		ycent = (r+.5)*yres/nrows;
		while (yread <= ycent+yrad) {
			if (yread < yres) {
				if (freadscan(scanin[yread%barsize],
						xres, in) < 0) {
					fprintf(stderr,
						"%s: truncated input (y=%d)\n",
						progname, yres-1-yread);
					quit(1);
				}
				pass2scan(scanin[yread%barsize], yread);
			}
			yread++;
		}
		if (obarsize > 0)
			scan2sync(r);
		for (c = 0; c < ncols; c++) {
			xcent = (c+.5)*xres/ncols;
			if (thresh > FTINY)
				dothresh(xcent, ycent, c, r);
			else if (rad > FTINY)
				dogauss(scanout[c], xcent, ycent, c, r);
			else
				dobox(scanout[c], xcent, ycent, c, r);
		}
		if (scanout != NULL && fwritescan(scanout, ncols, stdout) < 0) {
			fprintf(stderr, "%s: write error in pass2\n", progname);
			quit(1);
		}
	}
					/* skip leftover input */
	while (yread < yres) {
		if (freadscan(scanin[0], xres, in) < 0)
			break;
		yread++;
	}
	scan2flush();			/* flush output */
}


scan2init()			/* prepare scanline arrays */
{
	COLOR	ctmp;
	double	d;
	register int  i;

	xbrad = xres/ncols/2 + 1;
	ybrad = yres/nrows/2 + 1;
	if (rad > FTINY) {
		if (nrows >= yres && ncols >= xres)
			rad *= (y_r + x_c)/2.0;

		if (thresh > FTINY) {
			orad = CHECKRAD*THRESHRAD*rad + 1;
			xrad = orad/x_c + xbrad;
			yrad = orad/y_r + ybrad;
			obarsize = 2*orad + 1;
		} else {
			xrad = CHECKRAD*rad/x_c + 1;
			yrad = CHECKRAD*rad/y_r + 1;
		}
		initmask();		/* initialize filter table */
	} else {
		xrad = xbrad;
		yrad = ybrad;
	}
	barsize = 2*yrad + 1;
	scanin = (COLOR **)malloc(barsize*sizeof(COLOR *));
	if (scanin == NULL)
		goto memerr;
	for (i = 0; i < barsize; i++) {
		scanin[i] = (COLOR *)malloc(xres*sizeof(COLOR));
		if (scanin[i] == NULL)
			goto memerr;
	}
	if (obarsize > 0) {
		scoutbar = (COLOR **)malloc(obarsize*sizeof(COLOR *));
		greybar = (float **)malloc(obarsize*sizeof(float *));
		if ((scoutbar == NULL) | (greybar == NULL))
			goto memerr;
		for (i = 0; i < obarsize; i++) {
			scoutbar[i] = (COLOR *)malloc(ncols*sizeof(COLOR));
			greybar[i] = (float *)malloc(ncols*sizeof(float));
			if ((scoutbar[i] == NULL) | (greybar[i] == NULL))
				goto memerr;
		}
	} else {
		scanout = (COLOR *)malloc(ncols*sizeof(COLOR));
		if (scanout == NULL)
			goto memerr;
	}
					/* record pixel aspect ratio */
	if (!correctaspect) {
		d = order & YMAJOR ? x_c/y_r : y_r/x_c ;
		if (!FEQ(d,1.0))
			fputaspect(d, stdout);
	}
					/* record exposure */
	d = (*ourbright)(exposure);
	if (!FEQ(d,1.0))
		fputexpos(d, stdout);
					/* record color correction */
	copycolor(ctmp, exposure);
	scalecolor(ctmp, 1.0/d);
	if (!FEQ(colval(ctmp,RED),colval(ctmp,GRN)) ||
			!FEQ(colval(ctmp,GRN),colval(ctmp,BLU)))
		fputcolcor(ctmp, stdout);
	printf("\n");
					/* write out resolution */
	fputresolu(order, ncols, nrows, stdout);
	return;
memerr:
	fprintf(stderr, "%s: out of memory\n", progname);
	quit(1);
}


scan2sync(r)			/* synchronize grey averages and output scan */
int  r;
{
	static int  nextrow = 0;
	COLOR  ctmp;
	int  ybot;
	register int  c;
					/* average input scanlines */
	while (nextrow <= r+orad && nextrow < nrows) {
		ybot = (nextrow+.5)*yres/nrows;
		for (c = 0; c < ncols; c++) {
			dobox(ctmp, (int)((c+.5)*xres/ncols),ybot, c,nextrow);
			greybar[nextrow%obarsize][c] = (*ourbright)(ctmp);
		}
					/* and zero output scanline */
		memset((char *)scoutbar[nextrow%obarsize], '\0', ncols*sizeof(COLOR));
		nextrow++;
	}
					/* point to top scanline for output */
	if (r-orad >= 0)
		scanout = scoutbar[(r-orad)%obarsize];
	else
		scanout = NULL;
}


scan2flush()			/* flush output buffer */
{
	register int  r;

	for (r = nrows-orad; r < nrows; r++)
		if (fwritescan(scoutbar[r%obarsize], ncols, stdout) < 0)
			break;
	if (fflush(stdout) < 0) {
		fprintf(stderr, "%s: write error at end of pass2\n", progname);
		quit(1);
	}
}


void
quit(code)		/* remove temporary file and exit */
int  code;
{
	if (tfname != NULL)
		unlink(tfname);
	exit(code);
}
