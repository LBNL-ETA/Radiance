#ifndef lint
static const char	RCSid[] = "$Id: pcompos.c,v 2.21 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  pcompos.c - program to composite pictures.
 *
 *     6/30/87
 */

#include  <stdio.h>

#include  <math.h>

#ifdef MSDOS
#include  <fcntl.h>
#endif

#include  <time.h>

#include  "color.h"

#include  "resolu.h"

#define	 MAXFILE	64

#define  HASMIN		1
#define  HASMAX		2

					/* output picture size */
int  xsiz = 0;
int  ysiz = 0;
					/* input dimensions */
int  xmin = 0;
int  ymin = 0;
int  xmax = 0;
int  ymax = 0;

COLR  bgcolr = BLKCOLR;			/* background color */

int  labelht = 24;			/* label height */

int  checkthresh = 0;			/* check threshold value */

char  Command[] = "<Command>";
char  Label[] = "<Label>";

char  *progname;

struct {
	char  *name;			/* file or command name */
	FILE  *fp;			/* stream pointer */
	int  xres, yres;		/* picture size */
	int  xloc, yloc;		/* anchor point */
	int  flags;			/* HASMIN, HASMAX */
	COLR  thmin, thmax;		/* thresholds */
} input[MAXFILE];		/* our input files */

int  nfile;			/* number of files */

char  ourfmt[LPICFMT+1] = PICFMT;
int  wrongformat = 0;

FILE  *lblopen();
void  quit();


tabputs(s)			/* print line preceded by a tab */
char  *s;
{
	char  fmt[32];

	if (isheadid(s))
		return(0);
	if (formatval(fmt, s)) {
		if (globmatch(ourfmt, fmt)) {
			wrongformat = 0;
			strcpy(ourfmt, fmt);
		} else
			wrongformat = 1;
	} else {
		putc('\t', stdout);
		fputs(s, stdout);
	}
	return(0);
}


main(argc, argv)
int  argc;
char  *argv[];
{
	int  ncolumns = 0;
	int  autolabel = 0;
	int  curcol = 0, x0 = 0, curx = 0, cury = 0, spacing = 0;
	int  xsgn, ysgn;
	char  *thislabel;
	int  an;
#ifdef MSDOS
	extern int  _fmode;
	_fmode = O_BINARY;
	setmode(fileno(stdin), O_BINARY);
	setmode(fileno(stdout), O_BINARY);
#endif
	progname = argv[0];

	for (an = 1; an < argc && argv[an][0] == '-'; an++)
		switch (argv[an][1]) {
		case 'x':
			xsiz = atoi(argv[++an]);
			break;
		case 'y':
			ysiz = atoi(argv[++an]);
			break;
		case 'b':
			setcolr(bgcolr, atof(argv[an+1]),
					atof(argv[an+2]),
					atof(argv[an+3]));
			an += 3;
			break;
		case 'a':
			ncolumns = atoi(argv[++an]);
			break;
		case 's':
			spacing = atoi(argv[++an]);
			break;
		case 'o':
			curx = x0 = atoi(argv[++an]);
			cury = atoi(argv[++an]);
			break;
		case 'l':
			switch (argv[an][2]) {
			case 'a':
				autolabel++;
				break;
			case 'h':
				labelht = atoi(argv[++an]);
				break;
			case '\0':
				goto dofiles;
			default:
				goto userr;
			}
			break;
		case '\0':
		case 't':
			goto dofiles;
		default:
			goto userr;
		}
dofiles:
	newheader("RADIANCE", stdout);
	fputnow(stdout);
	for (nfile = 0; an < argc; nfile++) {
		if (nfile >= MAXFILE)
			goto toomany;
		thislabel = NULL;
		input[nfile].flags = 0;
		xsgn = ysgn = '-';
		while (an < argc && (argv[an][0] == '-' || argv[an][0] == '+'
				|| argv[an][0] == '=')) {
			switch (argv[an][1]) {
			case 't':
				checkthresh = 1;
				if (argv[an][0] == '-') {
					input[nfile].flags |= HASMIN;
					setcolr(input[nfile].thmin,
							atof(argv[an+1]),
							atof(argv[an+1]),
							atof(argv[an+1]));
				} else if (argv[an][0] == '+') {
					input[nfile].flags |= HASMAX;
					setcolr(input[nfile].thmax,
							atof(argv[an+1]),
							atof(argv[an+1]),
							atof(argv[an+1]));
				} else
					goto userr;
				an++;
				break;
			case 'l':
				if (strcmp(argv[an], "-l"))
					goto userr;
				thislabel = argv[++an];
				break;
			case '+':
			case '-':
			case '0':
				if (argv[an][0] != '=')
					goto userr;
				xsgn = argv[an][1];
				ysgn = argv[an][2];
				if (ysgn != '+' && ysgn != '-' && ysgn != '0')
					goto userr;
				break;
			case '\0':
				if (argv[an][0] == '-')
					goto getfile;
				goto userr;
			default:
				goto userr;
			}
			an++;
		}
getfile:
		if (argc-an < (ncolumns ? 1 : 3))
			goto userr;
		if (autolabel && thislabel == NULL)
			thislabel = argv[an];
		if (!strcmp(argv[an], "-")) {
			input[nfile].name = "<stdin>";
			input[nfile].fp = stdin;
		} else {
			if (argv[an][0] == '!') {
				input[nfile].name = Command;
				input[nfile].fp = popen(argv[an]+1, "r");
			} else {
				input[nfile].name = argv[an];
				input[nfile].fp = fopen(argv[an], "r");
			}
			if (input[nfile].fp == NULL) {
				perror(argv[an]);
				quit(1);
			}
		}
		an++;
						/* get header */
		printf("%s:\n", input[nfile].name);
		getheader(input[nfile].fp, tabputs, NULL);
		if (wrongformat) {
			fprintf(stderr, "%s: incompatible input format\n",
					input[nfile].name);
			quit(1);
		}
						/* get picture size */
		if (fgetresolu(&input[nfile].xres, &input[nfile].yres,
				input[nfile].fp) < 0) {
			fprintf(stderr, "%s: bad picture size\n",
					input[nfile].name);
			quit(1);
		}
		if (ncolumns > 0) {
			if (curcol >= ncolumns) {
				cury = ymax + spacing;
				curx = x0;
				curcol = 0;
			}
			input[nfile].xloc = curx;
			input[nfile].yloc = cury;
			curx += input[nfile].xres + spacing;
			curcol++;
		} else {
			input[nfile].xloc = atoi(argv[an++]);
			if (xsgn == '+')
				input[nfile].xloc -= input[nfile].xres;
			else if (xsgn == '0')
				input[nfile].xloc -= input[nfile].xres/2;
			input[nfile].yloc = atoi(argv[an++]);
			if (ysgn == '+')
				input[nfile].yloc -= input[nfile].yres;
			else if (ysgn == '0')
				input[nfile].yloc -= input[nfile].yres/2;
		}
		if (input[nfile].xloc < xmin)
			xmin = input[nfile].xloc;
		if (input[nfile].yloc < ymin)
			ymin = input[nfile].yloc;
		if (input[nfile].xloc+input[nfile].xres > xmax)
			xmax = input[nfile].xloc+input[nfile].xres;
		if (input[nfile].yloc+input[nfile].yres > ymax)
			ymax = input[nfile].yloc+input[nfile].yres;
		if (thislabel != NULL) {
			if (++nfile >= MAXFILE)
				goto toomany;
			input[nfile].name = Label;
			input[nfile].flags = 0;
			input[nfile].xres = input[nfile-1].xres;
			input[nfile].yres = labelht;
			if ((input[nfile].fp = lblopen(thislabel,
					&input[nfile].xres,
					&input[nfile].yres)) == NULL)
				goto labelerr;
			input[nfile].xloc = input[nfile-1].xloc;
			input[nfile].yloc = input[nfile-1].yloc +
					input[nfile-1].yres-input[nfile].yres;
		}
	}
	if (xsiz <= 0)
		xsiz = xmax;
	else if (xsiz > xmax)
		xmax = xsiz;
	if (ysiz <= 0)
		ysiz = ymax;
	else if (ysiz > ymax)
		ymax = ysiz;
					/* add new header info. */
	printargs(argc, argv, stdout);
	if (strcmp(ourfmt, PICFMT))
		fputformat(ourfmt, stdout);	/* print format if known */
	putchar('\n');
	fprtresolu(xsiz, ysiz, stdout);

	compos();
	
	quit(0);
userr:
	fprintf(stderr,
	"Usage: %s [-x xr][-y yr][-b r g b][-a n][-s p][-o x0 y0][-la][-lh h] ",
			progname);
	fprintf(stderr, "[-t min1][+t max1][-l lab][=SS] pic1 x1 y1 ..\n");
	quit(1);
toomany:
	fprintf(stderr, "%s: only %d files and labels allowed\n",
			progname, MAXFILE);
	quit(1);
labelerr:
	fprintf(stderr, "%s: error opening label\n", progname);
	quit(1);
}


compos()				/* composite pictures */
{
	COLR  *scanin, *scanout;
	int  y;
	register int  x, i;

	scanin = (COLR *)malloc((xmax-xmin)*sizeof(COLR));
	if (scanin == NULL)
		goto memerr;
	scanin -= xmin;
	if (checkthresh) {
		scanout = (COLR *)malloc(xsiz*sizeof(COLR));
		if (scanout == NULL)
			goto memerr;
	} else
		scanout = scanin;
	for (y = ymax-1; y >= 0; y--) {
		for (x = 0; x < xsiz; x++)
			copycolr(scanout[x], bgcolr);
		for (i = 0; i < nfile; i++) {
			if (input[i].yloc > y ||
					input[i].yloc+input[i].yres <= y)
				continue;
			if (freadcolrs(scanin+input[i].xloc,
					input[i].xres, input[i].fp) < 0) {
				fprintf(stderr, "%s: read error (y==%d)\n",
						input[i].name,
						y-input[i].yloc);
				quit(1);
			}
			if (y >= ysiz)
				continue;
			if (checkthresh) {
				x = input[i].xloc+input[i].xres;
				if (x > xsiz)
					x = xsiz;
				for (x--; x >= 0 && x >= input[i].xloc; x--) {
					if (input[i].flags & HASMIN &&
					cmpcolr(scanin[x], input[i].thmin) <= 0)
						continue;
					if (input[i].flags & HASMAX &&
					cmpcolr(scanin[x], input[i].thmax) >= 0)
						continue;
					copycolr(scanout[x], scanin[x]);
				}
			}
		}
		if (y >= ysiz)
			continue;
		if (fwritecolrs(scanout, xsiz, stdout) < 0) {
			perror(progname);
			quit(1);
		}
	}
	return;
memerr:
	perror(progname);
	quit(1);
}


int
cmpcolr(c1, c2)			/* compare two colr's (improvisation) */
register COLR  c1, c2;
{
	register int  i, j;

	j = 4;				/* check exponents first! */
	while (j--)
		if (i = c1[j] - c2[j])
			return(i);
	return(0);
}


FILE *
lblopen(s, xp, yp)		/* open pipe to label generator */
char  *s;
int  *xp, *yp;
{
	char  com[128];
	FILE  *fp;

	sprintf(com, "psign -s -.15 -a 2 -x %d -y %d '%.90s'", *xp, *yp, s);
	if ((fp = popen(com, "r")) == NULL)
		return(NULL);
	if (checkheader(fp, COLRFMT, NULL) < 0)
		goto err;
	if (fgetresolu(xp, yp, fp) < 0)
		goto err;
	return(fp);
err:
	pclose(fp);
	return(NULL);
}


void
quit(code)		/* exit gracefully */
int  code;
{
	register int  i;
				/* close input files */
	for (i = 0; i < nfile; i++)
		if (input[i].name == Command || input[i].name == Label)
			pclose(input[i].fp);
		else
			fclose(input[i].fp);
	exit(code);
}
