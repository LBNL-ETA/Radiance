#ifndef lint
static const char	RCSid[] = "$Id: cgraph.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  cgraph.c - routines for sending graphs to tty's.
 *
 *     Greg Ward
 *     7/7/86
 */


#include  <stdio.h>

#include  "mgvars.h"


#define  FHUGE  1e10

#define  MAXSIZE  10000		/* Maximum size in characters of output */

extern char  *progname;			/* argv[0] */

static char  outcarr[MAXSIZE];		/* character output array */

static double  xmin = XMAX, xmax = XMIN,	/* extrema */
		ymin = YMAX, ymax = YMIN;

static int  dwidth, dlength;		/* device width and length */

static int  nplottable;		/* number of plottable points */


cgraph(width, length)		/* do a character graph to stdout */
int  width, length;
{
	if (width * length > MAXSIZE) {
		fprintf(stderr, "%s: page too big\n", progname);
		quit(1);
	}
	dwidth = width;
	dlength = length;
	climits();			/* get min & max values */
	cplot();			/* do character plot */
}


climits()			/* get min & max values */
{
	int  i, cstretch();

	xmin = gparam[XMIN].flags & DEFINED ?
			varvalue(gparam[XMIN].name) :
			FHUGE ;
	xmax = gparam[XMAX].flags & DEFINED ?
			varvalue(gparam[XMAX].name) :
			-FHUGE ;
	ymin = gparam[YMIN].flags & DEFINED ?
			varvalue(gparam[YMIN].name) :
			FHUGE ;
	ymax = gparam[YMAX].flags & DEFINED ?
			varvalue(gparam[YMAX].name) :
			-FHUGE ;

	nplottable = 0;
	for (i = 0; i < MAXCUR; i++)
		mgcurve(i, cstretch);

	if (nplottable == 0) {
		fprintf(stderr, "%s: no plottable data\n", progname);
		quit(1);
	}
	printf("XMIN= %f  XMAX= %f  YMIN= %f  YMAX= %f\n",
			xmin, xmax, ymin, ymax);
}


cstretch(c, x, y)			/* stretch our boundaries */
int  c;
double  x, y;
{
	if (gparam[XMIN].flags & DEFINED &&
			x < xmin)
		return;
	if (gparam[XMAX].flags & DEFINED &&
			x > xmax)
		return;
	if (gparam[YMIN].flags & DEFINED &&
			y < ymin)
		return;
	if (gparam[YMAX].flags & DEFINED &&
			y > ymax)
		return;

	if (x < xmin)
		xmin = x;
	if (x > xmax)
		xmax = x;
	if (y < ymin)
		ymin = y;
	if (y > ymax)
		ymax = y;
		
	nplottable++;
}


cplot()				/* do character  plot */
{
	int  cpoint();
	int  i, j;
	register char  *op;

	for (op = outcarr+dlength*dwidth; op > outcarr; )
		*--op = ' ';

	for (i = 0; i < MAXCUR; i++)
		mgcurve(i, cpoint);

	for (i = 0; i < dlength; i++) {
		for (j = 0; j < dwidth; j++)
			putchar(*op++);
		putchar('\n');
	}
}


cpoint(c, x, y)			/* store a point */
int  c;
double  x, y;
{
	register int  ndx;

	if (x < xmin || x > xmax || y < ymin || y > ymax)
		return;

	ndx = (dlength-1)*(1.0 - (y - ymin)/(ymax - ymin)) + 0.5;
	ndx = dwidth*ndx + (dwidth-1)*(x-xmin)/(xmax-xmin) + 0.5;

	if (outcarr[ndx] == ' ')
		outcarr[ndx] = c+'A';
	else if (outcarr[ndx] > '1' && outcarr[ndx] < '9')
		outcarr[ndx]++;
	else if (outcarr[ndx] == '9')
		outcarr[ndx] = '*';
	else if (outcarr[ndx] != '*')
		outcarr[ndx] = '2';
}
