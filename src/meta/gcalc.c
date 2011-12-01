#ifndef lint
static const char	RCSid[] = "$Id: gcalc.c,v 1.4 2011/12/01 18:47:35 greg Exp $";
#endif
/*
 *  gcalc.c - routines to do calculations on graph files.
 *
 *     7/7/86
 *     Greg Ward Larson
 */

#include  <stdio.h>
#include  <string.h>

#include  "mgvars.h"


static double  xsum, xxsum, ysum, yysum, xysum;
static double  xmin, xmax, ymin, ymax;
static double  lastx, lasty, rsum;
static int  npts;

static void gcheader(char  *types);
static void calcpoint(int  c, double  x, double  y);
static void gcvalue(int  c, char  *types);

void
gcalc(			/* do a calculation */
	char  *types
)
{
	int  i;

	if (strchr(types, 'h') == NULL)
		gcheader(types);
	
	xmin = gparam[XMIN].flags & DEFINED ?
			varvalue(gparam[XMIN].name) :
			-1e10;
	xmax = gparam[XMAX].flags & DEFINED ?
			varvalue(gparam[XMAX].name) :
			1e10;

	for (i = 0; i < MAXCUR; i++) {
		xsum = xxsum = ysum = yysum = xysum = 0.0;
		rsum = 0.0;
		npts = 0;
		mgcurve(i, calcpoint);
		gcvalue(i, types);
	}
}


void
gcheader(			/* print header */
	register char  *types
)
{
	printf("__");
	while (*types)
		switch (*types++) {
		case 'n':
			printf("|_Label___________");
			break;
		case 'a':
			printf("|____Mean______S.D._");
			break;
		case 'm':
			printf("|_____Min_______Max_");
			break;
		case 'i':
			printf("|___Integ_");
			break;
		case 'l':
			printf("|____Slope_______Intcp____Corr_");
			break;
		default:
			break;
		}
	printf("\n");
}


void
gcvalue(		/* print the values for the curve */
	int  c,
	register char  *types
)
{
	double  d1, d2, d3, sqrt();

	if (npts < 1)
		return;

	printf("%c:", c+'A');
	while (*types)
		switch (*types++) {
		case 'n':
			printf("  %-16s", cparam[c][CLABEL].flags & DEFINED ?
					cparam[c][CLABEL].v.s : "");
			break;
		case 'a':
			d1 = sqrt((yysum - ysum*ysum/npts)/(npts-1));
			printf(" %9.2f %9.3f", ysum/npts, d1);
			break;
		case 'm':
			printf(" %9.2f %9.2f", ymin, ymax);
			break;
		case 'i':
			printf(" %9.2f", rsum);
			break;
		case 'l':
			d3 = xxsum - xsum*xsum/npts;
			d1 = (xysum - xsum*ysum/npts)/d3;
			d2 = (ysum - d1*xsum)/npts;
			d3 = d1*sqrt(d3/(yysum - ysum*ysum/npts));
			printf(" %10.6f %11.4f %7.3f", d1, d2, d3);
			break;
		default:
			break;
		}
	printf("\n");
}


void
calcpoint(		/* add a point to our calculation */
	int  c,
	double  x,
	double  y
)
{
	if (x < xmin || x > xmax)
		return;

	xsum += x;
	xxsum += x*x;
	ysum += y;
	yysum += y*y;
	xysum += x*y;

	if (npts) {
		if (y < ymin)
			ymin = y;
		else if (y > ymax)
			ymax = y;
	} else
		ymin = ymax = y;

	if (npts)
		rsum += ( y + lasty )*( x - lastx )/2.0;
	lastx = x;
	lasty = y;

	npts++;
}
