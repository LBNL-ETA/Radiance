#ifndef lint
static const char	RCSid[] = "$Id: pf2.c,v 2.6 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  pf2.c - routines used by pfilt.
 */

#include  <stdio.h>

#include  <stdlib.h>

#include  <math.h>

#include  "random.h"

#include  "color.h"

#define	 PI		3.14159265359

#define	 FTINY		(1e-6)

extern int  nrows, ncols;	/* number of rows and columns for output */

extern int  xres, yres;		/* x and y resolution */

extern int  avghot;		/* true means average in avgbrt spots */

extern double  hotlvl;		/* brightness considered "hot" */

extern int  npts;		/* # of points for stars */

extern double  spread;		/* spread for star points */

extern char  *progname;

extern COLOR  exposure;		/* exposure for frame */

extern double  (*ourbright)();	/* brightness calculation function */

#define	 AVGLVL		0.5	/* target mean brightness */

double	avgbrt;			/* average picture brightness */
long	npix;			/* # pixels in average */

typedef struct	hotpix {	/* structure for avgbrt pixels */
	struct hotpix  *next;	/* next in list */
	COLOR  val;		/* pixel color */
	short  x, y;		/* pixel position */
	float  slope;		/* random slope for diffraction */
}  HOTPIX;

HOTPIX	*head;			/* head of avgbrt pixel list */

double	sprdfact;		/* computed spread factor */


pass1init()			/* prepare for first pass */
{
	avgbrt = 0.0;
	npix = 0;
	head = NULL;
}


pass1default()			/* for single pass */
{
	avgbrt = AVGLVL;
	npix = 1;
	head = NULL;
}


pass1scan(scan, y)		/* process first pass scanline */
register COLOR	*scan;
int  y;
{
	double	cbrt;
	register int  x;
	register HOTPIX	 *hp;

	for (x = 0; x < xres; x++) {
	
		cbrt = (*ourbright)(scan[x]);

		if (cbrt <= 0)
			continue;

		if (avghot || cbrt < hotlvl) {
			avgbrt += cbrt;
			npix++;
		}
		if (npts && cbrt >= hotlvl) {
			hp = (HOTPIX *)malloc(sizeof(HOTPIX));
			if (hp == NULL) {
				fprintf(stderr, "%s: out of memory\n",
						progname);
				quit(1);
			}
			copycolor(hp->val, scan[x]);
			hp->x = x;
			hp->y = y;
			hp->slope = tan(PI*(0.5-(random()%npts+0.5)/npts));
			hp->next = head;
			head = hp;
		}
	}
}


pass2init()			/* prepare for final pass */
{
	if (!npix) {
		fprintf(stderr, "%s: picture too dark or too bright\n",
				progname);
		quit(1);
	}
	avgbrt /= (double)npix;

	scalecolor(exposure,  AVGLVL/avgbrt);
	
	sprdfact = spread / (hotlvl * (*ourbright)(exposure))
			* ((double)xres*xres + (double)yres*yres) / 4.0;
}


pass2scan(scan, y)		/* process final pass scanline */
register COLOR	*scan;
int  y;
{
	int  xmin, xmax;
	register int  x;
	register HOTPIX	 *hp;
	
	for (hp = head; hp != NULL; hp = hp->next) {
		if (hp->slope > FTINY) {
			xmin = (y - hp->y - 0.5)/hp->slope + hp->x;
			xmax = (y - hp->y + 0.5)/hp->slope + hp->x;
		} else if (hp->slope < -FTINY) {
			xmin = (y - hp->y + 0.5)/hp->slope + hp->x;
			xmax = (y - hp->y - 0.5)/hp->slope + hp->x;
		} else if (y == hp->y) {
			xmin = 0;
			xmax = xres-1;
		} else {
			xmin = 1;
			xmax = 0;
		}
		if (xmin < 0)
			xmin = 0;
		if (xmax >= xres)
			xmax = xres-1;
		for (x = xmin; x <= xmax; x++)
			starpoint(scan[x], x, y, hp);
	}
	for (x = 0; x < xres; x++)
		multcolor(scan[x], exposure);
}


starpoint(fcol, x, y, hp)		/* pixel is on the star's point */
COLOR  fcol;
int  x, y;
register HOTPIX	 *hp;
{
	COLOR  ctmp;
	double	d2;
	
	d2 = (double)(x - hp->x)*(x - hp->x) + (double)(y - hp->y)*(y - hp->y);
	if (d2 > sprdfact) {
		d2 = sprdfact / d2;
		if (d2 < FTINY)
			return;
		copycolor(ctmp, hp->val);
		scalecolor(ctmp, d2);
		addcolor(fcol, ctmp);
	} else if (d2 > FTINY) {
		addcolor(fcol, hp->val);
	}
}
