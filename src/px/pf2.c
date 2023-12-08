#ifndef lint
static const char	RCSid[] = "$Id: pf2.c,v 2.11 2023/12/08 17:56:26 greg Exp $";
#endif
/*
 *  pf2.c - routines used by pfilt.
 */

#include  "pfilt.h"
#include  "random.h"

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

static void starpoint(SCOLOR fcol, int  x, int  y, HOTPIX	 *hp);


void
pass1init(void)			/* prepare for first pass */
{
	avgbrt = 0.0;
	npix = 0;
	head = NULL;
}


void
pass1default(void)			/* for single pass */
{
	avgbrt = AVGLVL;
	npix = 1;
	head = NULL;
}


void
pass1scan(		/* process first pass scanline */
	COLORV	*scan,
	int  y
)
{
	double	cbrt;
	int  x;
	HOTPIX	 *hp;

	for (x = 0; x < xres; x++) {
	
		cbrt = (*ourbright)(scan+x*NCSAMP);

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
			scolor_color(hp->val, scan+x*NCSAMP);
			hp->x = x;
			hp->y = y;
			hp->slope = ttan(PI*(0.5-(irandom(npts)+0.5)/npts));
			hp->next = head;
			head = hp;
		}
	}
}


void
pass2init(void)			/* prepare for final pass */
{
	if (!npix) {
		fprintf(stderr, "%s: picture too dark or too bright\n",
				progname);
		quit(1);
	}
	avgbrt /= (double)npix;

	scalecolor(exposure,  AVGLVL/avgbrt);
	
	sprdfact = spread / (hotlvl * bright(exposure))
			* ((double)xres*xres + (double)yres*yres) / 4.0;
}


void
pass2scan(		/* process final pass scanline */
	COLORV	*scan,
	int  y
)
{
	int  xmin, xmax;
	int  x;
	HOTPIX	 *hp;
	
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
			starpoint(scan+x*NCSAMP, x, y, hp);
	}
	for (x = 0; x < xres; x++)
		smultcolor(scan+x*NCSAMP, exposure);
}


static void
starpoint(		/* pixel is on the star's point */
	SCOLOR  fcol,
	int  x,
	int  y,
	HOTPIX	 *hp
)
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
		saddcolor(fcol, ctmp);
	} else if (d2 > FTINY) {
		saddcolor(fcol, hp->val);
	}
}
