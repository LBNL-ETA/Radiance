/* RCSid: $Id: glare.h,v 2.3 2003/02/22 02:07:30 greg Exp $ */
/*
 * Common data structures for glare source finding routines
 */

#include "standard.h"
#include "view.h"
#include "color.h"
#include "setscan.h"

#define GLAREBR		7.0		/* glare source is this * avg. lum. */

#define SAMPDENS	75		/* default samples per unit in image */
#define TSAMPSTEP	10		/* sample step to compute threshold */

#define SEPS		1		/* sources this close ==> contig. */

#define SAMIN		.005		/* minimum solid angle for source */
#define MAXBUDDY	(4.*sqrt(SAMIN/PI))	/* max separation for pairing */

#define TOOSMALL(s)	((s)->brt*(s)->dom < threshold*SAMIN)

#define SABIG		.025		/* solid angle of splittable source */
#define LCORR		.12		/* linearity of splittable source */

extern VIEW	ourview;		/* our view */
extern VIEW	pictview;		/* picture view */
extern VIEW	leftview, rightview;	/* leftmost and rightmost views */

extern int	verbose;		/* verbose reporting */
extern char	*progname;		/* global argv[0] */

extern double	threshold;		/* threshold value for glare sources */

extern int	sampdens;		/* sample density */
extern ANGLE	glarang[];		/* glare calculation angles */
extern int	nglarangs;
extern double	maxtheta;		/* maximum glare angle (in radians) */
extern int	hsize;			/* horizontal size */

#define nglardirs	(2*nglarangs+1)
#define vsize		(sampdens-1)
#define hscale(v)	sqrt((double)(sampdens*sampdens - (v)*(v)))
#define hlim(v)		(int)(maxtheta*hscale(v))
#define h_theta(h,v)	(-(h)/hscale(v))

extern struct illum {
	float	theta;		/* glare direction */
	float	lcos, lsin;	/* cosine and sine to left view */
	float	rcos, rsin;	/* cosine and sine to right view */
	double	sum;		/* sum of indirect luminances */
	double	n;		/* number of values in sum */
} *indirect;		/* array of indirect illuminances */

struct srcspan {
	short	v;		/* vertical position */
	short	l, r;		/* left and right horizontal limits */
	float	brsum;		/* sum of brightnesses for this span */
	struct srcspan	*next;	/* next source span in list */
};

extern struct source {
	FVECT	dir;		/* source direction */
	double	dom;		/* solid angle of source */
	double	brt;		/* average source brightness */
	struct srcspan	*first;	/* first span for this source */
	struct source	*next;	/* next source in list */
} *donelist;			/* finished sources */

extern double	getviewpix();
extern double	pixsize();

extern long	npixinvw;	/* number of samples in view */
extern long	npixmiss;	/* number of samples missing */
