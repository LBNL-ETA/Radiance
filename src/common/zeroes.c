#ifndef lint
static const char	RCSid[] = "$Id: zeroes.c,v 2.4 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 *  zeroes.c - compute roots for various equations.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"


#include  <math.h>

#include  "fvect.h"


int
quadratic(r, a, b, c)		/* find real roots of quadratic equation */
double  *r;			/* roots in ascending order */
double  a, b, c; 
{
	double  disc;
	int  first;

	if (a < -FTINY)
		first = 1;
	else if (a > FTINY)
		first = 0;
	else if (fabs(b) > FTINY) {	/* solve linearly */
		r[0] = -c/b;
		return(1);
	} else
		return(0);		/* equation is c == 0 ! */
		
	b *= 0.5;			/* simplifies formula */
	
	disc = b*b - a*c;		/* discriminant */

	if (disc < -FTINY*FTINY)	/* no real roots */
		return(0);

	if (disc <= FTINY*FTINY) {	/* double root */
		r[0] = -b/a;
		return(1);
	}
	
	disc = sqrt(disc);

	r[first] = (-b - disc)/a;
	r[1-first] = (-b + disc)/a;

	return(2);
}
