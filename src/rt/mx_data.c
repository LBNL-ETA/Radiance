#ifndef lint
static const char	RCSid[] = "$Id: mx_data.c,v 2.7 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 *  mx_data.c - routine for stored mixtures.
 */

#include "copyright.h"

#include  "ray.h"

#include  "data.h"

#include  "func.h"

/*
 *  A stored mixture is specified:
 *
 *	modifier mixdata name
 *	6+ foremod backmod func dfname vfname v0 v1 .. xf
 *	0
 *	n A1 A2 ..
 *
 *  A picture mixture is specified as:
 *
 *	modifier mixpict name
 *	7+ foremod backmod func pfname vfname vx vy xf
 *	0
 *	n A1 A2 ..
 *
 *
 *  Vfname is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  Dfname is the data file.
 *  (Pfname is a picture file.)
 *  The dimensions of the data files and the number
 *  of variables must match.  The func is a single argument
 *  function in the case of mixdata (three argument in the case
 *  of mixpict), which returns the corrected data value given the
 *  interpolated value from the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 */


mx_data(m, r)			/* interpolate mixture data */
register OBJREC  *m;
RAY  *r;
{
	OBJECT	obj;
	double  coef;
	double  pt[MAXDIM];
	DATARRAY  *dp;
	OBJECT  mod[2];
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 6)
		objerror(m, USER, "bad # arguments");
	obj = objndx(m);
	for (i = 0; i < 2; i++)
		if (!strcmp(m->oargs.sarg[i], VOIDID))
			mod[i] = OVOID;
		else if ((mod[i] = lastmod(obj, m->oargs.sarg[i])) == OVOID) {
			sprintf(errmsg, "undefined modifier \"%s\"",
					m->oargs.sarg[i]);
			objerror(m, USER, errmsg);
		}
	dp = getdata(m->oargs.sarg[3]);
	i = (1 << dp->nd) - 1;
	mf = getfunc(m, 4, i<<5, 0);
	setfunc(m, r);
	errno = 0;
	for (i = 0; i < dp->nd; i++) {
		pt[i] = evalue(mf->ep[i]);
		if (errno)
			goto computerr;
	}
	coef = datavalue(dp, pt);
	errno = 0;
	coef = funvalue(m->oargs.sarg[2], 1, &coef);
	if (errno)
		goto computerr;
	if (raymixture(r, mod[0], mod[1], coef)) {
		if (m->omod != OVOID)
			objerror(m, USER, "inappropriate modifier");
		return(1);
	}
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


mx_pdata(m, r)			/* interpolate mixture picture */
register OBJREC  *m;
RAY  *r;
{
	OBJECT	obj;
	double	col[3], coef;
	double  pt[MAXDIM];
	DATARRAY  *dp;
	OBJECT  mod[2];
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 7)
		objerror(m, USER, "bad # arguments");
	obj = objndx(m);
	for (i = 0; i < 2; i++)
		if (!strcmp(m->oargs.sarg[i], VOIDID))
			mod[i] = OVOID;
		else if ((mod[i] = lastmod(obj, m->oargs.sarg[i])) == OVOID) {
			sprintf(errmsg, "undefined modifier \"%s\"",
					m->oargs.sarg[i]);
			objerror(m, USER, errmsg);
		}
	dp = getpict(m->oargs.sarg[3]);
	mf = getfunc(m, 4, 0x3<<5, 0);
	setfunc(m, r);
	errno = 0;
	pt[1] = evalue(mf->ep[0]);	/* y major ordering */
	pt[0] = evalue(mf->ep[1]);
	if (errno)
		goto computerr;
	for (i = 0; i < 3; i++)		/* get pixel from picture */
		col[i] = datavalue(dp+i, pt);
	errno = 0;			/* evaluate function on pixel */
	coef = funvalue(m->oargs.sarg[2], 3, col);
	if (errno)
		goto computerr;
	if (raymixture(r, mod[0], mod[1], coef)) {
		if (m->omod != OVOID)
			objerror(m, USER, "inappropriate modifier");
		return(1);
	}
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}
