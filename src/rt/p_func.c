#ifndef lint
static const char	RCSid[] = "$Id: p_func.c,v 2.8 2014/07/08 18:25:00 greg Exp $";
#endif
/*
 *  p_func.c - routine for procedural patterns.
 */

#include "copyright.h"

#include  "ray.h"
#include  "func.h"
#include  "rtotypes.h"

/*
 *	A procedural pattern can either be a brightness or a
 *  color function.  A brightness function is given as:
 *
 *	modifier brightfunc name
 *	2+ bvarname filename xf
 *	0
 *	n A1 A2 ..
 *
 *  A color function is given as:
 *
 *	modifier colorfunc name
 *	4+ rvarname gvarname bvarname filename xf
 *	0
 *	n A1 A2 ..
 *
 *  Filename is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 */


int
p_bfunc(			/* compute brightness pattern */
	OBJREC  *m,
	RAY  *r
)
{
	double  bval;
	MFUNC  *mf;

	if (m->oargs.nsargs < 2)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 1, 0x1, 0);
	setfunc(m, r);
	errno = 0;
	bval = evalue(mf->ep[0]);
	if (errno == EDOM || errno == ERANGE) {
		objerror(m, WARNING, "compute error");
		return(0);
	}
	scalecolor(r->pcol, bval);
	return(0);
}


int
p_cfunc(			/* compute color pattern */
	OBJREC  *m,
	RAY  *r
)
{
	COLOR  cval;
	MFUNC  *mf;

	if (m->oargs.nsargs < 4)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 3, 0x7, 0);
	setfunc(m, r);
	errno = 0;
	setcolor(cval, evalue(mf->ep[0]),
			evalue(mf->ep[1]),
			evalue(mf->ep[2]));
	if (errno == EDOM || errno == ERANGE) {
		objerror(m, WARNING, "compute error");
		return(0);
	}
	multcolor(r->pcol, cval);
	return(0);
}
