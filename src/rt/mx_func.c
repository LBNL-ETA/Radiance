#ifndef lint
static const char	RCSid[] = "$Id: mx_func.c,v 2.8 2004/03/30 16:13:01 schorsch Exp $";
#endif
/*
 *  mx_func.c - routine for mixture functions.
 */

#include "copyright.h"

#include  "ray.h"
#include  "func.h"
#include  "rtotypes.h"

/*
 *	A mixture function is specified:
 *
 *	modifier mixfunc name
 *	4+ foremod backmod varname vfname xf
 *	0
 *	n A1 A2 ..
 *
 *  Vfname is the name of the file where the variable definition
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 */


extern int
mx_func(			/* compute mixture function */
	register OBJREC  *m,
	RAY  *r
)
{
	OBJECT	obj;
	register int  i;
	double  coef;
	OBJECT  mod[2];
	register MFUNC  *mf;

	if (m->oargs.nsargs < 4)
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
	mf = getfunc(m, 3, 0x4, 0);
	setfunc(m, r);
	errno = 0;
	coef = evalue(mf->ep[0]);
	if (errno == EDOM || errno == ERANGE) {
		objerror(m, WARNING, "compute error");
		return(0);
	}
	if (raymixture(r, mod[0], mod[1], coef)) {
		if (m->omod != OVOID)
			objerror(m, USER, "inappropriate modifier");
		return(1);
	}
	return(0);
}
