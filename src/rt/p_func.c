#ifndef lint
static const char	RCSid[] = "$Id$";
#endif
/*
 *  p_func.c - routine for procedural patterns.
 */

#include "copyright.h"

#include  "ray.h"
#include  "func.h"
#include  "random.h"
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
 *  A spectral function is given as:
 *
 *	modifier specfunc name
 *	2+ sfunc filename xf
 *	0
 *	2+ nmA nmB A3 ..
 *
 *  Filename is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 *  For the "specfunc" primitive, sfunc(nm) is a function of wavelength
 *  and must be defined from nmA to nmB, and should average to 1 over
 *  its range.
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
	if ((errno == EDOM) | (errno == ERANGE)) {
		objerror(m, WARNING, "compute error");
		return(0);
	}
	scalescolor(r->pcol, bval);
	return(0);
}


int
p_cfunc(			/* compute color pattern */
	OBJREC  *m,
	RAY  *r
)
{
	SCOLOR  scval;
	MFUNC  *mf;

	if (m->oargs.nsargs < 4)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 3, 0x7, 0);
	setfunc(m, r);
	errno = 0;
	setscolor(scval, evalue(mf->ep[0]),
			evalue(mf->ep[1]),
			evalue(mf->ep[2]));
	if ((errno == EDOM) | (errno == ERANGE)) {
		objerror(m, WARNING, "compute error");
		return(0);
	}
	smultscolor(r->pcol, scval);
	return(0);
}


int
p_specfunc(			/* compute spectral pattern */
	OBJREC  *m,
	RAY  *r
)
{
	SCOLOR	scsamp;
	SCOLOR  scval;
	double	wl, wlmin, wlmax, wlstep;
	int	ns, i;

	if ((m->oargs.nsargs < 2) | (m->oargs.nfargs < 2))
		objerror(m, USER, "bad # arguments");
	if (m->oargs.farg[0] < m->oargs.farg[1]) {
		wlmin = m->oargs.farg[0];
		wlmax = m->oargs.farg[1];
	} else {
		wlmin = m->oargs.farg[1];
		wlmax = m->oargs.farg[0];
	}
	if (wlmin < WLPART[3]) wlmin = WLPART[3];
	if (wlmax > WLPART[0]) wlmax = WLPART[0];
	if (wlmin >= wlmax) {
		objerror(m, WARNING, "incompatible wavelength sampling");
		return(0);
	}
	wlstep = (wlmax - wlmin)/(double)MAXCSAMP;
	wl = 0.25/(double)NCSAMP*(WLPART[0] - WLPART[3]);
	if (wlstep < wl)		/* too fine to matter? */
		wlstep = wl;
	getfunc(m, 1, 0, 0);
	setfunc(m, r);
	errno = 0;
	ns = (wlmax - wlmin)/wlstep + .1;
	wl = wlmax - .5*wlstep;
	for (i = ns; i-- > 0; wl -= wlstep) {
		double	ws = wl + 0.9*(.5-frandom())*wlstep;
		scsamp[i] = funvalue(m->oargs.sarg[0], 1, &ws);
		if ((errno == EDOM) | (errno == ERANGE)) {
			objerror(m, WARNING, "compute error");
			return(0);
		}
	}
	convertscolorcol(scval, scsamp, ns, wlmin, wlmax);
	smultscolor(r->pcol, scval);
	return(0);
}
