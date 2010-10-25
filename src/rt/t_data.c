#ifndef lint
static const char	RCSid[] = "$Id: t_data.c,v 2.9 2010/10/25 22:57:45 greg Exp $";
#endif
/*
 *  t_data.c - routine for stored textures
 */

#include "copyright.h"

#include  "ray.h"
#include  "data.h"
#include  "func.h"
#include  "rtotypes.h"

/*
 *	A stored texture is specified as follows:
 *
 *	modifier texdata name
 *	8+ xfunc yfunc zfunc xdfname ydfname zdfname vfname v0 v1 .. xf
 *	0
 *	n A1 A2 .. An
 *
 *  Vfname is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The dfnames are the data file
 *  names.  The dimensions of the data files and the number
 *  of variables must match.  The funcs take three arguments to produce
 *  interpolated values from the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 */


extern int
t_data(			/* interpolate texture data */
	register OBJREC  *m,
	RAY  *r
)
{
	int  nv;
	FVECT  disp;
	double  dval[3], pt[MAXDIM];
	double  d;
	DATARRAY  *dp;
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 8)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[3]);
	i = (1 << (nv = dp->nd)) - 1;
	mf = getfunc(m, 6, i<<7, 1);
	setfunc(m, r);
	errno = 0;
	for (i = 0; i < nv; i++)
		pt[i] = evalue(mf->ep[i]);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	dval[0] = datavalue(dp, pt);
	for (i = 1; i < 3; i++) {
		dp = getdata(m->oargs.sarg[i+3]);
		if (dp->nd != nv)
			objerror(m, USER, "dimension error");
		dval[i] = datavalue(dp, pt);
	}
	errno = 0;
	for (i = 0; i < 3; i++)
		disp[i] = funvalue(m->oargs.sarg[i], 3, dval);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	if (mf->f != &unitxf)
		multv3(disp, disp, mf->f->xfm);
	if (r->rox != NULL) {
		multv3(disp, disp, r->rox->f.xfm);
		d = 1.0 / (mf->f->sca * r->rox->f.sca);
	} else
		d = 1.0 / mf->f->sca;
	VSUM(r->pert, r->pert, disp, d);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}
