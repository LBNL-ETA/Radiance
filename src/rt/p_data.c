#ifndef lint
static const char	RCSid[] = "$Id: p_data.c,v 2.7 2003/03/05 16:16:53 greg Exp $";
#endif
/*
 *  p_data.c - routine for stored patterns.
 */

#include "copyright.h"

#include  "ray.h"

#include  "data.h"

#include  "func.h"

/*
 *	A stored pattern can either be brightness or
 *  color data.  Brightness data is specified as:
 *
 *	modifier brightdata name
 *	4+ func dfname vfname v0 v1 .. xf
 *	0
 *	n A1 A2 ..
 *
 *  Color data is specified as:
 *
 *	modifier colordata name
 *	8+ rfunc gfunc bfunc rdfname gdfname bdfname vfname v0 v1 .. xf
 *	0
 *	n A1 A2 ..
 *
 *  Color picture data is specified as:
 *
 *	modifier colorpict name
 *	7+ rfunc gfunc bfunc pfname vfname vx vy xf
 *	0
 *	n A1 A2 ..
 *
 *  Vfname is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The dfnames are the data file
 *  names.  The dimensions of the data files and the number
 *  of variables must match.  The funcs take a single argument
 *  for brightdata, and three for colordata and colorpict to produce
 *  interpolated values from the file.  The xf is a transformation
 *  to get from the original coordinates to the current coordinates.
 */


p_bdata(m, r)			/* interpolate brightness data */
register OBJREC  *m;
RAY  *r;
{
	double  bval;
	double  pt[MAXDIM];
	DATARRAY  *dp;
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 4)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[1]);
	i = (1 << dp->nd) - 1;
	mf = getfunc(m, 2, i<<3, 0);
	setfunc(m, r);
	errno = 0;
	for (i = dp->nd; i-- > 0; ) {
		pt[i] = evalue(mf->ep[i]);
		if (errno == EDOM || errno == ERANGE)
			goto computerr;
	}
	bval = datavalue(dp, pt);
	errno = 0;
	bval = funvalue(m->oargs.sarg[0], 1, &bval);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	scalecolor(r->pcol, bval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


p_cdata(m, r)			/* interpolate color data */
register OBJREC  *m;
RAY  *r;
{
	double  col[3];
	COLOR  cval;
	double  pt[MAXDIM];
	int  nv;
	DATARRAY  *dp;
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 8)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[3]);
	i = (1 << (nv = dp->nd)) - 1;
	mf = getfunc(m, 6, i<<7, 0);
	setfunc(m, r);
	errno = 0;
	for (i = 0; i < nv; i++) {
		pt[i] = evalue(mf->ep[i]);
		if (errno == EDOM || errno == ERANGE)
			goto computerr;
	}
	col[0] = datavalue(dp, pt);
	for (i = 1; i < 3; i++) {
		dp = getdata(m->oargs.sarg[i+3]);
		if (dp->nd != nv)
			objerror(m, USER, "dimension error");
		col[i] = datavalue(dp, pt);
	}
	errno = 0;
	for (i = 0; i < 3; i++)
		if (fundefined(m->oargs.sarg[i]) < 3)
			colval(cval,i) = funvalue(m->oargs.sarg[i], 1, col+i);
		else
			colval(cval,i) = funvalue(m->oargs.sarg[i], 3, col);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	multcolor(r->pcol, cval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


p_pdata(m, r)			/* interpolate picture data */
register OBJREC  *m;
RAY  *r;
{
	double  col[3];
	COLOR  cval;
	double  pt[2];
	DATARRAY  *dp;
	register MFUNC  *mf;
	register int  i;

	if (m->oargs.nsargs < 7)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 4, 0x3<<5, 0);
	setfunc(m, r);
	errno = 0;
	pt[1] = evalue(mf->ep[0]);	/* y major ordering */
	pt[0] = evalue(mf->ep[1]);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	dp = getpict(m->oargs.sarg[3]);
	for (i = 0; i < 3; i++)
		col[i] = datavalue(dp+i, pt);
	errno = 0;
	for (i = 0; i < 3; i++)
		if (fundefined(m->oargs.sarg[i]) < 3)
			colval(cval,i) = funvalue(m->oargs.sarg[i], 1, col+i);
		else
			colval(cval,i) = funvalue(m->oargs.sarg[i], 3, col);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	multcolor(r->pcol, cval);
	return(0);

computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}
