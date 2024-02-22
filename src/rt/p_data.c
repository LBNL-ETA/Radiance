#ifndef lint
static const char	RCSid[] = "$Id: p_data.c,v 2.14 2024/02/22 20:11:54 greg Exp $";
#endif
/*
 *  p_data.c - routine for stored patterns.
 */

#include "copyright.h"

#include  "ray.h"
#include  "data.h"
#include  "func.h"
#include  "rtotypes.h"

/*
 *	A stored pattern can either be brightness,
 *  color, or spectral data.  Brightness data is specified as:
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
 *  A simple spectrum is specified as:
 *
 *	modifier spectrum name
 *	0
 *	0
 *	5+ nmA nmB s1 s2 s3 ..
 *
 *  A constant spectrum from a data file is given as:
 *
 *	modifier specfile name
 *	1 dfname
 *	0
 *	0
 *
 *  A spectral data file is given as:
 *
 *	modifier specdata name
 *	4+ sfunc dfname vfname v0 .. xf
 *	0
 *	n A1 A2 ..
 *
 *  A spectral image is given as:
 *
 *	modifier specpict name
 *	5+ sfunc sfname vfname vx vy xf
 *	0
 *	n A1 A2 ..
 *
 *  Vfname is the name of the file where the variable definitions
 *  can be found.  The list of real arguments can be accessed by
 *  definitions in the file.  The dfnames are the data file
 *  names.  The dimensions of the data files and the number
 *  of variables must match, except for specdata, which has a "hidden"
 *  last variable for the wavelength.  The funcs take a single argument
 *  for brightdata, three for colordata and colorpict, and two for
 *  specdata and specpict to modify interpolated values from the file.
 *  The xf is a transform spec to get from the original coordinates to
 *  the current coordinates.
 */


int
p_bdata(			/* interpolate brightness data */
	OBJREC  *m,
	RAY  *r
)
{
	double  bval;
	double  pt[MAXDDIM];
	DATARRAY  *dp;
	MFUNC  *mf;
	int  i;

	if (m->oargs.nsargs < 4)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[1]);
	i = (1 << dp->nd) - 1;
	mf = getfunc(m, 2, i<<3, 0);
	setfunc(m, r);
	errno = 0;
	for (i = dp->nd; i-- > 0; ) {
		pt[i] = evalue(mf->ep[i]);
		if ((errno == EDOM) | (errno == ERANGE))
			goto computerr;
	}
	bval = datavalue(dp, pt);
	errno = 0;
	bval = funvalue(m->oargs.sarg[0], 1, &bval);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	scalescolor(r->pcol, bval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


int
p_cdata(			/* interpolate color data */
	OBJREC  *m,
	RAY  *r
)
{
	double  col[3];
	COLOR  cval;
	double  pt[MAXDDIM];
	int  nv;
	DATARRAY  *dp;
	MFUNC  *mf;
	int  i;

	if (m->oargs.nsargs < 8)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[3]);
	i = (1 << (nv = dp->nd)) - 1;
	mf = getfunc(m, 6, i<<7, 0);
	setfunc(m, r);
	errno = 0;
	for (i = 0; i < nv; i++) {
		pt[i] = evalue(mf->ep[i]);
		if ((errno == EDOM) | (errno == ERANGE))
			goto computerr;
	}
	col[0] = datavalue(dp, pt);
	for (i = 1; i < 3; i++) {
		if (!strcmp(m->oargs.sarg[3+i], m->oargs.sarg[3])) {
			col[i] = col[0];	/* same data */
			continue;
		}
		dp = getdata(m->oargs.sarg[3+i]);
		if (dp->nd != nv)
			objerror(m, USER, "dimension error");
		col[i] = datavalue(dp, pt);
	}
	errno = 0;
	for (i = 0; i < 3; i++)
		if (i && fundefined(m->oargs.sarg[i]) < 3)
			colval(cval,i) = funvalue(m->oargs.sarg[i], 1, col+i);
		else
			colval(cval,i) = funvalue(m->oargs.sarg[i], 3, col);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	smultcolor(r->pcol, cval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


int
p_pdata(			/* interpolate picture data */
	OBJREC  *m,
	RAY  *r
)
{
	double  col[3];
	COLOR  cval;
	double  pt[2];
	DATARRAY  *dp;
	MFUNC  *mf;
	int  i;

	if (m->oargs.nsargs < 7)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 4, 0x3<<5, 0);
	setfunc(m, r);
	errno = 0;
	pt[1] = evalue(mf->ep[0]);	/* y major ordering */
	pt[0] = evalue(mf->ep[1]);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	dp = getpict(m->oargs.sarg[3]);
	for (i = 0; i < 3; i++)
		col[i] = datavalue(dp+i, pt);
	errno = 0;
	for (i = 0; i < 3; i++)
		if (i && fundefined(m->oargs.sarg[i]) < 3)
			colval(cval,i) = funvalue(m->oargs.sarg[i], 1, col+i);
		else
			colval(cval,i) = funvalue(m->oargs.sarg[i], 3, col);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	smultcolor(r->pcol, cval);
	return(0);

computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


int
p_spectrum(			/* simple constant spectrum */
	OBJREC  *m,
	RAY  *r
)
{
	COLORV	*scval;

	if ((scval = (COLORV *)m->os) == NULL) {
		COLORV	*sinp;
		double	hstep;
		int	i;
		if (m->oargs.nfargs < 5)
			objerror(m, USER, "bad # arguments");
		sinp = (COLORV *)malloc(sizeof(COLORV)*(m->oargs.nfargs-2));
		scval = (COLORV *)malloc(sizeof(COLORV)*NCSAMP);
		if ((sinp == NULL) | (scval == NULL))
			objerror(m, SYSTEM, "out of memory");
		for (i = m->oargs.nfargs-2; i--; )
			sinp[i] = (COLORV)m->oargs.farg[i+2];
		hstep = 0.5 * (m->oargs.farg[1] - m->oargs.farg[0]) /
				(m->oargs.nfargs-3.0);
		convertscolorcol(scval, sinp, m->oargs.nfargs-2,
				m->oargs.farg[0]-hstep, m->oargs.farg[1]+hstep);
		free(sinp);
		m->os = (void *)scval;
	}
	smultscolor(r->pcol, scval);
	return(0);
}


int
p_specfile(			/* constant spectrum from 1-D data file */
	OBJREC  *m,
	RAY  *r
)
{
	COLORV	*scval;

	if ((scval = (COLORV *)m->os) == NULL) {
		DATARRAY	*dp;
		COLORV		*sinp;
		double		step;
		int		i;
		if (m->oargs.nsargs != 1)
			objerror(m, USER, "bad # arguments");
		dp = getdata(m->oargs.sarg[0]);
		if (dp->nd != 1)
			objerror(m, USER, "data file must be 1-dimensional");

		sinp = (COLORV *)malloc(sizeof(COLORV)*dp->dim[0].ne);
		scval = (COLORV *)malloc(sizeof(COLORV)*NCSAMP);
		if ((sinp == NULL) | (scval == NULL))
			objerror(m, SYSTEM, "out of memory");
		step = dp->dim[0].siz / (dp->dim[0].ne - 1.0);
		for (i = dp->dim[0].ne; i-- > 0; ) {
			double	wl = dp->dim[0].org + i*step;
			sinp[i] = (COLORV)datavalue(dp, &wl);
		}
		convertscolorcol(scval, sinp, dp->dim[0].ne,
				dp->dim[0].org-.5*step,
				dp->dim[0].org+dp->dim[0].siz+.5*step);
		free(sinp);
		m->os = (void *)scval;
	}
	smultscolor(r->pcol, scval);
	return(0);
}


int
p_specdata(			/* varied spectrum from (N+1)-D file */
	OBJREC  *m,
	RAY  *r
)
{
	SCOLOR		scval;
	COLORV		*scdat;
	double		pt[MAXDDIM];
	DATARRAY	*dp;
	MFUNC		*mf;
	double		step;
	int		i;

	if (m->oargs.nsargs < 4)
		objerror(m, USER, "bad # arguments");
	dp = getdata(m->oargs.sarg[1]);
	if (dp->nd < 2)
		objerror(m, USER, "need at least 2-dimensional data");
	i = (1 << (dp->nd-1)) - 1;
	mf = getfunc(m, 2, i<<3, 0);
	setfunc(m, r);
	errno = 0;
	for (i = dp->nd-1; i-- > 0; ) {
		pt[i] = evalue(mf->ep[i]);
		if ((errno == EDOM) | (errno == ERANGE))
			goto computerr;
	}
	step = dp->dim[dp->nd-1].siz / (dp->dim[dp->nd-1].ne - 1.0);
	scdat = (COLORV *)malloc(sizeof(COLORV)*dp->dim[dp->nd-1].ne);
	if (scdat == NULL)
		objerror(m, SYSTEM, "out of memory");
	for (i = dp->dim[dp->nd-1].ne; i-- > 0; ) {
		double	bval[2];
		pt[dp->nd-1] = dp->dim[dp->nd-1].org + i*step;
		bval[0] = datavalue(dp, pt);
		bval[1] = pt[dp->nd-1];
		errno = 0;
		scdat[i] = funvalue(m->oargs.sarg[0], 2, bval);
		if ((errno == EDOM) | (errno == ERANGE))
			goto computerr;
	}
	convertscolorcol(scval, scdat, dp->dim[dp->nd-1].ne,
			dp->dim[dp->nd-1].org-.5*step,
			dp->dim[dp->nd-1].org+dp->dim[dp->nd-1].siz+.5*step);
	free(scdat);
	smultscolor(r->pcol, scval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}


int
p_specpict(			/* interpolate hyperspectral image data */
	OBJREC  *m,
	RAY  *r
)
{
	SCOLOR		scdat, scval;
	double		pt[3];
	DATARRAY	*dp;
	MFUNC		*mf;
	double		step;
	int 		i;

	if (m->oargs.nsargs < 5)
		objerror(m, USER, "bad # arguments");
	mf = getfunc(m, 2, 0x3<<3, 0);
	setfunc(m, r);
	errno = 0;
	pt[1] = evalue(mf->ep[0]);	/* y major ordering */
	pt[0] = evalue(mf->ep[1]);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	dp = getspec(m->oargs.sarg[1]);
	step = dp->dim[2].siz / (dp->dim[2].ne - 1.0);
	for (i = dp->dim[2].ne; i-- > 0; ) {
		double	bval[2];
		pt[2] = dp->dim[2].org + i*step;
		bval[0] = datavalue(dp, pt);
		bval[1] = pt[2];
		errno = 0;
		scdat[i] = funvalue(m->oargs.sarg[0], 2, bval);
		if ((errno == EDOM) | (errno == ERANGE))
			goto computerr;
	}
	convertscolorcol(scval, scdat, dp->dim[2].ne,
			dp->dim[2].org-.5*step,
			dp->dim[2].org+dp->dim[2].siz+.5*step);
	smultscolor(r->pcol, scval);
	return(0);
computerr:
	objerror(m, WARNING, "compute error");
	return(0);
}
