#ifndef lint
static const char RCSid[] = "$Id: mkillum4.c,v 2.2 2007/09/21 05:53:21 greg Exp $";
#endif
/*
 * Routines for handling BSDF data within mkillum
 */

#include "mkillum.h"
#include "paths.h"
#include "ezxml.h"


struct BSDF_data *
load_BSDF(		/* load BSDF data from file */
	char *fname
)
{
	char			*path;
	ezxml_t			fl, wld, wdb;
	struct BSDF_data	*dp;
	
	path = getpath(fname, getrlibpath(), R_OK);
	if (path == NULL) {
		sprintf(errmsg, "cannot find BSDF file \"%s\"", fname);
		error(WARNING, errmsg);
		return(NULL);
	}
	fl = ezxml_parse_file(path);
	if (fl == NULL) {
		sprintf(errmsg, "cannot open BSDF \"%s\"", path);
		error(WARNING, errmsg);
		return(NULL);
	}
	dp = (struct BSDF_data *)malloc(sizeof(struct BSDF_data));
	if (dp == NULL)
		goto memerr;
	for (wld = ezxml_child(fl, "WavelengthData");
				fl != NULL; fl = fl->next) {
		if (strcmp(ezxml_child(wld, "Wavelength")->txt, "Visible"))
			continue;
		wdb = ezxml_child(wld, "WavelengthDataBlock");
		if (wdb == NULL) continue;
		if (strcmp(ezxml_child(wdb, "WavelengthDataDirection")->txt,
					"Transmission Front"))
			continue;
	}
	/* etc... */
	ezxml_free(fl);
	return(dp);
memerr:
	error(SYSTEM, "out of memory in load_BSDF");
	return NULL;	/* pro forma return */
}


void
free_BSDF(		/* free BSDF data structure */
	struct BSDF_data *b
)
{
	if (b == NULL)
		return;
	free(b->inc_dir);
	free(b->inc_rad);
	free(b->out_dir);
	free(b->out_rad);
	free(b->bsdf);
	free(b);
}


void
r_BSDF_incvec(		/* compute random input vector at given location */
	FVECT v,
	struct BSDF_data *b,
	int i,
	double rv,
	MAT4 xm
)
{
	FVECT	pert;
	double	rad;
	int	j;
	
	getBSDF_incvec(v, b, i);
	rad = getBSDF_incrad(b, i);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	normalize(v);
}


void
r_BSDF_outvec(		/* compute random output vector at given location */
	FVECT v,
	struct BSDF_data *b,
	int o,
	double rv,
	MAT4 xm
)
{
	FVECT	pert;
	double	rad;
	int	j;
	
	getBSDF_outvec(v, b, o);
	rad = getBSDF_outrad(b, o);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	normalize(v);
}


#define  FEQ(a,b)	((a)-(b) <= 1e-7 && (b)-(a) <= 1e-7)

static int
addrot(			/* compute rotation (x,y,z) => (xp,yp,zp) */
	char *xfarg[],
	FVECT xp,
	FVECT yp,
	FVECT zp
)
{
	static char	bufs[3][16];
	int	bn = 0;
	char	**xfp = xfarg;
	double	theta;

	theta = atan2(yp[2], zp[2]);
	if (!FEQ(theta,0.0)) {
		*xfp++ = "-rx";
		sprintf(bufs[bn], "%f", theta*(180./PI));
		*xfp++ = bufs[bn++];
	}
	theta = asin(-xp[2]);
	if (!FEQ(theta,0.0)) {
		*xfp++ = "-ry";
		sprintf(bufs[bn], " %f", theta*(180./PI));
		*xfp++ = bufs[bn++];
	}
	theta = atan2(xp[1], xp[0]);
	if (!FEQ(theta,0.0)) {
		*xfp++ = "-rz";
		sprintf(bufs[bn], "%f", theta*(180./PI));
		*xfp++ = bufs[bn++];
	}
	*xfp = NULL;
	return(xfp - xfarg);
}


int
getBSDF_xfm(		/* compute transform for the given surface */
	MAT4 xm,
	FVECT nrm,
	UpDir ud
)
{
	char	*xfargs[7];
	XF	myxf;
	FVECT	updir, xdest, ydest;

	updir[0] = updir[1] = updir[2] = 0.;
	switch (ud) {
	case UDzneg:
		updir[2] = -1.;
		break;
	case UDyneg:
		updir[1] = -1.;
		break;
	case UDxneg:
		updir[0] = -1.;
		break;
	case UDxpos:
		updir[0] = 1.;
		break;
	case UDypos:
		updir[1] = 1.;
		break;
	case UDzpos:
		updir[2] = 1.;
		break;
	case UDunknown:
		error(WARNING, "unspecified up direction");
		return(0);
	}
	fcross(xdest, updir, nrm);
	if (normalize(xdest) == 0.0)
		return(0);
	fcross(ydest, nrm, xdest);
	xf(&myxf, addrot(xfargs, xdest, ydest, nrm), xfargs);
	copymat4(xm, myxf.xfm);
	return(1);
}


void
redistribute(		/* pass distarr ray sums through BSDF */
	struct BSDF_data *b,
	int nalt,
	int nazi,
	FVECT u,
	FVECT v,
	FVECT w,
	MAT4 xm
)
{
	MAT4	mymat;
	COLORV	*outarr;
	float	*inpcoef;
	COLORV	*cp, *csum;
	uint16	*distcnt;
	FVECT	dv;
	double	oom, wt;
	int	i, j, o;
	int	cnt;
	COLOR	col;
					/* allocate temporary memory */
	outarr = (COLORV *)malloc(b->nout * sizeof(COLOR));
	distcnt = (uint16 *)calloc(nalt*nazi, sizeof(uint16));
	inpcoef = (float *)malloc(b->ninc * sizeof(float));
	if ((outarr == NULL) | (distcnt == NULL) | (inpcoef == NULL))
		error(SYSTEM, "out of memory in redistribute");
					/* compose matrix */
	for (i = 3; i--; ) {
		mymat[i][0] = u[i];
		mymat[i][1] = v[i];
		mymat[i][2] = w[i];
		mymat[i][3] = 0.;
	}
	mymat[3][0] = mymat[3][1] = mymat[3][2] = 0.;
	mymat[3][3] = 1.;
	if (xm != NULL)
		multmat4(mymat, xm, mymat);
	for (i = 3; i--; ) {		/* make sure it's normalized */
		wt = 1./sqrt(	mymat[0][i]*mymat[0][i] +
				mymat[1][i]*mymat[1][i] +
				mymat[2][i]*mymat[2][i]	);
		for (j = 3; j--; )
			mymat[j][i] *= wt;
	}
					/* pass through BSDF */
	for (i = b->ninc; i--; ) {		/* get input coefficients */
		getBSDF_incvec(dv, b, i);
		multv3(dv, dv, mymat);
		wt = getBSDF_incrad(b, i);
		inpcoef[i] = PI*wt*wt * dv[2];	/* solid_angle*cosine(theta) */
	}
	for (o = b->nout; o--; ) {
		csum = &outarr[3*o];
		setcolor(csum, 0., 0., 0.);
		oom = getBSDF_outrad(b, o);
		oom *= oom * PI;
		for (i = b->ninc; i--; ) {
			wt = BSDF_data(b,i,o) * inpcoef[i] / oom;
			cp = &distarr[3*i];
			copycolor(col, cp);
			scalecolor(col, wt);
			addcolor(csum, col);
		}
		wt = 1./b->ninc;
		scalecolor(csum, wt);
	}
	free(inpcoef);
	newdist(nalt*nazi);		/* resample distribution */
	for (o = b->nout; o--; ) {
		getBSDF_outvec(dv, b, o);
		multv3(dv, dv, mymat);
		j = (.5 + atan2(dv[1],dv[0])*(.5/PI))*nazi + .5;
		if (j >= nazi) j = 0;
		i = (0.9999 - dv[2]*dv[2])*nalt;
		csum = &distarr[3*(i*nazi + j)];
		cp = &outarr[3*o];
		addcolor(csum, cp);
		++distcnt[i*nazi + j];
	}
	free(outarr);
					/* fill in missing bits */
	for (i = nalt; i--; )
	    for (j = nazi; j--; ) {
		int	ii, jj, alt, azi;
		if (distcnt[i*nazi + j])
			continue;
		csum = &distarr[3*(i*nazi + j)];
		setcolor(csum, 0., 0., 0.);
		cnt = 0;
		for (o = 0; !cnt; o++)
		    for (ii = -o; ii <= o; ii++) {
			alt = i + ii;
			if (alt < 0) continue;
			if (alt >= nalt) break;
			for (jj = -o; jj <= o; jj++) {
			    if (ii*ii + jj*jj != o*o)
				continue;
			    azi = j + jj;
			    if (azi >= nazi) azi -= nazi;
			    else if (azi < 0) azi += nazi;
			    if (!distcnt[alt*nazi + azi])
				continue;
			    cp = &distarr[3*(alt*nazi + azi)];
			    addcolor(csum, cp);
			    cnt += distcnt[alt*nazi + azi];
			}
		    }
		wt = 1./cnt;
		scalecolor(csum, wt);
	    }
					/* finish averages */
	for (i = nalt; i--; )
	    for (j = nazi; j--; ) {
		if ((cnt = distcnt[i*nazi + j]) <= 1)
			continue;
		csum = &distarr[3*(i*nazi + j)];
		wt = 1./cnt;
		scalecolor(csum, wt);
	    }
	free(distcnt);
}
