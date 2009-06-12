#ifndef lint
static const char RCSid[] = "$Id: bsdf.c,v 2.1 2009/06/12 17:37:37 greg Exp $";
#endif
/*
 * Routines for handling BSDF data
 */

#include "standard.h"
#include "bsdf.h"
#include "paths.h"
#include "ezxml.h"
#include <ctype.h>

#define MAXLATS		46		/* maximum number of latitudes */

/* BSDF angle specification */
typedef struct {
	char	name[64];		/* basis name */
	int	nangles;		/* total number of directions */
	struct {
		float	tmin;			/* starting theta */
		short	nphis;			/* number of phis (0 term) */
	}	lat[MAXLATS+1];		/* latitudes */
} ANGLE_BASIS;

#define	MAXABASES	3		/* limit on defined bases */

static ANGLE_BASIS	abase_list[MAXABASES] = {
	{
		"LBNL/Klems Full", 145,
		{ {-5., 1},
		{5., 8},
		{15., 16},
		{25., 20},
		{35., 24},
		{45., 24},
		{55., 24},
		{65., 16},
		{75., 12},
		{90., 0} }
	}, {
		"LBNL/Klems Half", 73,
		{ {-6.5, 1},
		{6.5, 8},
		{19.5, 12},
		{32.5, 16},
		{46.5, 20},
		{61.5, 12},
		{76.5, 4},
		{90., 0} }
	}, {
		"LBNL/Klems Quarter", 41,
		{ {-9., 1},
		{9., 8},
		{27., 12},
		{46., 12},
		{66., 8},
		{90., 0} }
	}
};

static int	nabases = 3;	/* current number of defined bases */


static int
ab_getvec(		/* get vector for this angle basis index */
	FVECT v,
	int ndx,
	void *p
)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	int	li;
	double	alt, azi, d;
	
	if ((ndx < 0) | (ndx >= ab->nangles))
		return(0);
	for (li = 0; ndx >= ab->lat[li].nphis; li++)
		ndx -= ab->lat[li].nphis;
	alt = PI/180.*0.5*(ab->lat[li].tmin + ab->lat[li+1].tmin);
	azi = 2.*PI*ndx/ab->lat[li].nphis;
	v[2] = d = cos(alt);
	d = sqrt(1. - d*d);	/* sin(alt) */
	v[0] = cos(azi)*d;
	v[1] = sin(azi)*d;
	return(1);
}


static int
ab_getndx(		/* get index corresponding to the given vector */
	FVECT v,
	void *p
)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	int	li, ndx;
	double	alt, azi, d;

	if ((v[2] < -1.0) | (v[2] > 1.0))
		return(-1);
	alt = 180.0/PI*acos(v[2]);
	azi = 180.0/PI*atan2(v[1], v[0]);
	if (azi < 0.0) azi += 360.0;
	for (li = 1; ab->lat[li].tmin <= alt; li++)
		if (!ab->lat[li].nphis)
			return(-1);
	--li;
	ndx = (int)((1./360.)*azi*ab->lat[li].nphis + 0.5);
	if (ndx >= ab->lat[li].nphis) ndx = 0;
	while (li--)
		ndx += ab->lat[li].nphis;
	return(ndx);
}


static double
ab_getohm(		/* get solid angle for this angle basis index */
	int ndx,
	void *p
)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	int	li;
	double	theta, theta1;
	
	if ((ndx < 0) | (ndx >= ab->nangles))
		return(0);
	for (li = 0; ndx >= ab->lat[li].nphis; li++)
		ndx -= ab->lat[li].nphis;
	theta1 = PI/180. * ab->lat[li+1].tmin;
	if (ab->lat[li].nphis == 1) {		/* special case */
		if (ab->lat[li].tmin > FTINY)
			error(USER, "unsupported BSDF coordinate system");
		return(2.*PI*(1. - cos(theta1)));
	}
	theta = PI/180. * ab->lat[li].tmin;
	return(2.*PI*(cos(theta) - cos(theta1))/(double)ab->lat[li].nphis);
}


static int
ab_getvecR(		/* get reverse vector for this angle basis index */
	FVECT v,
	int ndx,
	void *p
)
{
	if (!ab_getvec(v, ndx, p))
		return(0);

	v[0] = -v[0];
	v[1] = -v[1];
	v[2] = -v[2];

	return(1);
}


static int
ab_getndxR(		/* get index corresponding to the reverse vector */
	FVECT v,
	void *p
)
{
	FVECT  v2;
	
	v2[0] = -v[0];
	v2[1] = -v[1];
	v2[2] = -v[2];

	return ab_getndx(v2, p);
}


static void
load_bsdf_data(		/* load BSDF distribution for this wavelength */
	struct BSDF_data *dp,
	ezxml_t wdb
)
{
	char  *cbasis = ezxml_txt(ezxml_child(wdb,"ColumnAngleBasis"));
	char  *rbasis = ezxml_txt(ezxml_child(wdb,"RowAngleBasis"));
	char  *sdata;
	int  i;
	
	if ((cbasis == NULL) | (rbasis == NULL)) {
		error(WARNING, "missing column/row basis for BSDF");
		return;
	}
	/* XXX need to add routines for loading in foreign bases */
	for (i = nabases; i--; )
		if (!strcmp(cbasis, abase_list[i].name)) {
			dp->ninc = abase_list[i].nangles;
			dp->ib_priv = (void *)&abase_list[i];
			dp->ib_vec = ab_getvecR;
			dp->ib_ndx = ab_getndxR;
			dp->ib_ohm = ab_getohm;
			break;
		}
	if (i < 0) {
		sprintf(errmsg, "unsupported ColumnAngleBasis '%s'", cbasis);
		error(WARNING, errmsg);
		return;
	}
	for (i = nabases; i--; )
		if (!strcmp(rbasis, abase_list[i].name)) {
			dp->nout = abase_list[i].nangles;
			dp->ob_priv = (void *)&abase_list[i];
			dp->ob_vec = ab_getvec;
			dp->ob_ndx = ab_getndx;
			dp->ob_ohm = ab_getohm;
			break;
		}
	if (i < 0) {
		sprintf(errmsg, "unsupported RowAngleBasis '%s'", cbasis);
		error(WARNING, errmsg);
		return;
	}
				/* read BSDF data */
	sdata  = ezxml_txt(ezxml_child(wdb,"ScatteringData"));
	if (sdata == NULL) {
		error(WARNING, "missing BSDF ScatteringData");
		return;
	}
	dp->bsdf = (float *)malloc(sizeof(float)*dp->ninc*dp->nout);
	if (dp->bsdf == NULL)
		error(SYSTEM, "out of memory in load_bsdf_data");
	for (i = 0; i < dp->ninc*dp->nout; i++) {
		char  *sdnext = fskip(sdata);
		if (sdnext == NULL) {
			error(WARNING, "bad/missing BSDF ScatteringData");
			free(dp->bsdf); dp->bsdf = NULL;
			return;
		}
		while (*sdnext && isspace(*sdnext))
			sdnext++;
		if (*sdnext == ',') sdnext++;
		dp->bsdf[i] = atof(sdata);
		sdata = sdnext;
	}
	while (isspace(*sdata))
		sdata++;
	if (*sdata) {
		sprintf(errmsg, "%d extra characters after BSDF ScatteringData",
				strlen(sdata));
		error(WARNING, errmsg);
	}
}


static int
check_bsdf_data(	/* check that BSDF data is sane */
	struct BSDF_data *dp
)
{
	double *	omega_arr;
	double		dom, hemi_total;
	int		nneg;
	int		i, o;

	if (dp == NULL || dp->bsdf == NULL)
		return(0);
	omega_arr = (double *)calloc(dp->nout, sizeof(double));
	if (omega_arr == NULL)
		error(SYSTEM, "out of memory in check_bsdf_data");
	hemi_total = .0;
	for (o = dp->nout; o--; ) {
		FVECT	v;
		dom = getBSDF_outohm(dp,o);
		if (dom <= .0) {
			error(WARNING, "zero/negative solid angle");
			continue;
		}
		if (!getBSDF_outvec(v,dp,o) || v[2] < -FTINY) {
			error(WARNING, "illegal outgoing BSDF direction");
			free(omega_arr);
			return(0);
		}
		hemi_total += omega_arr[o] = dom*v[2];
	}
	if ((hemi_total > 1.02*PI) | (hemi_total < 0.98*PI)) {
		sprintf(errmsg, "outgoing BSDF hemisphere off by %.1f%%",
				100.*(hemi_total/PI - 1.));
		error(WARNING, errmsg);
	}
	dom = PI / hemi_total;		/* normalize solid angles */
	for (o = dp->nout; o--; )
		omega_arr[o] *= dom;
	nneg = 0;
	for (i = dp->ninc; i--; ) {
		hemi_total = .0;
		for (o = dp->nout; o--; ) {
			double	f = BSDF_value(dp,i,o);
			if (f > .0)
				hemi_total += f*omega_arr[o];
			else if (f < -FTINY)
				++nneg;
		}
		if (hemi_total > 1.02) {
			sprintf(errmsg, "BSDF direction passes %.1f%% of light",
					100.*hemi_total);
			error(WARNING, errmsg);
		}
	}
	free(omega_arr);
	if (nneg > 0) {
		sprintf(errmsg, "%d negative BSDF values", nneg);
		error(WARNING, errmsg);
		return(0);
	}
	return(1);
}

struct BSDF_data *
load_BSDF(		/* load BSDF data from file */
	char *fname
)
{
	char			*path;
	ezxml_t			fl, wtl, wld, wdb;
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
	if (ezxml_error(fl)[0]) {
		sprintf(errmsg, "BSDF \"%s\" %s", path, ezxml_error(fl));
		error(WARNING, errmsg);
		ezxml_free(fl);
		return(NULL);
	}
	if (strcmp(ezxml_name(fl), "WindowElement")) {
		sprintf(errmsg,
			"BSDF \"%s\": top level node not 'WindowElement'",
				path);
		error(WARNING, errmsg);
		ezxml_free(fl);
		return(NULL);
	}
	wtl = ezxml_child(ezxml_child(fl, "Optical"), "Layer");
	dp = (struct BSDF_data *)calloc(1, sizeof(struct BSDF_data));
	for (wld = ezxml_child(wtl, "WavelengthData");
				wld != NULL; wld = wld->next) {
		if (strcmp(ezxml_txt(ezxml_child(wld,"Wavelength")), "Visible"))
			continue;
		wdb = ezxml_child(wld, "WavelengthDataBlock");
		if (wdb == NULL) continue;
		if (strcmp(ezxml_txt(ezxml_child(wdb,"WavelengthDataDirection")),
					"Transmission Front"))
			continue;
		load_bsdf_data(dp, wdb);	/* load front BTDF */
		break;				/* ignore the rest */
	}
	ezxml_free(fl);				/* done with XML file */
	if (!check_bsdf_data(dp)) {
		sprintf(errmsg, "bad/missing BTDF data in \"%s\"", path);
		error(WARNING, errmsg);
		free_BSDF(dp);
		dp = NULL;
	}
	return(dp);
}


void
free_BSDF(		/* free BSDF data structure */
	struct BSDF_data *b
)
{
	if (b == NULL)
		return;
	if (b->bsdf != NULL)
		free(b->bsdf);
	free(b);
}


int
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
	
	if (!getBSDF_incvec(v, b, i))
		return(0);
	rad = sqrt(getBSDF_incohm(b, i) / PI);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	return(normalize(v) != 0.0);
}


int
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
	
	if (!getBSDF_outvec(v, b, o))
		return(0);
	rad = sqrt(getBSDF_outohm(b, o) / PI);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	return(normalize(v) != 0.0);
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

	if (yp[2]*yp[2] + zp[2]*zp[2] < 2.*FTINY*FTINY) {
		/* Special case for X' along Z-axis */
		theta = -atan2(yp[0], yp[1]);
		*xfp++ = "-ry";
		*xfp++ = xp[2] < 0.0 ? "90" : "-90";
		*xfp++ = "-rz";
		sprintf(bufs[bn], "%f", theta*(180./PI));
		*xfp++ = bufs[bn++];
		return(xfp - xfarg);
	}
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
getBSDF_xfm(		/* compute BSDF orient. -> world orient. transform */
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
