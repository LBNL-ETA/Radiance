#ifndef lint
static const char RCSid[] = "$Id: bsdf.c,v 2.12 2011/01/06 04:40:22 greg Exp $";
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

#define	MAXABASES	7		/* limit on defined bases */

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

#define  FEQ(a,b)	((a)-(b) <= 1e-6 && (b)-(a) <= 1e-6)

static int
fequal(double a, double b)
{
	if (b != .0)
		a = a/b - 1.;
	return((a <= 1e-6) & (a >= -1e-6));
}

// returns the name of the given tag
#ifdef ezxml_name
#undef ezxml_name
static char *
ezxml_name(ezxml_t xml)
{
	if (xml == NULL)
		return(NULL);
	return(xml->name);
}
#endif

// returns the given tag's character content or empty string if none
#ifdef ezxml_txt
#undef ezxml_txt
static char *
ezxml_txt(ezxml_t xml)
{
	if (xml == NULL)
		return("");
	return(xml->txt);
}
#endif


static int
ab_getvec(		/* get vector for this angle basis index */
	FVECT v,
	int ndx,
	void *p
)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	int	li;
	double	pol, azi, d;
	
	if ((ndx < 0) | (ndx >= ab->nangles))
		return(0);
	for (li = 0; ndx >= ab->lat[li].nphis; li++)
		ndx -= ab->lat[li].nphis;
	pol = PI/180.*0.5*(ab->lat[li].tmin + ab->lat[li+1].tmin);
	azi = 2.*PI*ndx/ab->lat[li].nphis;
	v[2] = d = cos(pol);
	d = sqrt(1. - d*d);	/* sin(pol) */
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
	double	pol, azi, d;

	if ((v[2] < -1.0) | (v[2] > 1.0))
		return(-1);
	pol = 180.0/PI*acos(v[2]);
	azi = 180.0/PI*atan2(v[1], v[0]);
	if (azi < 0.0) azi += 360.0;
	for (li = 1; ab->lat[li].tmin <= pol; li++)
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
load_angle_basis(	/* load custom BSDF angle basis */
	ezxml_t wab
)
{
	char	*abname = ezxml_txt(ezxml_child(wab, "AngleBasisName"));
	ezxml_t	wbb;
	int	i;
	
	if (!abname || !*abname)
		return;
	for (i = nabases; i--; )
		if (!strcasecmp(abname, abase_list[i].name))
			return;		/* assume it's the same */
	if (nabases >= MAXABASES)
		error(INTERNAL, "too many angle bases");
	strcpy(abase_list[nabases].name, abname);
	abase_list[nabases].nangles = 0;
	for (i = 0, wbb = ezxml_child(wab, "AngleBasisBlock");
			wbb != NULL; i++, wbb = wbb->next) {
		if (i >= MAXLATS)
			error(INTERNAL, "too many latitudes in custom basis");
		abase_list[nabases].lat[i+1].tmin = atof(ezxml_txt(
				ezxml_child(ezxml_child(wbb,
					"ThetaBounds"), "UpperTheta")));
		if (!i)
			abase_list[nabases].lat[i].tmin =
					-abase_list[nabases].lat[i+1].tmin;
		else if (!fequal(atof(ezxml_txt(ezxml_child(ezxml_child(wbb,
					"ThetaBounds"), "LowerTheta"))),
				abase_list[nabases].lat[i].tmin))
			error(WARNING, "theta values disagree in custom basis");
		abase_list[nabases].nangles +=
			abase_list[nabases].lat[i].nphis =
				atoi(ezxml_txt(ezxml_child(wbb, "nPhis")));
	}
	abase_list[nabases++].lat[i].nphis = 0;
}


static double
to_meters(		/* return factor to convert given unit to meters */
	const char *unit
)
{
	if (unit == NULL) return(1.);		/* safe assumption? */
	if (!strcasecmp(unit, "Meter")) return(1.);
	if (!strcasecmp(unit, "Foot")) return(.3048);
	if (!strcasecmp(unit, "Inch")) return(.0254);
	if (!strcasecmp(unit, "Centimeter")) return(.01);
	if (!strcasecmp(unit, "Millimeter")) return(.001);
	sprintf(errmsg, "unknown dimensional unit '%s'", unit);
	error(USER, errmsg);
}


static void
load_geometry(		/* load geometric dimensions and description (if any) */
	struct BSDF_data *dp,
	ezxml_t wdb
)
{
	ezxml_t		geom;
	double		cfact;
	const char	*fmt, *mgfstr;

	dp->dim[0] = dp->dim[1] = dp->dim[2] = 0;
	dp->mgf = NULL;
	if ((geom = ezxml_child(wdb, "Width")) != NULL)
		dp->dim[0] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((geom = ezxml_child(wdb, "Height")) != NULL)
		dp->dim[1] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((geom = ezxml_child(wdb, "Thickness")) != NULL)
		dp->dim[2] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((geom = ezxml_child(wdb, "Geometry")) == NULL ||
			(mgfstr = ezxml_txt(geom)) == NULL)
		return;
	if ((fmt = ezxml_attr(geom, "format")) != NULL &&
			strcasecmp(fmt, "MGF")) {
		sprintf(errmsg, "unrecognized geometry format '%s'", fmt);
		error(WARNING, errmsg);
		return;
	}
	cfact = to_meters(ezxml_attr(geom, "unit"));
	dp->mgf = (char *)malloc(strlen(mgfstr)+32);
	if (dp->mgf == NULL)
		error(SYSTEM, "out of memory in load_geometry");
	if (cfact < 0.99 || cfact > 1.01)
		sprintf(dp->mgf, "xf -s %.5f\n%s\nxf\n", cfact, mgfstr);
	else
		strcpy(dp->mgf, mgfstr);
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
	
	if ((!cbasis || !*cbasis) | (!rbasis || !*rbasis)) {
		error(WARNING, "missing column/row basis for BSDF");
		return;
	}
	for (i = nabases; i--; )
		if (!strcasecmp(cbasis, abase_list[i].name)) {
			dp->ninc = abase_list[i].nangles;
			dp->ib_priv = (void *)&abase_list[i];
			dp->ib_vec = ab_getvecR;
			dp->ib_ndx = ab_getndxR;
			dp->ib_ohm = ab_getohm;
			break;
		}
	if (i < 0) {
		sprintf(errmsg, "undefined ColumnAngleBasis '%s'", cbasis);
		error(WARNING, errmsg);
		return;
	}
	for (i = nabases; i--; )
		if (!strcasecmp(rbasis, abase_list[i].name)) {
			dp->nout = abase_list[i].nangles;
			dp->ob_priv = (void *)&abase_list[i];
			dp->ob_vec = ab_getvec;
			dp->ob_ndx = ab_getndx;
			dp->ob_ohm = ab_getohm;
			break;
		}
	if (i < 0) {
		sprintf(errmsg, "undefined RowAngleBasis '%s'", cbasis);
		error(WARNING, errmsg);
		return;
	}
				/* read BSDF data */
	sdata  = ezxml_txt(ezxml_child(wdb,"ScatteringData"));
	if (!sdata || !*sdata) {
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
				(int)strlen(sdata));
		error(WARNING, errmsg);
	}
}


static int
check_bsdf_data(	/* check that BSDF data is sane */
	struct BSDF_data *dp
)
{
	double		*omega_iarr, *omega_oarr;
	double		dom, contrib, hemi_total, full_total;
	int		nneg;
	FVECT		v;
	int		i, o;

	if (dp == NULL || dp->bsdf == NULL)
		return(0);
	omega_iarr = (double *)calloc(dp->ninc, sizeof(double));
	omega_oarr = (double *)calloc(dp->nout, sizeof(double));
	if ((omega_iarr == NULL) | (omega_oarr == NULL))
		error(SYSTEM, "out of memory in check_bsdf_data");
					/* incoming projected solid angles */
	hemi_total = .0;
	for (i = dp->ninc; i--; ) {
		dom = getBSDF_incohm(dp,i);
		if (dom <= .0) {
			error(WARNING, "zero/negative incoming solid angle");
			continue;
		}
		if (!getBSDF_incvec(v,dp,i) || v[2] > FTINY) {
			error(WARNING, "illegal incoming BSDF direction");
			free(omega_iarr); free(omega_oarr);
			return(0);
		}
		hemi_total += omega_iarr[i] = dom * -v[2];
	}
	if ((hemi_total > 1.02*PI) | (hemi_total < 0.98*PI)) {
		sprintf(errmsg, "incoming BSDF hemisphere off by %.1f%%",
				100.*(hemi_total/PI - 1.));
		error(WARNING, errmsg);
	}
	dom = PI / hemi_total;		/* fix normalization */
	for (i = dp->ninc; i--; )
		omega_iarr[i] *= dom;
					/* outgoing projected solid angles */
	hemi_total = .0;
	for (o = dp->nout; o--; ) {
		dom = getBSDF_outohm(dp,o);
		if (dom <= .0) {
			error(WARNING, "zero/negative outgoing solid angle");
			continue;
		}
		if (!getBSDF_outvec(v,dp,o) || v[2] < -FTINY) {
			error(WARNING, "illegal outgoing BSDF direction");
			free(omega_iarr); free(omega_oarr);
			return(0);
		}
		hemi_total += omega_oarr[o] = dom * v[2];
	}
	if ((hemi_total > 1.02*PI) | (hemi_total < 0.98*PI)) {
		sprintf(errmsg, "outgoing BSDF hemisphere off by %.1f%%",
				100.*(hemi_total/PI - 1.));
		error(WARNING, errmsg);
	}
	dom = PI / hemi_total;		/* fix normalization */
	for (o = dp->nout; o--; )
		omega_oarr[o] *= dom;
	nneg = 0;			/* check outgoing totals */
	for (i = 0; i < dp->ninc; i++) {
		hemi_total = .0;
		for (o = dp->nout; o--; ) {
			double	f = BSDF_value(dp,i,o);
			if (f >= .0)
				hemi_total += f*omega_oarr[o];
			else {
				nneg += (f < -FTINY);
				BSDF_value(dp,i,o) = .0f;
			}
		}
		if (hemi_total > 1.01) {
			sprintf(errmsg,
			"incoming BSDF direction %d passes %.1f%% of light",
					i, 100.*hemi_total);
			error(WARNING, errmsg);
		}
	}
	if (nneg) {
		sprintf(errmsg, "%d negative BSDF values (ignored)", nneg);
		error(WARNING, errmsg);
	}
	full_total = .0;		/* reverse roles and check again */
	for (o = 0; o < dp->nout; o++) {
		hemi_total = .0;
		for (i = dp->ninc; i--; )
			hemi_total += BSDF_value(dp,i,o) * omega_iarr[i];

		if (hemi_total > 1.01) {
			sprintf(errmsg,
			"outgoing BSDF direction %d collects %.1f%% of light",
					o, 100.*hemi_total);
			error(WARNING, errmsg);
		}
		full_total += hemi_total*omega_oarr[o];
	}
	full_total /= PI;
	if (full_total > 1.00001) {
		sprintf(errmsg, "BSDF transfers %.4f%% of light",
				100.*full_total);
		error(WARNING, errmsg);
	}
	free(omega_iarr); free(omega_oarr);
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
	if (strcasecmp(ezxml_txt(ezxml_child(ezxml_child(wtl,
			"DataDefinition"), "IncidentDataStructure")),
			"Columns")) {
		sprintf(errmsg,
			"BSDF \"%s\": unsupported IncidentDataStructure",
				path);
		error(WARNING, errmsg);
		ezxml_free(fl);
		return(NULL);
	}		
	load_angle_basis(ezxml_child(ezxml_child(wtl,
				"DataDefinition"), "AngleBasis"));
	dp = (struct BSDF_data *)calloc(1, sizeof(struct BSDF_data));
	load_geometry(dp, ezxml_child(wtl, "Material"));
	for (wld = ezxml_child(wtl, "WavelengthData");
				wld != NULL; wld = wld->next) {
		if (strcasecmp(ezxml_txt(ezxml_child(wld,"Wavelength")), "Visible"))
			continue;
		wdb = ezxml_child(wld, "WavelengthDataBlock");
		if (wdb == NULL) continue;
		if (strcasecmp(ezxml_txt(ezxml_child(wdb,"WavelengthDataDirection")),
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
	if (b->mgf != NULL)
		free(b->mgf);
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
	UpDir ud,
	char *xfbuf
)
{
	char	*xfargs[7];
	XF	myxf;
	FVECT	updir, xdest, ydest;
	int	i;

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
	if (xfbuf == NULL)
		return(1);
				/* return xf arguments as well */
	for (i = 0; xfargs[i] != NULL; i++) {
		*xfbuf++ = ' ';
		strcpy(xfbuf, xfargs[i]);
		while (*xfbuf) ++xfbuf;
	}
	return(1);
}
