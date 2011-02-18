#ifndef lint
static const char RCSid[] = "$Id: bsdf_m.c,v 3.2 2011/02/18 00:41:44 greg Exp $";
#endif
/*
 *  bsdf_m.c
 *  
 *  Definitions supporting BSDF matrices
 *
 *  Created by Greg Ward on 2/2/11.
 *  Copyright 2011 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include "ezxml.h"
#include "bsdf.h"
#include "bsdf_m.h"

#ifndef FTINY
#define FTINY		1e-6
#endif

/* Function return codes */
#define	RC_GOOD		1
#define RC_FAIL		0
#define RC_FORMERR	(-1)
#define RC_DATERR	(-2)
#define RC_UNSUPP	(-3)
#define RC_INTERR	(-4)
#define RC_MEMERR	(-5)

#define MAXLATS		46		/* maximum number of latitudes */

/* BSDF angle specification */
typedef struct {
	char	name[64];		/* basis name */
	int	nangles;		/* total number of directions */
	struct {
		float	tmin;			/* starting theta */
		int	nphis;			/* number of phis (0 term) */
	} lat[MAXLATS+1];		/* latitudes */
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

static int
fequal(double a, double b)
{
	if (b != .0)
		a = a/b - 1.;
	return (a <= 1e-6) & (a >= -1e-6);
}

/* returns the name of the given tag */
#ifdef ezxml_name
#undef ezxml_name
static char *
ezxml_name(ezxml_t xml)
{
	if (xml == NULL)
		return NULL;
	return xml->name;
}
#endif

/* returns the given tag's character content or empty string if none */
#ifdef ezxml_txt
#undef ezxml_txt
static char *
ezxml_txt(ezxml_t xml)
{
	if (xml == NULL)
		return "";
	return xml->txt;
}
#endif

/* Convert error to standard BSDF code */
static SDError
convert_errcode(int ec)
{
	switch (ec) {
	case RC_GOOD:
		return SDEnone;
	case RC_FORMERR:
		return SDEformat;
	case RC_DATERR:
		return SDEdata;
	case RC_UNSUPP:
		return SDEsupport;
	case RC_INTERR:
		return SDEinternal;
	case RC_MEMERR:
		return SDEmemory;
	}
	return SDEunknown;
}

/* Allocate a BSDF matrix of the given size */
static SDMat *
SDnewMatrix(int ni, int no)
{
	SDMat		*sm;

	if ((ni <= 0) | (no <= 0)) {
		strcpy(SDerrorDetail, "Empty BSDF matrix request");
		return NULL;
	}
	sm = (SDMat *)malloc(sizeof(SDMat) + (ni*no - 1)*sizeof(float));
	if (sm == NULL) {
		sprintf(SDerrorDetail, "Cannot allocate %dx%d BSDF matrix",
				ni, no);
		return NULL;
	}
	memset(sm, 0, sizeof(SDMat)-sizeof(float));
	sm->ninc = ni;
	sm->nout = no;
	
	return sm;
}

/* Free a BSDF matrix */
#define	SDfreeMatrix		free

/* get vector for this angle basis index */
static int
ab_getvec(FVECT v, int ndx, double randX, void *p)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	double		rx[2];
	int		li;
	double		pol, azi, d;
	
	if ((ndx < 0) | (ndx >= ab->nangles))
		return RC_FAIL;
	for (li = 0; ndx >= ab->lat[li].nphis; li++)
		ndx -= ab->lat[li].nphis;
	SDmultiSamp(rx, 2, randX);
	pol = M_PI/180.*( (1.-rx[0])*ab->lat[li].tmin +
				rx[0]*ab->lat[li+1].tmin );
	azi = 2.*M_PI*(ndx + rx[1] - .5)/ab->lat[li].nphis;
	v[2] = d = cos(pol);
	d = sqrt(1. - d*d);	/* sin(pol) */
	v[0] = cos(azi)*d;
	v[1] = sin(azi)*d;
	return RC_GOOD;
}

/* get index corresponding to the given vector */
static int
ab_getndx(const FVECT v, void *p)
{
	ANGLE_BASIS  *ab = (ANGLE_BASIS *)p;
	int	li, ndx;
	double	pol, azi, d;

	if (v == NULL)
		return -1;
	if ((v[2] < .0) | (v[2] > 1.0))
		return -1;
	pol = 180.0/M_PI*acos(v[2]);
	azi = 180.0/M_PI*atan2(v[1], v[0]);
	if (azi < 0.0) azi += 360.0;
	for (li = 1; ab->lat[li].tmin <= pol; li++)
		if (!ab->lat[li].nphis)
			return -1;
	--li;
	ndx = (int)((1./360.)*azi*ab->lat[li].nphis + 0.5);
	if (ndx >= ab->lat[li].nphis) ndx = 0;
	while (li--)
		ndx += ab->lat[li].nphis;
	return ndx;
}

/* compute square of real value */
static double sq(double x) { return x*x; }

/* get projected solid angle for this angle basis index */
static double
ab_getohm(int ndx, void *p)
{
	static int	last_li = -1;
	static double	last_ohm;
	ANGLE_BASIS	*ab = (ANGLE_BASIS *)p;
	int		li;
	double		theta, theta1;
	
	if ((ndx < 0) | (ndx >= ab->nangles))
		return -1.;
	for (li = 0; ndx >= ab->lat[li].nphis; li++)
		ndx -= ab->lat[li].nphis;
	if (li == last_li)			/* cached latitude? */
		return last_ohm;
	last_li = li;
	theta1 = M_PI/180. * ab->lat[li+1].tmin;
	if (ab->lat[li].nphis == 1)		/* special case */
		return last_ohm = M_PI*(1. - sq(cos(theta1)));
	theta = M_PI/180. * ab->lat[li].tmin;
	return last_ohm = M_PI*(sq(cos(theta)) - sq(cos(theta1))) /
				(double)ab->lat[li].nphis;
}

/* get reverse vector for this angle basis index */
static int
ab_getvecR(FVECT v, int ndx, double randX, void *p)
{
	int	na = (*(ANGLE_BASIS *)p).nangles;

	if (!ab_getvec(v, ndx, randX, p))
		return RC_FAIL;

	v[0] = -v[0];
	v[1] = -v[1];
	v[2] = -v[2];

	return RC_GOOD;
}

/* get index corresponding to the reverse vector */
static int
ab_getndxR(const FVECT v, void *p)
{
	FVECT  v2;
	
	v2[0] = -v[0];
	v2[1] = -v[1];
	v2[2] = -v[2];

	return ab_getndx(v2, p);
}

/* load custom BSDF angle basis */
static int
load_angle_basis(ezxml_t wab)
{
	char	*abname = ezxml_txt(ezxml_child(wab, "AngleBasisName"));
	ezxml_t	wbb;
	int	i;
	
	if (!abname || !*abname)
		return RC_FAIL;
	for (i = nabases; i--; )
		if (!strcasecmp(abname, abase_list[i].name))
			return RC_GOOD;	/* assume it's the same */
	if (nabases >= MAXABASES) {
		sprintf(SDerrorDetail, "Out of angle bases reading '%s'",
				abname);
		return RC_INTERR;
	}
	strcpy(abase_list[nabases].name, abname);
	abase_list[nabases].nangles = 0;
	for (i = 0, wbb = ezxml_child(wab, "AngleBasisBlock");
			wbb != NULL; i++, wbb = wbb->next) {
		if (i >= MAXLATS) {
			sprintf(SDerrorDetail, "Too many latitudes for '%s'",
				abname);
			return RC_INTERR;
		}
		abase_list[nabases].lat[i+1].tmin = atof(ezxml_txt(
				ezxml_child(ezxml_child(wbb,
					"ThetaBounds"), "UpperTheta")));
		if (!i)
			abase_list[nabases].lat[i].tmin =
					-abase_list[nabases].lat[i+1].tmin;
		else if (!fequal(atof(ezxml_txt(ezxml_child(ezxml_child(wbb,
					"ThetaBounds"), "LowerTheta"))),
				abase_list[nabases].lat[i].tmin)) {
			sprintf(SDerrorDetail, "Theta values disagree in '%s'",
				abname);
			return RC_DATERR;
		}
		abase_list[nabases].nangles +=
			abase_list[nabases].lat[i].nphis =
				atoi(ezxml_txt(ezxml_child(wbb, "nPhis")));
		if (abase_list[nabases].lat[i].nphis <= 0 ||
				(abase_list[nabases].lat[i].nphis == 1 &&
				abase_list[nabases].lat[i].tmin > FTINY)) {
			sprintf(SDerrorDetail, "Illegal phi count in '%s'",
				abname);
			return RC_DATERR;
		}
	}
	abase_list[nabases++].lat[i].nphis = 0;
	return RC_GOOD;
}

/* compute min. proj. solid angle and max. direct hemispherical scattering */
static int
get_extrema(SDSpectralDF *df)
{
	SDMat	*dp = (SDMat *)df->comp[0].dist;
	double	*ohma;
	int	i, o;
					/* initialize extrema */
	df->minProjSA = M_PI;
	df->maxHemi = .0;
	ohma = (double *)malloc(dp->nout*sizeof(double));
	if (ohma == NULL)
		return RC_MEMERR;
					/* get outgoing solid angles */
	for (o = dp->nout; o--; )
		if ((ohma[o] = mBSDF_outohm(dp,o)) < df->minProjSA)
			df->minProjSA = ohma[o];
					/* compute hemispherical sums */
	for (i = dp->ninc; i--; ) {
		double	hemi = .0;
		for (o = dp->nout; o--; )
			hemi += ohma[o] * mBSDF_value(dp, i, o);
		if (hemi > df->maxHemi)
			df->maxHemi = hemi;
	}
	free(ohma);
					/* need incoming solid angles, too? */
	if (dp->ninc < dp->nout || dp->ib_ohm != dp->ob_ohm ||
				dp->ib_priv != dp->ob_priv) {
		double	ohm;
		for (i = dp->ninc; i--; )
			if ((ohm = mBSDF_incohm(dp,i)) < df->minProjSA)
				df->minProjSA = ohm;
	}
	return (df->maxHemi <= 1.01);
}

/* skip integer in string */
static char *
i_skip(char  *s)
{
	while (isspace(*s))
		s++;
	if (*s == '-' || *s == '+')
		s++;
	if (!isdigit(*s))
		return(NULL);
	do
		s++;
	while (isdigit(*s));
	return(s);
}

/* skip float in string */
static char *
f_skip(char  *s)
{
	register char  *cp;

	while (isspace(*s))
		s++;
	if (*s == '-' || *s == '+')
		s++;
	cp = s;
	while (isdigit(*cp))
		cp++;
	if (*cp == '.') {
		cp++; s++;
		while (isdigit(*cp))
			cp++;
	}
	if (cp == s)
		return(NULL);
	if (*cp == 'e' || *cp == 'E')
		return(i_skip(cp+1));
	return(cp);
}

/* load BSDF distribution for this wavelength */
static int
load_bsdf_data(SDData *sd, ezxml_t wdb, int rowinc)
{
	char		*cbasis = ezxml_txt(ezxml_child(wdb,"ColumnAngleBasis"));
	char		*rbasis = ezxml_txt(ezxml_child(wdb,"RowAngleBasis"));
	SDSpectralDF	*df;
	SDMat		*dp;
	char		*sdata;
	int		inbi, outbi;
	int		i;

	if ((!cbasis || !*cbasis) | (!rbasis || !*rbasis)) {
		sprintf(SDerrorDetail, "Missing column/row basis for BSDF '%s'",
				sd->name);
		return RC_FORMERR;
	}
					/* allocate BSDF component */
	sdata = ezxml_txt(ezxml_child(wdb, "WavelengthDataDirection"));
	if (!strcasecmp(sdata, "Transmission Front")) {
		if (sd->tf != NULL)
			SDfreeSpectralDF(sd->tf);
		if ((sd->tf = SDnewSpectralDF(1)) == NULL)
			return RC_MEMERR;
		df = sd->tf;
	} else if (!strcasecmp(sdata, "Reflection Front")) {
		if (sd->rf != NULL)
			SDfreeSpectralDF(sd->rf);
		if ((sd->rf = SDnewSpectralDF(1)) == NULL)
			return RC_MEMERR;
		df = sd->rf;
	} else if (!strcasecmp(sdata, "Reflection Back")) {
		if (sd->rb != NULL)
			SDfreeSpectralDF(sd->rb);
		if ((sd->rb = SDnewSpectralDF(1)) == NULL)
			return RC_MEMERR;
		df = sd->rb;
	} else
		return RC_FAIL;
					/* get angle bases */
	sdata = ezxml_txt(ezxml_child(wdb,"ColumnAngleBasis"));
	if (!sdata || !*sdata) {
		sprintf(SDerrorDetail, "Missing column basis for BSDF '%s'",
				sd->name);
		return RC_FORMERR;
	}
	for (inbi = nabases; inbi--; )
		if (!strcasecmp(cbasis, abase_list[inbi].name))
			break;
	if (inbi < 0) {
		sprintf(SDerrorDetail, "Undefined ColumnAngleBasis '%s'",
				cbasis);
		return RC_FORMERR;
	}
	sdata = ezxml_txt(ezxml_child(wdb,"RowAngleBasis"));
	if (!sdata || !*sdata) {
		sprintf(SDerrorDetail, "Missing row basis for BSDF '%s'",
				sd->name);
		return RC_FORMERR;
	}
	for (outbi = nabases; outbi--; )
		if (!strcasecmp(rbasis, abase_list[outbi].name))
			break;
	if (outbi < 0) {
		sprintf(SDerrorDetail, "Undefined RowAngleBasis '%s'", cbasis);
		return RC_FORMERR;
	}
					/* allocate BSDF matrix */
	dp = SDnewMatrix(abase_list[inbi].nangles, abase_list[outbi].nangles);
	if (dp == NULL)
		return RC_MEMERR;
	dp->ib_priv = (void *)&abase_list[inbi];
	dp->ob_priv = (void *)&abase_list[outbi];
	if (df == sd->tf) {
		dp->ib_vec = ab_getvecR;
		dp->ib_ndx = ab_getndxR;
		dp->ob_vec = ab_getvec;
		dp->ob_ndx = ab_getndx;
	} else if (df == sd->rf) {
		dp->ib_vec = ab_getvec;
		dp->ib_ndx = ab_getndx;
		dp->ob_vec = ab_getvec;
		dp->ob_ndx = ab_getndx;
	} else /* df == sd->rb */ {
		dp->ib_vec = ab_getvecR;
		dp->ib_ndx = ab_getndxR;
		dp->ob_vec = ab_getvecR;
		dp->ob_ndx = ab_getndxR;
	}
	dp->ib_ohm = ab_getohm;
	dp->ob_ohm = ab_getohm;
	df->comp[0].cspec[0] = c_dfcolor; /* XXX monochrome for now */
	df->comp[0].dist = dp;
	df->comp[0].func = &SDhandleMtx;
					/* read BSDF data */
	sdata  = ezxml_txt(ezxml_child(wdb,"ScatteringData"));
	if (!sdata || !*sdata) {
		sprintf(SDerrorDetail, "Missing BSDF ScatteringData in '%s'",
				sd->name);
		return RC_FORMERR;
	}
	for (i = 0; i < dp->ninc*dp->nout; i++) {
		char  *sdnext = f_skip(sdata);
		if (sdnext == NULL) {
			sprintf(SDerrorDetail,
				"Bad/missing BSDF ScatteringData in '%s'",
					sd->name);
			return RC_FORMERR;
		}
		while (*sdnext && isspace(*sdnext))
			sdnext++;
		if (*sdnext == ',') sdnext++;
		if (rowinc) {
			int	r = i/dp->nout;
			int	c = i - c*dp->nout;
			mBSDF_value(dp,r,c) = atof(sdata);
		} else
			dp->bsdf[i] = atof(sdata);
		sdata = sdnext;
	}
	return get_extrema(df);
}

/* Subtract minimum (diffuse) scattering amount from BSDF */
static double
subtract_min(SDMat *sm)
{
	float	minv = sm->bsdf[0];
	int	n = sm->ninc*sm->nout;
	int	i;
	
	for (i = n; --i; )
		if (sm->bsdf[i] < minv)
			minv = sm->bsdf[i];
	for (i = n; i--; )
		sm->bsdf[i] -= minv;

	return minv*M_PI;		/* be sure to include multiplier */
}

/* Extract and separate diffuse portion of BSDF */
static void
extract_diffuse(SDValue *dv, SDSpectralDF *df)
{
	int	n;

	if (df == NULL || df->ncomp <= 0) {
		dv->spec = c_dfcolor;
		dv->cieY = .0;
		return;
	}
	dv->spec = df->comp[0].cspec[0];
	dv->cieY = subtract_min((SDMat *)df->comp[0].dist);
					/* in case of multiple components */
	for (n = df->ncomp; --n; ) {
		double	ymin = subtract_min((SDMat *)df->comp[n].dist);
		c_cmix(&dv->spec, dv->cieY, &dv->spec, ymin, &df->comp[n].cspec[0]);
		dv->cieY += ymin;
	}
	df->maxHemi -= dv->cieY;	/* correct minimum hemispherical */
	dv->spec.clock++;		/* make sure everything is set */
	c_ccvt(&dv->spec, C_CSXY+C_CSSPEC);
}

/* Load a BSDF matrix from an open XML file */
SDError
SDloadMtx(SDData *sd, ezxml_t fl)
{
	ezxml_t			wtl, wld, wdb;
	int			rowIn;
	struct BSDF_data	*dp;
	char			*txt;
	int			rval;
	
	if (strcmp(ezxml_name(fl), "WindowElement")) {
		sprintf(SDerrorDetail,
			"BSDF \"%s\": top level node not 'WindowElement'",
				sd->name);
		return SDEformat;
	}
	wtl = ezxml_child(ezxml_child(fl, "Optical"), "Layer");
	txt = ezxml_txt(ezxml_child(ezxml_child(wtl,
			"DataDefinition"), "IncidentDataStructure"));
	if (!strcasecmp(txt, "Rows"))
		rowIn = 1;
	else if (!strcasecmp(txt, "Columns"))
		rowIn = 0;
	else {
		sprintf(SDerrorDetail,
			"BSDF \"%s\": unsupported IncidentDataStructure",
				sd->name);
		return SDEsupport;
	}
				/* get angle basis */
	rval = load_angle_basis(ezxml_child(ezxml_child(wtl,
				"DataDefinition"), "AngleBasis"));
	if (rval < 0)
		goto err_return;
				/* load BSDF components */
	for (wld = ezxml_child(wtl, "WavelengthData");
				wld != NULL; wld = wld->next) {
		if (strcasecmp(ezxml_txt(ezxml_child(wld,"Wavelength")),
				"Visible"))
			continue;	/* just visible for now */
		for (wdb = ezxml_child(wld, "WavelengthDataBlock");
					wdb != NULL; wdb = wdb->next)
			if ((rval = load_bsdf_data(sd, wdb, rowIn)) < 0)
				goto err_return;
	}
				/* separate diffuse components */
	extract_diffuse(&sd->rLambFront, sd->rf);
	extract_diffuse(&sd->rLambBack, sd->rb);
	extract_diffuse(&sd->tLamb, sd->tf);
				/* return success */
	return SDEnone;
err_return:			/* jump here on failure */
	if (sd->rf != NULL) {
		SDfreeSpectralDF(sd->rf);
		sd->rf = NULL;
	}
	if (sd->rb != NULL) {
		SDfreeSpectralDF(sd->rb);
		sd->rb = NULL;
	}
	if (sd->tf != NULL) {
		SDfreeSpectralDF(sd->tf);
		sd->tf = NULL;
	}
	return convert_errcode(rval);
}

/* Get Matrix BSDF value */
static int
SDgetMtxBSDF(float coef[SDmaxCh], const FVECT outVec,
				const FVECT inVec, const void *dist)
{
	const SDMat	*dp = (const SDMat *)dist;
	int		i_ndx, o_ndx;
					/* get angle indices */
	i_ndx = mBSDF_incndx(dp, inVec);
	o_ndx = mBSDF_outndx(dp, outVec);
					/* try reciprocity if necessary */
	if ((i_ndx < 0) & (o_ndx < 0)) {
		i_ndx = mBSDF_incndx(dp, outVec);
		o_ndx = mBSDF_outndx(dp, inVec);
	}
	if ((i_ndx < 0) | (o_ndx < 0))
		return 0;		/* nothing from this component */
	coef[0] = mBSDF_value(dp, i_ndx, o_ndx);
	return 1;			/* XXX monochrome for now */
}

/* Query solid angle for vector */
static SDError
SDqueryMtxProjSA(double *psa, const FVECT vec, int qflags, const void *dist)
{
	const SDMat	*dp = (const SDMat *)dist;

	if (!(qflags & SDqueryInc+SDqueryOut))
		return SDEargument;
	if (qflags & SDqueryInc) {
		double	inc_psa = mBSDF_incohm(dp, mBSDF_incndx(dp, vec));
		if (inc_psa < .0)
			return SDEinternal;
		switch (qflags & SDqueryMin+SDqueryMax) {
		case SDqueryMax:
			if (inc_psa > psa[0])
				psa[0] = inc_psa;
			break;
		case SDqueryMin+SDqueryMax:
			if (inc_psa > psa[1])
				psa[1] = inc_psa;
			/* fall through */
		case SDqueryMin:
			if (inc_psa < psa[0])
				psa[0] = inc_psa;
			break;
		case 0:
			psa[0] = inc_psa;
			break;
		}
	}
	if (qflags & SDqueryOut) {
		double	out_psa = mBSDF_outohm(dp, mBSDF_outndx(dp, vec));
		if (out_psa < .0)
			return SDEinternal;
		switch (qflags & SDqueryMin+SDqueryMax) {
		case SDqueryMax:
			if (out_psa > psa[0])
				psa[0] = out_psa;
			break;
		case SDqueryMin+SDqueryMax:
			if (out_psa > psa[1])
				psa[1] = out_psa;
			/* fall through */
		case SDqueryMin:
			if (out_psa < psa[0])
				psa[0] = out_psa;
			break;
		case 0:
			psa[(qflags&SDqueryInc)!=0] = out_psa;
			break;
		}
	}
	return SDEnone;
}

/* Compute new cumulative distribution from BSDF */
static int
make_cdist(SDMatCDst *cd, const FVECT inVec, SDMat *dp, int rev)
{
	const unsigned	maxval = ~0;
	double		*cmtab, scale;
	int		o;

	cmtab = (double *)malloc((cd->calen+1)*sizeof(double));
	if (cmtab == NULL)
		return 0;
	cmtab[0] = .0;
	for (o = 0; o < cd->calen; o++) {
		if (rev)
			cmtab[o+1] = mBSDF_value(dp, o, cd->indx) *
					(*dp->ib_ohm)(o, dp->ib_priv);
		else
			cmtab[o+1] = mBSDF_value(dp, cd->indx, o) *
					(*dp->ob_ohm)(o, dp->ob_priv);
		cmtab[o+1] += cmtab[o];
	}
	cd->cTotal = cmtab[cd->calen];
	scale = (double)maxval / cd->cTotal;
	cd->carr[0] = 0;
	for (o = 1; o < cd->calen; o++)
		cd->carr[o] = scale*cmtab[o] + .5;
	cd->carr[cd->calen] = maxval;
	free(cmtab);
	return 1;
}

/* Get cumulative distribution for matrix BSDF */
static const SDCDst *
SDgetMtxCDist(const FVECT inVec, SDComponent *sdc)
{
	SDMat		*dp = (SDMat *)sdc->dist;
	int		reverse;
	SDMatCDst	myCD;
	SDMatCDst	*cd, *cdlast;

	if (dp == NULL)
		return NULL;
	memset(&myCD, 0, sizeof(myCD));
	myCD.indx = mBSDF_incndx(dp, inVec);
	if (myCD.indx >= 0) {
		myCD.ob_priv = dp->ob_priv;
		myCD.ob_vec = dp->ob_vec;
		myCD.calen = dp->nout;
		reverse = 0;
	} else {			/* try reciprocity */
		myCD.indx = mBSDF_outndx(dp, inVec);
		if (myCD.indx < 0)
			return NULL;
		myCD.ob_priv = dp->ib_priv;
		myCD.ob_vec = dp->ib_vec;
		myCD.calen = dp->ninc;
		reverse = 1;
	}
	cdlast = NULL;			/* check for it in cache list */
	for (cd = (SDMatCDst *)sdc->cdList;
				cd != NULL; cd = (SDMatCDst *)cd->next) {
		if (cd->indx == myCD.indx && (cd->calen == myCD.calen) &
					(cd->ob_priv == myCD.ob_priv) &
					(cd->ob_vec == myCD.ob_vec))
			break;
		cdlast = cd;
	}
	if (cd == NULL) {		/* need to allocate new entry */
		cd = (SDMatCDst *)malloc(sizeof(SDMatCDst) +
					myCD.calen*sizeof(myCD.carr[0]));
		if (cd == NULL)
			return NULL;
		*cd = myCD;		/* compute cumulative distribution */
		if (!make_cdist(cd, inVec, dp, reverse)) {
			free(cd);
			return NULL;
		}
		cdlast = cd;
	}
	if (cdlast != NULL) {		/* move entry to head of cache list */
		cdlast->next = cd->next;
		cd->next = sdc->cdList;
		sdc->cdList = (SDCDst *)cd;
	}
	return (SDCDst *)cd;		/* ready to go */
}

/* Sample cumulative distribution */
static SDError
SDsampMtxCDist(FVECT outVec, double randX, const SDCDst *cdp)
{
	const unsigned	maxval = ~0;
	const SDMatCDst	*mcd = (const SDMatCDst *)cdp;
	const unsigned	target = randX*maxval;
	int		i, iupper, ilower;
					/* binary search to find index */
	ilower = 0; iupper = mcd->calen;
	while ((i = (iupper + ilower) >> 1) != ilower)
		if ((long)target >= (long)mcd->carr[i])
			ilower = i;
		else
			iupper = i;
					/* localize random position */
	randX = (randX*maxval - mcd->carr[ilower]) /
			(double)(mcd->carr[iupper] - mcd->carr[ilower]);
					/* convert index to vector */
	if ((*mcd->ob_vec)(outVec, i, randX, mcd->ob_priv))
		return SDEnone;
	strcpy(SDerrorDetail, "BSDF sampling fault");
	return SDEinternal;
}

/* Fixed resolution BSDF methods */
SDFunc			SDhandleMtx = {
				&SDgetMtxBSDF,
				&SDqueryMtxProjSA,
				&SDgetMtxCDist,
				&SDsampMtxCDist,
				&SDfreeMatrix,
			};
