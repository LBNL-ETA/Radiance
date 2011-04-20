#ifndef lint
static const char RCSid[] = "$Id: bsdf.c,v 2.23 2011/04/20 14:44:05 greg Exp $";
#endif
/*
 *  bsdf.c
 *  
 *  Definitions for bidirectional scattering distribution functions.
 *
 *  Created by Greg Ward on 1/10/11.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ezxml.h"
#include "hilbert.h"
#include "bsdf.h"
#include "bsdf_m.h"
#include "bsdf_t.h"

/* English ASCII strings corresponding to ennumerated errors */
const char		*SDerrorEnglish[] = {
				"No error",
				"Memory error",
				"File input/output error",
				"File format error",
				"Illegal argument",
				"Invalid data",
				"Unsupported feature",
				"Internal program error",
				"Unknown error"
			};

/* Additional information on last error (ASCII English) */
char			SDerrorDetail[256];

/* Cache of loaded BSDFs */
struct SDCache_s	*SDcacheList = NULL;

/* Retain BSDFs in cache list */
int			SDretainSet = SDretainNone;

/* Report any error to the indicated stream (in English) */
SDError
SDreportEnglish(SDError ec, FILE *fp)
{
	if (!ec)
		return SDEnone;
	if ((ec < SDEnone) | (ec > SDEunknown)) {
		SDerrorDetail[0] = '\0';
		ec = SDEunknown;
	}
	if (fp == NULL)
		return ec;
	fputs(SDerrorEnglish[ec], fp);
	if (SDerrorDetail[0]) {
		fputs(": ", fp);
		fputs(SDerrorDetail, fp);
	}
	fputc('\n', fp);
	if (fp != stderr)
		fflush(fp);
	return ec;
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
	sprintf(SDerrorDetail, "Unknown dimensional unit '%s'", unit);
	return(-1.);
}

/* Load geometric dimensions and description (if any) */
static SDError
SDloadGeometry(SDData *sd, ezxml_t wdb)
{
	ezxml_t		geom;
	double		cfact;
	const char	*fmt, *mgfstr;

	if (wdb == NULL)		/* no geometry section? */
		return SDEnone;
	sd->dim[0] = sd->dim[1] = sd->dim[2] = .0;
	if ((geom = ezxml_child(wdb, "Width")) != NULL)
		sd->dim[0] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((geom = ezxml_child(wdb, "Height")) != NULL)
		sd->dim[1] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((geom = ezxml_child(wdb, "Thickness")) != NULL)
		sd->dim[2] = atof(ezxml_txt(geom)) *
				to_meters(ezxml_attr(geom, "unit"));
	if ((sd->dim[0] < 0) | (sd->dim[1] < 0) | (sd->dim[2] < 0)) {
		sprintf(SDerrorDetail, "Negative size in \"%s\"", sd->name);
		return SDEdata;
	}
	if ((geom = ezxml_child(wdb, "Geometry")) == NULL ||
			(mgfstr = ezxml_txt(geom)) == NULL)
		return SDEnone;
	if ((fmt = ezxml_attr(geom, "format")) != NULL &&
			strcasecmp(fmt, "MGF")) {
		sprintf(SDerrorDetail,
			"Unrecognized geometry format '%s' in \"%s\"",
					fmt, sd->name);
		return SDEsupport;
	}
	cfact = to_meters(ezxml_attr(geom, "unit"));
	sd->mgf = (char *)malloc(strlen(mgfstr)+32);
	if (sd->mgf == NULL) {
		strcpy(SDerrorDetail, "Out of memory in SDloadGeometry");
		return SDEmemory;
	}
	if (cfact < 0.99 || cfact > 1.01)
		sprintf(sd->mgf, "xf -s %.5f\n%s\nxf\n", cfact, mgfstr);
	else
		strcpy(sd->mgf, mgfstr);
	return SDEnone;
}

/* Load a BSDF struct from the given file (free first and keep name) */
SDError
SDloadFile(SDData *sd, const char *fname)
{
	SDError		lastErr;
	ezxml_t		fl, wtl;

	if ((sd == NULL) | (fname == NULL || !*fname))
		return SDEargument;
				/* free old data, keeping name */
	SDfreeBSDF(sd);
				/* parse XML file */
	fl = ezxml_parse_file(fname);
	if (fl == NULL) {
		sprintf(SDerrorDetail, "Cannot open BSDF \"%s\"", fname);
		return SDEfile;
	}
	if (ezxml_error(fl)[0]) {
		sprintf(SDerrorDetail, "BSDF \"%s\" %s", fname, ezxml_error(fl));
		ezxml_free(fl);
		return SDEformat;
	}
	if (strcmp(ezxml_name(fl), "WindowElement")) {
		sprintf(SDerrorDetail,
			"BSDF \"%s\": top level node not 'WindowElement'",
				sd->name);
		ezxml_free(fl);
		return SDEformat;
	}
	wtl = ezxml_child(ezxml_child(fl, "Optical"), "Layer");
	if (wtl == NULL) {
		sprintf(SDerrorDetail, "BSDF \"%s\": no optical layer'",
				sd->name);
		ezxml_free(fl);
		return SDEformat;
	}
				/* load geometry if present */
	lastErr = SDloadGeometry(sd, ezxml_child(wtl, "Material"));
	if (lastErr)
		return lastErr;
				/* try loading variable resolution data */
	lastErr = SDloadTre(sd, wtl);
				/* check our result */
	switch (lastErr) {
	case SDEformat:
	case SDEdata:
	case SDEsupport:	/* possibly we just tried the wrong format */
		lastErr = SDloadMtx(sd, wtl);
		break;
	default:		/* variable res. OK else serious error */
		break;
	}
				/* done with XML file */
	ezxml_free(fl);
	
	if (lastErr) {		/* was there a load error? */
		SDfreeBSDF(sd);
		return lastErr;
	}
				/* remove any insignificant components */
	if (sd->rf != NULL && sd->rf->maxHemi <= .001) {
		SDfreeSpectralDF(sd->rf); sd->rf = NULL;
	}
	if (sd->rb != NULL && sd->rb->maxHemi <= .001) {
		SDfreeSpectralDF(sd->rb); sd->rb = NULL;
	}
	if (sd->tf != NULL && sd->tf->maxHemi <= .001) {
		SDfreeSpectralDF(sd->tf); sd->tf = NULL;
	}
				/* return success */
	return SDEnone;
}

/* Allocate new spectral distribution function */
SDSpectralDF *
SDnewSpectralDF(int nc)
{
	SDSpectralDF	*df;
	
	if (nc <= 0) {
		strcpy(SDerrorDetail, "Zero component spectral DF request");
		return NULL;
	}
	df = (SDSpectralDF *)malloc(sizeof(SDSpectralDF) +
					(nc-1)*sizeof(SDComponent));
	if (df == NULL) {
		sprintf(SDerrorDetail,
				"Cannot allocate %d component spectral DF", nc);
		return NULL;
	}
	df->minProjSA = .0;
	df->maxHemi = .0;
	df->ncomp = nc;
	memset(df->comp, 0, nc*sizeof(SDComponent));
	return df;
}

/* Free cached cumulative distributions for BSDF component */
void
SDfreeCumulativeCache(SDSpectralDF *df)
{
	int	n;
	SDCDst	*cdp;

	if (df == NULL)
		return;
	for (n = df->ncomp; n-- > 0; )
		while ((cdp = df->comp[n].cdList) != NULL) {
			df->comp[n].cdList = cdp->next;
			free(cdp);
		}
}

/* Free a spectral distribution function */
void
SDfreeSpectralDF(SDSpectralDF *df)
{
	int	n;

	if (df == NULL)
		return;
	SDfreeCumulativeCache(df);
	for (n = df->ncomp; n-- > 0; )
		(*df->comp[n].func->freeSC)(df->comp[n].dist);
	free(df);
}

/* Shorten file path to useable BSDF name, removing suffix */
void
SDclipName(char *res, const char *fname)
{
	const char	*cp, *dot = NULL;
	
	for (cp = fname; *cp; cp++)
		if (*cp == '.')
			dot = cp;
	if ((dot == NULL) | (dot < fname+2))
		dot = cp;
	if (dot - fname >= SDnameLn)
		fname = dot - SDnameLn + 1;
	while (fname < dot)
		*res++ = *fname++;
	*res = '\0';
}

/* Initialize an unused BSDF struct (simply clears to zeroes) */
void
SDclearBSDF(SDData *sd, const char *fname)
{
	if (sd == NULL)
		return;
	memset(sd, 0, sizeof(SDData));
	if (fname == NULL)
		return;
	SDclipName(sd->name, fname);
}

/* Free data associated with BSDF struct */
void
SDfreeBSDF(SDData *sd)
{
	if (sd == NULL)
		return;
	if (sd->mgf != NULL) {
		free(sd->mgf);
		sd->mgf = NULL;
	}
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
	sd->rLambFront.cieY = .0;
	sd->rLambFront.spec.flags = 0;
	sd->rLambBack.cieY = .0;
	sd->rLambBack.spec.flags = 0;
	sd->tLamb.cieY = .0;
	sd->tLamb.spec.flags = 0;
}

/* Find writeable BSDF by name, or allocate new cache entry if absent */
SDData *
SDgetCache(const char *bname)
{
	struct SDCache_s	*sdl;
	char			sdnam[SDnameLn];

	if (bname == NULL)
		return NULL;

	SDclipName(sdnam, bname);
	for (sdl = SDcacheList; sdl != NULL; sdl = sdl->next)
		if (!strcmp(sdl->bsdf.name, sdnam)) {
			sdl->refcnt++;
			return &sdl->bsdf;
		}

	sdl = (struct SDCache_s *)calloc(1, sizeof(struct SDCache_s));
	if (sdl == NULL)
		return NULL;

	strcpy(sdl->bsdf.name, sdnam);
	sdl->next = SDcacheList;
	SDcacheList = sdl;

	sdl->refcnt = 1;
	return &sdl->bsdf;
}

/* Get loaded BSDF from cache (or load and cache it on first call) */
/* Report any problem to stderr and return NULL on failure */
const SDData *
SDcacheFile(const char *fname)
{
	SDData		*sd;
	SDError		ec;
	
	if (fname == NULL || !*fname)
		return NULL;
	SDerrorDetail[0] = '\0';
	if ((sd = SDgetCache(fname)) == NULL) {
		SDreportEnglish(SDEmemory, stderr);
		return NULL;
	}
	if (!SDisLoaded(sd) && (ec = SDloadFile(sd, fname))) {
		SDreportEnglish(ec, stderr);
		SDfreeCache(sd);
		return NULL;
	}
	return sd;
}

/* Free a BSDF from our cache (clear all if NULL) */
void
SDfreeCache(const SDData *sd)
{
	struct SDCache_s	*sdl, *sdLast = NULL;

	if (sd == NULL) {		/* free entire list */
		while ((sdl = SDcacheList) != NULL) {
			SDcacheList = sdl->next;
			SDfreeBSDF(&sdl->bsdf);
			free(sdl);
		}
		return;
	}
	for (sdl = SDcacheList; sdl != NULL; sdl = (sdLast=sdl)->next)
		if (&sdl->bsdf == sd)
			break;
	if (sdl == NULL || (sdl->refcnt -= (sdl->refcnt > 0)))
		return;			/* missing or still in use */
					/* keep unreferenced data? */
	if (SDisLoaded(sd) && SDretainSet) {
		if (SDretainSet == SDretainAll)
			return;		/* keep everything */
					/* else free cumulative data */
		SDfreeCumulativeCache(sd->rf);
		SDfreeCumulativeCache(sd->rb);
		SDfreeCumulativeCache(sd->tf);
		return;
	}
					/* remove from list and free */
	if (sdLast == NULL)
		SDcacheList = sdl->next;
	else
		sdLast->next = sdl->next;
	SDfreeBSDF(&sdl->bsdf);
	free(sdl);
}

/* Sample an individual BSDF component */
SDError
SDsampComponent(SDValue *sv, FVECT outVec, const FVECT inVec,
			double randX, SDComponent *sdc)
{
	float		coef[SDmaxCh];
	SDError		ec;
	const SDCDst	*cd;
	double		d;
	int		n;
					/* check arguments */
	if ((sv == NULL) | (outVec == NULL) | (inVec == NULL) | (sdc == NULL))
		return SDEargument;
					/* get cumulative distribution */
	cd = (*sdc->func->getCDist)(inVec, sdc);
	if (cd == NULL)
		return SDEmemory;
	if (cd->cTotal <= 1e-7) {	/* anything to sample? */
		sv->spec = c_dfcolor;
		sv->cieY = .0;
		memset(outVec, 0, 3*sizeof(double));
		return SDEnone;
	}
	sv->cieY = cd->cTotal;
					/* compute sample direction */
	ec = (*sdc->func->sampCDist)(outVec, randX, cd);
	if (ec)
		return ec;
					/* get BSDF color */
	n = (*sdc->func->getBSDFs)(coef, outVec, inVec, sdc->dist);
	if (n <= 0) {
		strcpy(SDerrorDetail, "BSDF sample value error");
		return SDEinternal;
	}
	sv->spec = sdc->cspec[0];
	d = coef[0];
	while (--n) {
		c_cmix(&sv->spec, d, &sv->spec, coef[n], &sdc->cspec[n]);
		d += coef[n];
	}
					/* make sure everything is set */
	c_ccvt(&sv->spec, C_CSXY+C_CSSPEC);
	return SDEnone;
}

#define	MS_MAXDIM	15

/* Convert 1-dimensional random variable to N-dimensional */
void
SDmultiSamp(double t[], int n, double randX)
{
	unsigned	nBits;
	double		scale;
	bitmask_t	ndx, coord[MS_MAXDIM];
	
	while (n > MS_MAXDIM)		/* punt for higher dimensions */
		t[--n] = rand()*(1./(RAND_MAX+.5));
	nBits = (8*sizeof(bitmask_t) - 1) / n;
	ndx = randX * (double)((bitmask_t)1 << (nBits*n));
					/* get coordinate on Hilbert curve */
	hilbert_i2c(n, nBits, ndx, coord);
					/* convert back to [0,1) range */
	scale = 1. / (double)((bitmask_t)1 << nBits);
	while (n--)
		t[n] = scale * ((double)coord[n] + rand()*(1./(RAND_MAX+.5)));
}

#undef MS_MAXDIM

/* Generate diffuse hemispherical sample */
static void
SDdiffuseSamp(FVECT outVec, int outFront, double randX)
{
					/* convert to position on hemisphere */
	SDmultiSamp(outVec, 2, randX);
	SDsquare2disk(outVec, outVec[0], outVec[1]);
	outVec[2] = 1. - outVec[0]*outVec[0] - outVec[1]*outVec[1];
	if (outVec[2] > 0)		/* a bit of paranoia */
		outVec[2] = sqrt(outVec[2]);
	if (!outFront)			/* going out back? */
		outVec[2] = -outVec[2];
}

/* Query projected solid angle coverage for non-diffuse BSDF direction */
SDError
SDsizeBSDF(double *projSA, const FVECT v1, const RREAL *v2,
				int qflags, const SDData *sd)
{
	SDSpectralDF	*rdf, *tdf;
	SDError		ec;
	int		i;
					/* check arguments */
	if ((projSA == NULL) | (v1 == NULL))
		return SDEargument;
					/* initialize extrema */
	switch (qflags) {
	case SDqueryMax:
		projSA[0] = .0;
		break;
	case SDqueryMin+SDqueryMax:
		projSA[1] = .0;
		/* fall through */
	case SDqueryMin:
		projSA[0] = 10.;
		break;
	case 0:
		return SDEargument;
	}
	if (v1[2] > 0)			/* front surface query? */
		rdf = sd->rf;
	else
		rdf = sd->rb;
	tdf = NULL;			/* transmitted component? */
	if (v2 != NULL && v1[2] > 0 ^ v2[2] > 0) {
		rdf = NULL;
		tdf = sd->tf;
	}
	ec = SDEdata;			/* run through components */
	for (i = (rdf==NULL) ? 0 : rdf->ncomp; i--; ) {
		ec = (*rdf->comp[i].func->queryProjSA)(projSA, v1, v2,
						qflags, rdf->comp[i].dist);
		if (ec)
			return ec;
	}
	for (i = (tdf==NULL) ? 0 : tdf->ncomp; i--; ) {
		ec = (*tdf->comp[i].func->queryProjSA)(projSA, v1, v2,
						qflags, tdf->comp[i].dist);
		if (ec)
			return ec;
	}
	if (ec) {			/* all diffuse? */
		projSA[0] = M_PI;
		if (qflags == SDqueryMin+SDqueryMax)
			projSA[1] = M_PI;
	}
	return SDEnone;
}

/* Return BSDF for the given incident and scattered ray vectors */
SDError
SDevalBSDF(SDValue *sv, const FVECT outVec, const FVECT inVec, const SDData *sd)
{
	int		inFront, outFront;
	SDSpectralDF	*sdf;
	float		coef[SDmaxCh];
	int		nch, i;
					/* check arguments */
	if ((sv == NULL) | (outVec == NULL) | (inVec == NULL) | (sd == NULL))
		return SDEargument;
					/* whose side are we on? */
	inFront = (inVec[2] > 0);
	outFront = (outVec[2] > 0);
					/* start with diffuse portion */
	if (inFront & outFront) {
		*sv = sd->rLambFront;
		sdf = sd->rf;
	} else if (!(inFront | outFront)) {
		*sv = sd->rLambBack;
		sdf = sd->rb;
	} else /* inFront ^ outFront */ {
		*sv = sd->tLamb;
		sdf = sd->tf;
	}
	sv->cieY *= 1./M_PI;
					/* add non-diffuse components */
	i = (sdf != NULL) ? sdf->ncomp : 0;
	while (i-- > 0) {
		nch = (*sdf->comp[i].func->getBSDFs)(coef, outVec, inVec,
							sdf->comp[i].dist);
		while (nch-- > 0) {
			c_cmix(&sv->spec, sv->cieY, &sv->spec,
					coef[nch], &sdf->comp[i].cspec[nch]);
			sv->cieY += coef[nch];
		}
	}
					/* make sure everything is set */
	c_ccvt(&sv->spec, C_CSXY+C_CSSPEC);
	return SDEnone;
}

/* Compute directional hemispherical scattering at this incident angle */
double
SDdirectHemi(const FVECT inVec, int sflags, const SDData *sd)
{
	double		hsum;
	SDSpectralDF	*rdf;
	const SDCDst	*cd;
	int		i;
					/* check arguments */
	if ((inVec == NULL) | (sd == NULL))
		return .0;
					/* gather diffuse components */
	if (inVec[2] > 0) {
		hsum = sd->rLambFront.cieY;
		rdf = sd->rf;
	} else /* !inFront */ {
		hsum = sd->rLambBack.cieY;
		rdf = sd->rb;
	}
	if ((sflags & SDsampDf+SDsampR) != SDsampDf+SDsampR)
		hsum = .0;
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT)
		hsum += sd->tLamb.cieY;
					/* gather non-diffuse components */
	i = ((sflags & SDsampSp+SDsampR) == SDsampSp+SDsampR &&
			rdf != NULL) ? rdf->ncomp : 0;
	while (i-- > 0) {		/* non-diffuse reflection */
		cd = (*rdf->comp[i].func->getCDist)(inVec, &rdf->comp[i]);
		if (cd != NULL)
			hsum += cd->cTotal;
	}
	i = ((sflags & SDsampSp+SDsampT) == SDsampSp+SDsampT &&
			sd->tf != NULL) ? sd->tf->ncomp : 0;
	while (i-- > 0) {		/* non-diffuse transmission */
		cd = (*sd->tf->comp[i].func->getCDist)(inVec, &sd->tf->comp[i]);
		if (cd != NULL)
			hsum += cd->cTotal;
	}
	return hsum;
}

/* Sample BSDF direction based on the given random variable */
SDError
SDsampBSDF(SDValue *sv, FVECT outVec, const FVECT inVec,
			double randX, int sflags, const SDData *sd)
{
	SDError		ec;
	int		inFront;
	SDSpectralDF	*rdf;
	double		rdiff;
	float		coef[SDmaxCh];
	int		i, j, n, nr;
	SDComponent	*sdc;
	const SDCDst	**cdarr = NULL;
					/* check arguments */
	if ((sv == NULL) | (outVec == NULL) | (inVec == NULL) | (sd == NULL) |
			(randX < 0) | (randX >= 1.))
		return SDEargument;
					/* whose side are we on? */
	inFront = (inVec[2] > 0);
					/* remember diffuse portions */
	if (inFront) {
		*sv = sd->rLambFront;
		rdf = sd->rf;
	} else /* !inFront */ {
		*sv = sd->rLambBack;
		rdf = sd->rb;
	}
	if ((sflags & SDsampDf+SDsampR) != SDsampDf+SDsampR)
		sv->cieY = .0;
	rdiff = sv->cieY;
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT)
		sv->cieY += sd->tLamb.cieY;
					/* gather non-diffuse components */
	i = nr = ((sflags & SDsampSp+SDsampR) == SDsampSp+SDsampR &&
			rdf != NULL) ? rdf->ncomp : 0;
	j = ((sflags & SDsampSp+SDsampT) == SDsampSp+SDsampT &&
			sd->tf != NULL) ? sd->tf->ncomp : 0;
	n = i + j;
	if (n > 0 && (cdarr = (const SDCDst **)malloc(n*sizeof(SDCDst *))) == NULL)
		return SDEmemory;
	while (j-- > 0) {		/* non-diffuse transmission */
		cdarr[i+j] = (*sd->tf->comp[j].func->getCDist)(inVec, &sd->tf->comp[j]);
		if (cdarr[i+j] == NULL) {
			free(cdarr);
			return SDEmemory;
		}
		sv->cieY += cdarr[i+j]->cTotal;
	}
	while (i-- > 0) {		/* non-diffuse reflection */
		cdarr[i] = (*rdf->comp[i].func->getCDist)(inVec, &rdf->comp[i]);
		if (cdarr[i] == NULL) {
			free(cdarr);
			return SDEmemory;
		}
		sv->cieY += cdarr[i]->cTotal;
	}
	if (sv->cieY <= 1e-7) {		/* anything to sample? */
		sv->cieY = .0;
		memset(outVec, 0, 3*sizeof(double));
		return SDEnone;
	}
					/* scale random variable */
	randX *= sv->cieY;
					/* diffuse reflection? */
	if (randX < rdiff) {
		SDdiffuseSamp(outVec, inFront, randX/rdiff);
		goto done;
	}
	randX -= rdiff;
					/* diffuse transmission? */
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT) {
		if (randX < sd->tLamb.cieY) {
			sv->spec = sd->tLamb.spec;
			SDdiffuseSamp(outVec, !inFront, randX/sd->tLamb.cieY);
			goto done;
		}
		randX -= sd->tLamb.cieY;
	}
					/* else one of cumulative dist. */
	for (i = 0; i < n && randX < cdarr[i]->cTotal; i++)
		randX -= cdarr[i]->cTotal;
	if (i >= n)
		return SDEinternal;
					/* compute sample direction */
	sdc = (i < nr) ? &rdf->comp[i] : &sd->tf->comp[i-nr];
	ec = (*sdc->func->sampCDist)(outVec, randX/cdarr[i]->cTotal, cdarr[i]);
	if (ec)
		return ec;
					/* compute color */
	j = (*sdc->func->getBSDFs)(coef, outVec, inVec, sdc->dist);
	if (j <= 0) {
		sprintf(SDerrorDetail, "BSDF \"%s\" sampling value error",
				sd->name);
		return SDEinternal;
	}
	sv->spec = sdc->cspec[0];
	rdiff = coef[0];
	while (--j) {
		c_cmix(&sv->spec, rdiff, &sv->spec, coef[j], &sdc->cspec[j]);
		rdiff += coef[j];
	}
done:
	if (cdarr != NULL)
		free(cdarr);
					/* make sure everything is set */
	c_ccvt(&sv->spec, C_CSXY+C_CSSPEC);
	return SDEnone;
}

/* Compute World->BSDF transform from surface normal and up (Y) vector */
SDError	
SDcompXform(RREAL vMtx[3][3], const FVECT sNrm, const FVECT uVec)
{
	if ((vMtx == NULL) | (sNrm == NULL) | (uVec == NULL))
		return SDEargument;
	VCOPY(vMtx[2], sNrm);
	if (normalize(vMtx[2]) == 0)
		return SDEargument;
	fcross(vMtx[0], uVec, vMtx[2]);
	if (normalize(vMtx[0]) == 0)
		return SDEargument;
	fcross(vMtx[1], vMtx[2], vMtx[0]);
	return SDEnone;
}

/* Compute inverse transform */
SDError	
SDinvXform(RREAL iMtx[3][3], RREAL vMtx[3][3])
{
	RREAL	mTmp[3][3];
	double	d;

	if ((iMtx == NULL) | (vMtx == NULL))
		return SDEargument;
					/* compute determinant */
	mTmp[0][0] = vMtx[2][2]*vMtx[1][1] - vMtx[2][1]*vMtx[1][2];
	mTmp[0][1] = vMtx[2][1]*vMtx[0][2] - vMtx[2][2]*vMtx[0][1];
	mTmp[0][2] = vMtx[1][2]*vMtx[0][1] - vMtx[1][1]*vMtx[0][2];
	d = vMtx[0][0]*mTmp[0][0] + vMtx[1][0]*mTmp[0][1] + vMtx[2][0]*mTmp[0][2];
	if (d == 0) {
		strcpy(SDerrorDetail, "Zero determinant in matrix inversion");
		return SDEargument;
	}
	d = 1./d;			/* invert matrix */
	mTmp[0][0] *= d; mTmp[0][1] *= d; mTmp[0][2] *= d;
	mTmp[1][0] = d*(vMtx[2][0]*vMtx[1][2] - vMtx[2][2]*vMtx[1][0]);
	mTmp[1][1] = d*(vMtx[2][2]*vMtx[0][0] - vMtx[2][0]*vMtx[0][2]);
	mTmp[1][2] = d*(vMtx[1][0]*vMtx[0][2] - vMtx[1][2]*vMtx[0][0]);
	mTmp[2][0] = d*(vMtx[2][1]*vMtx[1][0] - vMtx[2][0]*vMtx[1][1]);
	mTmp[2][1] = d*(vMtx[2][0]*vMtx[0][1] - vMtx[2][1]*vMtx[0][0]);
	mTmp[2][2] = d*(vMtx[1][1]*vMtx[0][0] - vMtx[1][0]*vMtx[0][1]);
	memcpy(iMtx, mTmp, sizeof(mTmp));
	return SDEnone;
}

/* Transform and normalize direction (column) vector */
SDError
SDmapDir(FVECT resVec, RREAL vMtx[3][3], const FVECT inpVec)
{
	FVECT	vTmp;

	if ((resVec == NULL) | (inpVec == NULL))
		return SDEargument;
	if (vMtx == NULL) {		/* assume they just want to normalize */
		if (resVec != inpVec)
			VCOPY(resVec, inpVec);
		return (normalize(resVec) > 0) ? SDEnone : SDEargument;
	}
	vTmp[0] = DOT(vMtx[0], inpVec);
	vTmp[1] = DOT(vMtx[1], inpVec);
	vTmp[2] = DOT(vMtx[2], inpVec);
	if (normalize(vTmp) == 0)
		return SDEargument;
	VCOPY(resVec, vTmp);
	return SDEnone;
}

/*################################################################*/
/*######### DEPRECATED ROUTINES AWAITING PERMANENT REMOVAL #######*/

/*
 * Routines for handling BSDF data
 */

#include "standard.h"
#include "paths.h"
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
	if (b != 0)
		a = a/b - 1.;
	return((a <= 1e-6) & (a >= -1e-6));
}

/* Returns the name of the given tag */
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

/* Returns the given tag's character content or empty string if none */
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
		sprintf(errmsg, "undefined RowAngleBasis '%s'", rbasis);
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
		if (dom <= 0) {
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
		if (dom <= 0) {
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
			if (f >= 0)
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
		if (strcasecmp(ezxml_txt(ezxml_child(wld,"Wavelength")),
				"Visible"))
			continue;
		for (wdb = ezxml_child(wld, "WavelengthDataBlock");
				wdb != NULL; wdb = wdb->next)
			if (!strcasecmp(ezxml_txt(ezxml_child(wdb,
					"WavelengthDataDirection")),
					"Transmission Front"))
				break;
		if (wdb != NULL) {		/* load front BTDF */
			load_bsdf_data(dp, wdb);
			break;			/* ignore the rest */
		}
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

/*######### END DEPRECATED ROUTINES #######*/
/*################################################################*/
