#ifndef lint
static const char RCSid[] = "$Id: bsdf.c,v 2.60 2021/04/28 00:59:10 greg Exp $";
#endif
/*
 *  bsdf.c
 *  
 *  Definitions for bidirectional scattering distribution functions.
 *
 *  Created by Greg Ward on 1/10/11.
 *
 */

#define	_USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
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

/* Pointer to error list in preferred language */
const char		**SDerrorList = SDerrorEnglish;

/* Additional information on last error (ASCII English) */
char			SDerrorDetail[256];

/* Empty distribution for getCDist() calls that fail for some reason */
const SDCDst		SDemptyCD;

/* Cache of loaded BSDFs */
struct SDCache_s	*SDcacheList = NULL;

/* Retain BSDFs in cache list? */
int			SDretainSet = SDretainNone;

/* Maximum cache size for any given BSDF? */
unsigned long		SDmaxCache = 0;		/* 0 == unlimited */

/* Report any error to the indicated stream */
SDError
SDreportError(SDError ec, FILE *fp)
{
	if (!ec)
		return SDEnone;
	if ((ec < SDEnone) | (ec > SDEunknown)) {
		SDerrorDetail[0] = '\0';
		ec = SDEunknown;
	}
	if (fp == NULL)
		return ec;
	fputs(SDerrorList[ec], fp);
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
SDloadGeometry(SDData *sd, ezxml_t wtl)
{
	ezxml_t		node, matl, geom;
	double		cfact;
	const char	*fmt = NULL, *mgfstr;

	SDerrorDetail[0] = '\0';
	sd->matn[0] = '\0'; sd->makr[0] = '\0';
	sd->dim[0] = sd->dim[1] = sd->dim[2] = 0;
	matl = ezxml_child(wtl, "Material");
	if (matl != NULL) {			/* get material info. */
		if ((node = ezxml_child(matl, "Name")) != NULL) {
			strncpy(sd->matn, ezxml_txt(node), SDnameLn);
			if (sd->matn[SDnameLn-1])
				strcpy(sd->matn+(SDnameLn-4), "...");
		}
		if ((node = ezxml_child(matl, "Manufacturer")) != NULL) {
			strncpy(sd->makr, ezxml_txt(node), SDnameLn);
			if (sd->makr[SDnameLn-1])
				strcpy(sd->makr+(SDnameLn-4), "...");
		}
		if ((node = ezxml_child(matl, "Width")) != NULL)
			sd->dim[0] = atof(ezxml_txt(node)) *
					to_meters(ezxml_attr(node, "unit"));
		if ((node = ezxml_child(matl, "Height")) != NULL)
			sd->dim[1] = atof(ezxml_txt(node)) *
					to_meters(ezxml_attr(node, "unit"));
		if ((node = ezxml_child(matl, "Thickness")) != NULL)
			sd->dim[2] = atof(ezxml_txt(node)) *
					to_meters(ezxml_attr(node, "unit"));
		if ((sd->dim[0] < 0) | (sd->dim[1] < 0) | (sd->dim[2] < 0)) {
			if (!SDerrorDetail[0])
				sprintf(SDerrorDetail, "Negative dimension in \"%s\"",
							sd->name);
			return SDEdata;
		}
	}
	sd->mgf = NULL;
	geom = ezxml_child(wtl, "Geometry");
	if (geom == NULL)			/* no actual geometry? */
		return SDEnone;
	fmt = ezxml_attr(geom, "format");
	if (fmt != NULL && strcasecmp(fmt, "MGF")) {
		sprintf(SDerrorDetail,
			"Unrecognized geometry format '%s' in \"%s\"",
					fmt, sd->name);
		return SDEsupport;
	}
	if ((node = ezxml_child(geom, "MGFblock")) == NULL ||
			(mgfstr = ezxml_txt(node)) == NULL)
		return SDEnone;
	while (isspace(*mgfstr))
		++mgfstr;
	if (!*mgfstr)
		return SDEnone;
	cfact = to_meters(ezxml_attr(node, "unit"));
	if (cfact <= 0)
		return SDEformat;
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
	wtl = ezxml_child(fl, "FileType");
	if (wtl != NULL && strcmp(ezxml_txt(wtl), "BSDF")) {
		sprintf(SDerrorDetail,
			"XML \"%s\": wrong FileType (must be 'BSDF')",
				sd->name);
		ezxml_free(fl);
		return SDEformat;
	}
	wtl = ezxml_child(ezxml_child(fl, "Optical"), "Layer");
	if (wtl == NULL) {
		sprintf(SDerrorDetail, "BSDF \"%s\": no optical layers",
				sd->name);
		ezxml_free(fl);
		return SDEformat;
	}
				/* load geometry if present */
	lastErr = SDloadGeometry(sd, wtl);
	if (lastErr) {
		ezxml_free(fl);
		return lastErr;
	}
				/* try loading variable resolution data */
	lastErr = SDloadTre(sd, wtl);
				/* check our result */
	if (lastErr == SDEsupport)	/* try matrix BSDF if not tree data */
		lastErr = SDloadMtx(sd, wtl);
		
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
	if (sd->tb != NULL && sd->tb->maxHemi <= .001) {
		SDfreeSpectralDF(sd->tb); sd->tb = NULL;
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

/* Add component(s) to spectral distribution function */
SDSpectralDF *
SDaddComponent(SDSpectralDF *odf, int nadd)
{
	SDSpectralDF	*df;

	if (odf == NULL)
		return SDnewSpectralDF(nadd);
	if (nadd <= 0)
		return odf;
	df = (SDSpectralDF *)realloc(odf, sizeof(SDSpectralDF) +
				(odf->ncomp+nadd-1)*sizeof(SDComponent));
	if (df == NULL) {
		sprintf(SDerrorDetail,
			"Cannot add %d component(s) to spectral DF", nadd);
		SDfreeSpectralDF(odf);
		return NULL;
	}
	memset(df->comp+df->ncomp, 0, nadd*sizeof(SDComponent));
	df->ncomp += nadd;
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
		if (df->comp[n].dist != NULL)
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
		else if (*cp == '/')
			dot = NULL;
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
	if (sd->tb != NULL) {
		SDfreeSpectralDF(sd->tb);
		sd->tb = NULL;
	}
	sd->rLambFront.cieY = .0;
	sd->rLambFront.spec.flags = 0;
	sd->rLambBack.cieY = .0;
	sd->rLambBack.spec.flags = 0;
	sd->tLambFront.cieY = .0;
	sd->tLambFront.spec.flags = 0;
	sd->tLambBack.cieY = .0;
	sd->tLambBack.spec.flags = 0;
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
	/* PLACE MUTEX LOCK HERE FOR THREAD-SAFE */
	if ((sd = SDgetCache(fname)) == NULL) {
		SDreportError(SDEmemory, stderr);
		return NULL;
	}
	if (!SDisLoaded(sd) && (ec = SDloadFile(sd, fname))) {
		SDreportError(ec, stderr);
		SDfreeCache(sd);
		sd = NULL;
	}
	/* END MUTEX LOCK */
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
		SDfreeCumulativeCache(sd->tb);
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
SDsampComponent(SDValue *sv, FVECT ioVec, double randX, SDComponent *sdc)
{
	float		coef[SDmaxCh];
	SDError		ec;
	FVECT		inVec;
	const SDCDst	*cd;
	double		d;
	int		n;
					/* check arguments */
	if ((sv == NULL) | (ioVec == NULL) | (sdc == NULL))
		return SDEargument;
					/* get cumulative distribution */
	VCOPY(inVec, ioVec);
	sv->cieY = 0;
	cd = (*sdc->func->getCDist)(inVec, sdc);
	if (cd != NULL)
		sv->cieY = cd->cTotal;
	if (sv->cieY <= 1e-6) {		/* nothing to sample? */
		sv->spec = c_dfcolor;
		memset(ioVec, 0, sizeof(FVECT));
		return SDEnone;
	}
					/* compute sample direction */
	ec = (*sdc->func->sampCDist)(ioVec, randX, cd);
	if (ec)
		return ec;
					/* get BSDF color */
	n = (*sdc->func->getBSDFs)(coef, ioVec, inVec, sdc);
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
	c_ccvt(&sv->spec, C_CSXY);	/* make sure (x,y) is set */
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

	if (n <= 0)			/* check corner cases */
		return;
	if (randX < 0) randX = 0;
	else if (randX >= 1.) randX = 0.999999999999999;
	if (n == 1) {
		t[0] = randX;
		return;
	}
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
	outVec[2] = sqrt(outVec[2]*(outVec[2]>0));
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
	if ((projSA == NULL) | (v1 == NULL) | (sd == NULL))
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
	if (v1[2] > 0) {		/* front surface query? */
		rdf = sd->rf;
		tdf = (sd->tf != NULL) ? sd->tf : sd->tb;
	} else {
		rdf = sd->rb;
		tdf = (sd->tb != NULL) ? sd->tb : sd->tf;
	}
	if (v2 != NULL) {		/* bidirectional? */
		if (v1[2] > 0 ^ v2[2] > 0)
			rdf = NULL;
		else
			tdf = NULL;
	}
	ec = SDEdata;			/* run through components */
	for (i = (rdf==NULL) ? 0 : rdf->ncomp; i--; ) {
		ec = (*rdf->comp[i].func->queryProjSA)(projSA, v1, v2,
						qflags, &rdf->comp[i]);
		if (ec)
			return ec;
	}
	for (i = (tdf==NULL) ? 0 : tdf->ncomp; i--; ) {
		ec = (*tdf->comp[i].func->queryProjSA)(projSA, v1, v2,
						qflags, &tdf->comp[i]);
		if (ec)
			return ec;
	}
	if (ec) {			/* all diffuse? */
		projSA[0] = M_PI;
		if (qflags == SDqueryMin+SDqueryMax)
			projSA[1] = M_PI;
	} else if (qflags == SDqueryMin+SDqueryMax && projSA[0] > projSA[1])
		projSA[0] = projSA[1];
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
	} else if (inFront) {
		*sv = sd->tLambFront;
		sdf = (sd->tf != NULL) ? sd->tf : sd->tb;
	} else /* outFront & !inFront */ {
		*sv = sd->tLambBack;
		sdf = (sd->tb != NULL) ? sd->tb : sd->tf;
	}
	sv->cieY *= 1./M_PI;
					/* add non-diffuse components */
	i = (sdf != NULL) ? sdf->ncomp : 0;
	while (i-- > 0) {
		nch = (*sdf->comp[i].func->getBSDFs)(coef, outVec, inVec,
							&sdf->comp[i]);
		while (nch-- > 0) {
			c_cmix(&sv->spec, sv->cieY, &sv->spec,
					coef[nch], &sdf->comp[i].cspec[nch]);
			sv->cieY += coef[nch];
		}
	}
	c_ccvt(&sv->spec, C_CSXY);	/* make sure (x,y) is set */
	return SDEnone;
}

/* Compute directional hemispherical scattering at this incident angle */
double
SDdirectHemi(const FVECT inVec, int sflags, const SDData *sd)
{
	double		hsum;
	SDSpectralDF	*rdf, *tdf;
	const SDCDst	*cd;
	int		i;
					/* check arguments */
	if ((inVec == NULL) | (sd == NULL))
		return .0;
					/* gather diffuse components */
	if (inVec[2] > 0) {
		hsum = sd->rLambFront.cieY;
		rdf = sd->rf;
		tdf = (sd->tf != NULL) ? sd->tf : sd->tb;
	} else /* !inFront */ {
		hsum = sd->rLambBack.cieY;
		rdf = sd->rb;
		tdf = (sd->tb != NULL) ? sd->tb : sd->tf;
	}
	if ((sflags & SDsampDf+SDsampR) != SDsampDf+SDsampR)
		hsum = .0;
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT)
		hsum += (inVec[2] > 0) ?
				sd->tLambFront.cieY : sd->tLambBack.cieY;
					/* gather non-diffuse components */
	i = (((sflags & SDsampSp+SDsampR) == SDsampSp+SDsampR) &
			(rdf != NULL)) ? rdf->ncomp : 0;
	while (i-- > 0) {		/* non-diffuse reflection */
		cd = (*rdf->comp[i].func->getCDist)(inVec, &rdf->comp[i]);
		if (cd != NULL)
			hsum += cd->cTotal;
	}
	i = (((sflags & SDsampSp+SDsampT) == SDsampSp+SDsampT) &
			(tdf != NULL)) ? tdf->ncomp : 0;
	while (i-- > 0) {		/* non-diffuse transmission */
		cd = (*tdf->comp[i].func->getCDist)(inVec, &tdf->comp[i]);
		if (cd != NULL)
			hsum += cd->cTotal;
	}
	return hsum;
}

/* Sample BSDF direction based on the given random variable */
SDError
SDsampBSDF(SDValue *sv, FVECT ioVec, double randX, int sflags, const SDData *sd)
{
	SDError		ec;
	FVECT		inVec;
	int		inFront;
	SDSpectralDF	*rdf, *tdf;
	double		rdiff;
	float		coef[SDmaxCh];
	int		i, j, n, nr;
	SDComponent	*sdc;
	const SDCDst	**cdarr = NULL;
					/* check arguments */
	if ((sv == NULL) | (ioVec == NULL) | (sd == NULL) |
			(randX < 0) | (randX >= 1.))
		return SDEargument;
					/* whose side are we on? */
	VCOPY(inVec, ioVec);
	inFront = (inVec[2] > 0);
					/* remember diffuse portions */
	if (inFront) {
		*sv = sd->rLambFront;
		rdf = sd->rf;
		tdf = (sd->tf != NULL) ? sd->tf : sd->tb;
	} else /* !inFront */ {
		*sv = sd->rLambBack;
		rdf = sd->rb;
		tdf = (sd->tb != NULL) ? sd->tb : sd->tf;
	}
	if ((sflags & SDsampDf+SDsampR) != SDsampDf+SDsampR)
		sv->cieY = .0;
	rdiff = sv->cieY;
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT)
		sv->cieY += inFront ? sd->tLambFront.cieY : sd->tLambBack.cieY;
					/* gather non-diffuse components */
	i = nr = (((sflags & SDsampSp+SDsampR) == SDsampSp+SDsampR) &
			(rdf != NULL)) ? rdf->ncomp : 0;
	j = (((sflags & SDsampSp+SDsampT) == SDsampSp+SDsampT) &
			(tdf != NULL)) ? tdf->ncomp : 0;
	n = i + j;
	if (n > 0 && (cdarr = (const SDCDst **)malloc(n*sizeof(SDCDst *))) == NULL)
		return SDEmemory;
	while (j-- > 0) {		/* non-diffuse transmission */
		cdarr[i+j] = (*tdf->comp[j].func->getCDist)(inVec, &tdf->comp[j]);
		if (cdarr[i+j] == NULL)
			cdarr[i+j] = &SDemptyCD;
		sv->cieY += cdarr[i+j]->cTotal;
	}
	while (i-- > 0) {		/* non-diffuse reflection */
		cdarr[i] = (*rdf->comp[i].func->getCDist)(inVec, &rdf->comp[i]);
		if (cdarr[i] == NULL)
			cdarr[i] = &SDemptyCD;
		sv->cieY += cdarr[i]->cTotal;
	}
	if (sv->cieY <= 1e-6) {		/* anything to sample? */
		sv->cieY = .0;
		memset(ioVec, 0, sizeof(FVECT));
		return SDEnone;
	}
					/* scale random variable */
	randX *= sv->cieY;
					/* diffuse reflection? */
	if (randX < rdiff) {
		SDdiffuseSamp(ioVec, inFront, randX/rdiff);
		goto done;
	}
	randX -= rdiff;
					/* diffuse transmission? */
	if ((sflags & SDsampDf+SDsampT) == SDsampDf+SDsampT) {
		const SDValue	*sdt = inFront ? &sd->tLambFront : &sd->tLambBack;
		if (randX < sdt->cieY) {
			sv->spec = sdt->spec;
			SDdiffuseSamp(ioVec, !inFront, randX/sdt->cieY);
			goto done;
		}
		randX -= sdt->cieY;
	}
					/* else one of cumulative dist. */
	for (i = 0; i < n && randX >= cdarr[i]->cTotal; i++)
		randX -= cdarr[i]->cTotal;
	if (i >= n)
		return SDEinternal;
					/* compute sample direction */
	sdc = (i < nr) ? &rdf->comp[i] : &tdf->comp[i-nr];
	ec = (*sdc->func->sampCDist)(ioVec, randX/cdarr[i]->cTotal, cdarr[i]);
	if (ec)
		return ec;
					/* compute color */
	j = (*sdc->func->getBSDFs)(coef, ioVec, inVec, sdc);
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
	c_ccvt(&sv->spec, C_CSXY);	/* make sure (x,y) is set */
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
