#ifndef lint
static const char RCSid[] = "$Id: loadbsdf.c,v 3.3 2011/02/19 23:42:09 greg Exp $";
#endif
/*
 * Simple interface for loading BSDF, Radiance-specific search
 */

#include "rtio.h"
#include "rterror.h"
#include "bsdf.h"
#include "paths.h"

/* Convert error from BSDF library */
char *
transSDError(SDError ec)
{
	if (!SDerrorDetail[0])
		return(strcpy(errmsg, SDerrorEnglish[ec]));

	sprintf(errmsg, "%s: %s", SDerrorEnglish[ec], SDerrorDetail);
	return(errmsg);
}

/* Make sure we're not over 100% scattering for this component */
static void
checkDF(const char *nm, double amt, const SDSpectralDF *dfp, const char *desc)
{
	if (dfp != NULL)
		amt += dfp->maxHemi;
	if (amt <= 1.01)
		return;
	sprintf(errmsg, "BSDF \"%s\" has %.1f%% %s", nm, amt*100., desc);
	error(WARNING, errmsg);
}

/* Load a BSDF file and perform some basic checks */
SDData *
loadBSDF(char *fname)
{
	SDData	*sd;
	SDError	ec;
	char	*pname;

	sd = SDgetCache(fname);			/* look up or allocate */
	if (sd == NULL)
		error(SYSTEM, "out of memory in loadBSDF");
	if (SDisLoaded(sd))			/* already in memory? */
		return(sd);
						/* else find and load it */
	pname = getpath(fname, getrlibpath(), R_OK);
	if (pname == NULL) {
		sprintf(errmsg, "cannot find BSDF file \"%s\"", fname);
		error(USER, errmsg);
	}
	ec = SDloadFile(sd, pname);
	if (ec)
		error(USER, transSDError(ec));
						/* simple checks */
	checkDF(sd->name, sd->rLambFront.cieY, sd->rf, "front reflection");
	checkDF(sd->name, sd->rLambBack.cieY, sd->rb, "rear reflection");
	checkDF(sd->name, sd->tLamb.cieY, sd->tf, "transmission");
#if 0
fprintf(stderr, "Loaded BSDF '%s' (file \"%s\")\n", sd->name, pname);
fprintf(stderr, "Front diffuse reflectance: %.1f%%\n", sd->rLambFront.cieY*100.);
fprintf(stderr, "Back diffuse reflectance: %.1f%%\n", sd->rLambBack.cieY*100.);
fprintf(stderr, "Diffuse transmittance: %.1f%%\n", sd->tLamb.cieY*100.);
if (sd->rf)
fprintf(stderr, "Maximum direct hemispherical front reflection: %.1f%%\n",
sd->rf->maxHemi*100.);
if (sd->rb)
fprintf(stderr, "Maximum direct hemispherical back reflection: %.1f%%\n",
sd->rb->maxHemi*100.);
if (sd->tf)
fprintf(stderr, "Maximum direct hemispherical transmission: %.1f%%\n",
sd->tf->maxHemi*100.);
#endif
	SDretainSet = SDretainAll;		/* keep data in core */
	return(sd);
}
