#ifndef lint
static const char RCSid[] = "$Id: loadbsdf.c,v 3.13 2021/03/27 17:50:18 greg Exp $";
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
		return(strcpy(errmsg, SDerrorList[ec]));

	sprintf(errmsg, "%s: %s", SDerrorList[ec], SDerrorDetail);
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
		error(SYSTEM, errmsg);
	}
	ec = SDloadFile(sd, pname);
	if (ec)
		error(USER, transSDError(ec));
						/* simple checks */
	checkDF(sd->name, sd->rLambFront.cieY, sd->rf, "front reflection");
	checkDF(sd->name, sd->rLambBack.cieY, sd->rb, "rear reflection");
	checkDF(sd->name, sd->tLambFront.cieY, sd->tf, "front transmission");
	checkDF(sd->name, sd->tLambBack.cieY, sd->tb, "back transmission");
#ifdef DEBUG
{
float	rgb[3];
fprintf(stderr, "Loaded BSDF '%s' (file \"%s\")\n", sd->name, pname);
ccy2rgb(&sd->rLambFront.spec, sd->rLambFront.cieY, rgb);
fprintf(stderr, "Front diffuse RGB: %.4f %.4f %.4f\n", rgb[0], rgb[1], rgb[2]);
ccy2rgb(&sd->rLambBack.spec, sd->rLambBack.cieY, rgb);
fprintf(stderr, "Back diffuse RGB: %.4f %.4f %.4f\n", rgb[0], rgb[1], rgb[2]);
ccy2rgb(&sd->tLambFront.spec, sd->tLamb.cieY, rgb);
fprintf(stderr, "Front diffuse RGB transmittance: %.4f %.4f %.4f\n", rgb[0], rgb[1], rgb[2]);
ccy2rgb(&sd->tLambBack.spec, sd->tLamb.cieY, rgb);
fprintf(stderr, "Back diffuse RGB transmittance: %.4f %.4f %.4f\n", rgb[0], rgb[1], rgb[2]);
if (sd->rf)
fprintf(stderr, "Maximum direct hemispherical front reflection: %.3f%%\n",
sd->rf->maxHemi*100.);
if (sd->rb)
fprintf(stderr, "Maximum direct hemispherical back reflection: %.3f%%\n",
sd->rb->maxHemi*100.);
if (sd->tf)
fprintf(stderr, "Maximum direct hemispherical front transmission: %.3f%%\n",
sd->tf->maxHemi*100.);
if (sd->tb)
fprintf(stderr, "Maximum direct hemispherical back transmission: %.3f%%\n",
sd->tb->maxHemi*100.);
}
#endif
	SDretainSet = SDretainAll;		/* keep data in core */
#ifdef  SMLMEM
	SDmaxCache = 5L*1024*1024;
#else
	SDmaxCache = 250L*1024*1024;
#endif
	return(sd);
}
