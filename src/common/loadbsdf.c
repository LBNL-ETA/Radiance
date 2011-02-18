#ifndef lint
static const char RCSid[] = "$Id: loadbsdf.c,v 3.1 2011/02/18 02:41:55 greg Exp $";
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
	static char	mymess[128];

	if (!SDerrorDetail[0])
		return(strcpy(mymess, SDerrorEnglish[ec]));
	sprintf(mymess, "%s: %s", SDerrorEnglish[ec], SDerrorDetail);
	return(mymess);
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
loadBSDF(char *name)
{
	SDData	*sd = SDgetCache(name);
	SDError	ec;
	char	*pname;

	if (sd == NULL)
		error(SYSTEM, "out of memory in loadBSDF");
	if (SDisLoaded(sd))
		return(sd);
	
	pname = getpath(name, getrlibpath(), R_OK);
	if (pname == NULL) {
		sprintf(errmsg, "cannot find BSDF file \"%s\"", name);
		error(USER, errmsg);
	}
	ec = SDloadFile(sd, pname);
	if (ec)
		error(USER, transSDError(ec));
						/* simple checks */
	checkDF(name, sd->rLambFront.cieY, sd->rf, "front reflection");
	checkDF(name, sd->rLambBack.cieY, sd->rb, "rear reflection");
	checkDF(name, sd->tLamb.cieY, sd->tf, "transmission");

	SDretainSet = SDretainAll;		/* keep data in core */
	return(sd);
}
