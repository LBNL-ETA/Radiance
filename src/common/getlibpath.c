#ifndef lint
static const char	RCSid[] = "$Id: getlibpath.c,v 2.6 2013/12/19 16:38:12 greg Exp $";
#endif
/*
 * Return Radiance library search path
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include <stdio.h>

#include "rtio.h"
#include "paths.h"

char *
getrlibpath()
{
	static char	*libpath = NULL;

	if (libpath == NULL)
		if ((libpath = getenv(ULIBVAR)) == NULL)
			libpath = DEFPATH;

	return(libpath);
}
