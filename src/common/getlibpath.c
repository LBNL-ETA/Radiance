#ifndef lint
static const char	RCSid[] = "$Id: getlibpath.c,v 2.4 2003/05/13 17:58:32 greg Exp $";
#endif
/*
 * Return Radiance library search path
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdio.h>

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
