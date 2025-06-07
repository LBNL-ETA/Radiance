#ifndef lint
static const char	RCSid[] = "$Id: getlibpath.c,v 2.7 2025/06/07 05:09:45 greg Exp $";
#endif
/*
 * Return Radiance library search path
 *
 *  External symbols declared in paths.h
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
