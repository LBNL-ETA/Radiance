#ifndef lint
static const char	RCSid[] = "$Id: getlibpath.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
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
getlibpath()
{
	static char	*libpath = NULL;

	if (libpath == NULL)
		if ((libpath = getenv(ULIBVAR)) == NULL)
			libpath = DEFPATH;

	return(libpath);
}
