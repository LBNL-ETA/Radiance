#ifndef lint
static const char	RCSid[] = "$Id: getlibpath.c,v 2.5 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Return Radiance library search path
 *
 *  External symbols declared in standard.h
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
