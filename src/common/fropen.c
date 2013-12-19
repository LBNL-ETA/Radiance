#ifndef lint
static const char	RCSid[] = "$Id: fropen.c,v 2.16 2013/12/19 16:38:12 greg Exp $";
#endif
/*
 * Find and open a Radiance library file.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include <stdio.h>

#include "rtio.h"
#include "paths.h"


FILE *
frlibopen(			/* find file and open for reading */
	char  *fname
)
{
	FILE  *fp;
	char  pname[PATH_MAX];
	char  *sp, *cp;

	if (fname == NULL)
		return(NULL);

	if (ISABS(fname) || fname[0] == '.')	/* absolute path */
		return(fopen(fname, "r"));
						/* check search path */
	sp = getrlibpath();
	do {
		cp = pname;
		while (*sp && (*cp = *sp++) != PATHSEP)
			cp++;
		if (cp > pname && !ISDIRSEP(cp[-1]))
			*cp++ = DIRSEP;
		strcpy(cp, fname);
		if ((fp = fopen(pname, "r")) != NULL)
			return(fp);			/* got it! */
	} while (*sp);
						/* not found */
	return(NULL);
}
