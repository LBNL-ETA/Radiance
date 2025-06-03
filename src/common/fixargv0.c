#ifndef lint
static const char	RCSid[] = "$Id: fixargv0.c,v 2.10 2025/06/03 21:31:51 greg Exp $";
#endif
/*
 * Fix argv[0] and assign global progname variable
 *
 *  External symbols declared in paths.h
 */

#include "paths.h"
#include <ctype.h>

char	*progname = NULL;	/* global argv[0] */

char *
fixargv0(char *av0)		/* extract command name from full path */
{
	char  *cp = av0;

	while (*cp) cp++;		/* start from end */

	while (cp-- > av0)
		switch (*cp) {
		CASEDIRSEP:			/* remove directory */
			av0 = cp+1;
			break;
#if defined(_WIN32) || defined(_WIN64)	/* only do for Windows: */
		case '.':			/* remove extension */
			*cp = '\0';
			continue;
		default:			/* convert to lower case */
			*cp = tolower(*cp);
			continue;
#endif
		}
	return(progname = av0);
}
