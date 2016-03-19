#ifndef lint
static const char	RCSid[] = "$Id: fixargv0.c,v 2.6 2016/03/19 00:44:46 schorsch Exp $";
#endif
/*
 * Fix argv[0] for DOS environments
 *
 *  External symbols declared in paths.h
 */

#include "copyright.h"

#include <ctype.h>
#include <string.h>

extern char *
fixargv0(av0)			/* extract command name from full path */
char  *av0;
{
	register char  *cp = av0;

	while (*cp) cp++;		/* start from end */
	while (cp-- > av0)
		switch (*cp) {		/* fix up command name */
		case '.':			/* remove extension */
			*cp = '\0';
			continue;
		case '\\':			/* remove directory */
			strcpy(av0, cp+1);
			break;
		default:			/* convert to lower case */
			*cp = tolower(*cp);
			continue;
		}
	return(av0);
}


