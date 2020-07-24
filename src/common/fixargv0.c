#ifndef lint
static const char	RCSid[] = "$Id: fixargv0.c,v 2.9 2020/07/24 16:58:16 greg Exp $";
#endif
/*
 * Fix argv[0] for DOS environments
 *
 *  External symbols declared in paths.h
 */

#include <ctype.h>
#include <string.h>

char *
fixargv0(char *av0)		/* extract command name from full path */
{
	char  *cp = av0, *end;

	while (*cp) cp++;		/* start from end */
	end = cp;
	while (cp-- > av0)
		switch (*cp) {		/* fix up command name */
		case '.':			/* remove extension */
			*cp = '\0';
			end = cp;
			continue;
		case '\\':			/* remove directory */
		case '/':
			/* make sure the original pointer remains the same */
			memmove(av0, cp+1, end-cp);
			return(av0);
		default:			/* convert to lower case */
			*cp = tolower(*cp);
			continue;
		}
	return(av0);
}
