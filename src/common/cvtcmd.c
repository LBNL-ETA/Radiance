#ifndef lint
static const char RCSid[] = "$Id: cvtcmd.c,v 3.3 2019/07/21 16:14:52 greg Exp $";
#endif
/*
 *  Convert a set of arguments into a command line for pipe() or system()
 *
 *  Defines in paths.h
 */

#include "paths.h"

/* Check if any of the characters in str2 are found in str1 */
int
matchany(const char *str1, const char *str2)
{
	while (*str1) {
		const char	*cp = str2;
		while (*cp)
			if (*cp++ == *str1)
				return(*str1);
		++str1;
	}
	return(0);
}

/* Convert a set of arguments into a command line, being careful of specials */
char *
convert_commandline(char *cmd, const int len, char *av[])
{
	char	*cp;

	for (cp = cmd; *av != NULL; av++) {
		const int	n = strlen(*av);
		if (cp+n >= cmd+(len-3))
			return(NULL);
		if (matchany(*av, SPECIALS)) {
			const int	quote =
#ifdef ALTQUOT
				strchr(*av, QUOTCHAR) ? ALTQUOT :
#endif
					QUOTCHAR;
			*cp++ = quote;
			strcpy(cp, *av);
			cp += n;
			*cp++ = quote;
		} else {
			strcpy(cp, *av);
			cp += n;
		}
		*cp++ = ' ';
	}
	if (cp <= cmd)
		return(NULL);
	*--cp = '\0';
	return(cmd);
}
