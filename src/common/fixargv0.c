#ifndef lint
static const char	RCSid[] = "$Id: fixargv0.c,v 2.11 2025/06/07 05:09:45 greg Exp $";
#endif
/*
 * Fix argv[0] and assign global progname variable
 *
 *  External symbols declared in paths.h
 */

#include "rtio.h"
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


void
printargs(		/* print command arguments to a file */
	int  ac,
	char  **av,
	FILE  *fp
)
{
	if (ac <= 0) return;

	if (progname == NULL)
		fixargv0(av[0]);

	if (progname >= av[0] && progname - av[0] < strlen(av[0]))
		fputword(progname, fp);
	else
		fputword(av[0], fp);
	while (--ac > 0) {
		fputc(' ', fp);
		fputword(*++av, fp);
	}
	fputc('\n', fp);
}
