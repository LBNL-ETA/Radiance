#ifndef lint
static const char	RCSid[] = "$Id: fputword.c,v 3.2 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * Read white space separated words from stream
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdio.h>

#include <ctype.h>


void
fputword(s, fp)			/* put (quoted) word to file stream */
char  *s;
FILE  *fp;
{
	register char	*cp;
					/* check if quoting needed */
	for (cp = s; *cp; cp++)
		if (isspace(*cp)) {
			int	quote;
			if (index(s, '"'))
				quote = '\'';
			else
				quote = '"';
			fputc(quote, fp);
			fputs(s, fp);
			fputc(quote, fp);
			return;
		}
					/* output sans quotes */
	fputs(s, fp);
}
