#ifndef lint
static const char	RCSid[] = "$Id: fputword.c,v 3.4 2003/06/27 06:53:21 greg Exp $";
#endif
/*
 * Read white space separated words from stream
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include "rtio.h"

#include <ctype.h>


void
fputword(s, fp)			/* put (quoted) word to file stream */
char  *s;
FILE  *fp;
{
	int		hasspace = 0;
	int		quote = 0;
	register char	*cp;
					/* check if quoting needed */
	for (cp = s; *cp; cp++)
		if (isspace(*cp))
			hasspace++;
		else if (*cp == '"')
			quote = '\'';
		else if (*cp == '\'')
			quote = '"';

	if (hasspace || quote) {	/* output with quotes */
		if (!quote) quote = '"';
		fputc(quote, fp);
		fputs(s, fp);
		fputc(quote, fp);
	} else				/* output sans quotes */
		fputs(s, fp);
}
