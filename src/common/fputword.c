#ifndef lint
static const char	RCSid[] = "$Id: fputword.c,v 3.11 2024/08/03 15:32:59 greg Exp $";
#endif
/*
 * Write word to stream, quoting as necessary
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include <ctype.h>

#include "rtio.h"

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
		else if ((cp > s) & (*cp == '"') && cp[1])
			quote = '\'';
		else if ((cp > s) & (*cp == '\'') && cp[1])
			quote = '"';

	if (!*s | hasspace | quote) {	/* output with quotes */
		if (!quote) quote = '"';
		fputc(quote, fp);
		fputs(s, fp);
		fputc(quote, fp);
	} else				/* output sans quotes */
		fputs(s, fp);
}
