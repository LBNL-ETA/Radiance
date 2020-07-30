#ifndef lint
static const char	RCSid[] = "$Id: fgetline.c,v 2.11 2020/07/30 17:25:01 greg Exp $";
#endif
/*
 * fgetline.c - read line with escaped newlines.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"


char *
fgetline(		/* read in line with escapes, elide final newline */
	char  *s,
	int  n,
	FILE  *fp
)
{
	char  *cp = s;
	int  c = EOF;

	while (--n > 0 && (c = getc(fp)) != EOF) {
		if (c == '\r' && (c = getc(fp)) != '\n') {
			ungetc(c, fp);		/* must be Apple file */
			c = '\n';
		}
		if (c == '\n' && (cp == s || cp[-1] != '\\'))
			break;
		*cp++ = c;
	}
	if ((cp == s) & (c == EOF))
		return(NULL);
	*cp = '\0';
#if defined(_WIN32) || defined(_WIN64)
	for (cp = s; *cp; cp++)			/* remove escaped newlines */
		if (cp[0] == '\\' && cp[1] == '\n') {
			*cp++ = ' '; *cp = ' ';
		}
#endif
	return(s);
}
