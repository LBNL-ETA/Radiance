#ifndef lint
static const char	RCSid[] = "$Id: fgetline.c,v 2.4 2003/06/27 06:53:21 greg Exp $";
#endif
/*
 * fgetline.c - read line with escaped newlines.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"


char *
fgetline(s, n, fp)	/* read in line with escapes, elide final newline */
char  *s;
int  n;
register FILE  *fp;
{
	register char  *cp = s;
	register int  c = EOF;

	while (--n > 0 && (c = getc(fp)) != EOF) {
		if (c == '\r')
			continue;
		if (c == '\n' && (cp == s || cp[-1] != '\\'))
			break;
		*cp++ = c;
	}
	if (cp == s && c == EOF)
		return(NULL);
	*cp = '\0';
	return(s);
}
