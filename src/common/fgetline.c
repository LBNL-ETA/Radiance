#ifndef lint
static const char	RCSid[] = "$Id: fgetline.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * fgetline.c - read line with escaped newlines.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include  <stdio.h>


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
