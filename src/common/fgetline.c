#ifndef lint
static const char	RCSid[] = "$Id: fgetline.c,v 2.8 2014/07/09 22:46:51 greg Exp $";
#endif
/*
 * fgetline.c - read line with escaped newlines.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"

#ifdef getc_unlocked		/* avoid horrendous overhead of flockfile */
#undef getc
#define getc    getc_unlocked
#endif


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
	return(s);
}
