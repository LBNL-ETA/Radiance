#ifndef lint
static const char	RCSid[] = "$Id: fgetword.c,v 2.5 2003/07/27 22:12:01 schorsch Exp $";
#endif
/*
 * Read white space separated words from stream
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"

#include  <ctype.h>


char *
fgetword(s, n, fp)		/* get (quoted) word up to n-1 characters */
char  *s;
int  n;
register FILE  *fp;
{
	int  quote = '\0';
	register char  *cp;
	register int  c;
					/* skip initial white space */
	do
		c = getc(fp);
	while (isspace(c));
					/* check for quote */
	if ((c == '"') | (c == '\'')) {
		quote = c;
		c = getc(fp);
	}
					/* check for end of file */
	if (c == EOF)
		return(NULL);
					/* get actual word */
	cp = s;
	do {
		if (--n <= 0)		/* check length limit */
			break;
		*cp++ = c;
		c = getc(fp);
	} while (c != EOF && !(quote ? c==quote : isspace(c)));
	*cp = '\0';
	if ((c != EOF) & (!quote))	/* replace space */
		ungetc(c, fp);
	return(s);
}
