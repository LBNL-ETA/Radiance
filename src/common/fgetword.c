#ifndef lint
static const char	RCSid[] = "$Id: fgetword.c,v 2.8 2017/08/10 19:10:44 greg Exp $";
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
fgetword(			/* get (quoted) word up to n-1 characters */
	char  *s,
	int  n,
	FILE  *fp
)
{
	int  quote = '\0';
	char  *cp;
	int  c;
					/* sanity checks */
	if ((s == NULL) | (n < 2))
		return(NULL);
					/* skip initial white space */
	do
		c = getc(fp);
	while (isspace(c));
					/* check for quote */
	if ((c == '"') | (c == '\'')) {
		quote = c;
		c = getc(fp);
	}
					/* get actual word */
	cp = s;
	while (c != EOF && !(quote ? c==quote : isspace(c))) {
		if (--n <= 0)		/* check length limit */
			break;
		*cp++ = c;
		c = getc(fp);
	}
	*cp = '\0';
	if (c == EOF) {			/* hit end-of-file? */
		if (cp == s)		/* and got nothing? */
			return(NULL);
	} else if (!quote)		/* replace white character */
		ungetc(c, fp);
	return(s);
}
