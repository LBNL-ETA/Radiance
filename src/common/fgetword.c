#ifndef lint
static const char	RCSid[] = "$Id: fgetword.c,v 2.9 2023/02/10 04:32:19 greg Exp $";
#endif
/*
 * Read white space separated words from stream
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"

#include  <ctype.h>


#define isquote(c)	(((c) == '"') | ((c) == '\''))


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
	if (isquote(c)) {
		quote = c;
		c = getc(fp);
	}
	cp = s;				/* get actual word */
	while (c != EOF) {
		if (c == quote)		/* end quote? */
			quote = '\0';
		else if (!quote && isquote(c))
			quote = c;	/* started new quote */
		else {
			if (!quote && isspace(c))
				break;	/* end of word */
			if (--n <= 0)
				break;	/* hit length limit */
			*cp++ = c;
		}
		c = getc(fp);		/* get next character */
	}
	*cp = '\0';
	if ((c == EOF) & (cp == s))	/* hit end-of-file? */
		return(NULL);
	return(s);
}
