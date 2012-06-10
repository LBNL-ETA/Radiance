#ifndef lint
static const char	RCSid[] = "$Id: badarg.c,v 2.7 2012/06/10 17:50:29 greg Exp $";
#endif
/*
 * Check argument list against format string.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include <ctype.h>

#include "rtio.h"


int
badarg(		/* check argument list */
int	ac,
char	**av,
char	*fl
)
{
	int	i;
	char	*s;

	if (fl == NULL)
		fl = "";		/* no arguments? */
	for (i = 1; *fl; i++,av++,fl++) {
		if (i > ac || *av == NULL)
			return(-1);
		s = *av;
		switch (*fl) {
		case 's':		/* string */
			while (isspace(*s))
				++s;
			if (*s == '\0')
				return(i);
			while (*s)
				if (!isascii(*s++))
					return(i);
			break;
		case 'i':		/* integer */
			if (!isintd(s, " \t\r\n"))
				return(i);
			break;
		case 'f':		/* float */
			if (!isfltd(s, " \t\r\n"))
				return(i);
			break;
		default:		/* bad call! */
			return(-1);
		}
	}
	return(0);		/* all's well */
}
