#ifndef lint
static const char	RCSid[] = "$Id: badarg.c,v 2.4 2003/06/07 12:50:20 schorsch Exp $";
#endif
/*
 * Check argument list against format string.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <ctype.h>

#include "standard.h"


int
badarg(ac, av, fl)		/* check argument list */
int	ac;
register char	**av;
register char	*fl;
{
	register int  i;

	if (fl == NULL)
		fl = "";		/* no arguments? */
	for (i = 1; *fl; i++,av++,fl++) {
		if (i > ac || *av == NULL)
			return(-1);
		switch (*fl) {
		case 's':		/* string */
			if (**av == '\0' || isspace(**av))
				return(i);
			break;
		case 'i':		/* integer */
			if (!isintd(*av, " \t\r\n"))
				return(i);
			break;
		case 'f':		/* float */
			if (!isfltd(*av, " \t\r\n"))
				return(i);
			break;
		default:		/* bad call! */
			return(-1);
		}
	}
	return(0);		/* all's well */
}
