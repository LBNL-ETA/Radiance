#ifndef lint
static const char	RCSid[] = "$Id: wputs.c,v 3.4 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Default warning output function.
 */

#include "copyright.h"

int	nowarn = 0;		/* don't print warnings? */

void
wputs(s)
char	*s;
{
	if (!nowarn)
		eputs(s);
}
