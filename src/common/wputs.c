#ifndef lint
static const char	RCSid[] = "$Id: wputs.c,v 3.5 2003/06/07 12:50:21 schorsch Exp $";
#endif
/*
 * Default warning output function.
 */

#include "copyright.h"

#include "standard.h"

int	nowarn = 0;		/* don't print warnings? */

void
wputs(s)
char	*s;
{
	if (!nowarn)
		eputs(s);
}
