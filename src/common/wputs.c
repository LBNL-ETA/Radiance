#ifndef lint
static const char	RCSid[] = "$Id: wputs.c,v 3.6 2003/07/30 10:11:06 schorsch Exp $";
#endif
/*
 * Default warning output function.
 */

#include "copyright.h"

#include "standard.h"

int	nowarn = 0;		/* don't print warnings? */

void
wputs(char *s)
{
	if (!nowarn)
		eputs(s);
}
