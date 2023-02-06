#ifndef lint
static const char	RCSid[] = "$Id: wputs.c,v 3.7 2023/02/06 22:40:21 greg Exp $";
#endif
/*
 * Default warning output function.
 */

#include "copyright.h"

#include <stdio.h>

#include "rterror.h"

int	nowarn = 0;		/* don't print warnings? */

void
wputs(const char *s)
{
	if (!nowarn)
		eputs(s);
}
