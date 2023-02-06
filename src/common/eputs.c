#ifndef lint
static const char	RCSid[] = "$Id: eputs.c,v 2.5 2023/02/06 22:40:21 greg Exp $";
#endif
/*
 * Default error output function.
 */

#include "copyright.h"

#include <stdio.h>

#include "rterror.h"

void
eputs(const char *s)			/* error message */
{
	fputs(s, stderr);
}
