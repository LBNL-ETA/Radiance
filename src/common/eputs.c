#ifndef lint
static const char	RCSid[] = "$Id: eputs.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * Default error output function.
 */

#include "copyright.h"

#include <stdio.h>

void
eputs(s)			/* error message */
char  *s;
{
	fputs(s, stderr);
}
