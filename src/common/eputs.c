#ifndef lint
static const char	RCSid[] = "$Id: eputs.c,v 2.4 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Default error output function.
 */

#include "copyright.h"

#include <stdio.h>

#include "rterror.h"

void
eputs(s)			/* error message */
char  *s;
{
	fputs(s, stderr);
}
