#ifndef lint
static const char	RCSid[] = "$Id: quit.c,v 2.5 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Default program quit routine.
 */

#include "copyright.h"

#include <stdlib.h>

#include "rterror.h"

void
quit(code)			/* quit program */
int  code;
{
	exit(code);
}
