#ifndef lint
static const char	RCSid[] = "$Id: quit.c,v 2.4 2003/06/07 12:50:20 schorsch Exp $";
#endif
/*
 * Default program quit routine.
 */

#include "copyright.h"

#include <stdlib.h>

void
quit(code)			/* quit program */
int  code;
{
	exit(code);
}
