#ifndef lint
static const char	RCSid[] = "$Id: quit.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * Default program quit routine.
 */

#include "copyright.h"

void
quit(code)			/* quit program */
int  code;
{
	exit(code);
}
