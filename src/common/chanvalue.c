#ifndef lint
static const char	RCSid[] = "$Id: chanvalue.c,v 3.2 2003/06/07 12:50:20 schorsch Exp $";
#endif
/*
 * Dummy definition of chanvalue() for calcomp routines.
 */

#include "standard.h"

double
chanvalue(int n)
{
	eputs("Call to unsupported chanvalue routine\n");
	quit(1);
	return 0.0; /* pro forma return */
}
