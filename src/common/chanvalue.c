#ifndef lint
static const char	RCSid[] = "$Id: chanvalue.c,v 3.3 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Dummy definition of chanvalue() for calcomp routines.
 */

#include "rterror.h"
#include "calcomp.h"

double
chanvalue(int n)
{
	eputs("Call to unsupported chanvalue routine\n");
	quit(1);
	return 0.0; /* pro forma return */
}
