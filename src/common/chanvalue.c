#ifndef lint
static const char	RCSid[] = "$Id: chanvalue.c,v 3.1 2003/02/22 02:07:21 greg Exp $";
#endif
/*
 * Dummy definition of chanvalue() for calcomp routines.
 */

double
chanvalue(int n)
{
	eputs("Call to unsupported chanvalue routine\n");
	quit(1);
}
