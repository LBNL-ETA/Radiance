#ifndef lint
static const char	RCSid[] = "$Id: biggerlib.c,v 3.5 2003/08/04 22:37:53 greg Exp $";
#endif
/*
 *  biggerlib.c - functions for an even bigger library.
 *
 *   10/2/86
 */

#include <stdio.h>
#include <math.h>

#include "calcomp.h"

double  argument(int);
static double  l_j0(char *), l_j1(char *), l_jn(char *);
static double  l_y0(char *), l_y1(char *), l_yn(char *);
static double  l_erf(char *), l_erfc(char *);


void
biggerlib()			/* expand the library */
{
				/* the Bessel functions */
    funset("j0", 1, ':', l_j0);
    funset("j1", 1, ':', l_j1);
    funset("jn", 2, ':', l_jn);
    funset("y0", 1, ':', l_y0);
    funset("y1", 1, ':', l_y1);
    funset("yn", 2, ':', l_yn);
    funset("erf", 1, ':', l_erf);
    funset("erfc", 1, ':', l_erfc);
}


static double
l_j0(char *nm)
{
    return(j0(argument(1)));
}


static double
l_j1(char *nm)
{
    return(j1(argument(1)));
}


static double
l_jn(char *nm)
{
    return(jn((int)(argument(1)+.5), argument(2)));
}


static double
l_y0(char *nm)
{
    return(y0(argument(1)));
}


static double
l_y1(char *nm)
{
    return(y1(argument(1)));
}


static double
l_yn(char *nm)
{
    return(yn((int)(argument(1)+.5), argument(2)));
}


static double
l_erf(char *nm)
{
    extern double  erf();

    return(erf(argument(1)));
}


static double
l_erfc(char *nm)
{
    extern double  erfc();

    return(erfc(argument(1)));
}
