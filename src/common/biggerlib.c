#ifndef lint
static const char	RCSid[] = "$Id: biggerlib.c,v 3.1 2003/02/22 02:07:21 greg Exp $";
#endif
/*
 *  biggerlib.c - functions for an even bigger library.
 *
 *   10/2/86
 */

#include <math.h>

double  argument();
static double  l_j0(), l_j1(), l_jn(), l_y0(), l_y1(), l_yn();
static double  l_erf(), l_erfc();


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
l_j0()
{
    return(j0(argument(1)));
}


static double
l_j1()
{
    return(j1(argument(1)));
}


static double
l_jn()
{
    return(jn((int)(argument(1)+.5), argument(2)));
}


static double
l_y0()
{
    return(y0(argument(1)));
}


static double
l_y1()
{
    return(y1(argument(1)));
}


static double
l_yn()
{
    return(yn((int)(argument(1)+.5), argument(2)));
}


static double
l_erf()
{
    extern double  erf();

    return(erf(argument(1)));
}


static double
l_erfc()
{
    extern double  erfc();

    return(erfc(argument(1)));
}
