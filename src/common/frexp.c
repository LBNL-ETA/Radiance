#ifndef lint
static const char	RCSid[] = "$Id: frexp.c,v 1.4 2003/11/12 18:06:54 greg Exp $";
#endif

double
frexp(x, ip)		/* call it paranoia, I've seen the lib version */
double  x;
int  *ip;
{
	int  neg;
	register int  i;

	if (neg = (x < 0.0))
		x = -x;
	else if (x == 0.0) {
		*ip = 0;
		return(0.0);
	}
	if (x < 0.5)
		for (i = 0; x < 0.5; i--)
			x *= 2.0;
	else
		for (i = 0; x >= 1.0; i++)
			x *= 0.5;
	*ip = i;
	return (neg ? -x : x);
}
