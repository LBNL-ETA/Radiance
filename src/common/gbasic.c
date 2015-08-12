#include "gbasic.h"


int		gb_epsorder(g3Float a,g3Float b)
{
	a = a - b;
	return gb_signum(a);
}

int		gb_epseq(g3Float a,g3Float b)
{
	return (fabs(a - b) < GB_EPSILON);
}

int		gb_inrange(g3Float a,g3Float rb,g3Float re)
{
	if (((a >= rb) && (a <= re)) || ((a <= rb) && (a >= re)))
		return 1;
	return 0;
}

int		gb_signum(g3Float a)
{
	return ((gb_epseq(a,0.0)) ? 0 : ((a < 0.0) ? -1 : 1));
}


g3Float	gb_max(g3Float a,g3Float b)
{
	return (a < b) ? b : a;
}

g3Float	gb_min(g3Float a,g3Float b)
{
	return (a < b) ? a : b;
}

int		gb_getroots(g3Float* r1,g3Float* r2,g3Float a,g3Float b,g3Float c)
{
	g3Float d;
	if (gb_epseq(a,0.0)) {
		if (gb_epseq(b,0.0))
			return 0;
		*r1 = *r2 = -c/b;
		return 1;
	}
	d = b*b - 4*a*c;
	if (d < -GB_EPSILON) 
		return 0;
	if (gb_epseq(d,0.0)) {
		*r1 = *r2 = 0.5*b/a;
		return 1;
	}
	d = sqrt(d);
	*r1 = 0.5*(b + d)/a;
	*r2 = 0.5*(b - d)/a;
	return 2;
}
