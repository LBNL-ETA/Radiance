#ifndef lint
static const char	RCSid[] = "$Id: point.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include "local4014.h"
#include "lib4014.h"

extern void
point(
	int xi,
	int yi
)
{
	move(xi,yi);
	cont(xi,yi);
}
