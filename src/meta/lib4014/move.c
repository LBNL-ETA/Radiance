#ifndef lint
static const char	RCSid[] = "$Id: move.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include "local4014.h"
#include "lib4014.h"

extern void
move(
	int xi,
	int yi
)
{
	putch(035);
	cont(xi,yi);
}
