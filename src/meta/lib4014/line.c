#ifndef lint
static const char	RCSid[] = "$Id: line.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include "local4014.h"
#include "lib4014.h"

extern void
line(
	int x0,
	int y0,
	int x1,
	int y1
)
{
	move(x0,y0);
	cont(x1,y1);
}
