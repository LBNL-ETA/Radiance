#ifndef lint
static const char	RCSid[] = "$Id: space.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include "local4014.h"
#include "lib4014.h"

extern void
space(
	int x0,
	int y0,
	int x1,
	int y1
)
{
	botx = 0.;
	boty = 0.;
	obotx = x0;
	oboty = y0;
	if(scaleflag)
		return;
	scalex = 3120./(x1-x0);
	scaley = 3120./(y1-y0);
}
