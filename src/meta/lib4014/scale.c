#ifndef lint
static const char	RCSid[] = "$Id: scale.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include "local4014.h"
#include "lib4014.h"

extern void
scale(
	char i,
	float x,
	float y;
)
{
	switch(i) {
	default:
		return;
	case 'c':
		x *= 2.54;
		y *= 2.54;
	case 'i':
		x /= 200;
		y /= 200;
	case 'u':
		scalex = 1/x;
		scaley = 1/y;
	}
	scaleflag = 1;
}
