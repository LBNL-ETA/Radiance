#ifndef lint
static const char	RCSid[] = "$Id: erase.c,v 1.3 2003/11/15 02:13:37 schorsch Exp $";
#endif

#include <stdio.h>

#include "platform.h"
#include "local4014.h"
#include "lib4014.h"

extern void
erase(void)
{
	putch(033);
	putch(014);
	fflush(stdout);
	ohiy= -1;
	ohix = -1;
	oextra = -1;
	oloy = -1;
	sleep(2);
}
