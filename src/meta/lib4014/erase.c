#ifndef lint
static const char	RCSid[] = "$Id: erase.c,v 1.2 2003/10/27 10:28:59 schorsch Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)erase.c	4.1 (Berkeley) 6/27/83";
#endif

#include <stdio.h>

#include "platform.h"

extern int ohiy;
extern int ohix;
extern int oloy;
extern int oextra;
erase(){
		putch(033);
		putch(014);
		fflush(stdout);
		ohiy= -1;
		ohix = -1;
		oextra = -1;
		oloy = -1;
		sleep(2);
}
