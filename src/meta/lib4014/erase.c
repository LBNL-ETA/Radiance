#ifndef lint
static const char	RCSid[] = "$Id: erase.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)erase.c	4.1 (Berkeley) 6/27/83";
#endif

#include <stdio.h>

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
