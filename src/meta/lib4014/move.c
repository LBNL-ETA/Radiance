#ifndef lint
static const char	RCSid[] = "$Id: move.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)move.c	4.1 (Berkeley) 6/27/83";
#endif

move(xi,yi){
	putch(035);
	cont(xi,yi);
}
