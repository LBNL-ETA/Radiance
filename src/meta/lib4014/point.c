#ifndef lint
static const char	RCSid[] = "$Id: point.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)point.c	4.1 (Berkeley) 6/27/83";
#endif

point(xi,yi){
	move(xi,yi);
	cont(xi,yi);
}
