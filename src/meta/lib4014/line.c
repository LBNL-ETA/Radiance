#ifndef lint
static const char	RCSid[] = "$Id: line.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)line.c	4.1 (Berkeley) 6/27/83";
#endif

line(x0,y0,x1,y1){
	move(x0,y0);
	cont(x1,y1);
}
