#ifndef lint
static const char	RCSid[] = "$Id: circle.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)circle.c	4.1 (Berkeley) 6/27/83";
#endif

circle(x,y,r){
	arc(x,y,x+r,y,x+r,y);
}
