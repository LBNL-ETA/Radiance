#ifndef lint
static const char	RCSid[] = "$Id: primout.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  Routines for primitive output
 *
 *      1/10/85
 */


#include  "meta.h"


FILE  *pout = NULL;		/* the primitive output stream */


plseg(a0, xstart, ystart, xend, yend)		/* plot line segment */

int	a0, xstart, ystart, xend, yend;

{
    PRIMITIVE	p;
    int		reverse;

    if (xstart < xend) {
	p.xy[XMN] = xstart;
	p.xy[XMX] = xend;
	reverse = FALSE;
    } else {
	p.xy[XMN] = xend;
	p.xy[XMX] = xstart;
	reverse = TRUE;
    }

    if (ystart < yend) {
	p.xy[YMN] = ystart;
        p.xy[YMX] = yend;
    } else {
	p.xy[YMN] = yend;
	p.xy[YMX] = ystart;
	reverse = ystart > yend && !reverse;
    }

    p.com = PLSEG;
    p.arg0 = (reverse << 6) | a0;
    p.args = NULL;

    writep(&p, pout);

}






pprim(co, a0, xmin, ymin, xmax, ymax, s)	/* print primitive */

int	co, a0, xmin, ymin, xmax, ymax;
char	*s;

{
    PRIMITIVE	p;

    p.com = co;
    p.arg0 = a0;
    p.xy[XMN] = xmin;
    p.xy[YMN] = ymin;
    p.xy[XMX] = xmax;
    p.xy[YMX] = ymax;
    p.args = s;

    writep(&p, pout);

}




pglob(co, a0, s)			/* print global */

int  co, a0;
char  *s;

{
    PRIMITIVE  p;
    
    p.com = co;
    p.arg0 = a0;
    p.xy[XMN] = p.xy[YMN] = p.xy[XMX] = p.xy[YMX] = -1;
    p.args = s;

    writep(&p, pout);

}
