#ifndef lint
static const char	RCSid[] = "$Id: primout.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  Routines for primitive output
 *
 *      1/10/85
 */


#include  "meta.h"


FILE  *pout = NULL;		/* the primitive output stream */


void
plseg(		/* plot line segment */
	int	a0,
	int xstart,
	int ystart,
	int xend,
	int yend
)

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


void
pprim(	/* print primitive */
	int	co,
	int a0,
	int xmin,
	int ymin,
	int xmax,
	int ymax,
	char	*s
)

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



void
pglob(			/* print global */
	int  co,
	int  a0,
	char  *s
)
{
    PRIMITIVE  p;
    
    p.com = co;
    p.arg0 = a0;
    p.xy[XMN] = p.xy[YMN] = p.xy[XMX] = p.xy[YMX] = -1;
    p.args = s;

    writep(&p, pout);

}
