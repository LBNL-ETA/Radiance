#ifndef lint
static const char	RCSid[] = "$Id: mqdraw.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *   Interface to MacIntosh QuickDraw routines
 *
 *      6/5/85
 */


#include  "meta.h"

#include  "macplot.h"


static short  curpat = -1;	/* current line drawing pattern */

static short  curpsiz = -1;	/* current pen size */

static short  curpmod = -1;	/* current pen mode */

static int  xpos = -1,		/* current position */
	    ypos = -1;



printstr(p)		/* output a string */

register PRIMITIVE  *p;

{
    char  *ctop(), *ptoc();

    MoveTo(mapx(p->xy[XMN]), mapy(p->xy[YMN]));
    DrawString(ctop(p->args));
    ptoc(p->args);
    xpos = -1;
    ypos = -1;

}





plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
    static short  right = FALSE;
    int  x1, x2, y1, y2;
    short  pp, ps;

    pp = (p->arg0 >> 4) & 03;
    if (p->arg0 & 0100 && pp != 0)
	pp += 03;
    if (pp != curpat) {
	PenPat(macpat[pp]);
	curpat = pp;
    }

    ps = WIDTH((p->arg0 >> 2) & 03);
    if (ps != curpsiz) {
	PenSize(CONV(ps, dxsize) + 1, CONV(ps, dysize) + 1);
	curpsiz = ps;
    }

    if (curpmod != patOr) {
	PenMode(patOr);
	curpmod = patOr;
    }

    x1 = mapx(p->xy[XMN]);
    x2 = mapx(p->xy[XMX]);

    if (p->arg0 & 0100) {
	y1 = mapy(p->xy[YMX]);
	y2 = mapy(p->xy[YMN]);
    } else {
	y1 = mapy(p->xy[YMN]);
	y2 = mapy(p->xy[YMX]);
    }

    if (x1 == xpos && y1 == ypos) {
	LineTo(x2, y2);
	xpos = x2;
	ypos = y2;
    } else if (x2 == xpos && y2 == ypos) {
	LineTo(x1, y1);
	xpos = x1;
	ypos = y1;
    } else if (right = !right) {
	MoveTo(x1, y1);
	LineTo(x2, y2);
	xpos = x2;
	ypos = y2;
    } else {
	MoveTo(x2, y2);
	LineTo(x1, y1);
	xpos = x1;
	ypos = y1;
    }

}




setfill(a0)			/* set filling mode */

register int  a0;

{
    short  pp, pm;

    pp = pati[(a0 >> 2) & 03];
    if (pp != curpat) {
	PenPat(macpat[pp]);
	curpat = pp;
    }

    if (a0 & 0100)
	pm = patXor;
    else
	pm = patOr;
    if (pm != curpmod) {
	PenMode(pm);
	curpmod = pm;
    }

}




fillrect(p)			/* fill a rectangle */

register PRIMITIVE  *p;

{
    static Rect  r;

    r.left = mapx(p->xy[XMN]);
    r.right = mapx(p->xy[XMX]);
    r.top = mapy(p->xy[YMX]);
    r.bottom = mapy(p->xy[YMN]);

    setfill(p->arg0);
    PaintRect(&r);

}



filltri(p)			/* fill a triangle */

register PRIMITIVE  *p;

{
    PolyHandle  polyh;
    int  x[4], y[4];
    int  skipv;
    register int  i;

    polyh = OpenPoly();

    x[0] = x[1] = mapx(p->xy[XMN]);
    x[2] = x[3] = mapx(p->xy[XMX]);
    y[1] = y[2] = mapy(p->xy[YMN]);
    y[0] = y[3] = mapy(p->xy[YMX]);

    skipv = (p->arg0 >> 4) & 03;
    if (skipv == 3)
	MoveTo(x[2], y[2]);
    else
	MoveTo(x[3], y[3]);

    for (i = 0; i < 4; i++)
	if (i != skipv)
	    LineTo(x[i], y[i]);
    
    ClosePoly();
    setfill(p->arg0);
    PaintPoly(polyh);
    KillPoly(polyh);
    xpos = -1;
    ypos = -1;

}




xform(xp, yp, p)		/* transform a point according to p */

register int  *xp, *yp;
register PRIMITIVE  *p;

{
    int  x, y;

    switch (p->arg0 & 060) {
	case 0:			/* right */
	    x = *xp;
	    y = *yp;
	    break;
	case 020:		/* up */
	    x = (XYSIZE-1) - *yp;
	    y = *xp;
	    break;
	case 040:		/* left */
	    x = (XYSIZE-1) - *xp;
	    y = (XYSIZE-1) - *yp;
	    break;
	case 060:		/* down */
	    x = *yp;
	    y = (XYSIZE-1) - *xp;
	    break;
    }
    
    *xp = CONV(x, p->xy[XMX] - p->xy[XMN]) + p->xy[XMN];
    *yp = CONV(y, p->xy[YMX] - p->xy[YMN]) + p->xy[YMN];

}



fillpoly(p)			/* fill a polygon */

register PRIMITIVE  *p;

{
    int  x0, y0, curx, cury;
    PolyHandle  polyh;
    char  *nextscan();
    register char  *s;

    polyh = OpenPoly();

    if ((s = nextscan(nextscan(p->args, "%d", &x0), "%d", &y0)) == NULL)
        error(USER, "illegal polygon spec in fillpoly");
    
    xform(&x0, &y0, p);
    x0 = mapx(x0); y0 = mapy(y0);
    MoveTo(x0, y0);
    
    while ((s = nextscan(nextscan(s, "%d", &curx), "%d", &cury)) != NULL) {
        xform(&curx, &cury, p);
        curx = mapx(curx); cury = mapy(cury);
        LineTo(curx, cury);
    }
    LineTo(x0, y0);
    
    ClosePoly();

    if (p->arg0 & 0100) {		/* draw border */
        if (curpat != 0) {
            PenPat(macpat[0]);
            curpat = 0;
        }
        if (curpsiz != 1) {
            PenSize(1, 1);
            curpsiz = 1;
        }
	if (curpmod != patOr) {
	    PenMode(patOr);
	    curpmod = patOr;
	}
	FramePoly(polyh);
    }
    
    setfill(p->arg0 & 077);
    PaintPoly(polyh);
    KillPoly(polyh);
    xpos = -1;
    ypos = -1;

}
