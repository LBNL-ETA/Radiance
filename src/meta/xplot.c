#ifndef lint
static const char	RCSid[] = "$Id: xplot.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *   X window plotting functions
 *
 *      2/26/86
 */


#include  "meta.h"

#include  "plot.h"

#undef  TRUE

#undef  FALSE

#include  <X/Xlib.h>


#define  BORWIDTH  5

#define  mapx(x)  CONV(x,dxsize)
#define  mapy(y)  CONV((XYSIZE-1)-(y),dysize)

#define  fill(xmin,ymin,xmax,ymax,pm) \
		XTileSet(wind,xmin,ymin,(xmax)-(xmin)+1,(ymax)-(ymin)+1,pm)

#define  MAXVERT  128

#define  MAXCOLOR	65535

				/* X window fonts */
static struct {
	char  *name;		/* font name */
	Font  id;		/* font id */
} font[] = {
	{"9x15", 0},
	{"8x13", 0},
	{"6x13", 0},
	{"vtdwidth", 0},
	{"accordbfx", 0}
};

					/* map to our matrix string types */
static int  fontmap[16] = {0,3,0,3,1,3,1,3,2,4,2,4,2,4,2,4};

				/* colors */
int  pixel[4];

int  dxsize, dysize;		/* window size */

Window  wind;			/* our window */

extern Pixmap  pixpat();



init(name, geom)		/* initialize window */
char  *name;
char  *geom;
{
    WindowInfo  Info;
    OpaqueFrame  mainframe;
    char  defgeom[32];
    Color  cdef[4];
    int  dummy;

    if (XOpenDisplay("") == NULL)
	error(SYSTEM, "can't open display");
    mainframe.bdrwidth = BORWIDTH;
    mainframe.border = BlackPixmap;
    mainframe.background = WhitePixmap;
    dxsize = DisplayWidth() - 2*BORWIDTH-4;
    dysize = DisplayHeight() - 2*BORWIDTH-26;
    if (dxsize > dysize)
	dxsize = dysize;
    else
	dysize = dxsize;
    sprintf(defgeom, "=%dx%d+2+25", dxsize, dysize);
    wind = XCreate("Metafile display", progname, geom,
		defgeom, &mainframe, 64, 64);
    if (wind == 0)
	error(SYSTEM, "can't create window");
    XStoreName(wind, name);
    XMapWindow(wind);
    XSelectInput(wind, ButtonPressed|ExposeWindow|ExposeRegion);
    dxsize = mainframe.width;
    dysize = mainframe.height;
    if (DisplayPlanes() < 2 || XGetColorCells(0, 4, 0, &dummy, pixel) == 0) {
	pixel[0] = pixel[1] = pixel[2] = pixel[3] = BlackPixel;
    } else {
	cdef[0].red = cdef[0].green = cdef[0].blue = 0;
	cdef[0].pixel = pixel[0];
	cdef[1].red = MAXCOLOR; cdef[1].green = cdef[1].blue = 0;
	cdef[1].pixel = pixel[1];
	cdef[2].red = cdef[2].blue = 0; cdef[2].green = MAXCOLOR;
	cdef[2].pixel = pixel[2];
	cdef[3].red = cdef[3].green = 0; cdef[3].blue = MAXCOLOR;
	cdef[3].pixel = pixel[3];
	XStoreColors(4, cdef);
    }
}



endpage()		/* end of this graph */
{
    XExposeEvent  evnt;

    XFeep(0);
    for ( ; ; ) {
	XNextEvent(&evnt);
	if (evnt.type == ExposeWindow) {
	    dxsize = evnt.width;
	    dysize = evnt.height;
	    if (dxsize > dysize)
		dxsize = dysize;
	    else
		dysize = dxsize;
	    replay(0, 0, XYSIZE, XYSIZE);
	} else if (evnt.type == ExposeRegion) {
	    replay((int)((long)XYSIZE*evnt.x/dxsize),
			(int)((long)XYSIZE*(dysize-evnt.y-evnt.height)/dysize),
			(int)((long)XYSIZE*(evnt.x+evnt.width)/dxsize),
			(int)((long)XYSIZE*(dysize-evnt.y)/dysize));
	} else
	    break;
    }
    XClear(wind);

}



printstr(p)		/* output a string */

register PRIMITIVE  *p;

{
    int  fn;

    fn = fontmap[(p->arg0 >> 2) & 017];		/* get font number */
						/* set font */
    if (font[fn].id == 0) {
	font[fn].id = XGetFont(font[fn].name);
	if (font[fn].id == 0) {
	    sprintf(errmsg, "can't open font \"%s\"", font[fn].name);
	    error(SYSTEM, errmsg);
	}
    }

    XTextPad(wind, mapx(p->xy[XMN]), mapy(p->xy[YMN]), p->args, strlen(p->args),
		font[fn].id, 0, 0, pixel[p->arg0&03], WhitePixel,
		GXcopy, AllPlanes);

}





plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
    int  x1, y1, x2, y2;
    int  pp, ps, pw, ph;
    Pixmap  pm;

    pp = (p->arg0 >> 4) & 03;
    if (p->arg0 & 0100 && pp != 0)
	pp += 03;

    pm = pixpat(pp, p->arg0&03);

    ps = WIDTH((p->arg0 >> 2) & 03);
    pw = CONV((ps)/2, dxsize);
    ph = CONV((ps)/2, dysize);

    x1 = mapx(p->xy[XMN]);
    x2 = mapx(p->xy[XMX]);
    if (p->arg0 & 0100) {		/* reverse slope */
	y1 = mapy(p->xy[YMX]);
	y2 = mapy(p->xy[YMN]);
    } else {
	y1 = mapy(p->xy[YMN]);
	y2 = mapy(p->xy[YMX]);
    }

    if (pm != BlackPixmap || pw+ph > 0)
	paintline(x1,y1,x2-x1,y2-y1,pw,ph,pm,0L,0L,-1);
    else
	XLine(wind,x1,y1,x2,y2,2*pw+1,2*ph+1,BlackPixel,GXcopy,AllPlanes);


}



/*
 *  This routine paints a line with calls to fill().  The line can
 *    start and end at arbitrary points on a longer line segment.
 */

paintline(x, y, run, rise, hrad, vrad, pm, run2, rise2, n)

register int  x, y;
int  run, rise;
int  hrad, vrad;
Pixmap  pm;
long  run2, rise2;
int  n;

{
    int  xstep, ystep;

    if (run >= 0)
	xstep = 1;
    else {
	xstep = -1;
	run = -run;
    }
    if (rise >= 0)
	ystep = 1;
    else {
	ystep = -1;
	rise = -rise;
    }
    if (n < 0)
	n = max(run, rise);

    if (run > rise)
	while (n >= 0)
	    if (run2 >= rise2) {
		fill(x-hrad, y-vrad, x+hrad, y+vrad, pm);
		n--;
		x += xstep;
		rise2 += rise;
	    } else {
		y += ystep;
		run2 += run;
	    }
    else
	while (n >= 0)
	    if (rise2 >= run2) {
		fill(x-hrad, y-vrad, x+hrad, y+vrad, pm);
		n--;
		y += ystep;
		run2 += run;
	    } else {
		x += xstep;
		rise2 += rise;
	    }

}




Pixmap
pixpat(pp, pc)			/* return Pixmap for pattern number */

int  pp, pc;

{
    static int  curpat = -1;		/* current pattern */
    static int  curcol = -1;		/* current color */
    static Pixmap  curpix = 0;		/* current Pixmap id */
    short  bitval[16];
    Bitmap  bm;
    register int  i, j;

    if (pp == curpat && pc == curcol)	/* already set */
	return(curpix);

    if (pp == 0 && pc == 0)		/* all black */
	return(BlackPixmap);

    if (curpix != 0)			/* free old Pixmap */
	XFreePixmap(curpix);

    if (pp == 0)			/* constant Pixmap */
	curpix = XMakeTile(pixel[pc]);
    else {
#if  (PATSIZE != 16)
	error(SYSTEM, "bad code segment in pixpat");
#else
	for (i = 0; i < 16; i++) {
	    bitval[15-i] = 0;
	    for (j = 0; j < 16; j++)
		bitval[15-i] |= ((pattern[pp][i>>3][j] >> (i&07)) & 01) << j;
	}
	bm = XStoreBitmap(16, 16, bitval);
	curpix = XMakePixmap(bm, pixel[pc], WhitePixel);
	XFreeBitmap(bm);
#endif
    }
    curpat = pp;
    curcol = pc;
    return(curpix);

}



fillrect(p)			/* fill a rectangle */

register PRIMITIVE  *p;

{
    int  left, right, top, bottom;

    left = mapx(p->xy[XMN]);
    bottom = mapy(p->xy[YMN]);
    right = mapx(p->xy[XMX]);
    top = mapy(p->xy[YMX]);

    fill(left, top, right, bottom, pixpat(pati[(p->arg0>>2)&03],p->arg0&03));

}



#ifdef  nyet

filltri(p)			/* fill a triangle */

register PRIMITIVE  *p;

{
    Vertex  v[4];
    int  left, right, top, bottom;

    left = mapx(p->xy[XMN]);
    right = mapx(p->xy[XMX]);
    top = mapy(p->xy[YMX]);
    bottom = mapy(p->xy[YMN]);

    switch (p->arg0 & 060) {
    case 0:			/* right (& down) */
	v[0].x = left; v[0].y = bottom;
	v[1].x = right; v[1].y = bottom;
	v[2].x = right; v[2].y = top;
	break;
    case 020:			/* up */
	v[0].x = right; v[0].y = bottom;
	v[1].x = right; v[1].y = top;
	v[2].x = left; v[2].y = top;
	break;
    case 040:			/* left */
	v[0].x = right; v[0].y = top;
	v[1].x = left; v[1].y = top;
	v[2].x = left; v[2].y = bottom;
	break;
    case 060:			/* down */
	v[0].x = left; v[0].y = top;
	v[1].x = left; v[1].y = bottom;
	v[2].x = right; v[2].y = bottom;
	break;
    }
    v[3].x = v[0].x; v[3].y = v[0].y;
    v[0].flags = v[1].flags = v[2].flags = v[3].flags = 0;
    
    XDrawTiled(wind, v, 4, pixpat(pati[(p->arg0>>2)&03],p->arg0&03),
		GXcopy, AllPlanes);

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

#else

filltri() {}
fillpoly() {}

#endif
