#ifndef lint
static const char	RCSid[] = "$Id: x11plot.c,v 1.4 2007/11/21 18:51:05 greg Exp $";
#endif
/*
 *   X window plotting functions
 *
 *      2/26/86
 */


#include  "meta.h"

#include  "plot.h"

#include "string.h"

#undef  TRUE

#undef  FALSE

#include  <X11/Xlib.h>

extern void	replay(int xmin, int ymin, int xmax, int ymax);

#define  BORWIDTH  5

#define  BlackPix BlackPixel(dpy,0 )
#define  WhitePix WhitePixel(dpy,0 )

#define  mapx(x)  CONV(x,dxsize)
#define  mapy(y)  CONV((XYSIZE-1)-(y),dysize)

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

#define	LONG_DASHED_LIST_LENGTH		2
#define	SHORT_DASHED_LIST_LENGTH	2
#define	DOTTED_LIST_LENGTH		2

int dash_list_length[] = {
	0,
	LONG_DASHED_LIST_LENGTH,
	SHORT_DASHED_LIST_LENGTH,
	DOTTED_LIST_LENGTH,
	};

unsigned char long_dashed[LONG_DASHED_LIST_LENGTH] =	{4,7};
unsigned char short_dashed[SHORT_DASHED_LIST_LENGTH] =	{4,3};
unsigned char dotted[DOTTED_LIST_LENGTH] =		{5,2};

unsigned char *Linetypes[] = {
	0,
	long_dashed,
	short_dashed,
	dotted,
	};

Display	*dpy;
Window	wind;			/* our window */
Colormap cmap;
GC	gc;
Font	curfont;		/* current font */
int	curwidth;		/* current width of lines */
int	curfill;		/* current fill style (FillTiled or FillSolid) */
int	curcol;			/* current color */
int	curlinetype;		/* current line style */
XGCValues gcval;
int	pixel[4];
int	dxsize, dysize;		/* window size */
int	debug = False;		/* use XSynchronize if true */

static void
adjustsize()
{
	if (dxsize > dysize)
		dxsize = dysize;
	else
		dysize = dxsize;
}


void
init(name, geom)		/* initialize window */
char  *name;
char  *geom;
{
    char  defgeom[32];
    XColor  cdef,dum;
    XEvent  evnt;

    curfont = curfill = curcol = -1;
    curlinetype = 0;

    if ((dpy=XOpenDisplay("")) == NULL)
	error(SYSTEM, "can't open display");
    
    /* for debugging so we don't have to flush always */
    if (debug) 
	(void) XSynchronize(dpy, True);

    dxsize = DisplayWidth(dpy,0) - 2*BORWIDTH-4;
    dysize = DisplayHeight(dpy,0) - 2*BORWIDTH-100;
    adjustsize();

    sprintf(defgeom, "=%dx%d+2+25", dxsize, dysize);
    /* XUseGeometry(dpy,0,geom,defgeom,BORWIDTH,100,100,100,100,
		&xoff,&yoff,&dxsize,&dysize); */
    gc = DefaultGC(dpy,0);			/* get default gc */
    cmap = DefaultColormap(dpy,0);		/* and colormap */

    wind = XCreateSimpleWindow(dpy,DefaultRootWindow(dpy),0,0,dxsize,dysize,
		BORWIDTH,BlackPix,WhitePix);
    if (wind == 0)
	error(SYSTEM, "can't create window");
    XStoreName(dpy, wind, name);
    XMapRaised(dpy, wind);
    XSelectInput(dpy, wind, StructureNotifyMask | ButtonPressMask | ExposureMask);
    if (DisplayPlanes(dpy,0) < 2) {	/* less than 2 color planes, use black */
	    pixel[0] = pixel[1] = pixel[2] = pixel[3] = BlackPix;
    } else {
	if (XAllocNamedColor(dpy, cmap, "black", &dum, &cdef)==0)
		error(USER,"cannot allocate black!");
	pixel[0] = cdef.pixel;
	if (XAllocNamedColor(dpy, cmap, "red", &dum, &cdef)==0)
		{
		error(WARNING,"cannot allocate red");
		cdef.pixel = BlackPix;
		}
	pixel[1] = cdef.pixel;
	if (XAllocNamedColor(dpy, cmap, "green", &dum, &cdef)==0)
		{
		error(WARNING,"cannot allocate green");
		cdef.pixel = BlackPix;
		}
	pixel[2] = cdef.pixel;
	if (XAllocNamedColor(dpy, cmap, "blue", &dum, &cdef)==0)
		{
		error(WARNING,"cannot allocate blue");
		cdef.pixel = BlackPix;
		}
	pixel[3] = cdef.pixel;
    }

    while (1)
	{
	XNextEvent(dpy, &evnt);
	if (evnt.type == ConfigureNotify) /* wait for first ConfigureNotify */
		break;
	}
    dxsize = evnt.xconfigure.width;
    dysize = evnt.xconfigure.height;
    adjustsize();
    while (1)
	{
	XNextEvent(dpy, &evnt);
	if (evnt.type == Expose)	/* wait for first Expose */
		break;
	}
}


void
endpage()		/* end of this graph */
{
    XEvent  evnt;
    int    quit=False;

    XFlush(dpy);
    XBell(dpy, 0);
    while ( !quit ) {
	XNextEvent(dpy, &evnt);
	switch (evnt.type) {
	  case ConfigureNotify:
		dxsize = evnt.xconfigure.width;
		dysize = evnt.xconfigure.height;
		adjustsize();
		break;
	  case Expose:
		replay((int)((long)XYSIZE*evnt.xexpose.x/dxsize),
		   (int)((long)XYSIZE*(dysize-evnt.xexpose.y-evnt.xexpose.height)/dysize),
		   (int)((long)XYSIZE*(evnt.xexpose.x+evnt.xexpose.width)/dxsize),
		   (int)((long)XYSIZE*(dysize-evnt.xexpose.y)/dysize));
		break;
	  case ButtonPress:
		quit = True;
		break;
	} /* switch */
    } /* for */
    XClearWindow(dpy, wind);

}



void
printstr(p)		/* output a string */

register PRIMITIVE  *p;

{
    int  fn;
    int  col;

    fn = fontmap[(p->arg0 >> 2) & 017];		/* get font number */
						/* set font */
    if (font[fn].id == 0) 
	{
	font[fn].id = XLoadFont(dpy, font[fn].name);
	if (font[fn].id == 0) 
	    {
	    sprintf(errmsg, "can't open font \"%s\"", font[fn].name);
	    error(SYSTEM, errmsg);
	    }
	}
    col = p->arg0 & 03;
    if (curcol != col)
	{
	curcol = col;
	XSetForeground(dpy, gc, pixel[col]);
	}
    if (curfont != font[fn].id) 
	{
	curfont = font[fn].id;
	XSetFont(dpy, gc, curfont);
	}

    XDrawImageString(dpy, wind, gc, mapx(p->xy[XMN]), mapy(p->xy[YMN]),
		p->args, strlen(p->args));
}





void
plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
    int  x1, y1, x2, y2;
    int  linetype, ps, pw;
    int col;

    linetype = (p->arg0 >> 4) & 03;	/* line style (solid, dashed, etc) */

    col = p->arg0 & 03;			/* color */

    ps = WIDTH((p->arg0 >> 2) & 03);
    pw = CONV((ps)/2, dxsize);

    x1 = mapx(p->xy[XMN]);
    x2 = mapx(p->xy[XMX]);
    if (p->arg0 & 0100) {		/* reverse slope */
	y1 = mapy(p->xy[YMX]);
	y2 = mapy(p->xy[YMN]);
    } else {
	y1 = mapy(p->xy[YMN]);
	y2 = mapy(p->xy[YMX]);
    }

    if (curcol != col || curwidth != pw)
	{
	curcol = col;
	gcval.foreground = pixel[col];
	curwidth = pw;
	gcval.line_width = pw*2+1;	/* convert to thickness in pixels */
	gcval.join_style = JoinRound;
	XChangeGC(dpy, gc, GCJoinStyle|GCLineWidth|GCForeground, &gcval);
	}
    if (curlinetype != linetype)
	{
	curlinetype = linetype;
	if (linetype==0)
		gcval.line_style = LineSolid;
	else
		{
		gcval.line_style = LineOnOffDash;
		XSetDashes(dpy, gc, 0, Linetypes[linetype], dash_list_length[linetype]);
		}
	XChangeGC(dpy, gc, GCLineStyle, &gcval);
	}
    XDrawLine(dpy,wind,gc,x1,y1,x2,y2);
}


void
pXFlush()
{
	XFlush(dpy);
}


#ifdef  nyet

static void
fill(xmin,ymin,xmax,ymax,pm)
int xmin,ymin,xmax,ymax;
Pixmap pm;
{
	if (pm != 0 && curpat != pm)
		{
		XSetTile(dpy, gc, pm);
		curpat = pm;
		}
	XFillRectangle(dpy, wind, gc, xmin, ymin, xmax-xmin+1, ymax-ymin+1);
}


void
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




void
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
    if (curfill != FillTiled)
	    {
	    XSetFillStyle(dpy, gc, FillTiled);
	    curfill = FillTiled;
	    XSetTile(dpy, gc, pixpat(pati[(p->arg0>>2)&03],p->arg0&03));
	    }
    XDrawFIlled(dpy, wind, gc, v, 4);
}


vpod
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


void
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

void filltri(PRIMITIVE *p) {}
void fillpoly(PRIMITIVE *p) {}
void fillrect(PRIMITIVE *p) {}

#endif
