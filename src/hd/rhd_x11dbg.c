#ifndef lint
static const char	RCSid[] = "$Id: rhd_x11dbg.c,v 3.2 2003/07/21 22:30:18 schorsch Exp $";
#endif
/*
 * X11 driver for holodeck display.
 * Based on rview driver.
 */

#include "standard.h"
#include "rhd_qtree.h"

#include  <X11/Xlib.h>
#include  <X11/cursorfont.h>
#include  <X11/Xutil.h>

#include  "x11icon.h"

#ifndef FEQ
#define FEQ(a,b)	((a)-(b) <= FTINY && (a)-(b) >= -FTINY)
#endif

#define	CTRL(c)		((c)-'@')

#define GAMMA		2.2		/* default gamma correction */

#define MOVPCT		10		/* percent distance to move */
#define MOVDIR(b)	((b)==Button1 ? 1 : (b)==Button2 ? 0 : -1)

#define MINWIDTH	480		/* minimum graphics window width */
#define MINHEIGHT	400		/* minimum graphics window height */

#define VIEWDIST	356		/* assumed viewing distance (mm) */

#define BORWIDTH	5		/* border width */

#define  ourscreen	DefaultScreen(ourdisplay)
#define  ourroot	RootWindow(ourdisplay,ourscreen)
#define  ourmask	(StructureNotifyMask|ExposureMask|KeyPressMask|\
			ButtonPressMask|ButtonReleaseMask)

#define  levptr(etype)	((etype *)&currentevent)

struct driver	odev;			/* global device driver structure */

static XEvent  currentevent;		/* current event */

static int  ncolors = 0;		/* color table size */
static int  mapped = 0;			/* window is mapped? */
static unsigned long  *pixval = NULL;	/* allocated pixels */
static unsigned long  ourblack=0, ourwhite=1;

static Display  *ourdisplay = NULL;	/* our display */
static XVisualInfo  ourvinfo;		/* our visual information */
static Window  gwind = 0;		/* our graphics window */
static GC  ourgc = 0;			/* our graphics context for drawing */
static Colormap ourmap = 0;		/* our color map */

static double	pwidth, pheight;	/* pixel dimensions (mm) */

static int	inpresflags;		/* input result flags */

static int	headlocked = 0;		/* lock vertical motion */

static int  getpixels(), xnewcolr(), freepixels(), resizewindow(),
		getevent(), getkey(), moveview(), getmove(), fixwindow();
static unsigned long  true_pixel();


static int
mytmflags()			/* figure out tone mapping flags */
{
	extern char	*progname;
	register char	*cp, *tail;
					/* find basic name */
	for (cp = tail = progname; *cp; cp++)
		if (*cp == '/')
			tail = cp+1;
	for (cp = tail; *cp && *cp != '.'; cp++)
		;
	if (cp-tail == 3 && !strncmp(tail, "x11", 3))
		return(TM_F_CAMERA);
	if (cp-tail == 4 && !strncmp(tail, "x11h", 4))
		return(TM_F_HUMAN);
	error(USER, "illegal driver name");
}


dev_open(id)			/* initialize X11 driver */
char  *id;
{
	extern char  *getenv();
	char  *gv;
	double	gamval = GAMMA;
	int  nplanes;
	XSetWindowAttributes	ourwinattr;
	XWMHints  ourxwmhints;
	XSizeHints	oursizhints;
					/* set quadtree globals */
	qtMinNodesiz = 2;
					/* open display server */
	ourdisplay = XOpenDisplay(NULL);
	if (ourdisplay == NULL)
		error(USER, "cannot open X-windows; DISPLAY variable set?\n");
					/* find a usable visual */
	nplanes = DisplayPlanes(ourdisplay, ourscreen);
	if (XMatchVisualInfo(ourdisplay,ourscreen,
				nplanes>12?nplanes:24,TrueColor,&ourvinfo) ||
			XMatchVisualInfo(ourdisplay,ourscreen,
				nplanes>12?nplanes:24,DirectColor,&ourvinfo)) {
		ourblack = 0;
		ourwhite = ourvinfo.red_mask |
				ourvinfo.green_mask |
				ourvinfo.blue_mask ;
	} else {
		if (nplanes < 4)
			error(INTERNAL, "not enough colors\n");
		if (!XMatchVisualInfo(ourdisplay,ourscreen,
					nplanes,PseudoColor,&ourvinfo) &&
				!XMatchVisualInfo(ourdisplay,ourscreen,
					nplanes,GrayScale,&ourvinfo))
			error(INTERNAL, "unsupported visual type\n");
		ourblack = BlackPixel(ourdisplay,ourscreen);
		ourwhite = WhitePixel(ourdisplay,ourscreen);
	}
					/* set gamma and tone mapping */
	if ((gv = XGetDefault(ourdisplay, "radiance", "gamma")) != NULL
			|| (gv = getenv("DISPLAY_GAMMA")) != NULL)
		gamval = atof(gv);
	if (tmInit(mytmflags(), stdprims, gamval) == NULL)
		error(SYSTEM, "not enough memory in dev_open");
					/* open window */
	ourwinattr.background_pixel = ourblack;
	ourwinattr.border_pixel = ourblack;
	ourwinattr.event_mask = ourmask;
					/* this is stupid */
	ourwinattr.colormap = XCreateColormap(ourdisplay, ourroot,
				ourvinfo.visual, AllocNone);
	gwind = XCreateWindow(ourdisplay, ourroot, 0, 0,
		DisplayWidth(ourdisplay,ourscreen)-2*BORWIDTH,
		DisplayHeight(ourdisplay,ourscreen)-2*BORWIDTH,
		BORWIDTH, ourvinfo.depth, InputOutput, ourvinfo.visual,
		CWBackPixel|CWBorderPixel|CWColormap|CWEventMask, &ourwinattr);
	if (gwind == 0)
		error(SYSTEM, "cannot create window\n");
   	XStoreName(ourdisplay, gwind, id);
					/* get graphics context */
	ourgc = XCreateGC(ourdisplay, gwind, 0, NULL);
					/* set window manager hints */
	ourxwmhints.flags = InputHint|IconPixmapHint;
	ourxwmhints.input = True;
	ourxwmhints.icon_pixmap = XCreateBitmapFromData(ourdisplay,
			gwind, x11icon_bits, x11icon_width, x11icon_height);
	XSetWMHints(ourdisplay, gwind, &ourxwmhints);
	oursizhints.min_width = MINWIDTH;
	oursizhints.min_height = MINHEIGHT;
	oursizhints.flags = PMinSize;
	XSetNormalHints(ourdisplay, gwind, &oursizhints);
					/* map the window and get its size */
	XMapWindow(ourdisplay, gwind);
	dev_input();
					/* allocate our leaf pile */
	if (!qtAllocLeaves(DisplayWidth(ourdisplay,ourscreen) *
			DisplayHeight(ourdisplay,ourscreen) /
			(qtMinNodesiz*qtMinNodesiz)))
		error(SYSTEM, "insufficient memory for leaf storage");

					/* figure out sensible view */
	pwidth = (double)DisplayWidthMM(ourdisplay, ourscreen) /
			DisplayWidth(ourdisplay, ourscreen);
	pheight = (double)DisplayHeightMM(ourdisplay, ourscreen) /
			DisplayHeight(ourdisplay, ourscreen);
	odev.v = stdview;
	odev.name = id;
	odev.v.type = VT_PER;
	odev.v.horiz = 2.*180./PI * atan(0.5/VIEWDIST*pwidth*odev.hres);
	odev.v.vert = 2.*180./PI * atan(0.5/VIEWDIST*pheight*odev.vres);
	odev.ifd = ConnectionNumber(ourdisplay);
}


dev_close()			/* close our display */
{
	freepixels();
	XFreeGC(ourdisplay, ourgc);
	XDestroyWindow(ourdisplay, gwind);
	gwind = 0;
	ourgc = 0;
	XCloseDisplay(ourdisplay);
	ourdisplay = NULL;
	qtFreeLeaves();
	tmDone(NULL);
	odev.v.type = 0;
	odev.hres = odev.vres = 0;
	odev.ifd = -1;
}


dev_view(nv)			/* assign new driver view */
VIEW	*nv;
{
	if (nv != &odev.v) {
		if (!FEQ(nv->horiz,odev.v.horiz) ||	/* resize window? */
				!FEQ(nv->vert,odev.v.vert)) {
			odev.hres = 2.*VIEWDIST/pwidth *
					tan(PI/180./2.*nv->horiz);
			odev.vres = 2.*VIEWDIST/pheight *
					tan(PI/180./2.*nv->vert);
			XResizeWindow(ourdisplay, gwind, odev.hres, odev.vres);
		}
		odev.v = *nv;
	}
	qtReplant();
}


dev_auxcom(cmd, args)		/* process an auxiliary command */
char	*cmd, *args;
{
	sprintf(errmsg, "%s: unknown command", cmd);
	error(COMMAND, errmsg);
}


VIEW *
dev_auxview(n, hvres)		/* return nth auxiliary view */
int	n;
int	hvres[2];
{
	if (n)
		return(NULL);
	hvres[0] = odev.hres; hvres[1] = odev.vres;
	return(&odev.v);
}


int
dev_input()			/* get X11 input */
{
	inpresflags = 0;

	do
		getevent();

	while (XQLength(ourdisplay) > 0);

	return(inpresflags);
}


dev_paintr(rgb, xmin, ymin, xmax, ymax)		/* fill a rectangle */
BYTE	rgb[3];
int  xmin, ymin, xmax, ymax;
{
	unsigned long  pixel;

	if (!mapped)
		return;
	if (ncolors > 0)
		pixel = pixval[get_pixel(rgb, xnewcolr)];
	else
		pixel = true_pixel(rgb);
	XSetForeground(ourdisplay, ourgc, pixel);
/*	XFillRectangle(ourdisplay, gwind, 
		ourgc, xmin, odev.vres-ymax, xmax-xmin, ymax-ymin);
*/
	xmin = (xmin+xmax)>>1;
	ymin = (ymin+ymax)>>1;
	XDrawPoint(ourdisplay, gwind, ourgc, xmin, odev.vres-ymin);
}


int
dev_flush()			/* flush output */
{
	qtUpdate();
	return(XPending(ourdisplay));
}


static
xnewcolr(ndx, r, g, b)		/* enter a color into hardware table */
int  ndx;
int  r, g, b;
{
	XColor  xcolor;

	xcolor.pixel = pixval[ndx];
	xcolor.red = r << 8;
	xcolor.green = g << 8;
	xcolor.blue = b << 8;
	xcolor.flags = DoRed|DoGreen|DoBlue;

	XStoreColor(ourdisplay, ourmap, &xcolor);
}


static int
getpixels()				/* get the color map */
{
	XColor  thiscolor;
	register int  i, j;

	if (ncolors > 0)
		return(ncolors);
	if (ourvinfo.visual == DefaultVisual(ourdisplay,ourscreen)) {
		ourmap = DefaultColormap(ourdisplay,ourscreen);
		goto loop;
	}
newmap:
	ourmap = XCreateColormap(ourdisplay,gwind,ourvinfo.visual,AllocNone);
loop:
	for (ncolors = ourvinfo.colormap_size;
			ncolors > ourvinfo.colormap_size/3;
			ncolors = ncolors*.937) {
		pixval = (unsigned long *)malloc(ncolors*sizeof(unsigned long));
		if (pixval == NULL)
			return(ncolors = 0);
		if (XAllocColorCells(ourdisplay,ourmap,0,NULL,0,pixval,ncolors))
			break;
		free((void *)pixval);
		pixval = NULL;
	}
	if (pixval == NULL) {
		if (ourmap == DefaultColormap(ourdisplay,ourscreen))
			goto newmap;		/* try it with our map */
		else
			return(ncolors = 0);	/* failed */
	}
	if (ourmap != DefaultColormap(ourdisplay,ourscreen))
		for (i = 0; i < ncolors; i++) {	/* reset black and white */
			if (pixval[i] != ourblack && pixval[i] != ourwhite)
				continue;
			thiscolor.pixel = pixval[i];
			thiscolor.flags = DoRed|DoGreen|DoBlue;
			XQueryColor(ourdisplay,
					DefaultColormap(ourdisplay,ourscreen),
					&thiscolor);
			XStoreColor(ourdisplay, ourmap, &thiscolor);
			for (j = i; j+1 < ncolors; j++)
				pixval[j] = pixval[j+1];
			ncolors--;
			i--;
		}
	XSetWindowColormap(ourdisplay, gwind, ourmap);
	return(ncolors);
}


static
freepixels()				/* free our pixels */
{
	if (ncolors == 0)
		return;
	XFreeColors(ourdisplay,ourmap,pixval,ncolors,0L);
	free((void *)pixval);
	pixval = NULL;
	ncolors = 0;
	if (ourmap != DefaultColormap(ourdisplay,ourscreen))
		XFreeColormap(ourdisplay, ourmap);
	ourmap = 0;
}


static unsigned long
true_pixel(rgb)			/* return true pixel value for color */
register BYTE	rgb[3];
{
	register unsigned long  rval;

	rval = ourvinfo.red_mask*rgb[RED]/255 & ourvinfo.red_mask;
	rval |= ourvinfo.green_mask*rgb[GRN]/255 & ourvinfo.green_mask;
	rval |= ourvinfo.blue_mask*rgb[BLU]/255 & ourvinfo.blue_mask;
	return(rval);
}


static
getevent()			/* get next event */
{
	XNextEvent(ourdisplay, levptr(XEvent));
	switch (levptr(XEvent)->type) {
	case ConfigureNotify:
		resizewindow(levptr(XConfigureEvent));
		break;
	case UnmapNotify:
		mapped = 0;
		freepixels();
		break;
	case MapNotify:
		if (ourvinfo.class == PseudoColor ||
				ourvinfo.class == GrayScale) {
			if (getpixels() == 0)
				error(SYSTEM, "cannot allocate colors\n");
			new_ctab(ncolors);
		}
		mapped = 1;
		break;
	case Expose:
		fixwindow(levptr(XExposeEvent));
		break;
	case KeyPress:
		getkey(levptr(XKeyPressedEvent));
		break;
	case ButtonPress:
		getmove(levptr(XButtonPressedEvent));
		break;
	}
}


static
ilclip(dp, wp)			/* clip world coordinates to device */
int	dp[2][2];
FVECT	wp[2];
{
	static FVECT	vmin = {0.,0.,0.}, vmax = {1.,1.,FHUGE};
	FVECT	ip[2];
				/* not exactly right, but who cares? */
	viewloc(ip[0], &odev.v, wp[0]);
	viewloc(ip[1], &odev.v, wp[1]);
	if (!clip(ip[0], ip[1], vmin, vmax))
		return(0);
	dp[0][0] = ip[0][0]*odev.hres;
	dp[0][1] = ip[0][1]*odev.vres;
	dp[1][0] = ip[1][0]*odev.hres;
	dp[1][1] = ip[1][1]*odev.vres;
	return(1);
}


static
draw3dline(wp)			/* draw 3d line in world coordinates */
FVECT	wp[2];
{
	int	dp[2][2];

	if (!ilclip(dp, wp))
		return;
	XDrawLine(ourdisplay, gwind, ourgc,
			dp[0][0], odev.vres-1 - dp[0][1],
			dp[1][0], odev.vres-1 - dp[1][1]);
}


static
draw_grids()			/* draw holodeck section grids */
{
	static BYTE	gridrgb[3] = {0x0, 0xff, 0xff};
	unsigned long  pixel;

	if (!mapped || odev.v.type != VT_PER)
		return;
	if (ncolors > 0)
		pixel = pixval[get_pixel(gridrgb, xnewcolr)];
	else
		pixel = true_pixel(gridrgb);
	XSetForeground(ourdisplay, ourgc, pixel);
					/* draw each grid line */
	gridlines(draw3dline);
}


static
moveview(dx, dy, move)		/* move our view */
int	dx, dy, move;
{
	VIEW	nv;
	double	d;
	register int	i, li;
				/* start with old view */
	nv = odev.v;
				/* change view direction */
	if (move) {
		if ((li = qtFindLeaf(dx, dy)) < 0)
			return(0);	/* not on window */
		for (i = 0; i < 3; i++)
			nv.vdir[i] = qtL.wp[li][i] - nv.vp[i];
	} else {
		if (viewray(nv.vp, nv.vdir, &odev.v,
				(dx+.5)/odev.hres, (dy+.5)/odev.vres) < -FTINY)
			return(0);	/* outside view */
	}
				/* move viewpoint */
	if (move > 0)
		for (i = 0; i < 3; i++)
			nv.vp[i] += MOVPCT/100. * nv.vdir[i];
	else if (move < 0)
		for (i = 0; i < 3; i++)
			nv.vp[i] -= MOVPCT/100. * nv.vdir[i];
	if (move && headlocked) {
		d = 0;		/* bring head back to same height */
		for (i = 0; i < 3; i++)
			d += odev.v.vup[i] * (odev.v.vp[i] - nv.vp[i]);
		for (i = 0; i < 3; i++)
			nv.vp[i] += d * odev.v.vup[i];
	}
	if (setview(&nv) != NULL)
		return(0);	/* illegal view */
	dev_view(&nv);
	inpresflags |= DEV_NEWVIEW;
	return(1);
}


static
getmove(ebut)				/* get view change */
XButtonPressedEvent	*ebut;
{
	int	whichbutton = ebut->button;
	int	oldnodesiz = qtMinNodesiz;
	Window	rootw, childw;
	int	rootx, rooty, wx, wy;
	unsigned int	statemask;

	qtMinNodesiz = 16;		/* for quicker update */
	XNoOp(ourdisplay);

	while (!XCheckMaskEvent(ourdisplay,
			ButtonReleaseMask, levptr(XEvent))) {

		if (!XQueryPointer(ourdisplay, gwind, &rootw, &childw,
				&rootx, &rooty, &wx, &wy, &statemask))
			break;		/* on another screen */

		if (!moveview(wx, odev.vres-1-wy, MOVDIR(whichbutton))) {
			sleep(1);
			continue;
		}
		XClearWindow(ourdisplay, gwind);
		qtUpdate();
		draw_grids();
	}
	if (!(inpresflags & DEV_NEWVIEW)) {	/* do final motion */
		whichbutton = levptr(XButtonReleasedEvent)->button;
		wx = levptr(XButtonReleasedEvent)->x;
		wy = levptr(XButtonReleasedEvent)->y;
		moveview(wx, odev.vres-1-wy, MOVDIR(whichbutton));
	}
	dev_flush();

	qtMinNodesiz = oldnodesiz;	/* restore quadtree resolution */
}


static
getkey(ekey)				/* get input key */
register XKeyPressedEvent  *ekey;
{
	int  n;
	char	buf[8];

	n = XLookupString(ekey, buf, sizeof(buf), NULL, NULL);
	if (n != 1)
		return;
	switch (buf[0]) {
	case 'h':			/* turn on height motion lock */
		headlocked = 1;
		return;
	case 'H':			/* turn off height motion lock */
		headlocked = 0;
		return;
	case CTRL('S'):
	case 'p':			/* pause computation */
		inpresflags |= DEV_WAIT;
		return;
	case 'v':			/* spit out view */
		inpresflags |= DEV_PUTVIEW;
		return;
	case CTRL('Q'):
	case '\n':
	case '\r':			/* resume computation */
		inpresflags |= DEV_RESUME;
		return;
	case CTRL('R'):			/* redraw screen */
		if (ncolors > 0)
			new_ctab(ncolors);
		qtRedraw(0, 0, odev.hres, odev.vres);
		return;
	case CTRL('L'):			/* refresh from server */
		if (inpresflags & DEV_REDRAW)
			return;
		XClearWindow(ourdisplay, gwind);
		draw_grids();
		XFlush(ourdisplay);
		qtCompost(100);			/* unload the old tree */
		if (ncolors > 0)
			new_ctab(ncolors);
		inpresflags |= DEV_REDRAW;	/* resend values from server */
		return;
	case CTRL('D'):
	case 'Q':
	case 'q':			/* quit the program */
		inpresflags |= DEV_SHUTDOWN;
		return;
	default:
		XBell(ourdisplay, 0);
		return;
	}
}


static
fixwindow(eexp)				/* repair damage to window */
register XExposeEvent  *eexp;
{
	if (odev.hres == 0 || odev.vres == 0) {	/* first exposure */
		odev.hres = eexp->width;
		odev.vres = eexp->height;
		inpresflags |= DEV_NEWSIZE;
	}
	qtRedraw(eexp->x, odev.vres - eexp->y - eexp->height,
			eexp->x + eexp->width, odev.vres - eexp->y);
}


static
resizewindow(ersz)			/* resize window */
register XConfigureEvent  *ersz;
{
	if (ersz->width == odev.hres && ersz->height == odev.vres)
		return;

	odev.hres = ersz->width;
	odev.vres = ersz->height;

	odev.v.horiz = 2.*180./PI * atan(0.5/VIEWDIST*pwidth*odev.hres);
	odev.v.vert = 2.*180./PI * atan(0.5/VIEWDIST*pheight*odev.vres);

	inpresflags |= DEV_NEWSIZE|DEV_NEWVIEW;
}
