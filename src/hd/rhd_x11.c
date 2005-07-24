#ifndef lint
static const char	RCSid[] = "$Id: rhd_x11.c,v 3.41 2005/07/24 19:53:08 greg Exp $";
#endif
/*
 * X11 driver for holodeck display.
 * Based on rview driver.
 */

#include  <stdlib.h>
#include  <stdio.h>
#include  <X11/Xlib.h>
#include  <X11/cursorfont.h>
#include  <X11/Xutil.h>
#include  <time.h>

#include "platform.h"
#include "rtmath.h"
#include "rterror.h"
#include "plocate.h"
#include "rhdisp.h"
#include "rhd_qtree.h"
#include "x11icon.h"

#ifndef RAYQLEN
#define RAYQLEN		50000		/* max. rays to queue before flush */
#endif

#ifndef FEQ
#define FEQ(a,b)	((a)-(b) <= FTINY && (a)-(b) >= -FTINY)
#endif

#define GAMMA		2.2		/* default gamma correction */

#define FRAMESTATE(s)	(((s)&(ShiftMask|ControlMask))==(ShiftMask|ControlMask))

#define MOVPCT		7		/* percent distance to move /frame */
#define MOVDIR(b)	((b)==Button1 ? 1 : (b)==Button2 ? 0 : -1)
#define MOVDEG		(-5)		/* degrees to orbit CW/down /frame */
#define MOVORB(s)	((s)&ShiftMask ? 1 : (s)&ControlMask ? -1 : 0)

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

TMstruct	*tmGlobal;		/* global tone-mapping structure */

char odev_args[64];			/* command arguments */

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


static int mytmflags(void);
static void xnewcolr(int  ndx, int r, int g, int b);
static int getpixels(void);
static void freepixels(void);
static unsigned long true_pixel(register BYTE rgb[3]);
static void getevent(void);
static int ilclip(int dp[2][2], FVECT wp[2]);
static void draw3dline(FVECT wp[2]);
static void draw_grids(void);
static int moveview(int dx, int dy, int mov, int orb);
static void getframe(XButtonPressedEvent *ebut);
static void waitabit(void);
static void getmove(XButtonPressedEvent *ebut);
static void getkey(register XKeyPressedEvent *ekey);
static void fixwindow(register XExposeEvent *eexp);
static void resizewindow(register XConfigureEvent *ersz);


static int
mytmflags(void)			/* figure out tone mapping flags */
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
		return(TM_F_CAMERA|TM_F_NOSTDERR);
	if (cp-tail == 4 && !strncmp(tail, "x11h", 4))
		return(TM_F_HUMAN|TM_F_NOSTDERR);
	error(USER, "illegal driver name");
	return 0; /* pro forma return */
}


extern void
dev_open(			/* initialize X11 driver */
	char  *id
)
{
	static RGBPRIMS	myprims = STDPRIMS;
	char  *ev;
	double	gamval = GAMMA;
	RGBPRIMP	dpri = stdprims;
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
	if ((ev = XGetDefault(ourdisplay, "radiance", "gamma")) != NULL
			|| (ev = getenv("DISPLAY_GAMMA")) != NULL)
		gamval = atof(ev);
	if ((ev = getenv("DISPLAY_PRIMARIES")) != NULL &&
			sscanf(ev, "%f %f %f %f %f %f %f %f",
				&myprims[RED][CIEX],&myprims[RED][CIEY],
				&myprims[GRN][CIEX],&myprims[GRN][CIEY],
				&myprims[BLU][CIEX],&myprims[BLU][CIEY],
				&myprims[WHT][CIEX],&myprims[WHT][CIEY]) >= 6)
		dpri = myprims;
	tmGlobal = tmInit(mytmflags(), dpri, gamval);
	if (tmGlobal == NULL)
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
	ourxwmhints.icon_pixmap = XCreateBitmapFromData(ourdisplay, gwind,
			(char *)x11icon_bits, x11icon_width, x11icon_height);
	XSetWMHints(ourdisplay, gwind, &ourxwmhints);
	oursizhints.min_width = MINWIDTH;
	oursizhints.min_height = MINHEIGHT;
	oursizhints.flags = PMinSize;
	XSetNormalHints(ourdisplay, gwind, &oursizhints);
					/* figure out sensible view */
	pwidth = (double)DisplayWidthMM(ourdisplay, ourscreen) /
			DisplayWidth(ourdisplay, ourscreen);
	pheight = (double)DisplayHeightMM(ourdisplay, ourscreen) /
			DisplayHeight(ourdisplay, ourscreen);
	odev.v = stdview;
	odev.v.type = VT_PER;
					/* map the window and get its size */
	XMapWindow(ourdisplay, gwind);
	dev_input();			/* sets size and view angles */
					/* allocate our leaf pile */
	if (!qtAllocLeaves(DisplayWidth(ourdisplay,ourscreen) *
			DisplayHeight(ourdisplay,ourscreen) * 3 /
			(qtMinNodesiz*qtMinNodesiz*2)))
		error(SYSTEM, "insufficient memory for leaf storage");
	odev.name = id;
	odev.ifd = ConnectionNumber(ourdisplay);
}


extern void
dev_close(void)			/* close our display */
{
	freepixels();
	XFreeGC(ourdisplay, ourgc);
	XDestroyWindow(ourdisplay, gwind);
	gwind = 0;
	ourgc = 0;
	XCloseDisplay(ourdisplay);
	ourdisplay = NULL;
	qtFreeLeaves();
	tmDone(tmGlobal);
	odev.v.type = 0;
	odev.hres = odev.vres = 0;
	odev.ifd = -1;
}


extern void
dev_clear(void)			/* clear our quadtree */
{
	qtCompost(100);
	if (ncolors > 0)
		new_ctab(ncolors);
	rayqleft = 0;			/* hold off update */
}


extern int
dev_view(			/* assign new driver view */
	VIEW	*nv
)
{
	if (nv->type == VT_PAR ||		/* check view legality */
			nv->horiz > 160. || nv->vert > 160.) {
		error(COMMAND, "illegal view type/angle");
		nv->type = VT_PER;
		nv->horiz = odev.v.horiz;
		nv->vert = odev.v.vert;
		return(0);
	}
	if (nv->vfore > FTINY) {
		error(COMMAND, "cannot handle fore clipping");
		nv->vfore = 0.;
		return(0);
	}
	if (nv != &odev.v) {
		if (!FEQ(nv->horiz,odev.v.horiz) ||	/* resize window? */
				!FEQ(nv->vert,odev.v.vert)) {
			int	dw = DisplayWidth(ourdisplay,ourscreen);
			int	dh = DisplayHeight(ourdisplay,ourscreen);

			dw -= 25;	/* for window frame */
			dh -= 50;
			odev.hres = 2.*VIEWDIST/pwidth *
					tan(PI/180./2.*nv->horiz);
			odev.vres = 2.*VIEWDIST/pheight *
					tan(PI/180./2.*nv->vert);
			if (odev.hres > dw) {
				odev.vres = dw * odev.vres / odev.hres;
				odev.hres = dw;
			}
			if (odev.vres > dh) {
				odev.hres = dh * odev.hres / odev.vres;
				odev.vres = dh;
			}
			XResizeWindow(ourdisplay, gwind, odev.hres, odev.vres);
			dev_input();	/* wait for resize event */
		}
		odev.v = *nv;
	}
	qtReplant();
	return(1);
}


extern void
dev_section(		/* add octree for geometry rendering */
	char	*ofn
)
{
	/* unimplemented */
}


extern void
dev_auxcom(		/* process an auxiliary command */
	char	*cmd,
	char	*args
)
{
	sprintf(errmsg, "%s: unknown command", cmd);
	error(COMMAND, errmsg);
}


extern VIEW *
dev_auxview(		/* return nth auxiliary view */
	int	n,
	int	hvres[2]
)
{
	if (n)
		return(NULL);
	hvres[0] = odev.hres; hvres[1] = odev.vres;
	return(&odev.v);
}


extern int
dev_input(void)			/* get X11 input */
{
	inpresflags = 0;

	do
		getevent();

	while (XPending(ourdisplay) > 0);

	odev.inpready = 0;

	return(inpresflags);
}


extern void
dev_paintr(		/* fill a rectangle */
	BYTE	rgb[3],
	int  xmin,
	int  ymin,
	int  xmax,
	int  ymax
)
{
	unsigned long  pixel;

	if (!mapped)
		return;
	if (ncolors > 0)
		pixel = pixval[get_pixel(rgb, xnewcolr)];
	else
		pixel = true_pixel(rgb);
	XSetForeground(ourdisplay, ourgc, pixel);
	XFillRectangle(ourdisplay, gwind, 
		ourgc, xmin, odev.vres-ymax, xmax-xmin, ymax-ymin);
}


extern int
dev_flush(void)			/* flush output */
{
	qtUpdate();
	rayqleft = RAYQLEN;
	return(odev.inpready = XPending(ourdisplay));
}


static void
xnewcolr(		/* enter a color into hardware table */
	int  ndx,
	int r,
	int g,
	int b
)
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
getpixels(void)				/* get the color map */
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


static void
freepixels(void)				/* free our pixels */
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
true_pixel(			/* return true pixel value for color */
	register BYTE	rgb[3]
)
{
	register unsigned long  rval;

	rval = ourvinfo.red_mask*rgb[RED]/255 & ourvinfo.red_mask;
	rval |= ourvinfo.green_mask*rgb[GRN]/255 & ourvinfo.green_mask;
	rval |= ourvinfo.blue_mask*rgb[BLU]/255 & ourvinfo.blue_mask;
	return(rval);
}


static void
getevent(void)			/* get next event */
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
		if (FRAMESTATE(levptr(XButtonPressedEvent)->state))
			getframe(levptr(XButtonPressedEvent));
		else
			getmove(levptr(XButtonPressedEvent));
		break;
	}
}


static int
ilclip(			/* clip world coordinates to device */
	int	dp[2][2],
	FVECT	wp[2]
)
{
	static FVECT	vmin = {0.,0.,0.}, vmax = {1.-FTINY,1.-FTINY,FHUGE};
	FVECT	wpc[2], ip[2];
	double	d, d0, d1;
				/* check for points behind view */
	d = DOT(odev.v.vp, odev.v.vdir);
	d0 = DOT(wp[0], odev.v.vdir) - d;
	d1 = DOT(wp[1], odev.v.vdir) - d;
				/* work on copy of world points */
	if (d0 <= d1) {
		VCOPY(wpc[0], wp[0]); VCOPY(wpc[1], wp[1]);
	} else {
		d = d0; d0 = d1; d1 = d;
		VCOPY(wpc[1], wp[0]); VCOPY(wpc[0], wp[1]);
	}
	if (d0 <= FTINY) {
		if (d1 <= FTINY) return(0);
		VSUB(wpc[0], wpc[0], wpc[1]);
		d = .99*d1/(d1-d0);
		VSUM(wpc[0], wpc[1], wpc[0], d);
	}
				/* get view coordinates and clip to window */
	viewloc(ip[0], &odev.v, wpc[0]);
	viewloc(ip[1], &odev.v, wpc[1]);
	if (!clip(ip[0], ip[1], vmin, vmax))
		return(0);
	dp[0][0] = ip[0][0]*odev.hres;
	dp[0][1] = ip[0][1]*odev.vres;
	dp[1][0] = ip[1][0]*odev.hres;
	dp[1][1] = ip[1][1]*odev.vres;
	return(1);
}


static void
draw3dline(			/* draw 3d line in world coordinates */
	FVECT	wp[2]
)
{
	int	dp[2][2];

	if (!ilclip(dp, wp))
		return;
	XDrawLine(ourdisplay, gwind, ourgc,
			dp[0][0], odev.vres-1 - dp[0][1],
			dp[1][0], odev.vres-1 - dp[1][1]);
}


static void
draw_grids(void)			/* draw holodeck section grids */
{
	static BYTE	gridrgb[3] = {0x0, 0xff, 0xff};
	unsigned long  pixel;

	if (ncolors > 0)
		pixel = pixval[get_pixel(gridrgb, xnewcolr)];
	else
		pixel = true_pixel(gridrgb);
	XSetForeground(ourdisplay, ourgc, pixel);
					/* draw each grid line */
	gridlines(draw3dline);
}


static int
moveview(	/* move our view */
	int	dx,
	int	dy,
	int	mov,
	int	orb
)
{
	VIEW	nv;
	FVECT	odir, v1;
	double	d;
	register int	li;
				/* start with old view */
	nv = odev.v;
				/* change view direction */
	if (mov | orb) {
		if ((li = qtFindLeaf(dx, dy)) < 0)
			return(0);	/* not on window */
		VSUM(odir, qtL.wp[li], nv.vp, -1.);
	} else {
		if (viewray(nv.vp, nv.vdir, &odev.v,
				(dx+.5)/odev.hres, (dy+.5)/odev.vres) < -FTINY)
			return(0);	/* outside view */
	}
	if (orb && mov) {		/* orbit left/right */
		spinvector(odir, odir, nv.vup, d=MOVDEG*PI/180.*mov);
		VSUM(nv.vp, qtL.wp[li], odir, -1.);
		spinvector(nv.vdir, nv.vdir, nv.vup, d);
	} else if (orb) {		/* orbit up/down */
		fcross(v1, odir, nv.vup);
		if (normalize(v1) == 0.)
			return(0);
		spinvector(odir, odir, v1, d=MOVDEG*PI/180.*orb);
		VSUM(nv.vp, qtL.wp[li], odir, -1.);
		spinvector(nv.vdir, nv.vdir, v1, d);
	} else if (mov) {		/* move forward/backward */
		d = MOVPCT/100. * mov;
		VSUM(nv.vp, nv.vp, odir, d);
	}
	if (!mov ^ !orb && headlocked) {	/* restore head height */
		VSUM(v1, odev.v.vp, nv.vp, -1.);
		d = DOT(v1, odev.v.vup);
		VSUM(nv.vp, nv.vp, odev.v.vup, d);
	}
	if (setview(&nv) != NULL)
		return(0);	/* illegal view */
	dev_view(&nv);
	inpresflags |= DFL(DC_SETVIEW);
	return(1);
}


static void
getframe(				/* get focus frame */
	XButtonPressedEvent	*ebut
)
{
	int	startx = ebut->x, starty = ebut->y;
	int	endx, endy, midx, midy;
	FVECT	v1;
	int	li;
					/* get mouse drag */
	XMaskEvent(ourdisplay, ButtonReleaseMask, levptr(XEvent));
	endx = levptr(XButtonReleasedEvent)->x;
	endy = levptr(XButtonReleasedEvent)->y;
	midx = (startx + endx) >> 1;
	midy = (starty + endy) >> 1;
					/* set focus distance */
	if ((li = qtFindLeaf(midx, midy)) < 0)
		return;			/* not on window */
	VCOPY(v1, qtL.wp[li]);
	odev.v.vdist = sqrt(dist2(odev.v.vp, v1));
					/* set frame for rendering */
	if ((endx == startx) | (endy == starty))
		return;
	if (endx < startx) {register int c = endx; endx = startx; startx = c;}
	if (endy < starty) {register int c = endy; endy = starty; starty = c;}
	sprintf(odev_args, "%.3f %.3f %.3f %.3f",
			(startx+.5)/odev.hres, 1.-(endy+.5)/odev.vres,
			(endx+.5)/odev.hres, 1.-(starty+.5)/odev.vres);
	inpresflags |= DFL(DC_FOCUS);
}


static void
waitabit(void)				/* pause a moment */
{
	struct timespec	ts;
	ts.tv_sec = 0;
	ts.tv_nsec = 100000000;
	nanosleep(&ts, NULL);
}


static void
getmove(				/* get view change */
	XButtonPressedEvent	*ebut
)
{
	int	movdir = MOVDIR(ebut->button);
	int	movorb = MOVORB(ebut->state);
	int	oldnodesiz = qtMinNodesiz;
	Window	rootw, childw;
	int	rootx, rooty, wx, wy;
	unsigned int	statemask;

	qtMinNodesiz = 16;		/* for quicker update */
	XNoOp(ourdisplay);

	while (!XCheckMaskEvent(ourdisplay,
			ButtonReleaseMask, levptr(XEvent))) {
		waitabit();
		if (!XQueryPointer(ourdisplay, gwind, &rootw, &childw,
				&rootx, &rooty, &wx, &wy, &statemask))
			break;		/* on another screen */

		if (!moveview(wx, odev.vres-1-wy, movdir, movorb)) {
			sleep(1);
			continue;
		}
		XClearWindow(ourdisplay, gwind);
		qtUpdate();
		draw_grids();
	}
	if (!(inpresflags & DFL(DC_SETVIEW))) {	/* do final motion */
		movdir = MOVDIR(levptr(XButtonReleasedEvent)->button);
		wx = levptr(XButtonReleasedEvent)->x;
		wy = levptr(XButtonReleasedEvent)->y;
		moveview(wx, odev.vres-1-wy, movdir, movorb);
	}
	dev_flush();

	qtMinNodesiz = oldnodesiz;	/* restore quadtree resolution */
}


static void
getkey(				/* get input key */
	register XKeyPressedEvent  *ekey
)
{
	Window	rootw, childw;
	int	rootx, rooty, wx, wy;
	unsigned int	statemask;
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
	case 'l':			/* retrieve last view */
		inpresflags |= DFL(DC_LASTVIEW);
		return;
	case 'f':			/* frame view position */
		if (!XQueryPointer(ourdisplay, gwind, &rootw, &childw,
				&rootx, &rooty, &wx, &wy, &statemask))
			return;		/* on another screen */
		sprintf(odev_args, "%.4f %.4f", (wx+.5)/odev.hres,
				1.-(wy+.5)/odev.vres);
		inpresflags |= DFL(DC_FOCUS);
		return;
	case 'F':			/* unfocus */
		odev_args[0] = '\0';
		inpresflags |= DFL(DC_FOCUS);
		return;
	case 'p':			/* pause computation */
		inpresflags |= DFL(DC_PAUSE);
		return;
	case 'v':			/* spit out view */
		inpresflags |= DFL(DC_GETVIEW);
		return;
	case '\n':
	case '\r':			/* resume computation */
		inpresflags |= DFL(DC_RESUME);
		return;
	case CTRL('R'):			/* redraw screen */
		if (ncolors > 0)
			new_ctab(ncolors);
		qtRedraw(0, 0, odev.hres, odev.vres);
		return;
	case CTRL('L'):			/* refresh from server */
		if (inpresflags & DFL(DC_REDRAW))
			return;
		XClearWindow(ourdisplay, gwind);
		draw_grids();
		XFlush(ourdisplay);
		qtCompost(100);			/* unload the old tree */
		if (ncolors > 0)
			new_ctab(ncolors);
		inpresflags |= DFL(DC_REDRAW);	/* resend values from server */
		rayqleft = 0;			/* hold off update */
		return;
	case 'K':			/* kill rtrace process(es) */
		inpresflags |= DFL(DC_KILL);
		break;
	case 'R':			/* restart rtrace */
		inpresflags |= DFL(DC_RESTART);
		break;
	case 'C':			/* clobber holodeck */
		inpresflags |= DFL(DC_CLOBBER);
		break;
	case 'q':			/* quit the program */
		inpresflags |= DFL(DC_QUIT);
		return;
	default:
		XBell(ourdisplay, 0);
		return;
	}
}


static void
fixwindow(				/* repair damage to window */
	register XExposeEvent  *eexp
)
{
	if (odev.hres == 0 || odev.vres == 0) {	/* first exposure */
		odev.hres = eexp->width;
		odev.vres = eexp->height;
	}
	qtRedraw(eexp->x, odev.vres - eexp->y - eexp->height,
			eexp->x + eexp->width, odev.vres - eexp->y);
}


static void
resizewindow(			/* resize window */
	register XConfigureEvent  *ersz
)
{
	if (ersz->width == odev.hres && ersz->height == odev.vres)
		return;

	odev.hres = ersz->width;
	odev.vres = ersz->height;

	odev.v.horiz = 2.*180./PI * atan(0.5/VIEWDIST*pwidth*odev.hres);
	odev.v.vert = 2.*180./PI * atan(0.5/VIEWDIST*pheight*odev.vres);

	inpresflags |= DFL(DC_SETVIEW);
}
