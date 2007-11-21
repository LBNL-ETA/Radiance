#ifndef lint
static const char	RCSid[] = "$Id: x11.c,v 2.33 2007/11/21 16:42:03 greg Exp $";
#endif
/*
 *  x11.c - driver for X-windows version 11
 */

#include "copyright.h"

#include  "standard.h"
#include  <sys/ioctl.h>
#if  !defined(FNDELAY) && defined(O_NONBLOCK)
#define  FNDELAY  O_NONBLOCK
#endif

#include  <X11/Xlib.h>
#include  <X11/cursorfont.h>
#include  <X11/Xutil.h>

#include  "platform.h"
#include  "color.h"
#include  "driver.h"
#include  "x11twind.h"
#include  "x11icon.h"

#define GAMMA		2.2		/* default exponent correction */

#define MINWIDTH	(32*COMCW)	/* minimum graphics window width */
#define MINHEIGHT	(MINWIDTH/2)	/* minimum graphics window height */

#define BORWIDTH	5		/* border width */
#define COMHEIGHT	(COMLH*COMCH)	/* command line height (pixels) */

#define COMFN		"8x13"		/* command line font name */
#define COMLH		3		/* number of command lines */
#define COMCW		8		/* approx. character width (pixels) */
#define COMCH		14		/* approx. character height (pixels) */

#define  ourscreen	DefaultScreen(ourdisplay)
#define  ourroot	RootWindow(ourdisplay,ourscreen)

#define  levptr(etype)	((etype *)&currentevent)

static XEvent  currentevent;		/* current event */

static int  ncolors = 0;		/* color table size */
static int  mapped = 0;			/* window is mapped? */
static unsigned long  *pixval = NULL;	/* allocated pixels */
static unsigned long  ourblack=0, ourwhite=1;

static Display  *ourdisplay = NULL;	/* our display */

static XVisualInfo  ourvinfo;		/* our visual information */

static Window  gwind = 0;		/* our graphics window */

static Cursor  pickcursor = 0;		/* cursor used for picking */

static int  gwidth, gheight;		/* graphics window size */

static int  comheight;			/* desired comline height */
static TEXTWIND  *comline = NULL;	/* our command line */

static char  c_queue[64];		/* input queue */
static int  c_first = 0;		/* first character in queue */
static int  c_last = 0;			/* last character in queue */

static GC  ourgc = 0;			/* our graphics context for drawing */

static Colormap ourmap = 0;		/* our color map */

#define IC_X11		0
#define IC_IOCTL	1
#define IC_READ		2

static int  inpcheck;			/* whence to check input */

static void x11_errout(char  *msg);

static dr_closef_t x11_close;
static dr_clearf_t x11_clear;
static dr_paintrf_t x11_paintr;
static dr_getcurf_t x11_getcur;
static dr_comoutf_t x11_comout;
static dr_cominf_t x11_comin;
static dr_flushf_t x11_flush;

static dr_cominf_t std_comin;
static dr_comoutf_t std_comout;

static struct driver  x11_driver = {
	x11_close, x11_clear, x11_paintr, x11_getcur,
	NULL, NULL, x11_flush, 1.0
};

static dr_getchf_t x11_getc;

static void freepixels(void);
static int getpixels(void);
static dr_newcolrf_t xnewcolr;
static unsigned long true_pixel(COLOR  col);
static void getevent(void);
static void getkey(XKeyPressedEvent  *ekey);
static void fixwindow(XExposeEvent  *eexp);
static void resizewindow(XConfigureEvent  *ersz);

extern dr_initf_t x11_init; /* XXX this should be in a seperate header file */


extern struct driver *
x11_init(		/* initialize driver */
	char  *name,
	char  *id
)
{
	char  *gv;
	int  nplanes;
	XSetWindowAttributes	ourwinattr;
	XWMHints  ourxwmhints;
	XSizeHints	oursizhints;
					/* open display server */
	ourdisplay = XOpenDisplay(NULL);
	if (ourdisplay == NULL) {
		eputs("cannot open X-windows; DISPLAY variable set?\n");
		return(NULL);
	}
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
		if (nplanes < 4) {
			eputs("not enough colors\n");
			return(NULL);
		}
		if (!XMatchVisualInfo(ourdisplay,ourscreen,
					nplanes,PseudoColor,&ourvinfo) &&
				!XMatchVisualInfo(ourdisplay,ourscreen,
					nplanes,GrayScale,&ourvinfo)) {
			eputs("unsupported visual type\n");
			return(NULL);
		}
		ourblack = BlackPixel(ourdisplay,ourscreen);
		ourwhite = WhitePixel(ourdisplay,ourscreen);
	}
					/* set gamma */
	if ((gv = XGetDefault(ourdisplay, "radiance", "gamma")) != NULL
			|| (gv = getenv("DISPLAY_GAMMA")) != NULL)
		make_gmap(atof(gv));
	else
		make_gmap(GAMMA);
					/* X11 command line or no? */
	if (!strcmp(name, "x11"))
		comheight = COMHEIGHT;
	else /* "x11d" */ {
		comheight = 0;
#ifndef  FNDELAY
		eputs("warning: x11d driver not fully functional on this machine\n");
#endif
	}
					/* open window */
	ourwinattr.background_pixel = ourblack;
	ourwinattr.border_pixel = ourblack;
					/* this is stupid */
	ourwinattr.colormap = XCreateColormap(ourdisplay, ourroot,
				ourvinfo.visual, AllocNone);
	gwind = XCreateWindow(ourdisplay, ourroot, 0, 0,
		DisplayWidth(ourdisplay,ourscreen)-2*BORWIDTH,
		DisplayHeight(ourdisplay,ourscreen)-2*BORWIDTH,
		BORWIDTH, ourvinfo.depth, InputOutput, ourvinfo.visual,
		CWBackPixel|CWBorderPixel|CWColormap, &ourwinattr);
	if (gwind == 0) {
		eputs("cannot create window\n");
		return(NULL);
	}
   	XStoreName(ourdisplay, gwind, id);
	/* create a cursor */
	pickcursor = XCreateFontCursor(ourdisplay, XC_diamond_cross);
	ourgc = XCreateGC(ourdisplay, gwind, 0, NULL);
	ourxwmhints.flags = InputHint|IconPixmapHint;
	ourxwmhints.input = True;
	ourxwmhints.icon_pixmap = XCreateBitmapFromData(ourdisplay,
			gwind, x11icon_bits, x11icon_width, x11icon_height);
	XSetWMHints(ourdisplay, gwind, &ourxwmhints);
	oursizhints.min_width = MINWIDTH;
	oursizhints.min_height = MINHEIGHT+comheight;
	oursizhints.flags = PMinSize;
	XSetNormalHints(ourdisplay, gwind, &oursizhints);
	XSelectInput(ourdisplay, gwind, ExposureMask);
	XMapWindow(ourdisplay, gwind);
	XWindowEvent(ourdisplay, gwind, ExposureMask, levptr(XEvent));
	gwidth = levptr(XExposeEvent)->width;
	gheight = levptr(XExposeEvent)->height - comheight;
	x11_driver.xsiz = gwidth < MINWIDTH ? MINWIDTH : gwidth;
	x11_driver.ysiz = gheight < MINHEIGHT ? MINHEIGHT : gheight;
	x11_driver.inpready = 0;
	mapped = 1;
					/* set i/o vectors */
	if (comheight) {
		x11_driver.comin = x11_comin;
		x11_driver.comout = x11_comout;
		erract[COMMAND].pf = x11_comout;
		if (erract[WARNING].pf != NULL)
			erract[WARNING].pf = x11_errout;
		inpcheck = IC_X11;
	} else {
		x11_driver.comin = std_comin;
		x11_driver.comout = std_comout;
		erract[COMMAND].pf = std_comout;
		inpcheck = IC_IOCTL;
	}
	return(&x11_driver);
}


static void
x11_close(void)			/* close our display */
{
	erract[COMMAND].pf = NULL;		/* reset error vectors */
	if (erract[WARNING].pf != NULL)
		erract[WARNING].pf = wputs;
	if (ourdisplay == NULL)
		return;
	if (comline != NULL) {
		xt_close(comline);
		comline = NULL;
	}
	freepixels();
	XFreeGC(ourdisplay, ourgc);
	XDestroyWindow(ourdisplay, gwind);
	gwind = 0;
	ourgc = 0;
	XFreeCursor(ourdisplay, pickcursor);
	XCloseDisplay(ourdisplay);
	ourdisplay = NULL;
}


static void
x11_clear(			/* clear our display */
	int  xres,
	int  yres
)
{
						/* check limits */
	if (xres < MINWIDTH)
		xres = MINWIDTH;
	if (yres < MINHEIGHT)
		yres = MINHEIGHT;
						/* resize window */
	if (xres != gwidth || yres != gheight) {
		XSelectInput(ourdisplay, gwind, 0);
		XResizeWindow(ourdisplay, gwind, xres, yres+comheight);
		gwidth = xres;
		gheight = yres;
		XFlush(ourdisplay);
		sleep(2);			/* wait for window manager */
		XSync(ourdisplay, 1);		/* discard input */
	}
	XClearWindow(ourdisplay, gwind);
						/* reinitialize color table */
	if (ourvinfo.class == PseudoColor || ourvinfo.class == GrayScale) {
		if (getpixels() == 0)
			eputs("cannot allocate colors\n");
		else
			new_ctab(ncolors);
	}
						/* get new command line */
	if (comline != NULL)
		xt_close(comline);
	if (comheight) {
		comline = xt_open(ourdisplay, gwind, 0, gheight, gwidth,
				comheight, 0, ourblack, ourwhite, COMFN);
		if (comline == NULL) {
			eputs("cannot open command line window\n");
			quit(1);
		}
		XSelectInput(ourdisplay, comline->w, ExposureMask);
						/* remove earmuffs */
		XSelectInput(ourdisplay, gwind,
		StructureNotifyMask|ExposureMask|KeyPressMask|ButtonPressMask);
	} else					/* remove earmuffs */
		XSelectInput(ourdisplay, gwind,
			StructureNotifyMask|ExposureMask|ButtonPressMask);
}


static void
x11_paintr(		/* fill a rectangle */
	COLOR  col,
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
		pixel = pixval[get_pixel(col, xnewcolr)];
	else
		pixel = true_pixel(col);
	XSetForeground(ourdisplay, ourgc, pixel);
	XFillRectangle(ourdisplay, gwind, 
		ourgc, xmin, gheight-ymax, xmax-xmin, ymax-ymin);
}


static void
x11_flush(void)			/* flush output */
{
	char	buf[256];
	int	n;
						/* check for input */
	XNoOp(ourdisplay);
	n = XPending(ourdisplay);			/* from X server */
	while (n-- > 0)
		getevent();
#ifdef FNDELAY
	if (inpcheck == IC_IOCTL) {			/* from stdin */
#ifdef FIONREAD
		if (ioctl(fileno(stdin), FIONREAD, &n) < 0) {
#else
		if (1) {
#endif
			if (fcntl(fileno(stdin), F_SETFL, FNDELAY) < 0) {
				eputs("cannot change input mode\n");
				quit(1);
			}
			inpcheck = IC_READ;
		} else
			x11_driver.inpready += n;
	}
	if (inpcheck == IC_READ) {
		n = read(fileno(stdin), buf, sizeof(buf)-1);
		if (n > 0) {
			buf[n] = '\0';
			tocombuf(buf, &x11_driver);
		}
	}
#endif
}


static void
x11_comin(		/* read in a command line */
	char  *inp,
	char  *prompt
)
{
	if (prompt != NULL) {
		x11_flush();		/* make sure we get everything */
		if (fromcombuf(inp, &x11_driver))
			return;
		xt_puts(prompt, comline);
	}
	xt_cursor(comline, TBLKCURS);
	editline(inp, x11_getc, x11_comout);
	xt_cursor(comline, TNOCURS);
}


static void
x11_comout(		/* output a string to command line */
	char  *outp
)
{
	if (comline == NULL || outp == NULL || !outp[0])
		return;
	xt_puts(outp, comline);
	if (outp[strlen(outp)-1] == '\n')
		XFlush(ourdisplay);
}


static void
x11_errout(			/* output an error message */
	char  *msg
)
{
	eputs(msg);		/* send to stderr also! */
	x11_comout(msg);
}


static void
std_comin(		/* read in command line from stdin */
	char  *inp,
	char  *prompt
)
{
	if (prompt != NULL) {
		if (fromcombuf(inp, &x11_driver))
			return;
		if (!x11_driver.inpready)
			std_comout(prompt);
	}
#ifdef FNDELAY
	if (inpcheck == IC_READ) {	/* turn off FNDELAY */
		if (fcntl(fileno(stdin), F_SETFL, 0) < 0) {
			eputs("cannot change input mode\n");
			quit(1);
		}
		inpcheck = IC_IOCTL;
	}
#endif
	if (gets(inp) == NULL) {
		strcpy(inp, "quit");
		return;
	}
	x11_driver.inpready -= strlen(inp) + 1;
	if (x11_driver.inpready < 0)
		x11_driver.inpready = 0;
}


static void
std_comout(		/* write out string to stdout */
	char	*outp
)
{
	fputs(outp, stdout);
	fflush(stdout);
}


static int
x11_getcur(		/* get cursor position */
	int  *xp,
	int  *yp
)
{
	while (XGrabPointer(ourdisplay, gwind, True, ButtonPressMask,
			GrabModeAsync, GrabModeAsync, None, pickcursor,
			CurrentTime) != GrabSuccess)
		sleep(2);

	do
		getevent();
	while (c_last <= c_first && levptr(XEvent)->type != ButtonPress);
	*xp = levptr(XButtonPressedEvent)->x;
	*yp = gheight-1 - levptr(XButtonPressedEvent)->y;
	XUngrabPointer(ourdisplay, CurrentTime);
	XFlush(ourdisplay);				/* insure release */
	if (c_last > c_first)			/* key pressed */
		return(x11_getc());
						/* button pressed */
	if (levptr(XButtonPressedEvent)->button == Button1)
		return(MB1);
	if (levptr(XButtonPressedEvent)->button == Button2)
		return(MB2);
	if (levptr(XButtonPressedEvent)->button == Button3)
		return(MB3);
	return(ABORT);
}


static void
xnewcolr(		/* enter a color into hardware table */
	int  ndx,
	int  r,
	int  g,
	int  b
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
	COLOR  col
)
{
	unsigned long  rval;
	BYTE  rgb[3];

	map_color(rgb, col);
	rval = ourvinfo.red_mask*rgb[RED]/255 & ourvinfo.red_mask;
	rval |= ourvinfo.green_mask*rgb[GRN]/255 & ourvinfo.green_mask;
	rval |= ourvinfo.blue_mask*rgb[BLU]/255 & ourvinfo.blue_mask;
	return(rval);
}


static int
x11_getc(void)			/* get a command character */
{
	while (c_last <= c_first) {
		c_first = c_last = 0;		/* reset */
		getevent();			/* wait for key */
	}
	x11_driver.inpready--;
	return(c_queue[c_first++]);
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
				eputs("cannot allocate colors\n");
			else
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
		break;
	}
}


static void
getkey(				/* get input key */
	register XKeyPressedEvent  *ekey
)
{
	register int  n;

	n = XLookupString(ekey, c_queue+c_last, sizeof(c_queue)-c_last,
				NULL, NULL);
	c_last += n;
	x11_driver.inpready += n;
}


static void
fixwindow(				/* repair damage to window */
	register XExposeEvent  *eexp
)
{
	char  buf[80];

	if (eexp->window == gwind) {
		sprintf(buf, "repaint %d %d %d %d\n",
			eexp->x, gheight - eexp->y - eexp->height,
			eexp->x + eexp->width, gheight - eexp->y);
		tocombuf(buf, &x11_driver);
	} else if (eexp->window == comline->w) {
		if (eexp->count == 0)
			xt_redraw(comline);
	}
}


static void
resizewindow(			/* resize window */
	register XConfigureEvent  *ersz
)
{
	if (ersz->width == gwidth && ersz->height-comheight == gheight)
		return;

	gwidth = ersz->width;
	gheight = ersz->height-comheight;
	x11_driver.xsiz = gwidth < MINWIDTH ? MINWIDTH : gwidth;
	x11_driver.ysiz = gheight < MINHEIGHT ? MINHEIGHT : gheight;

	tocombuf("new\n", &x11_driver);
}
