#ifndef lint
static const char RCSid[] = "$Id";
#endif
/*
 *  x11image.c - driver for X-windows
 *
 *     3/1/90
 */

/*
 *  Modified for X11
 *
 *  January 1990
 *
 *  Anat Grynberg  and Greg Ward
 */


#include  "standard.h"

#include  <string.h>
#include  <signal.h>
#include  <unistd.h>
#include  <sys/types.h>
#include  <wait.h>
#include  <X11/Xlib.h>
#include  <X11/cursorfont.h>
#include  <X11/Xutil.h>
#include  <X11/Xatom.h>

#include  "color.h"
#include  "tonemap.h"
#include  "clrtab.h"
#include  "view.h"
#include  "x11raster.h"
#include  "random.h"

#define  FONTNAME	"8x13"		/* text font we'll use */

#define  CTRL(c)	((c)-'@')

#define  BORWIDTH	5		/* border width */

#define  ICONSIZ	(8*10)		/* maximum icon dimension (even 8) */

#define  FIXWEIGHT	20		/* weight to add for fixation points */

#define  ourscreen	DefaultScreen(thedisplay)
#define  ourroot	RootWindow(thedisplay,ourscreen)

#define  revline(x0,y0,x1,y1)	XDrawLine(thedisplay,wind,revgc,x0,y0,x1,y1)

#define  redraw(x,y,w,h) patch_raster(wind,(x)-xoff,(y)-yoff,x,y,w,h,ourras)

double  gamcor = 2.2;			/* gamma correction */
char  *gamstr = NULL;			/* gamma value override */

int  dither = 1;			/* dither colors? */
int  fast = 0;				/* keep picture in Pixmap? */

char	*dispname = NULL;		/* our display name */

Window  wind = 0;			/* our output window */
unsigned long  ourblack=0, ourwhite=1;	/* black and white for this visual */
int  maxcolors = 0;			/* maximum colors */
int  greyscale = 0;			/* in grey */

int  scale = 0;				/* scalefactor; power of two */

int  xoff = 0;				/* x image offset */
int  yoff = 0;				/* y image offset */

int  parent = 0;			/* number of children, -1 if child */
int  sequential = 0;			/* display images in sequence */

char  *tout = "od";			/* output of 't' command */
int  tinterv = 0;			/* interval between mouse reports */

int  tmflags = TM_F_LINEAR;		/* tone mapping flags */

VIEW  ourview = STDVIEW;		/* image view parameters */
int  gotview = 0;			/* got parameters from file */

COLR  *scanline;			/* scan line buffer */
TMbright  *lscan;			/* encoded luminance scanline */
BYTE  *cscan;				/* encoded chroma scanline */
BYTE  *pscan;				/* compute pixel scanline */

RESOLU  inpres;				/* input resolution and ordering */
int  xmax, ymax;			/* picture dimensions */
int  width, height;			/* window size */
char  *fname = NULL;			/* input file name */
FILE  *fin = NULL;			/* input file */
long  *scanpos = NULL;			/* scan line positions in file */
int  cury = 0;				/* current scan location */

double  exposure = 1.0;			/* exposure compensation used */

int  wrongformat = 0;			/* input in another format? */

GC	ourgc;				/* standard graphics context */
GC	revgc;				/* graphics context with GXinvert */

int		*ourrank;		/* our visual class ranking */
XVisualInfo	ourvis;			/* our visual */
XRASTER		*ourras;		/* our stored image */
unsigned char	*ourdata;		/* our image data */

struct {
	int  xmin, ymin, xsiz, ysiz;
}  bbox = {0, 0, 0, 0};			/* current bbox */

char  *geometry = NULL;			/* geometry specification */

char  icondata[ICONSIZ*ICONSIZ/8];	/* icon bitmap data */
int  iconwidth = 0, iconheight = 0;

char  *progname;

char  errmsg[128];

BYTE  clrtab[256][3];			/* global color map */


Display  *thedisplay;
Atom  closedownAtom, wmProtocolsAtom;

int  sigrecv;

void  onsig(int i) { sigrecv++; }

typedef void doboxf_t(COLR *scn, int n, void *);

static gethfunc headline;
static void init(int argc, char **argv);
static void quiterr(char *err);
static int viscmp(XVisualInfo *v1, XVisualInfo *v2);
static void getbestvis(void);
static void getras(void);
static void getevent(void);
static int traceray(int xpos, int ypos);
static int docom(XKeyPressedEvent *ekey);
static void moveimage(XButtonPressedEvent *ebut);
static void getbox(XButtonPressedEvent *ebut);
static void trackrays(XButtonPressedEvent *ebut);
static void revbox(int x0, int y0, int x1, int y1);
static doboxf_t colavg;
static doboxf_t addfix;
static int avgbox(COLOR cavg);
static int dobox(doboxf_t *f, void *p);
static void make_tonemap(void);
static void tmap_colrs(COLR *scn, int len);
static void getmono(void);
static void add2icon(int y, COLR *scan);
static void getfull(void);
static void getgrey(void);
static void getmapped(void);
static void scale_rcolors(XRASTER *xr, double sf);
static int getscan(int y);


int
main(int  argc, char  *argv[])
{
	int  i;
	int  pid;
	
	progname = argv[0];
	fin = stdin;

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'c':			/* number of colors */
				maxcolors = atoi(argv[++i]);
				break;
			case 'b':			/* greyscale only */
				greyscale = !greyscale;
				break;
			case 'm':			/* monochrome */
				greyscale = 1;
				maxcolors = 2;
				break;
			case 'd':			/* display or dither */
				if (argv[i][2] == 'i')
					dispname = argv[++i];
				else
					dither = !dither;
				break;
			case 'f':			/* save pixmap */
				fast = !fast;
				break;
			case 's':			/* one at a time */
				sequential = !sequential;
				break;
			case 'o':			/* 't' output */
				tout = argv[i]+2;
				break;
			case 't':			/* msec interval */
				tinterv = atoi(argv[++i]);
				break;
			case 'e':			/* exposure comp. */
				i++;
				if (argv[i][0] == 'a') {
					tmflags = TM_F_CAMERA;
					break;
				}
				if (argv[i][0] == 'h') {
					tmflags = TM_F_HUMAN;
					break;
				}
				if (argv[i][0] != '+' && argv[i][0] != '-')
					goto userr;
				scale = atoi(argv[i]);
				break;
			case 'g':			/* gamma comp. */
				if (argv[i][2] == 'e')
					geometry = argv[++i];
				else
					gamstr = argv[++i];
				break;
			default:
				goto userr;
			}
		else if (argv[i][0] == '=')
			geometry = argv[i];
		else
			break;

	if (i > argc)
		goto userr;
	while (i < argc-1) {
		sigrecv = 0;
		signal(SIGCONT, onsig);
		if ((pid=fork()) == 0) {	/* a child for each picture */
			parent = -1;
			break;
		}
		if (pid < 0)
			quiterr("fork failed");
		parent++;
		while (!sigrecv)
			pause();	/* wait for wake-up call */
		i++;
	}
	if (i < argc) {			/* open picture file */
		fname = argv[i];
		fin = fopen(fname, "r");
		if (fin == NULL)
			quiterr("cannot open picture file");
	}
				/* get header */
	getheader(fin, headline, NULL);
				/* get picture dimensions */
	if (wrongformat || !fgetsresolu(&inpres, fin))
		quiterr("bad picture format");
	xmax = scanlen(&inpres);
	ymax = numscans(&inpres);
				/* set view parameters */
	if (gotview && setview(&ourview) != NULL)
		gotview = 0;
	if ((scanline = (COLR *)malloc(xmax*sizeof(COLR))) == NULL)
		quiterr("out of memory");

	init(argc, argv);			/* get file and open window */

	for ( ; ; )
		getevent();		/* main loop */
userr:
	fprintf(stderr,
"Usage: %s [-di disp][[-ge] spec][-b][-m][-d][-f][-c nclrs][-e spec][-g gamcor][-s][-ospec][-t intvl] pic ..\n",
			progname);
	exit(1);
}


static int
headline(		/* get relevant info from header */
	char	*s,
	void	*p
)
{
	char  fmt[32];

	if (isexpos(s))
		exposure *= exposval(s);
	else if (formatval(fmt, s))
		wrongformat = strcmp(fmt, COLRFMT);
	else if (isview(s) && sscanview(&ourview, s) > 0)
		gotview++;
	return(0);
}


static void
init(			/* get data and open window */
	int argc,
	char **argv
)
{
	XSetWindowAttributes	ourwinattr;
	XClassHint	xclshints;
	XWMHints	xwmhints;
	XSizeHints	xszhints;
	XTextProperty	windowName, iconName;
	XGCValues	xgcv;
	char	*name;
	register int	i;
	
	if (fname != NULL) {
		scanpos = (long *)malloc(ymax*sizeof(long));
		if (scanpos == NULL)
			quiterr("out of memory");
		for (i = 0; i < ymax; i++)
			scanpos[i] = -1;
		name = fname;
	} else
		name = progname;
				/* remove directory prefix from name */
	for (i = strlen(name); i-- > 0; )
		if (name[i] == '/')
			break;
	name += i+1;
	if ((thedisplay = XOpenDisplay(dispname)) == NULL)
		quiterr("cannot open display");
				/* set gamma value */
	if (gamstr == NULL)		/* get it from the X server */
		gamstr = XGetDefault(thedisplay, "radiance", "gamma");
	if (gamstr == NULL)		/* get it from the environment */
		gamstr = getenv("DISPLAY_GAMMA");
	if (gamstr != NULL)
		gamcor = atof(gamstr);
				/* get best visual for default screen */
	getbestvis();
				/* store image */
	getras();
				/* get size and position */
	xszhints.flags = 0;
	xszhints.width = xmax; xszhints.height = ymax;
	if (geometry != NULL) {
		i = XParseGeometry(geometry, &xszhints.x, &xszhints.y,
				(unsigned *)&xszhints.width,
				(unsigned *)&xszhints.height);
		if ((i&(WidthValue|HeightValue)) == (WidthValue|HeightValue))
			xszhints.flags |= USSize;
		else
			xszhints.flags |= PSize;
		if ((i&(XValue|YValue)) == (XValue|YValue)) {
			xszhints.flags |= USPosition;
			if (i & XNegative)
				xszhints.x += DisplayWidth(thedisplay,
				ourscreen)-1-xszhints.width-2*BORWIDTH;
			if (i & YNegative)
				xszhints.y += DisplayHeight(thedisplay,
				ourscreen)-1-xszhints.height-2*BORWIDTH;
		}
	}
	/* open window */
	i = CWEventMask|CWCursor|CWBackPixel|CWBorderPixel;
	ourwinattr.border_pixel = ourwhite;
	ourwinattr.background_pixel = ourblack;
	if (ourvis.visual != DefaultVisual(thedisplay,ourscreen)) {
		ourwinattr.colormap = newcmap(thedisplay, ourscreen, ourvis.visual);
		i |= CWColormap;
	}
	ourwinattr.event_mask = ExposureMask|KeyPressMask|ButtonPressMask|
			ButtonReleaseMask|ButtonMotionMask|StructureNotifyMask;
	ourwinattr.cursor = XCreateFontCursor(thedisplay, XC_diamond_cross);
	wind = XCreateWindow(thedisplay, ourroot, xszhints.x, xszhints.y,
			xszhints.width, xszhints.height, BORWIDTH,
			ourvis.depth, InputOutput, ourvis.visual,
			i, &ourwinattr);
	if (wind == 0)
		quiterr("cannot create window");
	width = xmax;
	height = ymax;
	/* prepare graphics drawing context */
	if ((xgcv.font = XLoadFont(thedisplay, FONTNAME)) == 0)
		quiterr("cannot get font");
	xgcv.foreground = ourblack;
	xgcv.background = ourwhite;
	ourgc = XCreateGC(thedisplay, wind, GCForeground|GCBackground|
			GCFont, &xgcv);
	xgcv.function = GXinvert;
	revgc = XCreateGC(thedisplay, wind, GCForeground|GCBackground|
			GCFunction, &xgcv);

	/* set up the window manager */
	xwmhints.flags = InputHint|IconPixmapHint;
	xwmhints.input = True;
	xwmhints.icon_pixmap = XCreateBitmapFromData(thedisplay,
			wind, icondata, iconwidth, iconheight);

   	windowName.encoding = iconName.encoding = XA_STRING;
	windowName.format = iconName.format = 8;
	windowName.value = (u_char *)name;
	windowName.nitems = strlen(windowName.value);
	iconName.value = (u_char *)name;
	iconName.nitems = strlen(windowName.value);

	xclshints.res_name = NULL;
	xclshints.res_class = "Ximage";
	XSetWMProperties(thedisplay, wind, &windowName, &iconName,
			argv, argc, &xszhints, &xwmhints, &xclshints);
	closedownAtom = XInternAtom(thedisplay, "WM_DELETE_WINDOW", False);
	wmProtocolsAtom = XInternAtom(thedisplay, "WM_PROTOCOLS", False);
	XSetWMProtocols(thedisplay, wind, &closedownAtom, 1);

	XMapWindow(thedisplay, wind);
} /* end of init */


static void
quiterr(		/* print message and exit */
	char  *err
)
{
	register int  es;
	int  cs;

	if ( (es = err != NULL) )
		fprintf(stderr, "%s: %s: %s\n", progname, 
				fname==NULL?"<stdin>":fname, err);
	if (thedisplay != NULL)
		XCloseDisplay(thedisplay);
	if ((parent < 0) & (sigrecv == 0))
		kill(getppid(), SIGCONT);
	while (parent > 0 && wait(&cs) != -1) {	/* wait for any children */
		if (es == 0)
			es = cs>>8 & 0xff;
		parent--;
	}
	exit(es);
}


static int
viscmp(		/* compare visual to see which is better, descending */
	register XVisualInfo	*v1,
	register XVisualInfo	*v2
)
{
	int	bad1 = 0, bad2 = 0;
	register int  *rp;

	if (v1->class == v2->class) {
		if (v1->class == TrueColor || v1->class == DirectColor) {
					/* prefer 24-bit to 32-bit */
			if (v1->depth == 24 && v2->depth == 32)
				return(-1);
			if (v1->depth == 32 && v2->depth == 24)
				return(1);
					/* go for maximum depth otherwise */
			return(v2->depth - v1->depth);
		}
					/* don't be too greedy */
		if (maxcolors <= 1<<v1->depth && maxcolors <= 1<<v2->depth)
			return(v1->depth - v2->depth);
		return(v2->depth - v1->depth);
	}
					/* prefer Pseudo when < 15-bit */
	if ((v1->class == TrueColor || v1->class == DirectColor) &&
			v1->depth < 15)
		bad1 = 1;
	if ((v2->class == TrueColor || v2->class == DirectColor) &&
			v2->depth < 15)
		bad2 = -1;
	if (bad1 | bad2)
		return(bad1+bad2);
					/* otherwise, use class ranking */
	for (rp = ourrank; *rp != -1; rp++) {
		if (v1->class == *rp)
			return(-1);
		if (v2->class == *rp)
			return(1);
	}
	return(0);
}


static void
getbestvis(void)			/* get the best visual for this screen */
{
#ifdef DEBUG
static char  vistype[][12] = {
		"StaticGray",
		"GrayScale",
		"StaticColor",
		"PseudoColor",
		"TrueColor",
		"DirectColor"
};
#endif
	static int	rankings[3][6] = {
		{TrueColor,DirectColor,PseudoColor,GrayScale,StaticGray,-1},
		{PseudoColor,GrayScale,StaticGray,-1},
		{PseudoColor,GrayScale,StaticGray,-1}
	};
	XVisualInfo	*xvi;
	int	vismatched;
	register int	i, j;

	if (greyscale) {
		ourrank = rankings[2];
		if (maxcolors < 2) maxcolors = 256;
	} else if (maxcolors >= 2 && maxcolors <= 256)
		ourrank = rankings[1];
	else {
		ourrank = rankings[0];
		maxcolors = 256;
	}
					/* find best visual */
	ourvis.screen = ourscreen;
	xvi = XGetVisualInfo(thedisplay,VisualScreenMask,&ourvis,&vismatched);
	if (xvi == NULL)
		quiterr("no visuals for this screen!");
#ifdef DEBUG
	fprintf(stderr, "Supported visuals:\n");
	for (i = 0; i < vismatched; i++)
		fprintf(stderr, "\ttype %s, depth %d\n",
				vistype[xvi[i].class], xvi[i].depth);
#endif
	for (i = 0, j = 1; j < vismatched; j++)
		if (viscmp(&xvi[i],&xvi[j]) > 0)
			i = j;
					/* compare to least acceptable */
	for (j = 0; ourrank[j++] != -1; )
		;
	ourvis.class = ourrank[--j];
	ourvis.depth = 1;
	if (viscmp(&xvi[i],&ourvis) > 0)
		quiterr("inadequate visuals on this screen");
					/* OK, we'll use it */
	ourvis = xvi[i];
#ifdef DEBUG
	fprintf(stderr, "Selected visual type %s, depth %d\n",
			vistype[ourvis.class], ourvis.depth);
#endif
					/* make appropriate adjustments */
	if (ourvis.class == GrayScale || ourvis.class == StaticGray)
		greyscale = 1;
	if (ourvis.depth <= 8 && ourvis.colormap_size < maxcolors)
		maxcolors = ourvis.colormap_size;
	if (ourvis.class == StaticGray) {
		ourblack = 0;
		ourwhite = 255;
	} else if (ourvis.class == PseudoColor) {
		ourblack = BlackPixel(thedisplay,ourscreen);
		ourwhite = WhitePixel(thedisplay,ourscreen);
		if ((ourblack|ourwhite) & ~255L) {
			ourblack = 0;
			ourwhite = 1;
		}
		if (maxcolors > 4)
			maxcolors -= 2;
	} else {
		ourblack = 0;
		ourwhite = ourvis.red_mask|ourvis.green_mask|ourvis.blue_mask;
	}
	XFree((char *)xvi);
}


static void
getras(void)				/* get raster file */
{
	if (maxcolors <= 2) {		/* monochrome */
		ourdata = (unsigned char *)malloc(ymax*((xmax+7)/8));
		if (ourdata == NULL)
			goto fail;
		ourras = make_raster(thedisplay, &ourvis, 1, ourdata,
				xmax, ymax, 8);
		if (ourras == NULL)
			goto fail;
		getmono();
	} else if ((ourvis.class == TrueColor) | (ourvis.class == DirectColor)) {
		int  datsiz = ourvis.depth>16 ? sizeof(int32) : sizeof(int16);
		ourdata = (unsigned char *)malloc(datsiz*xmax*ymax);
		if (ourdata == NULL)
			goto fail;
		ourras = make_raster(thedisplay, &ourvis, datsiz*8,
				ourdata, xmax, ymax, datsiz*8);
		if (ourras == NULL)
			goto fail;
		getfull();
	} else {
		ourdata = (unsigned char *)malloc(xmax*ymax);
		if (ourdata == NULL)
			goto fail;
		ourras = make_raster(thedisplay, &ourvis, 8, ourdata,
				xmax, ymax, 8);
		if (ourras == NULL)
			goto fail;
		if (greyscale)
			getgrey();
		else
			getmapped();
		if (ourvis.class != StaticGray && !init_rcolors(ourras,clrtab))
			goto fail;
	}
	return;
fail:
	quiterr("could not create raster image");
}


static void
getevent(void)				/* process the next event */
{
	XEvent xev;

	XNextEvent(thedisplay, &xev);
	switch ((int)xev.type) {
	case KeyPress:
		docom(&xev.xkey);
		break;
	case ConfigureNotify:
		width = xev.xconfigure.width;
		height = xev.xconfigure.height;
		break;
	case MapNotify:
		map_rcolors(ourras, wind);
		if (fast)
			make_rpixmap(ourras, wind);
		if ((!sequential) & (parent < 0) & (sigrecv == 0)) {
			kill(getppid(), SIGCONT);
			sigrecv--;
		}
		break;
	case UnmapNotify:
		if (!fast)
			unmap_rcolors(ourras);
		break;
	case Expose:
		redraw(xev.xexpose.x, xev.xexpose.y,
				xev.xexpose.width, xev.xexpose.height);
		break;
	case ButtonPress:
		if (xev.xbutton.state & (ShiftMask|ControlMask))
			moveimage(&xev.xbutton);
		else
			switch (xev.xbutton.button) {
			case Button1:
				getbox(&xev.xbutton);
				break;
			case Button2:
				traceray(xev.xbutton.x, xev.xbutton.y);
				break;
			case Button3:
				trackrays(&xev.xbutton);
				break;
			}
		break;
	case ClientMessage:
		if ((xev.xclient.message_type == wmProtocolsAtom) &&
				(xev.xclient.data.l[0] == closedownAtom))
			quiterr(NULL);
		break;
	}
}


static int
traceray(			/* print requested pixel data */
	int  xpos,
	int  ypos
)
{
	RREAL  hv[2];
	FVECT  rorg, rdir;
	COLOR  cval;
	register char  *cp;

	bbox.xmin = xpos; bbox.xsiz = 1;
	bbox.ymin = ypos; bbox.ysiz = 1;
	avgbox(cval);
	scalecolor(cval, 1./exposure);
	pix2loc(hv, &inpres, xpos-xoff, ypos-yoff);
	if (!gotview || viewray(rorg, rdir, &ourview, hv[0], hv[1]) < 0)
		rorg[0] = rorg[1] = rorg[2] =
		rdir[0] = rdir[1] = rdir[2] = 0.;

	for (cp = tout; *cp; cp++)	/* print what they asked for */
		switch (*cp) {
		case 'o':			/* origin */
			printf("%e %e %e ", rorg[0], rorg[1], rorg[2]);
			break;
		case 'd':			/* direction */
			printf("%e %e %e ", rdir[0], rdir[1], rdir[2]);
			break;
		case 'v':			/* radiance value */
			printf("%e %e %e ", colval(cval,RED),
					colval(cval,GRN), colval(cval,BLU));
			break;
		case 'l':			/* luminance */
			printf("%e ", luminance(cval));
			break;
		case 'p':			/* pixel position */
			printf("%d %d ", (int)(hv[0]*inpres.xr),
					(int)(hv[1]*inpres.yr));
			break;
		}
	putchar('\n');
	fflush(stdout);
	return(0);
}


static int
docom(				/* execute command */
	XKeyPressedEvent  *ekey
)
{
	char  buf[80];
	COLOR  cval;
	XColor  cvx;
	int  com, n;
	double  comp;
	RREAL  hv[2];

	n = XLookupString(ekey, buf, sizeof(buf), NULL, NULL); 
	if (n == 0)
		return(0);
	com = buf[0];
	switch (com) {			/* interpret command */
	case 'q':
	case 'Q':
	case CTRL('D'):				/* quit */
		quiterr(NULL);
	case '\n':
	case '\r':
	case 'l':
	case 'c':				/* value */
		if (!avgbox(cval))
			return(-1);
		switch (com) {
		case '\n':
		case '\r':				/* radiance */
			sprintf(buf, "%.3f", intens(cval)/exposure);
			break;
		case 'l':				/* luminance */
			sprintf(buf, "%.1fL", luminance(cval)/exposure);
			break;
		case 'c':				/* color */
			comp = pow(2.0, (double)scale);
			sprintf(buf, "(%.2f,%.2f,%.2f)",
					colval(cval,RED)*comp,
					colval(cval,GRN)*comp,
					colval(cval,BLU)*comp);
			break;
		}
		XDrawImageString(thedisplay, wind, ourgc,
				bbox.xmin, bbox.ymin+bbox.ysiz, buf, strlen(buf)); 
		return(0);
	case 'i':				/* identify (contour) */
		if (ourras->pixels == NULL)
			return(-1);
		n = ourdata[ekey->x-xoff+xmax*(ekey->y-yoff)];
		n = ourras->pmap[n];
		cvx.pixel = ourras->cdefs[n].pixel;
		cvx.red = random() & 65535;
		cvx.green = random() & 65535;
		cvx.blue = random() & 65535;
		cvx.flags = DoRed|DoGreen|DoBlue;
		XStoreColor(thedisplay, ourras->cmap, &cvx);
		return(0);
	case 'p':				/* position */
		pix2loc(hv, &inpres, ekey->x-xoff, ekey->y-yoff);
		sprintf(buf, "(%d,%d)", (int)(hv[0]*inpres.xr),
				(int)(hv[1]*inpres.yr));
		XDrawImageString(thedisplay, wind, ourgc, ekey->x, ekey->y,
					buf, strlen(buf));
		return(0);
	case 't':				/* trace */
		return(traceray(ekey->x, ekey->y));
	case 'a':				/* auto exposure */
		if (fname == NULL)
			return(-1);
		tmflags = TM_F_CAMERA;
		strcpy(buf, "auto exposure...");
		goto remap;
	case 'h':				/* human response */
		if (fname == NULL)
			return(-1);
		tmflags = TM_F_HUMAN;
		strcpy(buf, "human exposure...");
		goto remap;
	case '=':				/* adjust exposure */
	case '@':				/* adaptation level */
		if (!avgbox(cval))
			return(-1);
		comp = bright(cval);
		if (comp < 1e-20) {
			XBell(thedisplay, 0);
			return(-1);
		}
		if (com == '@')
			comp = 106./exposure/
			pow(1.219+pow(comp*WHTEFFICACY/exposure,.4),2.5);
		else
			comp = .5/comp;
		comp = log(comp)/.69315 - scale;
		n = comp < 0 ? comp-.5 : comp+.5 ;	/* round */
		if (tmflags != TM_F_LINEAR)
			tmflags = TM_F_LINEAR;	/* turn off tone mapping */
		else {
			if (n == 0)		/* else check if any change */
				return(0);
			scale_rcolors(ourras, pow(2.0, (double)n));
		}
		scale += n;
		sprintf(buf, "%+d", scale);
	remap:
		XDrawImageString(thedisplay, wind, ourgc,
				bbox.xmin, bbox.ymin+bbox.ysiz, buf, strlen(buf));
		XFlush(thedisplay);
		/* free(ourdata); 	This is done in XDestroyImage()! */
		free_raster(ourras);
		getras();
	/* fall through */
	case CTRL('R'):				/* redraw */
	case CTRL('L'):
		unmap_rcolors(ourras);
		XClearWindow(thedisplay, wind);
		map_rcolors(ourras, wind);
		if (fast)
			make_rpixmap(ourras, wind);
		redraw(0, 0, width, height);
		return(0);
	case 'f':				/* turn on fast redraw */
		fast = 1;
		make_rpixmap(ourras, wind);
		return(0);
	case 'F':				/* turn off fast redraw */
		fast = 0;
		free_rpixmap(ourras);
		return(0);
	case '0':				/* recenter origin */
		if ((xoff == 0) & (yoff == 0))
			return(0);
		xoff = yoff = 0;
		XClearWindow(thedisplay, wind);
		redraw(0, 0, width, height);
		return(0);
	case ' ':				/* clear */
		redraw(bbox.xmin, bbox.ymin, bbox.xsiz, bbox.ysiz);
		return(0);
	default:
		XBell(thedisplay, 0);
		return(-1);
	}
}


static void
moveimage(				/* shift the image */
	XButtonPressedEvent  *ebut
)
{
	XEvent	e;
	int	mxo, myo;

	XMaskEvent(thedisplay, ButtonReleaseMask|ButtonMotionMask, &e);
	while (e.type == MotionNotify) {
		mxo = e.xmotion.x;
		myo = e.xmotion.y;
		revline(ebut->x, ebut->y, mxo, myo);
		revbox(xoff+mxo-ebut->x, yoff+myo-ebut->y,
				xoff+mxo-ebut->x+xmax, yoff+myo-ebut->y+ymax);
		XMaskEvent(thedisplay,ButtonReleaseMask|ButtonMotionMask,&e);
		revline(ebut->x, ebut->y, mxo, myo);
		revbox(xoff+mxo-ebut->x, yoff+myo-ebut->y,
				xoff+mxo-ebut->x+xmax, yoff+myo-ebut->y+ymax);
	}
	xoff += e.xbutton.x - ebut->x;
	yoff += e.xbutton.y - ebut->y;
	XClearWindow(thedisplay, wind);
	redraw(0, 0, width, height);
}


static void
getbox(				/* get new bbox */
	XButtonPressedEvent  *ebut
)
{
	XEvent	e;

	XMaskEvent(thedisplay, ButtonReleaseMask|ButtonMotionMask, &e);
	while (e.type == MotionNotify) {
		revbox(ebut->x, ebut->y, bbox.xmin = e.xmotion.x, bbox.ymin = e.xmotion.y);
		XMaskEvent(thedisplay,ButtonReleaseMask|ButtonMotionMask,&e);
		revbox(ebut->x, ebut->y, bbox.xmin, bbox.ymin);
	}
	bbox.xmin = e.xbutton.x<0 ? 0 : (e.xbutton.x>=width ? width-1 : e.xbutton.x);
	bbox.ymin = e.xbutton.y<0 ? 0 : (e.xbutton.y>=height ? height-1 : e.xbutton.y);
	if (bbox.xmin > ebut->x) {
		bbox.xsiz = bbox.xmin - ebut->x + 1;
		bbox.xmin = ebut->x;
	} else {
		bbox.xsiz = ebut->x - bbox.xmin + 1;
	}
	if (bbox.ymin > ebut->y) {
		bbox.ysiz = bbox.ymin - ebut->y + 1;
		bbox.ymin = ebut->y;
	} else {
		bbox.ysiz = ebut->y - bbox.ymin + 1;
	}
}


static void
trackrays(				/* trace rays as mouse moves */
	XButtonPressedEvent  *ebut
)
{
	XEvent	e;
	unsigned long	lastrept;

	traceray(ebut->x, ebut->y);
	lastrept = ebut->time;
	XMaskEvent(thedisplay, ButtonReleaseMask|ButtonMotionMask, &e);
	while (e.type == MotionNotify) {
		if (e.xmotion.time >= lastrept + tinterv) {
			traceray(e.xmotion.x, e.xmotion.y);
			lastrept = e.xmotion.time;
		}
		XMaskEvent(thedisplay,ButtonReleaseMask|ButtonMotionMask,&e);
	}
}


static void
revbox(			/* draw bbox with reversed lines */
	int  x0,
	int  y0,
	int  x1,
	int  y1
)
{
	revline(x0, y0, x1, y0);
	revline(x0, y1, x1, y1);
	revline(x0, y0, x0, y1);
	revline(x1, y0, x1, y1);
}


static void
colavg(
	register COLR	*scn,
	register int	n,
	void *cavg
)
{
	COLOR	col;

	while (n--) {
		colr_color(col, *scn++);
		addcolor((COLORV*)cavg, col);
	}
}


static int
avgbox(				/* average color over current bbox */
	COLOR	cavg
)
{
	double	d;
	register int	rval;

	setcolor(cavg, 0., 0., 0.);
	rval = dobox(colavg, (char *)cavg);
	if (rval > 0) {
		d = 1./rval;
		scalecolor(cavg, d);
	}
	return(rval);
}


static int
dobox(				/* run function over bbox */
	//void	(*f)(),			/* function to call for each subscan */
	doboxf_t *f,			/* function to call for each subscan */
	void	*p			/* pointer to private data */
)
{
	int  left, right, top, bottom;
	int  y;

	left = bbox.xmin - xoff;
	right = left + bbox.xsiz;
	if (left < 0)
		left = 0;
	if (right > xmax)
		right = xmax;
	if (left >= right)
		return(0);
	top = bbox.ymin - yoff;
	bottom = top + bbox.ysiz;
	if (top < 0)
		top = 0;
	if (bottom > ymax)
		bottom = ymax;
	if (top >= bottom)
		return(0);
	for (y = top; y < bottom; y++) {
		if (getscan(y) == -1)
			return(-1);
		(*f)(scanline+left, right-left, p);
	}
	return((right-left)*(bottom-top));
}


static void
addfix(			/* add fixation points to histogram */
	COLR	*scn,
	int	n,
	void	*p  /*NOTUSED*/
)
{
	if (tmCvColrs(lscan, TM_NOCHROM, scn, n))
		goto tmerr;
	if (tmAddHisto(lscan, n, FIXWEIGHT))
		goto tmerr;
	return;
tmerr:
	quiterr("tone mapping error");
}


static void
make_tonemap(void)			/* initialize tone mapping */
{
	int  flags, y;

	if (tmflags != TM_F_LINEAR && fname == NULL) {
		fprintf(stderr, "%s: cannot adjust tone of standard input\n",
				progname);
		tmflags = TM_F_LINEAR;
	}
	if (tmflags == TM_F_LINEAR) {	/* linear with clamping */
		setcolrcor(pow, 1.0/gamcor);
		return;
	}
	flags = tmflags;		/* histogram adjustment */
	if (greyscale) flags |= TM_F_BW;
	if (tmTop != NULL) {		/* reuse old histogram if one */
		tmDone(tmTop);
		tmTop->flags = flags;
	} else {			/* else initialize */
		if ((lscan = (TMbright *)malloc(xmax*sizeof(TMbright))) == NULL)
			goto memerr;
		if (greyscale) {
			cscan = TM_NOCHROM;
			if ((pscan = (BYTE *)malloc(sizeof(BYTE)*xmax)) == NULL)
				goto memerr;
		} else if ((pscan=cscan = (BYTE *)malloc(3*sizeof(BYTE)*xmax))
				== NULL)
			goto memerr;
						/* initialize tm library */
		if (tmInit(flags, stdprims, gamcor) == NULL)
			goto memerr;
		if (tmSetSpace(stdprims, WHTEFFICACY/exposure))
			goto tmerr;
						/* compute picture histogram */
		for (y = 0; y < ymax; y++) {
			getscan(y);
			if (tmCvColrs(lscan, TM_NOCHROM, scanline, xmax))
				goto tmerr;
			if (tmAddHisto(lscan, xmax, 1))
				goto tmerr;
		}
	}
	tmDup();			/* add fixations to duplicate map */
	dobox(addfix, NULL);
					/* (re)compute tone mapping */
	if (tmComputeMapping(gamcor, 0., 0.))
		goto tmerr;
	return;
memerr:
	quiterr("out of memory in make_tonemap");
tmerr:
	quiterr("tone mapping error");
}


static void
tmap_colrs(		/* apply tone mapping to scanline */
	register COLR  *scn,
	int  len
)
{
	register BYTE  *ps;

	if (tmflags == TM_F_LINEAR) {
		if (scale)
			shiftcolrs(scn, len, scale);
		colrs_gambs(scn, len);
		return;
	}
	if (len > xmax)
		quiterr("code error 1 in tmap_colrs");
	if (tmCvColrs(lscan, cscan, scn, len))
		goto tmerr;
	if (tmMapPixels(pscan, lscan, cscan, len))
		goto tmerr;
	ps = pscan;
	if (greyscale)
		while (len--) {
			scn[0][RED] = scn[0][GRN] = scn[0][BLU] = *ps++;
			scn[0][EXP] = COLXS;
			scn++;
		}
	else
		while (len--) {
			scn[0][RED] = *ps++;
			scn[0][GRN] = *ps++;
			scn[0][BLU] = *ps++;
			scn[0][EXP] = COLXS;
			scn++;
		}
	return;
tmerr:
	quiterr("tone mapping error");
}


static void
getmono(void)			/* get monochrome data */
{
	register unsigned char	*dp;
	register int	x, err;
	int	y, errp;
	short	*cerr;

	if ((cerr = (short *)calloc(xmax,sizeof(short))) == NULL)
		quiterr("out of memory in getmono");
	dp = ourdata - 1;
	for (y = 0; y < ymax; y++) {
		getscan(y);
		add2icon(y, scanline);
		normcolrs(scanline, xmax, scale);
		err = 0;
		for (x = 0; x < xmax; x++) {
			if (!(x&7))
				*++dp = 0;
			errp = err;
			err += normbright(scanline[x]) + cerr[x];
			if (err > 127)
				err -= 255;
			else
				*dp |= 1<<(7-(x&07));
			err /= 3;
			cerr[x] = err + errp;
		}
	}
	free((void *)cerr);
}


static void
add2icon(		/* add a scanline to our icon data */
	int  y,
	COLR  *scan
)
{
	static short  cerr[ICONSIZ];
	static int  ynext;
	static char  *dp;
	COLR  clr;
	register int  err;
	register int	x, ti;
	int  errp;

	if (iconheight == 0) {		/* initialize */
		if (xmax <= ICONSIZ && ymax <= ICONSIZ) {
			iconwidth = xmax;
			iconheight = ymax;
		} else if (xmax > ymax) {
			iconwidth = ICONSIZ;
			iconheight = ICONSIZ*ymax/xmax;
			if (iconheight < 1)
				iconheight = 1;
		} else {
			iconwidth = ICONSIZ*xmax/ymax;
			if (iconwidth < 1)
				iconwidth = 1;
			iconheight = ICONSIZ;
		}
		ynext = 0;
		dp = icondata - 1;
	}
	if (y < ynext*ymax/iconheight)	/* skip this one */
		return;
	err = 0;
	for (x = 0; x < iconwidth; x++) {
		if (!(x&7))
			*++dp = 0;
		errp = err;
		ti = x*xmax/iconwidth;
		copycolr(clr, scan[ti]);
		normcolrs(&clr, 1, scale);
		err += normbright(clr) + cerr[x];
		if (err > 127)
			err -= 255;
		else
			*dp |= 1<<(x&07);
		err /= 3;
		cerr[x] = err + errp;
	}
	ynext++;
}


static void
getfull(void)			/* get full (24-bit) data */
{
	int	y;
	register uint32	*dp;
	register uint16	*dph;
	register int	x;
					/* initialize tone mapping */
	make_tonemap();
					/* read and convert file */
	dp = (uint32 *)ourdata;
	dph = (uint16 *)ourdata;
	for (y = 0; y < ymax; y++) {
		getscan(y);
		add2icon(y, scanline);
		tmap_colrs(scanline, xmax);
		switch (ourras->image->blue_mask) {
		case 0xff:		/* 24-bit RGB */
			for (x = 0; x < xmax; x++)
				*dp++ =	(uint32)scanline[x][RED] << 16 |
					(uint32)scanline[x][GRN] << 8 |
					(uint32)scanline[x][BLU] ;
			break;
		case 0xff0000:		/* 24-bit BGR */
			for (x = 0; x < xmax; x++)
				*dp++ =	(uint32)scanline[x][RED] |
					(uint32)scanline[x][GRN] << 8 |
					(uint32)scanline[x][BLU] << 16 ;
			break;
#if 0
		case 0x1f:		/* 15-bit RGB */
			for (x = 0; x < xmax; x++)
				*dph++ =	(scanline[x][RED] << 7 & 0x7c00) |
					(scanline[x][GRN] << 2 & 0x3e0) |
					(unsigned)scanline[x][BLU] >> 3 ;
			break;
		case 0x7c00:		/* 15-bit BGR */
			for (x = 0; x < xmax; x++)
				*dph++ =	(unsigned)scanline[x][RED] >> 3 |
					(scanline[x][GRN] << 2 & 0x3e0) |
					(scanline[x][BLU] << 7 & 0x7c00) ;
			break;
#endif
		default:		/* unknown */
			if (ourvis.depth > 16)
				for (x = 0; x < xmax; x++) {
					*dp = ourras->image->red_mask &
						ourras->image->red_mask*scanline[x][RED]/255;
					*dp |= ourras->image->green_mask &
						ourras->image->green_mask*scanline[x][GRN]/255;
					*dp++ |= ourras->image->blue_mask &
						ourras->image->blue_mask*scanline[x][BLU]/255;
				}
			else
				for (x = 0; x < xmax; x++) {
					*dph = ourras->image->red_mask &
						ourras->image->red_mask*scanline[x][RED]/255;
					*dph |= ourras->image->green_mask &
						ourras->image->green_mask*scanline[x][GRN]/255;
					*dph++ |= ourras->image->blue_mask &
						ourras->image->blue_mask*scanline[x][BLU]/255;
				}
			break;
		}
	}
}


static void
getgrey(void)			/* get greyscale data */
{
	int	y;
	register unsigned char	*dp;
	register int	x;
					/* initialize tone mapping */
	make_tonemap();
					/* read and convert file */
	dp = ourdata;
	for (y = 0; y < ymax; y++) {
		getscan(y);
		add2icon(y, scanline);
		tmap_colrs(scanline, xmax);
		if (maxcolors < 256)
			for (x = 0; x < xmax; x++)
				*dp++ =	((int32)scanline[x][GRN] *
					maxcolors + maxcolors/2) >> 8;
		else
			for (x = 0; x < xmax; x++)
				*dp++ =	scanline[x][GRN];
	}
	for (x = 0; x < maxcolors; x++)
		clrtab[x][RED] = clrtab[x][GRN] =
			clrtab[x][BLU] = ((int32)x*256 + 128)/maxcolors;
}


static void
getmapped(void)			/* get color-mapped data */
{
	int	y;
					/* make sure we can do it first */
	if (fname == NULL)
		quiterr("cannot map colors from standard input");
					/* initialize tone mapping */
	make_tonemap();
					/* make histogram */
	if (new_histo((int32)xmax*ymax) == -1)
		quiterr("cannot initialize histogram");
	for (y = 0; y < ymax; y++) {
		if (getscan(y) < 0)
			break;
		add2icon(y, scanline);
		tmap_colrs(scanline, xmax);
		cnt_colrs(scanline, xmax);
	}
					/* map pixels */
	if (!new_clrtab(maxcolors))
		quiterr("cannot create color map");
	for (y = 0; y < ymax; y++) {
		getscan(y);
		tmap_colrs(scanline, xmax);
		if (dither)
			dith_colrs(ourdata+y*xmax, scanline, xmax);
		else
			map_colrs(ourdata+y*xmax, scanline, xmax);
	}
}


static void
scale_rcolors(			/* scale color map */
	register XRASTER	*xr,
	double	sf
)
{
	register int	i;
	long	maxv;

	if (xr->pixels == NULL)
		return;

	sf = pow(sf, 1.0/gamcor);
	maxv = 65535/sf;

	for (i = xr->ncolors; i--; ) {
		xr->cdefs[i].red = xr->cdefs[i].red > maxv ?
				65535 :
				xr->cdefs[i].red * sf;
		xr->cdefs[i].green = xr->cdefs[i].green > maxv ?
				65535 :
				xr->cdefs[i].green * sf;
		xr->cdefs[i].blue = xr->cdefs[i].blue > maxv ?
				65535 :
				xr->cdefs[i].blue * sf;
	}
	XStoreColors(thedisplay, xr->cmap, xr->cdefs, xr->ncolors);
}


static int
getscan(
	int  y
)
{
	static int  trunced = -1;		/* truncated file? */
skipit:
	if (trunced >= 0 && y >= trunced) {
		memset(scanline, '\0', xmax*sizeof(COLR));
		return(-1);
	}
	if (y != cury) {
		if (scanpos == NULL || scanpos[y] == -1)
			return(-1);
		if (fseek(fin, scanpos[y], 0) == -1)
			quiterr("fseek error");
		cury = y;
	} else if (scanpos != NULL && scanpos[y] == -1)
		scanpos[y] = ftell(fin);

	if (freadcolrs(scanline, xmax, fin) < 0) {
		fprintf(stderr, "%s: %s: unfinished picture\n",
				progname, fname==NULL?"<stdin>":fname);
		trunced = y;
		goto skipit;
	}
	cury++;
	return(0);
}
