#ifndef lint
static const char	RCSid[] = "$Id: ximage.c,v 2.9 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  ximage.c - driver for X-windows
 *
 *     4/30/87
 *     3/3/88		Added calls to xraster & Paul Heckbert's ciq routines
 */

#include  "standard.h"

#include  <X/Xlib.h>
#include  <X/cursors/bcross.cursor>
#include  <X/cursors/bcross_mask.cursor>

#include  <sys/types.h>

#include  <ctype.h>

#include  <time.h>

#include  "color.h"

#include  "resolu.h"

#include  "xraster.h"

#include  "view.h"

#include  "pic.h"

#include  "random.h"

#define  controlshift(e)	(((XButtonEvent *)(e))->detail & (ShiftMask|ControlMask))

#define  FONTNAME	"9x15"		/* text font we'll use */

#define  CTRL(c)	((c)-'@')

#define  BORWIDTH	5		/* border width */
#define  BARHEIGHT	25		/* menu bar size */

double  gamcor = 2.2;			/* gamma correction */

XRASTER  *ourras = NULL;		/* our stored raster image */

int  dither = 1;			/* dither colors? */
int  fast = 0;				/* keep picture in Pixmap? */

Window  wind = 0;			/* our output window */
Font  fontid;				/* our font */

int  maxcolors = 0;			/* maximum colors */
int  greyscale = 0;			/* in grey */

int  scale = 0;				/* scalefactor; power of two */

int  xoff = 0;				/* x image offset */
int  yoff = 0;				/* y image offset */

VIEW  ourview = STDVIEW;		/* image view parameters */
int  gotview = 0;			/* got parameters from file */

COLR  *scanline;			/* scan line buffer */

int  xmax, ymax;			/* picture resolution */
int  width, height;			/* window size */
char  *fname = NULL;			/* input file name */
FILE  *fin = stdin;			/* input file */
long  *scanpos = NULL;			/* scan line positions in file */
int  cury = 0;				/* current scan location */

double  exposure = 1.0;			/* exposure compensation used */

int  wrongformat = 0;			/* input in another format */

struct {
	int  xmin, ymin, xsiz, ysiz;
}  box = {0, 0, 0, 0};			/* current box */

char  *geometry = NULL;			/* geometry specification */

char  *progname;

char  errmsg[128];


main(argc, argv)
int  argc;
char  *argv[];
{
	extern char  *getenv();
	char  *gv;
	int  headline();
	int  i;
	
	progname = argv[0];
	if ((gv = getenv("DISPLAY_GAMMA")) != NULL)
		gamcor = atof(gv);

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'c':
				maxcolors = atoi(argv[++i]);
				break;
			case 'b':
				greyscale = !greyscale;
				break;
			case 'm':
				maxcolors = 2;
				break;
			case 'd':
				dither = !dither;
				break;
			case 'f':
				fast = !fast;
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				scale = atoi(argv[++i]);
				break;
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			default:
				goto userr;
			}
		else if (argv[i][0] == '=')
			geometry = argv[i];
		else
			break;

	if (argc-i == 1) {
		fname = argv[i];
		fin = fopen(fname, "r");
		if (fin == NULL) {
			sprintf(errmsg, "can't open file \"%s\"", fname);
			quiterr(errmsg);
		}
	}
				/* get header */
	getheader(fin, headline, NULL);
				/* get picture dimensions */
	if (wrongformat || fgetresolu(&xmax, &ymax, fin) != (YMAJOR|YDECR))
		quiterr("bad picture size");
				/* set view parameters */
	if (gotview && setview(&ourview) != NULL)
		gotview = 0;
	if ((scanline = (COLR *)malloc(xmax*sizeof(COLR))) == NULL)
		quiterr("out of memory");

	init();			/* get file and open window */

	for ( ; ; )
		getevent();		/* main loop */
userr:
	fprintf(stderr,
	"Usage: %s [=geometry][-b][-m][-d][-f][-c ncolors][-e +/-stops] file\n",
			progname);
	quit(1);
}


int
headline(s)		/* get relevant info from header */
char  *s;
{
	char  fmt[32];

	if (isexpos(s))
		exposure *= exposval(s);
	else if (isformat(s)) {
		formatval(fmt, s);
		wrongformat = strcmp(fmt, COLRFMT);
	} else if (isview(s) && sscanview(&ourview, s) > 0)
		gotview++;
	return(0);
}


init()			/* get data and open window */
{
	register int  i;
	OpaqueFrame  mainframe;
	char  defgeom[32];
	
	if (fname != NULL) {
		scanpos = (long *)malloc(ymax*sizeof(long));
		if (scanpos == NULL)
			goto memerr;
		for (i = 0; i < ymax; i++)
			scanpos[i] = -1;
	}
	if (XOpenDisplay(NULL) == NULL)
		quiterr("can't open display; DISPLAY variable set?");
#ifdef notdef
	if (DisplayWidth() - 2*BORWIDTH < xmax ||
			DisplayHeight() - 2*BORWIDTH - BARHEIGHT < ymax)
		quiterr("resolution mismatch");
#endif
	if (maxcolors == 0) {		/* get number of available colors */
		maxcolors = 1<<DisplayPlanes();
		if (maxcolors > 4) maxcolors -= 2;
	}
				/* store image */
	getras();

	mainframe.bdrwidth = BORWIDTH;
	mainframe.border = WhitePixmap;
	mainframe.background = BlackPixmap;
	sprintf(defgeom, "=%dx%d+0+22", xmax, ymax);
	wind = XCreate("Radiance Image", progname, geometry, defgeom,
			&mainframe, 16, 16);
	if (wind == 0)
		quiterr("can't create window");
	width = mainframe.width;
	height = mainframe.height;
	fontid = XGetFont(FONTNAME);
	if (fontid == 0)
		quiterr("can't get font");
	XStoreName(wind, fname == NULL ? progname : fname);
	XDefineCursor(wind, XCreateCursor(bcross_width, bcross_height,
			bcross_bits, bcross_mask_bits,
			bcross_x_hot, bcross_y_hot,
			BlackPixel, WhitePixel, GXcopy));
	XSelectInput(wind, ButtonPressed|ButtonReleased|UnmapWindow
			|RightDownMotion|MiddleDownMotion|LeftDownMotion
			|KeyPressed|ExposeWindow|ExposeRegion);
	XMapWindow(wind);
	return;
memerr:
	quiterr("out of memory");
}


quiterr(err)		/* print message and exit */
char  *err;
{
	if (err != NULL)
		fprintf(stderr, "%s: %s\n", progname, err);

	exit(err == NULL ? 0 : 1);
}


void
eputs(s)
char	*s;
{
	fputs(s, stderr);
}


void
quit(code)
int  code;
{
	exit(code);
}


getras()				/* get raster file */
{
	colormap	ourmap;

	ourras = (XRASTER *)calloc(1, sizeof(XRASTER));
	if (ourras == NULL)
		goto memerr;
	ourras->width = xmax;
	ourras->height = ymax;
	if (maxcolors <= 2) {		/* monochrome */
		ourras->data.m = (unsigned short *)malloc(BitmapSize(xmax,ymax));
		if (ourras->data.m == NULL)
			goto memerr;
		getmono();
	} else {
		ourras->data.bz = (unsigned char *)malloc(BZPixmapSize(xmax,ymax));
		if (ourras->data.bz == NULL)
			goto memerr;
		if (greyscale)
			biq(dither,maxcolors,1,ourmap);
		else
			ciq(dither,maxcolors,1,ourmap);
		if (init_rcolors(ourras, ourmap) == 0)
			goto memerr;
	}
	return;
memerr:
	quiterr("out of memory");
}


getevent()				/* process the next event */
{
	WindowInfo  info;
	XEvent  e;

	XNextEvent(&e);
	switch (e.type) {
	case KeyPressed:
		docom(&e);
		break;
	case ExposeWindow:
		XQueryWindow(wind, &info);	/* in case resized */
		width = info.width;
		height = info.height;
	/* fall through */
	case ExposeRegion:
		fixwindow(&e);
		break;
	case UnmapWindow:
		unmap_rcolors(ourras);
		break;
	case ButtonPressed:
		if (controlshift(&e))
			moveimage(&e);
		else
			getbox(&e);
		break;
	}
}


fixwindow(eexp)				/* fix new or stepped-on window */
XExposeEvent  *eexp;
{
	redraw(eexp->x, eexp->y, eexp->width, eexp->height);
}


redraw(x, y, w, h)			/* redraw section of window */
int  x, y;
int  w, h;
{
	if (ourras->ncolors && map_rcolors(ourras) == NULL) {
		fprintf(stderr, "%s: cannot allocate colors\n", progname);
		return(-1);
	}
	if (fast)
		make_rpixmap(ourras);
	return(patch_raster(wind,x-xoff,y-yoff,x,y,w,h,ourras) ? 0 : -1);
}


docom(ekey)					/* execute command */
XKeyEvent  *ekey;
{
	char  buf[80];
	COLOR  cval;
	Color  cvx;
	char  *cp;
	int  n;
	double  comp;
	FVECT  rorg, rdir;

	cp = XLookupMapping(ekey, &n);
	if (n == 0)
		return(0);
	switch (*cp) {			/* interpret command */
	case 'q':
	case CTRL('D'):				/* quit */
		quit(0);
	case '\n':
	case '\r':
	case 'l':
	case 'c':				/* value */
		if (avgbox(cval) == -1)
			return(-1);
		switch (*cp) {
		case '\n':
		case '\r':				/* radiance */
			sprintf(buf, "%.3f", intens(cval)/exposure);
			break;
		case 'l':				/* luminance */
			sprintf(buf, "%.0fL", luminance(cval)/exposure);
			break;
		case 'c':				/* color */
			comp = pow(2.0, (double)scale);
			sprintf(buf, "(%.2f,%.2f,%.2f)",
					colval(cval,RED)*comp,
					colval(cval,GRN)*comp,
					colval(cval,BLU)*comp);
			break;
		}
		XText(wind, box.xmin, box.ymin, buf, strlen(buf),
				fontid, BlackPixel, WhitePixel);
		return(0);
	case 'i':				/* identify (contour) */
		if (ourras->pixels == NULL)
			return(-1);
		n = ourras->data.bz[ekey->x-xoff+BZPixmapSize(xmax,ekey->y-yoff)];
		n = ourras->pmap[n];
		cvx.pixel = ourras->cdefs[n].pixel;
		cvx.red = random() & 65535;
		cvx.green = random() & 65535;
		cvx.blue = random() & 65535;
		XStoreColor(&cvx);
		return(0);
	case 'p':				/* position */
		sprintf(buf, "(%d,%d)", ekey->x-xoff, ymax-1-ekey->y+yoff);
		XText(wind, ekey->x, ekey->y, buf, strlen(buf),
				fontid, BlackPixel, WhitePixel);
		return(0);
	case 't':				/* trace */
		if (!gotview) {
			XFeep(0);
			return(-1);
		}
		if (viewray(rorg, rdir, &ourview, (ekey->x-xoff+.5)/xmax,
				(ymax-1-ekey->y+yoff+.5)/ymax) < 0)
			return(-1);
		printf("%e %e %e ", rorg[0], rorg[1], rorg[2]);
		printf("%e %e %e\n", rdir[0], rdir[1], rdir[2]);
		fflush(stdout);
		return(0);
	case '=':				/* adjust exposure */
		if (avgbox(cval) == -1)
			return(-1);
		n = log(.5/bright(cval))/.69315 - scale;	/* truncate */
		if (n == 0)
			return(0);
		scale_rcolors(ourras, pow(2.0, (double)n));
		scale += n;
		sprintf(buf, "%+d", scale);
		XText(wind, box.xmin, box.ymin, buf, strlen(buf),
				fontid, BlackPixel, WhitePixel);
		XFlush();
		free_raster(ourras);
		getras();
	/* fall through */
	case CTRL('R'):				/* redraw */
	case CTRL('L'):
		unmap_rcolors(ourras);
		XClear(wind);
		return(redraw(0, 0, width, height));
	case ' ':				/* clear */
		return(redraw(box.xmin, box.ymin, box.xsiz, box.ysiz));
	default:
		XFeep(0);
		return(-1);
	}
}


moveimage(ep)				/* shift the image */
XButtonEvent  *ep;
{
	XButtonEvent  eb;

	XMaskEvent(ButtonReleased, &eb);
	xoff += eb.x - ep->x;
	yoff += eb.y - ep->y;
	XClear(wind);
	return(redraw(0, 0, width, height));
}


getbox(ebut)				/* get new box */
XButtonEvent  *ebut;
{
	union {
		XEvent  e;
		XButtonEvent  b;
		XMouseMovedEvent  m;
	}  e;

	XMaskEvent(ButtonReleased|MouseMoved, &e.e);
	while (e.e.type == MouseMoved) {
		revbox(ebut->x, ebut->y, box.xmin = e.m.x, box.ymin = e.m.y);
		XMaskEvent(ButtonReleased|MouseMoved, &e.e);
		revbox(ebut->x, ebut->y, box.xmin, box.ymin);
	}
	box.xmin = e.b.x<0 ? 0 : (e.b.x>=width ? width-1 : e.b.x);
	box.ymin = e.b.y<0 ? 0 : (e.b.y>=height ? height-1 : e.b.y);
	if (box.xmin > ebut->x) {
		box.xsiz = box.xmin - ebut->x + 1;
		box.xmin = ebut->x;
	} else {
		box.xsiz = ebut->x - box.xmin + 1;
	}
	if (box.ymin > ebut->y) {
		box.ysiz = box.ymin - ebut->y + 1;
		box.ymin = ebut->y;
	} else {
		box.ysiz = ebut->y - box.ymin + 1;
	}
}


revbox(x0, y0, x1, y1)			/* draw box with reversed lines */
int  x0, y0, x1, y1;
{
	XLine(wind, x0, y0, x1, y0, 1, 1, 0, GXinvert, AllPlanes);
	XLine(wind, x0, y1, x1, y1, 1, 1, 0, GXinvert, AllPlanes);
	XLine(wind, x0, y0, x0, y1, 1, 1, 0, GXinvert, AllPlanes);
	XLine(wind, x1, y0, x1, y1, 1, 1, 0, GXinvert, AllPlanes);
}


avgbox(clr)				/* average color over current box */
COLOR  clr;
{
	int  left, right, top, bottom;
	int  y;
	double  d;
	COLOR  ctmp;
	register int  x;

	setcolor(clr, 0.0, 0.0, 0.0);
	left = box.xmin - xoff;
	right = left + box.xsiz;
	if (left < 0)
		left = 0;
	if (right > xmax)
		right = xmax;
	if (left >= right)
		return(-1);
	top = box.ymin - yoff;
	bottom = top + box.ysiz;
	if (top < 0)
		top = 0;
	if (bottom > ymax)
		bottom = ymax;
	if (top >= bottom)
		return(-1);
	for (y = top; y < bottom; y++) {
		if (getscan(y) == -1)
			return(-1);
		for (x = left; x < right; x++) {
			colr_color(ctmp, scanline[x]);
			addcolor(clr, ctmp);
		}
	}
	d = 1.0/((right-left)*(bottom-top));
	scalecolor(clr, d);
	return(0);
}


getmono()			/* get monochrome data */
{
	register unsigned short	*dp;
	register int	x, err;
	int	y, errp;
	rgbpixel	*inl;
	short	*cerr;

	if ((inl = (rgbpixel *)malloc(xmax*sizeof(rgbpixel))) == NULL
			|| (cerr = (short *)calloc(xmax,sizeof(short))) == NULL)
		quiterr("out of memory in getmono");
	dp = ourras->data.m - 1;
	for (y = 0; y < ymax; y++) {
		picreadline3(y, inl);
		err = 0;
		for (x = 0; x < xmax; x++) {
			if (!(x&0xf))
				*++dp = 0;
			errp = err;
			err += rgb_bright(&inl[x]) + cerr[x];
			if (err > 127)
				err -= 255;
			else
				*dp |= 1<<(x&0xf);
			err /= 3;
			cerr[x] = err + errp;
		}
	}
	free((void *)inl);
	free((void *)cerr);
}


init_rcolors(xr, cmap)				/* (re)assign color values */
register XRASTER	*xr;
colormap	cmap;
{
	register int	i;
	register unsigned char	*p;

	xr->pmap = (int *)malloc(256*sizeof(int));
	if (xr->pmap == NULL)
		return(0);
	xr->cdefs = (Color *)malloc(256*sizeof(Color));
	if (xr->cdefs == NULL)
		return(0);
	for (i = 0; i < 256; i++)
		xr->pmap[i] = -1;
	xr->ncolors = 0;
	for (p = xr->data.bz, i = BZPixmapSize(xr->width,xr->height); i--; p++)
		if (xr->pmap[*p] == -1) {
			xr->cdefs[xr->ncolors].red = cmap[0][*p] << 8;
			xr->cdefs[xr->ncolors].green = cmap[1][*p] << 8;
			xr->cdefs[xr->ncolors].blue = cmap[2][*p] << 8;
			xr->cdefs[xr->ncolors].pixel = *p;
			xr->pmap[*p] = xr->ncolors++;
		}
	xr->cdefs = (Color *)realloc((char *)xr->cdefs, xr->ncolors*sizeof(Color));
	if (xr->cdefs == NULL)
		return(0);
	return(1);
}


scale_rcolors(xr, sf)			/* scale color map */
register XRASTER	*xr;
double	sf;
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
	XStoreColors(xr->ncolors, xr->cdefs);
}


getscan(y)
int  y;
{
	if (y != cury) {
		if (scanpos == NULL || scanpos[y] == -1)
			return(-1);
		if (fseek(fin, scanpos[y], 0) == -1)
			quiterr("fseek error");
		cury = y;
	} else if (scanpos != NULL && scanpos[y] == -1)
		scanpos[y] = ftell(fin);

	if (freadcolrs(scanline, xmax, fin) < 0)
		quiterr("read error");

	cury++;
	return(0);
}


picreadline3(y, l3)			/* read in 3-byte scanline */
int  y;
register rgbpixel  *l3;
{
	register int	i;
							/* read scanline */
	if (getscan(y) < 0)
		quiterr("cannot seek for picreadline");
							/* convert scanline */
	normcolrs(scanline, xmax, scale);
	for (i = 0; i < xmax; i++) {
		l3[i].r = scanline[i][RED];
		l3[i].g = scanline[i][GRN];
		l3[i].b = scanline[i][BLU];
	}
}


picwriteline(y, l)		/* add 8-bit scanline to image */
int  y;
pixel  *l;
{
	bcopy((char *)l, (char *)ourras->data.bz+BZPixmapSize(xmax,y), BZPixmapSize(xmax,1));
}


picreadcm(map)			/* do gamma correction */
colormap  map;
{
	extern double  pow();
	register int  i, val;

	for (i = 0; i < 256; i++) {
		val = pow((i+0.5)/256.0, 1.0/gamcor) * 256.0;
		map[0][i] = map[1][i] = map[2][i] = val;
	}
}


picwritecm(map)			/* handled elsewhere */
colormap  map;
{
#ifdef DEBUG
	register int i;

	for (i = 0; i < 256; i++)
		printf("%d %d %d\n", map[0][i],map[1][i],map[2][i]);
#endif
}
