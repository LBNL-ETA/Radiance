#ifndef lint
static const char	RCSid[] = "$Id: vgaimage.c,v 2.10 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  vgaimage.c - driver for VGA board under DOS
 */

#include  "standard.h"
#include  <graph.h>
#include  "color.h"
#include  "random.h"
#include  "resolu.h"
#include <dos.h>
#include <i86.h>

#define	 M_RDOWN	0x8
#define	 M_RUP		0x10
#define	 M_LDOWN	0x2
#define	 M_LUP		0x4
#define	 M_MOTION	0x1

int crad;
int mouse_event = 0;
int mouse_xpos = -1;
int mouse_ypos = -1;

#define	 hide_cursor()	move_cursor(-1,-1)
#define	 show_cursor()	move_cursor(mouse_xpos,mouse_ypos)

#define	 CTRL(c)	((c)-'@')

#define	 MAXWIDTH	1024
#define	 MAXHEIGHT	768

short  ourblack = 0; ourwhite = 1;

double	gamcor = 2.2;			/* gamma correction */

int  dither = 1;			/* dither colors? */

int  maxcolors = 0;			/* maximum colors */
int  minpix = 0;			/* minimum pixel value */
int  greyscale = 0;			/* in grey */

int  scale = 0;				/* scalefactor; power of two */

COLR  scanline[MAXWIDTH];		/* scan line buffer */

int  xmax, ymax;			/* picture dimensions */
FILE  *fin = stdin;			/* input file */
long  scanpos[MAXHEIGHT];		/* scan line positions in file */
int  cury = 0;				/* current scan location */

double	exposure = 1.0;			/* exposure compensation used */

int  wrongformat = 0;			/* input in another format? */

struct {
	int  xmin, ymin, xsiz, ysiz;
}  box = {0, 0, 0, 0};			/* current box */

int  initialized = 0;
int  cheight, cwidth;

#define postext(x,y)	_settextposition(1+(y)/cheight,1+(x)/cwidth)

char  *progname;

char  errmsg[128];

extern BYTE  clrtab[256][3];		/* global color map */

extern long  ftell();


main(argc, argv)
int  argc;
char  *argv[];
{
	extern char  *getenv(), *fixargv0();
	char  *gv;
	int  headline();
	int  i;

	progname = argv[0] = fixargv0(argv[0]);
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
		else
			break;

	if (i == argc-1) {		/* open picture */
		fin = fopen(argv[i], "r");
		if (fin == NULL) {
			sprintf(errmsg, "cannot open file \"%s\"", argv[i]);
			quiterr(errmsg);
		}
	} else if (i != argc)
		goto userr;
	setmode(fileno(fin), O_BINARY);
				/* get header */
	getheader(fin, headline, NULL);
				/* get picture dimensions */
	if (wrongformat || fgetresolu(&xmax, &ymax, fin) < 0)
		quiterr("bad picture format");
	if (xmax > MAXWIDTH | ymax > MAXHEIGHT)
		quiterr("input picture too large for VGA");

	init();			/* initialize and display */

	while (docommand())	/* loop on command */
		;
	quiterr(NULL);
userr:
	fprintf(stderr,
	"Usage: %s [-b][-m][-d][-c ncolors][-e +/-stops] file\n", progname);
	exit(1);
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
	}
	return(0);
}


init()			/* initialize and load display */
{
	static struct {
		short	mode;
		short	xsiz, ysiz;
	} video[] = {
		{_MRES256COLOR, 320, 200},
		{_VRES256COLOR, 640, 480},
		{_SVRES256COLOR, 800, 600},
		{_XRES256COLOR, 1024, 768},
		{-1, 0, 0}
	};
	struct videoconfig	config;
	register int	i;
				/* pick a card... */
	for (i = 0; video[i].mode != -1; i++)
		if (video[i].xsiz >= xmax && video[i].ysiz >= ymax)
			break;
	if (video[i].mode == -1)
		quiterr("input picture too large");
	if (_setvideomode(video[i].mode) == 0)
		quiterr("inadequate display card for picture");
	ms_init();
	initialized = 1;
	_getvideoconfig(&config);
	if (maxcolors == 0 | maxcolors > config.numcolors)
		maxcolors = config.numcolors-2;
	if (maxcolors <= config.numcolors-2)
		minpix = 2;
	else
		ourwhite = maxcolors-1;
	_settextcolor(ourwhite);
	cheight = config.numypixels/config.numtextrows;
	cwidth = config.numxpixels/config.numtextcols;
				/* clear scan position array */
	for (i = 0; i < ymax; i++)
		scanpos[i] = -1;
				/* display image */
	if (greyscale)
		greyimage();
	else
		mappedimage();
}


quiterr(err)		/* print message and exit */
char  *err;
{
	if (initialized) {
		ms_done();
		_setvideomode(_DEFAULTMODE);
	}
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


int
docommand()			/* execute command */
{
	char  buf[64];
	COLOR  cval;
	int  com;
	double	comp;

	while (!kbhit())
		watch_mouse();
	com = getch();
	switch (com) {			/* interpret command */
	case 'q':
	case CTRL('Z'):				/* quit */
		return(0);
	case '\n':
	case '\r':
	case 'l':
	case 'c':				/* value */
		if (avgbox(cval) == -1)
			return(-1);
		switch (com) {
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
		postext(box.xmin+box.xsiz/2, box.ymin+box.ysiz/2);
		hide_cursor();
		_outtext(buf);
		show_cursor();
		return(1);
	default:
		return(-1);
	}
}


watch_mouse()			/* look after mousie */
{
	int	a_x, a_y, l_x, l_y;

	if (mouse_event & M_MOTION)
		move_cursor(mouse_xpos, mouse_ypos);
	if (!(mouse_event & M_LDOWN))
		return;
	l_x = a_x = mouse_xpos; l_y = a_y = mouse_ypos;
	hide_cursor();
	revbox(a_x, a_y, l_x, l_y);		/* show box */
	do {
		mouse_event = 0;
		while (!mouse_event)
			;
		if (mouse_event & M_MOTION) {
			revbox(a_x, a_y, l_x, l_y);
			revbox(a_x, a_y, l_x=mouse_xpos, l_y=mouse_ypos);
		}
	} while (!(mouse_event & M_LUP));
	revbox(a_x, a_y, l_x, l_y);		/* hide box */
	show_cursor();
	box.xmin = mouse_xpos;
	box.ymin = mouse_ypos;
	if (box.xmin > a_x) {
		box.xsiz = box.xmin - a_x + 1;
		box.xmin = a_x;
	} else {
		box.xsiz = a_x - box.xmin + 1;
	}
	if (box.ymin > a_y) {
		box.ysiz = box.ymin - a_y + 1;
		box.ymin = a_y;
	} else {
		box.ysiz = a_y - box.ymin + 1;
	}
	mouse_event = 0;
}


revbox(x0, y0, x1, y1)			/* draw box with reversed lines */
int  x0, y0, x1, y1;
{
	_setplotaction(_GXOR);
	_setcolor(255);
	_moveto(x0, y0);
	_lineto(x1, y0);
	_lineto(x1, y1);
	_lineto(x0, y1);
	_lineto(x0, y0);
}


int
avgbox(clr)				/* average color over current box */
COLOR  clr;
{
	static COLOR  lc;
	static int  ll, lr, lt, lb;
	int  left, right, top, bottom;
	int  y;
	double	d;
	COLOR  ctmp;
	register int  x;

	setcolor(clr, 0.0, 0.0, 0.0);
	left = box.xmin;
	right = left + box.xsiz;
	if (left < 0)
		left = 0;
	if (right > xmax)
		right = xmax;
	if (left >= right)
		return(-1);
	top = box.ymin;
	bottom = top + box.ysiz;
	if (top < 0)
		top = 0;
	if (bottom > ymax)
		bottom = ymax;
	if (top >= bottom)
		return(-1);
	if (left == ll && right == lr && top == lt && bottom == lb) {
		copycolor(clr, lc);
		return(0);
	}
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
	ll = left; lr = right; lt = top; lb = bottom;
	copycolor(lc, clr);
	return(0);
}


setpalette()			/* set our palette using clrtab */
{
	long	cvals[256];
	register int	i;

	cvals[ourblack] = _BLACK;
	cvals[ourwhite] = _BRIGHTWHITE;
	for (i = 0; i < maxcolors; i++)
		cvals[i+minpix] = (long)clrtab[i][BLU]<<14 & 0x3f0000L |
				  clrtab[i][GRN]<<6 & 0x3f00 |
				  clrtab[i][RED]>>2;
	_remapallpalette(cvals);
}


greyimage()			/* display greyscale image */
{
	short	thiscolor, lastcolor = -1;
	int	y;
	register int	x;
					/* set gamma correction */
	setcolrgam(gamcor);
					/* set up color map */
	for (x = 0; x < maxcolors; x++)
		clrtab[x][RED] = clrtab[x][GRN] = clrtab[x][BLU] =
				((long)x*256 + 128)/maxcolors;
	setpalette();
	_setplotaction(_GPSET);
					/* read and display file */
	for (y = 0; y < ymax; y++) {
		getscan(y);
		if (scale)
			shiftcolrs(scanline, xmax, scale);
		for (x = 0; x < xmax; x++)
			scanline[x][GRN] = normbright(scanline[x]);
		colrs_gambs(scanline, xmax);
		if (maxcolors < 256)
			for (x = 0; x < xmax; x++) {
				thiscolor = ((long)scanline[x][GRN] *
						maxcolors + maxcolors/2) / 256;
				if (thiscolor != lastcolor)
					_setcolor((lastcolor=thiscolor)+minpix);
				_setpixel(x, y);
			}
		else
			for (x = 0; x < xmax; x++) {
				thiscolor = scanline[x][GRN];
				if (thiscolor != lastcolor)
					_setcolor((lastcolor=thiscolor)+minpix);
				_setpixel(x, y);
			}
	}
}


mappedimage()			/* display color-mapped image */
{
	BYTE	bscan[MAXWIDTH];
	int	y;
	register int	x;
					/* set gamma correction */
	setcolrgam(gamcor);
					/* make histogram */
	_outtext("Quantizing image -- Please wait...");
	new_histo();
	for (y = 0; y < ymax; y++) {
		if (getscan(y) < 0)
			quiterr("seek error in getmapped");
		if (scale)
			shiftcolrs(scanline, xmax, scale);
		colrs_gambs(scanline, xmax);
		cnt_colrs(scanline, xmax);
	}
					/* map pixels */
	if (!new_clrtab(maxcolors))
		quiterr("cannot create color map");
	setpalette();
	_setplotaction(_GPSET);
					/* display image */
	for (y = 0; y < ymax; y++) {
		if (getscan(y) < 0)
			quiterr("seek error in getmapped");
		if (scale)
			shiftcolrs(scanline, xmax, scale);
		colrs_gambs(scanline, xmax);
		if (dither)
			dith_colrs(bscan, scanline, xmax);
		else
			map_colrs(bscan, scanline, xmax);
		for (x = 0; x < xmax; x++) {
			if (x==0 || bscan[x] != bscan[x-1])
				_setcolor(bscan[x]+minpix);
			_setpixel(x, y);
		}
	}
}


getscan(y)
int  y;
{
	if (y != cury) {
		if (scanpos[y] == -1)
			return(-1);
		if (fseek(fin, scanpos[y], 0) == -1)
			quiterr("fseek error");
		cury = y;
	} else if (fin != stdin && scanpos[y] == -1)
		scanpos[y] = ftell(fin);

	if (freadcolrs(scanline, xmax, fin) < 0)
		quiterr("picture read error");

	cury++;
	return(0);
}


/*
 * Microsoft Mouse handling routines
 */


#pragma off (check_stack)
void _loadds far mouse_handler (int max, int mcx, int mdx)
{
#pragma aux mouse_handler parm [EAX] [ECX] [EDX]
	mouse_event = max;
	mouse_xpos = mcx;
	mouse_ypos = mdx;
}
#pragma on (check_stack)


void
move_cursor(newx, newy)		/* move cursor to new position */
int  newx, newy;
{
	static char  *imp = NULL;
	static int  curx = -1, cury = -1;
#define xcmin		(curx-crad<0 ? 0 : curx-crad)
#define ycmin		(cury-crad<0 ? 0 : cury-crad)
#define xcmax		(curx+crad>=xmax ? xmax-1 : curx+crad)
#define ycmax		(cury+crad>=ymax ? ymax-1 : cury+crad)

	if (newx == curx & newy == cury)
		return;
	if (imp == NULL &&
		(imp = bmalloc(_imagesize(0,0,2*crad+1,2*crad+1))) == NULL) {
		quiterr("out of memory in move_cursor");
	}
	if (curx >= 0 & cury >= 0)	/* clear old cursor */
		_putimage(xcmin, ycmin, imp, _GPSET);
					/* record new position */
	curx = newx; cury = newy;
	if (curx < 0 | cury < 0)
		return;		/* no cursor */
					/* save under new cursor */
	_getimage(xcmin, ycmin, xcmax, ycmax, imp);
					/* draw new cursor */
	_setplotaction(_GPSET);
	_setcolor(ourwhite);
	_rectangle(_GFILLINTERIOR, xcmin, cury-1, xcmax, cury+1);
	_rectangle(_GFILLINTERIOR, curx-1, ycmin, curx+1, ycmax);
	_setcolor(ourblack);
	_moveto(xcmin+1, cury);
	_lineto(xcmax-1, cury);
	_moveto(curx, ycmin+1);
	_lineto(curx, ycmax-1);
#undef xcmin
#undef ycmin
#undef xcmax
#undef ycmax
}


int
ms_init()
{
    struct SREGS sregs;
    union REGS inregs, outregs;
    int far *ptr;
    int (far *function_ptr)();

    segread(&sregs);

    /* check for mouse driver */

    inregs.w.ax = 0;
    int386 (0x33, &inregs, &outregs);
    if( outregs.w.ax != -1 ) {
	return(0);
    }

    crad = ymax/40;

    /* set screen limits */

    inregs.w.ax = 0x7;		/* horizontal resolution */
    inregs.w.cx = 0;
    inregs.w.dx = xmax-1;
    int386x( 0x33, &inregs, &outregs, &sregs );
    inregs.w.ax = 0x8;		/* vertical resolution */
    inregs.w.cx = 0;
    inregs.w.dx = ymax-1;
    int386x( 0x33, &inregs, &outregs, &sregs );

    /* install watcher */

    inregs.w.ax = 0xC;
    inregs.w.cx = M_LDOWN | M_LUP | M_MOTION;
    function_ptr = mouse_handler;
    inregs.x.edx = FP_OFF( function_ptr );
    sregs.es	 = FP_SEG( function_ptr );
    int386x( 0x33, &inregs, &outregs, &sregs );

    return(1);
}

ms_done()
{
    union REGS inregs, outregs;

    /* uninstall watcher */

    inregs.w.ax = 0;
    int386 (0x33, &inregs, &outregs);
}
