#ifndef lint
static const char	RCSid[] = "$Id: xglaresrc.c,v 2.7 2003/11/10 16:52:24 greg Exp $";
#endif
/*
 *  Circle sources in a displayed image.
 *
 *	18 Mar 1991	Greg Ward
 */

#include "standard.h"
#include "view.h"
#include "paths.h"
#include <signal.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#define XIM		"ximage"

#define NSEG		30		/* number of segments per circle */

#define  FONTNAME	"8x13"		/* text font we'll use */

float	col[3] = {1.,0.,0.};		/* color */

VIEW	ourview = STDVIEW;		/* view for picture */
RESOLU	pres;				/* picture resolution */

char	*progname;			/* program name */

Display	*theDisplay = NULL;		/* connection to server */

#define rwind		RootWindow(theDisplay,ourScreen)
#define ourScreen	DefaultScreen(theDisplay)

GC	vecGC, strGC;
Window	gwind;


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*windowname = NULL;
	FILE	*fp;

	progname = *argv++; argc--;
	while (argc > 0 && argv[0][0] == '-') {
		switch (argv[0][1]) {
		case 'n':
			windowname = *++argv;
			argc--;
			break;
		case 'c':
			col[0] = atof(*++argv);
			col[1] = atof(*++argv);
			col[2] = atof(*++argv);
			argc -= 3;
			break;
		}
		argv++; argc--;
	}
	if (argc < 1 || argc > 2) {
		fprintf(stderr,
		"Usage: %s [-n windowname][-c color] picture [glaresrc]\n",
				progname);
		exit(1);
	}
	init(argv[0], windowname);
	if (argc < 2)
		fp = stdin;
	else if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "%s: cannot open \"%s\"\n", progname, argv[1]);
		exit(1);
	}
	circle_sources(fp);
	exit(0);
}


init(pname, wname)		/* get view and find window */
char	*pname, *wname;
{
	extern Window	xfindwind();
	XWindowAttributes	wa;
	XColor	xc;
	XGCValues	gcv;
	register int	i;
					/* get the viewing parameters */
	if (viewfile(pname, &ourview, &pres) <= 0 ||
			setview(&ourview) != NULL) {
		fprintf(stderr, "%s: cannot get view from \"%s\"\n",
				progname, pname);
		exit(1);
	}
					/* open the display */
	if ((theDisplay = XOpenDisplay(NULL)) == NULL) {
		fprintf(stderr,
		"%s: cannot open display; DISPLAY variable set?\n",
				progname);
		exit(1);
	}
					/* find our window */
	if (wname == NULL) {
				/* remove directory prefix from name */
		for (i = strlen(pname); i-- > 0; )
			if (pname[i] == '/')
				break;
		wname = pname+i+1;
		i = 0;
	} else
		i = 1;
	gwind = xfindwind(theDisplay, rwind, wname, 4);
	if (gwind == None) {
		if (i) {
			fprintf(stderr, "%s: cannot find \"%s\" window\n",
					progname, wname);
			exit(2);
		}
					/* start ximage */
		if (vfork() == 0) {
			execlp(XIM, XIM, "-c", "256", pname, 0);
			perror(XIM);
			fprintf(stderr, "%s: cannot start %s\n",
					progname, XIM);
			kill(getppid(), SIGPIPE);
			_exit(1);
		}
		do
			sleep(8);
		while ((gwind=xfindwind(theDisplay,rwind,wname,4)) == None);
	} else
		XMapRaised(theDisplay, gwind);
	do {
		XGetWindowAttributes(theDisplay, gwind, &wa);
		sleep(6);
	} while (wa.map_state != IsViewable);
	if (wa.width != scanlen(&pres) || wa.height != numscans(&pres)) {
		fprintf(stderr,
		"%s: warning -- window seems to be the wrong size!\n",
				progname);
		if (pres.rt & YMAJOR) {
			pres.xr = wa.width;
			pres.yr = wa.height;
		} else {
			pres.xr = wa.height;
			pres.yr = wa.width;
		}
	}
					/* set graphics context */
	gcv.font = XLoadFont(theDisplay, FONTNAME);
	if (gcv.font == 0) {
		fprintf(stderr, "%s: cannot load font \"%s\"\n",
				progname, FONTNAME);
		exit(1);
	}
	xc.red = col[0] >= 1.0 ? 65535 : (unsigned)(65536*col[0]);
	xc.green = col[1] >= 1.0 ? 65535 : (unsigned)(65536*col[1]);
	xc.blue = col[2] >= 1.0 ? 65535 : (unsigned)(65536*col[2]);
	xc.flags = DoRed|DoGreen|DoBlue;
	gcv.background = xc.green >= 32768 ?
			BlackPixel(theDisplay,DefaultScreen(theDisplay)) :
			WhitePixel(theDisplay,DefaultScreen(theDisplay)) ;
	if (XAllocColor(theDisplay, wa.colormap, &xc)) {
		gcv.foreground = xc.pixel;
		vecGC = XCreateGC(theDisplay,gwind,
				GCForeground|GCBackground|GCFont,&gcv);
		strGC = vecGC;
	} else {
		gcv.function = GXinvert;
		vecGC = XCreateGC(theDisplay,gwind,GCFunction,&gcv);
		gcv.foreground = xc.green < 32768 ?
			BlackPixel(theDisplay,DefaultScreen(theDisplay)) :
			WhitePixel(theDisplay,DefaultScreen(theDisplay)) ;
		strGC = XCreateGC(theDisplay,gwind,
				GCForeground|GCBackground|GCFont,&gcv);
	}
}


circle_sources(fp)		/* circle sources listed in fp */
FILE	*fp;
{
	char	linbuf[256];
	int	reading = 0;
	FVECT	dir;
	double	dom, lum;

	while (fgets(linbuf, sizeof(linbuf), fp) != NULL)
		if (reading) {
			if (!strncmp(linbuf, "END", 3)) {
				XFlush(theDisplay);
				return;
			}
			if (sscanf(linbuf, "%lf %lf %lf %lf %lf",
					&dir[0], &dir[1], &dir[2],
					&dom, &lum) != 5)
				break;
			circle(dir, dom);
			value(dir, lum);
		} else if (!strcmp(linbuf, "BEGIN glare source\n"))
			reading++;

	fprintf(stderr, "%s: error reading glare sources\n", progname);
	exit(1);
}


circle(dir, dom)		/* indicate a solid angle on image */
FVECT	dir;
double	dom;
{
	FVECT	start, cur;
	XPoint	pt[NSEG+1];
	FVECT	pp;
	int	ip[2];
	register int	i;

	fcross(cur, dir, ourview.vup);
	if (normalize(cur) == 0.0)
		goto fail;
	spinvector(start, dir, cur, acos(1.-dom/(2.*PI)));
	for (i = 0; i <= NSEG; i++) {
		spinvector(cur, start, dir, 2.*PI*i/NSEG);
		cur[0] += ourview.vp[0];
		cur[1] += ourview.vp[1];
		cur[2] += ourview.vp[2];
		viewloc(pp, &ourview, cur);
		if (pp[2] <= 0.0)
			goto fail;
		loc2pix(ip, &pres, pp[0], pp[1]);
		pt[i].x = ip[0];
		pt[i].y = ip[1];
	}
	XDrawLines(theDisplay, gwind, vecGC, pt, NSEG+1, CoordModeOrigin);
	return;
fail:
	fprintf(stderr, "%s: cannot draw source at (%f,%f,%f)\n",
			progname, dir[0], dir[1], dir[2]);
}


value(dir, v)			/* print value on image */
FVECT	dir;
double	v;
{
	FVECT	pos;
	FVECT	pp;
	int	ip[2];
	char	buf[32];

	pos[0] = ourview.vp[0] + dir[0];
	pos[1] = ourview.vp[1] + dir[1];
	pos[2] = ourview.vp[2] + dir[2];
	viewloc(pp, &ourview, pos);
	if (pp[2] <= 0.0)
		return;
	loc2pix(ip, &pres, pp[0], pp[1]);
	sprintf(buf, "%.0f", v);
	XDrawImageString(theDisplay, gwind, strGC,
			ip[0], ip[1], buf, strlen(buf)); 
}
