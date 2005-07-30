#ifndef lint
static const char	RCSid[] = "$Id: xshowtrace.c,v 2.12 2005/07/30 22:02:29 greg Exp $";
#endif
/*
 *  Display an image and watch the rays get traced.
 *
 *	9/21/90	Greg Ward
 */

#include <X11/Xlib.h>

#include "standard.h"
#include "paths.h"
#include "view.h"

#define MAXDEPTH	32		/* ridiculous ray tree depth */

#ifdef  SMLFLT
#define  sscanvec(s,v)	(sscanf(s,"%f %f %f",v,v+1,v+2)==3)
#else
#define  sscanvec(s,v)	(sscanf(s,"%lf %lf %lf",v,v+1,v+2)==3)
#endif

char	rtcom[64] = "rtrace -h- -otp -fa -x 1";
char	xicom[] = "ximage -c 256";

VIEW	ourview = STDVIEW;		/* view for picture */
RESOLU	ourres;				/* picture resolution */

char	*progname;			/* program name */

char	*picture;			/* picture name */

FILE	*pin;				/* input stream */

Display	*theDisplay = NULL;		/* connection to server */

struct node {				/* ray tree node */
	int	ipt[2];
	struct node	*sister;
	struct node	*daughter;
};

#define newnode()	(struct node *)calloc(1, sizeof(struct node))

int	slow = 0;		/* slow trace? */

void mainloop(void);
static void freetree(struct node	*tp);
static void tracerays(struct node	*tp);
static int strtoipt(int	ipt[2], char	*str);
static void setvec(int	ipt[2]);
static void vector(int	ip1[2], int	ip2[2]);


int
main(		/* takes both the octree and the image */
	int	argc,
	char	*argv[]
)
{
	int	i;
	char	combuf[PATH_MAX];

	progname = argv[0];
	for (i = 1; i < argc-2; i++)
		if (!strcmp(argv[i], "-s"))
			slow++;
		else if (!strcmp(argv[i], "-T"))
			strcat(rtcom, " -oTp");
		else
			break;
	if (i > argc-2) {
		fprintf(stderr, "Usage: %s [-s][-T] [rtrace args] octree picture\n",
				progname);
		exit(1);
	}
	picture = argv[argc-1];
					/* get the viewing parameters */
	if (viewfile(picture, &ourview, &ourres) <= 0 ||
			setview(&ourview) != NULL) {
		fprintf(stderr, "%s: cannot get view from \"%s\"\n",
				progname, picture);
		exit(1);
	}
					/* open the display */
	if ((theDisplay = XOpenDisplay(NULL)) == NULL) {
		fprintf(stderr,
		"%s: cannot open display; DISPLAY variable set?\n",
				progname);
		exit(1);
	}
					/* build input command */
	sprintf(combuf, "%s \"%s\" | %s", xicom, picture, rtcom);
	for ( ; i < argc-1; i++) {
		strcat(combuf, " ");
		strcat(combuf, argv[i]);
	}
					/* start the damn thing */
	if ((pin = popen(combuf, "r")) == NULL)
		exit(1);
					/* loop on input */
	mainloop();
					/* close pipe and exit */
	pclose(pin);
	exit(0);
}


void
mainloop(void)				/* get and process input */
{
	static struct node	*sis[MAXDEPTH];
	register struct node	*newp;
	char	buf[128];
	int	level;
	register int	i;

	level = 0;
	while (fgets(buf, sizeof(buf), pin) != NULL) {
		if ((newp = newnode()) == NULL) {
			fprintf(stderr, "%s: memory error\n", progname);
			return;
		}
		for (i = 0; buf[i] == '\t'; i++)
			;
		if (strtoipt(newp->ipt, buf+i) < 0) {
			fprintf(stderr, "%s: bad read\n", progname);
			return;
		}
		newp->sister = sis[i];
		sis[i] = newp;
		if (i < level) {
			newp->daughter = sis[level];
			sis[level] = NULL;
		}
		level = i;
		if (i == 0) {
			setvec(sis[0]->ipt);
			tracerays(sis[0]);
			freetree(sis[0]);
			sis[0] = NULL;
			if (!slow)
				XFlush(theDisplay);
		}
	}
}


static void
freetree(				/* free a trace tree */
	struct node	*tp
)
{
	register struct node	*kid, *k2;

	for (kid = tp->daughter; kid != NULL; kid = k2) {
		k2 = kid->sister;
		freetree(kid);
	}
	free((void *)tp);
}


static void
tracerays(				/* trace a ray tree */
	struct node	*tp
)
{
	register struct node	*kid;

	for (kid = tp->daughter; kid != NULL; kid = kid->sister) {
		vector(tp->ipt, kid->ipt);
		tracerays(kid);
	}
}


static int
strtoipt(		/* convert string x y z to image point */
	int	ipt[2],
	char	*str
)
{
	FVECT	im_pt, pt;

	if (!sscanvec(str, pt))
		return(-1);
	if (DOT(pt,pt) <= FTINY)		/* origin is really infinity */
		ipt[0] = ipt[1] = -1;			/* null vector */
	else {
		viewloc(im_pt, &ourview, pt);
		loc2pix(ipt, &ourres, im_pt[0], im_pt[1]);
	}
	return(0);
}


#define rwind		RootWindow(theDisplay,ourScreen)
#define ourScreen	DefaultScreen(theDisplay)

GC	vecGC = 0;
Window	gwind = 0;
int	xoff, yoff;


static void
setvec(			/* set up vector drawing for pick */
	int	ipt[2]
)
{
	extern Window	xfindwind();
	XWindowAttributes	wa;
	XColor	xc;
	XGCValues	gcv;
	int	rx, ry, wx, wy;
	Window	rw, cw;
	unsigned int	pm;
					/* compute pointer location */
	if (gwind == 0) {
		register char	*wn;
		for (wn = picture; *wn; wn++);
		while (wn > picture && wn[-1] != '/') wn--;
		if ((gwind = xfindwind(theDisplay, rwind, wn, 4)) == 0) {
			fprintf(stderr, "%s: cannot find display window!\n",
					progname);
			exit(1);
		}
	}
	XQueryPointer(theDisplay, gwind, &rw, &cw, &rx, &ry, &wx, &wy, &pm);
	xoff = wx - ipt[0];
	yoff = wy - ipt[1];
					/* set graphics context */
	if (vecGC == 0) {
		XGetWindowAttributes(theDisplay, gwind, &wa);
		xc.red = 65535; xc.green = 0; xc.blue = 0;
		xc.flags = DoRed|DoGreen|DoBlue;
		if (XAllocColor(theDisplay, wa.colormap, &xc)) {
			gcv.foreground = xc.pixel;
			vecGC = XCreateGC(theDisplay,gwind,GCForeground,&gcv);
		} else {
			gcv.function = GXinvert;
			vecGC = XCreateGC(theDisplay,gwind,GCFunction,&gcv);
		}
	}
}


static void
vector(		/* draw a vector */
	int	ip1[2],
	int	ip2[2]
)
{
	if (ip2[0] == -1 && ip2[1] == -1)
		return;			/* null vector */
	XDrawLine(theDisplay, gwind, vecGC,
			ip1[0]+xoff, ip1[1]+yoff,
			ip2[0]+xoff, ip2[1]+yoff);
	if (slow) {
		XFlush(theDisplay);
		sleep(1);
	}
}
