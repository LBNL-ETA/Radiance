#ifndef lint
static const char	RCSid[] = "$Id: xdebugamb.c,v 3.1 2003/02/22 02:07:28 greg Exp $";
#endif
---
> static char SCCSid[] = "@(#)ambient.c 2.29 5/2/95 LBL";
68c68
< #define tracktime	(shm_boundary == NULL || ambfp == NULL)
---
> #define tracktime	0
249a250
> printf("check %e %e %e\n", av->pos[0], av->pos[1], av->pos[2]);
289a291
> printf("use %e %e %e\n", av->pos[0], av->pos[1], av->pos[2]);
508a511
> printf("insert %e %e %e\n", av->pos[0], av->pos[1], av->pos[2]);
 */
/*
 *  Display an image and watch the rays get traced.
 *
 *	9/21/90	Greg Ward
 */

#include "standard.h"
#include "view.h"
#include  <time.h>
#include "resolu.h"
#include <ctype.h>
#include <X11/Xlib.h>

#define MAXDEPTH	32		/* ridiculous ray tree depth */

#ifdef  SMLFLT
#define  sscanvec(s,v)	(sscanf(s,"%f %f %f",v,v+1,v+2)==3)
#else
#define  sscanvec(s,v)	(sscanf(s,"%lf %lf %lf",v,v+1,v+2)==3)
#endif

char	rtcom[] = "rtrace -h- -otp -fa -x 1";
char	xicom[] = "ximage -c 128";

VIEW	ourview = STDVIEW;		/* view for picture */
RESOLU	ourres;				/* picture resolution */

char	*progname;			/* program name */

char	*picture;			/* picture name */
float	*zbuffer;

#define zval(x,y)	(zbuffer[(y)*scanlen(&ourres)+(x)])

FILE	*pin;				/* input stream */

Display	*theDisplay = NULL;		/* connection to server */

struct node {				/* ray tree node */
	int	ipt[2];
	struct node	*sister;
	struct node	*daughter;
};

#define newnode()	(struct node *)calloc(1, sizeof(struct node))

int	slow = 0;		/* slow trace? */


main(argc, argv)		/* takes both the octree and the image */
int	argc;
char	*argv[];
{
	int	i;
	char	combuf[256];

	progname = argv[0];
	for (i = 1; i < argc-3; i++)
		if (!strcmp(argv[i], "-s"))
			slow++;
		else
			break;
	if (i > argc-3) {
		fprintf(stderr, "Usage: %s [-s] [rtrace args] octree picture zbuffer\n",
				progname);
		exit(1);
	}
	picture = argv[argc-2];
					/* get the viewing parameters */
	if (viewfile(picture, &ourview, &ourres) <= 0 ||
			setview(&ourview) != NULL) {
		fprintf(stderr, "%s: cannot get view from \"%s\"\n",
				progname, picture);
		exit(1);
	}
	zbload(argv[argc-1]);
					/* open the display */
	if ((theDisplay = XOpenDisplay(NULL)) == NULL) {
		fprintf(stderr,
		"%s: cannot open display; DISPLAY variable set?\n",
				progname);
		exit(1);
	}
					/* build input command */
	sprintf(combuf, "%s %s | %s", xicom, picture, rtcom);
	for ( ; i < argc-2; i++) {
		strcat(combuf, " ");
		strcat(combuf, argv[i]);
	}
					/* start the damn thing */
	if ((pin = popen(combuf, "r")) == NULL)
		exit(1);
	sleep(20);
					/* loop on input */
	mainloop();
					/* close pipe and exit */
	pclose(pin);
	exit(0);
}


zbload(fname)
char *fname;
{
	int fd;
	zbuffer = (float *)malloc(ourres.xr*ourres.yr*sizeof(float));
	fd = open(fname,0);
	if (fd<0) exit(1);
	read(fd,zbuffer,ourres.xr*ourres.yr*sizeof(float));
	close(fd);
}


mainloop()				/* get and process input */
{
	static struct node	*sis[MAXDEPTH];
	register struct node	*newp;
	char	buf[128];
	int	level;
	register int	i;

	level = 0;
	while (fgets(buf, sizeof(buf), pin) != NULL) {
		if (isalpha(buf[0])) {
			i = 0;
			while (isalpha(buf[++i]))
				;
			buf[i++] = '\0';
			ambval(buf, buf+i);
			continue;
		}
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


freetree(tp)				/* free a trace tree */
struct node	*tp;
{
	register struct node	*kid;

	for (kid = tp->daughter; kid != NULL; kid = kid->sister)
		freetree(kid);
	free((void *)tp);
}


tracerays(tp)				/* trace a ray tree */
struct node	*tp;
{
	register struct node	*kid;

	for (kid = tp->daughter; kid != NULL; kid = kid->sister) {
		vector(tp->ipt, kid->ipt);
		tracerays(kid);
	}
}


strtoipt(ipt, str)		/* convert string x y z to image point */
int	ipt[2];
char	*str;
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

GC	insGC = 0;
GC	accGC = 0;
GC	useGC = 0;
GC	vecGC = 0;
Window	gwind = 0;
int	xoff=0, yoff=0;


setvec(ipt)			/* set up vector drawing for pick */
int	*ipt;
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
	if (ipt != NULL) {
		XQueryPointer(theDisplay, gwind, &rw, &cw, &rx, &ry, &wx, &wy, &pm);
		xoff = wx - ipt[0];
		yoff = wy - ipt[1];
	}
					/* set graphics contexts */
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
		xc.red = 0; xc.green = 0; xc.blue = 65535;
		xc.flags = DoRed|DoGreen|DoBlue;
		XAllocColor(theDisplay, wa.colormap, &xc);
		gcv.foreground = xc.pixel;
		insGC = XCreateGC(theDisplay,gwind,GCForeground,&gcv);
		xc.red = 32000; xc.green = 0; xc.blue = 10000;
		xc.flags = DoRed|DoGreen|DoBlue;
		XAllocColor(theDisplay, wa.colormap, &xc);
		gcv.foreground = xc.pixel;
		accGC = XCreateGC(theDisplay,gwind,GCForeground,&gcv);
		xc.red = 0; xc.green = 65000; xc.blue = 0;
		xc.flags = DoRed|DoGreen|DoBlue;
		XAllocColor(theDisplay, wa.colormap, &xc);
		gcv.foreground = xc.pixel;
		useGC = XCreateGC(theDisplay,gwind,GCForeground,&gcv);
	}
}


vector(ip1, ip2)		/* draw a vector */
int	ip1[2], ip2[2];
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


ambval(typ, str)
char	*typ, *str;
{
	int	ipt[2];
	FVECT	im_pt, pt;

	if (!sscanvec(str, pt))
		return(-1);
	viewloc(im_pt, &ourview, pt);
	loc2pix(ipt, &ourres, im_pt[0], im_pt[1]);
	if (ipt[0] < 0 || ipt[0] > scanlen(&ourres) ||
			ipt[1] < 0 || ipt[1] > numscans(&ourres))
		return;
	if (im_pt[2] > zval(ipt[0],ipt[1])+.05)
		return;				/* obstructed */
	setvec(NULL);
	if (!strcmp(typ, "insert"))
		XDrawRectangle(theDisplay, gwind, insGC,
				ipt[0]+xoff, ipt[1]+yoff, 2, 2);
	else if (!strcmp(typ, "check"))
		XDrawRectangle(theDisplay, gwind, accGC,
				ipt[0]+xoff, ipt[1]+yoff, 2, 2);
	else if (!strcmp(typ, "use"))
		XDrawRectangle(theDisplay, gwind, useGC,
				ipt[0]+xoff, ipt[1]+yoff, 2, 2);
	XFlush(theDisplay);
}
