#ifndef lint
static const char	RCSid[] = "$Id: rhd_glx.c,v 3.29 2003/02/22 02:07:24 greg Exp $";
#endif
/*
 * OpenGL GLX driver for holodeck display.
 * Based on old GLX driver using cones.
 *
 * Define symbol STEREO for stereo viewing.
 * Define symbol DOBJ for display object viewing.
 */

#ifdef NOSTEREO
#ifdef STEREO
#undef STEREO
#else
#undef NOSTEREO
#endif
#endif

#include "standard.h"

#include <GL/glx.h>
#include <GL/glu.h>
#ifdef STEREO
#include <X11/extensions/SGIStereo.h>
#endif

#include "rhd_sample.h"
#ifdef DOBJ
#include "rhdobj.h"
#endif

#include "x11icon.h"

#ifndef RAYQLEN
#define RAYQLEN		1024		/* max. rays to queue before flush */
#endif

#ifndef FEQ
#define FEQ(a,b)	((a)-(b) <= FTINY && (a)-(b) >= -FTINY)
#endif

#define GAMMA		1.4		/* default gamma correction */

#define FRAMESTATE(s)	(((s)&(ShiftMask|ControlMask))==(ShiftMask|ControlMask))

#define MOVPCT		7		/* percent distance to move /frame */
#define MOVDIR(b)	((b)==Button1 ? 1 : (b)==Button2 ? 0 : -1)
#define MOVDEG		(-5)		/* degrees to orbit CW/down /frame */
#define MOVORB(s)	((s)&ShiftMask ? 1 : (s)&ControlMask ? -1 : 0)

#ifndef TARGETFPS
#define TARGETFPS	2.0		/* target frames/sec during motion */
#endif

#define MINWIDTH	480		/* minimum graphics window width */
#define MINHEIGHT	400		/* minimum graphics window height */

#define VIEWDIST	356		/* assumed viewing distance (mm) */

#define BORWIDTH	5		/* border width */

#define setstereobuf(bid)	(glXWaitGL(), \
				XSGISetStereoBuffer(ourdisplay, gwind, bid), \
				glXWaitX())

#define  ourscreen	DefaultScreen(ourdisplay)
#define  ourroot	RootWindow(ourdisplay,ourscreen)
#define  ourmask	(StructureNotifyMask|ExposureMask|KeyPressMask|\
			ButtonPressMask|ButtonReleaseMask)

#define  levptr(etype)	((etype *)&currentevent)

struct driver	odev;			/* global device driver structure */

char odev_args[64];			/* command arguments */

#ifdef STEREO
static VIEW	vwright;		/* right eye view */
#endif

static int	rayqleft = 0;		/* rays left to queue before flush */

static XEvent  currentevent;		/* current event */

static int  mapped = 0;			/* window is mapped? */
static unsigned long  ourblack=0, ourwhite=~0;

static Display  *ourdisplay = NULL;	/* our display */
static XVisualInfo  *ourvinf;		/* our visual information */
static Window  gwind = 0;		/* our graphics window */
static GLXContext	gctx;		/* our GLX context */

static double	pwidth, pheight;	/* pixel dimensions (mm) */

static double	mindpth, maxdpth;	/* min. and max. depth */

double	dev_zmin, dev_zmax;		/* fore and aft clipping plane dist. */

static int	inpresflags;		/* input result flags */

static int	headlocked = 0;		/* lock vertical motion */

static int  resizewindow(), getevent(), getkey(), moveview(), wipeclean(),
		setglpersp(), getframe(), getmove(), fixwindow(), mytmflags();

#ifdef STEREO
static int  pushright(), popright();
#endif
static double	getdistance();
#define mapdepth(d) ((d)> 0.9995 ? FHUGE: dev_zmin/ \
		     (1.-(d)*(1.-dev_zmin/dev_zmax)))

extern time_t	time();


dev_open(id)			/* initialize GLX driver */
char  *id;
{
	extern char	*getenv();
	static RGBPRIMS	myprims = STDPRIMS;
	static int	atlBest[] = {GLX_RGBA, GLX_RED_SIZE,8,
				GLX_GREEN_SIZE,8, GLX_BLUE_SIZE,8,
				GLX_DEPTH_SIZE,15, None};
	char	*ev;
	double	gamval = GAMMA;
	RGBPRIMP	dpri = stdprims;
	XSetWindowAttributes	ourwinattr;
	XWMHints	ourxwmhints;
	XSizeHints	oursizhints;
					/* check for unsupported stereo */
#ifdef NOSTEREO
	error(INTERNAL, "stereo display driver unavailable");
#endif
					/* open display server */
	ourdisplay = XOpenDisplay(NULL);
	if (ourdisplay == NULL)
		error(USER, "cannot open X-windows; DISPLAY variable set?\n");
#ifdef STEREO
	switch (XSGIQueryStereoMode(ourdisplay, ourroot)) {
	case STEREO_TOP:
	case STEREO_BOTTOM:
		break;
	case STEREO_OFF:
		error(USER,
	"wrong video mode: run \"/usr/gfx/setmon -n STR_TOP\" first");
	case X_STEREO_UNSUPPORTED:
		error(USER, "stereo mode not supported on this screen");
	default:
		error(INTERNAL, "unknown stereo mode");
	}
#endif
					/* find a usable visual */
	ourvinf = glXChooseVisual(ourdisplay, ourscreen, atlBest);
	if (ourvinf == NULL)
		error(USER, "no suitable visuals available");
					/* get a context */
	gctx = glXCreateContext(ourdisplay, ourvinf, NULL, GL_TRUE);
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
	if (tmInit(mytmflags(), dpri, gamval) == NULL)
		error(SYSTEM, "not enough memory in dev_open");
					/* open window */
	ourwinattr.background_pixel = ourblack;
	ourwinattr.border_pixel = ourblack;
	ourwinattr.event_mask = ourmask;
					/* this is stupid */
	ourwinattr.colormap = XCreateColormap(ourdisplay, ourroot,
				ourvinf->visual, AllocNone);
	gwind = XCreateWindow(ourdisplay, ourroot, 0, 0,
		DisplayWidth(ourdisplay,ourscreen)-2*BORWIDTH,
#ifdef STEREO
		(DisplayHeight(ourdisplay,ourscreen)-2*BORWIDTH)/2,
#else
		DisplayHeight(ourdisplay,ourscreen)-2*BORWIDTH,
#endif
		BORWIDTH, ourvinf->depth, InputOutput, ourvinf->visual,
		CWBackPixel|CWBorderPixel|CWColormap|CWEventMask, &ourwinattr);
	if (gwind == 0)
		error(SYSTEM, "cannot create window\n");
   	XStoreName(ourdisplay, gwind, id);
					/* set window manager hints */
	ourxwmhints.flags = InputHint|IconPixmapHint;
	ourxwmhints.input = True;
	ourxwmhints.icon_pixmap = XCreateBitmapFromData(ourdisplay,
			gwind, x11icon_bits, x11icon_width, x11icon_height);
	XSetWMHints(ourdisplay, gwind, &ourxwmhints);
	oursizhints.min_width = MINWIDTH;
#ifdef STEREO
	oursizhints.min_height = MINHEIGHT/2;
	oursizhints.max_width = DisplayWidth(ourdisplay,ourscreen)-2*BORWIDTH;
	oursizhints.max_height = (DisplayHeight(ourdisplay,ourscreen) -
					2*BORWIDTH)/2;
	oursizhints.flags = PMinSize|PMaxSize;
#else
	oursizhints.min_height = MINHEIGHT;
	oursizhints.flags = PMinSize;
#endif
	XSetNormalHints(ourdisplay, gwind, &oursizhints);
					/* set GLX context */
	glXMakeCurrent(ourdisplay, gwind, gctx);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_DITHER);
	glDisable(GL_CULL_FACE);
					/* figure out sensible view */
	pwidth = (double)DisplayWidthMM(ourdisplay, ourscreen) /
			DisplayWidth(ourdisplay, ourscreen);
	pheight = (double)DisplayHeightMM(ourdisplay, ourscreen) /
			DisplayHeight(ourdisplay, ourscreen);
#ifdef STEREO
	pheight *= 2.;
	setstereobuf(STEREO_BUFFER_LEFT);
#endif
	checkglerr("setting rendering parameters");
	copystruct(&odev.v, &stdview);
	odev.v.type = VT_PER;
					/* map the window */
	XMapWindow(ourdisplay, gwind);
	dev_input();			/* sets size and view angles */
					/* allocate our samples */
	if (!smInit(DisplayWidth(ourdisplay,ourscreen) *
			DisplayHeight(ourdisplay,ourscreen) / 10))
		error(SYSTEM, "insufficient memory for value storage");
	mindpth = FHUGE; maxdpth = FTINY;
	odev.name = id;
	odev.ifd = ConnectionNumber(ourdisplay);
}


dev_close()			/* close our display and free resources */
{
	smInit(0);
#ifdef DOBJ
	dobj_cleanup();
#endif
	glXMakeCurrent(ourdisplay, None, NULL);
	glXDestroyContext(ourdisplay, gctx);
	XDestroyWindow(ourdisplay, gwind);
	gwind = 0;
	XCloseDisplay(ourdisplay);
	ourdisplay = NULL;
	tmDone(NULL);
	odev.v.type = 0;
	odev.hres = odev.vres = 0;
	odev.ifd = -1;
}


dev_clear()			/* clear our representation */
{
	smInit(rsL.max_samp);
	wipeclean(1);
	rayqleft = 0;			/* hold off update */
}


int
dev_view(nv)			/* assign new driver view */
register VIEW	*nv;
{
	double	d;

	if (nv->type != VT_PER ||		/* check view legality */
			nv->horiz > 160. || nv->vert > 160.) {
		error(COMMAND, "illegal view type/angle");
		nv->type = odev.v.type;
		nv->horiz = odev.v.horiz;
		nv->vert = odev.v.vert;
		return(0);
	}
	if (nv != &odev.v) {
						/* resize window? */
		if (!FEQ(nv->horiz,odev.v.horiz) ||
				!FEQ(nv->vert,odev.v.vert)) {
			int	dw = DisplayWidth(ourdisplay,ourscreen);
			int	dh = DisplayHeight(ourdisplay,ourscreen);

			dw -= 25;	/* for window frame */
			dh -= 50;
#ifdef STEREO
			dh /= 2;
#endif
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
			dev_input();	/* get resize event */
		}
		copystruct(&odev.v, nv);	/* setview() already called */
	}
#ifdef STEREO
	copystruct(&vwright, nv);
	d = eyesepdist / sqrt(nv->hn2);
	VSUM(vwright.vp, nv->vp, nv->hvec, d);
	/* setview(&vwright);	-- Unnecessary */
#endif
	wipeclean(0);
	return(1);
}


dev_section(gfn, pfn)		/* add octree for geometry rendering */
char	*gfn, *pfn;
{
	/* unimplemented */
}


dev_auxcom(cmd, args)		/* process an auxiliary command */
char	*cmd, *args;
{
#ifdef DOBJ
	if (dobj_command(cmd, args) >= 0)
		return;
#endif
	sprintf(errmsg, "%s: unknown command", cmd);
	error(COMMAND, errmsg);
}


VIEW *
dev_auxview(n, hvres)		/* return nth auxiliary view */
int	n;
int	hvres[2];
{
	hvres[0] = odev.hres; hvres[1] = odev.vres;
	if (n == 0)
		return(&odev.v);
#ifdef STEREO
	if (n == 1)
		return(&vwright);
#endif
	return(NULL);
}


int
dev_input()			/* get X11 input */
{
	inpresflags = 0;

	do
		getevent();

	while (XPending(ourdisplay) > 0);

	odev.inpready = 0;

	return(inpresflags);
}


dev_value(c, d, p)		/* add a pixel value to our mesh */
COLR	c;
FVECT	d, p;
{
	double	depth;
#ifdef DOBJ
	if (dobj_lightsamp != NULL) {	/* in light source sampling */
		(*dobj_lightsamp)(c, d, p);
		return;
	}
#endif
	if (p != NULL) {		/* add depth to our range */
		depth = (p[0] - odev.v.vp[0])*d[0] +
			(p[1] - odev.v.vp[1])*d[1] +
			(p[2] - odev.v.vp[2])*d[2];
		if (depth > FTINY) {
			if (depth < mindpth)
				mindpth = depth;
			if (depth > maxdpth)
				maxdpth = depth;
		}
	}
	smNewSamp(c, d, p);		/* add to display representation */
	if (!--rayqleft)
		dev_flush();		/* flush output */
}


int
dev_flush()			/* flush output */
{
	if (mapped) {
#ifdef STEREO
		pushright();			/* update right eye */
		smUpdate(&vwright, MAXQUALITY);
#ifdef DOBJ
		dobj_render();			/* usually in foreground */
#endif
		popright();			/* update left eye */
#endif
		smUpdate(&odev.v, MAXQUALITY);
		checkglerr("rendering mesh");
#ifdef DOBJ
		dobj_render();
#endif
		glFlush();			/* flush OGL */
	}
	rayqleft = RAYQLEN;
					/* flush X11 and return # pending */
	return(odev.inpready = XPending(ourdisplay));
}


checkglerr(where)		/* check for GL or GLU error */
char	*where;
{
	register GLenum	errcode;

	while ((errcode = glGetError()) != GL_NO_ERROR) {
		sprintf(errmsg, "OpenGL error %s: %s",
				where, gluErrorString(errcode));
		error(WARNING, errmsg);
	}
}


#ifdef STEREO
static
pushright()			/* push on right view */
{
	double	d;

	setstereobuf(STEREO_BUFFER_RIGHT);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	d = -eyesepdist / sqrt(odev.v.hn2);
	glTranslated(d*odev.v.hvec[0], d*odev.v.hvec[1], d*odev.v.hvec[2]);
	checkglerr("setting right view");
}


static
popright()			/* pop off right view */
{
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	setstereobuf(STEREO_BUFFER_LEFT);
}
#endif


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
	if (cp > tail && cp[-1] == 'h')
		return(TM_F_HUMAN|TM_F_NOSTDERR);
	else
		return(TM_F_CAMERA|TM_F_NOSTDERR);
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
		break;
	case MapNotify:
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


static
draw3dline(wp)			/* draw 3d line in world coordinates */
register FVECT	wp[2];
{
	glVertex3d(wp[0][0], wp[0][1], wp[0][2]);
	glVertex3d(wp[1][0], wp[1][1], wp[1][2]);
}


static
draw_grids(fore)		/* draw holodeck section grids */
int	fore;
{
	if (fore)
		glColor4ub(0, 255, 255, 0);
	else
		glColor4ub(0, 0, 0, 0);
	glBegin(GL_LINES);		/* draw each grid line */
	gridlines(draw3dline);
	glEnd();
	checkglerr("drawing grid lines");
}

static double
getdistance(dx, dy, direc)	/* distance from fore plane along view ray */
int	dx, dy;
FVECT	direc;
{
	GLfloat	gldepth;
	double	dist;

	if (dx<0 | dx>=odev.hres | dy<0 | dy>=odev.vres)
		return(FHUGE);
	glReadPixels(dx,dy, 1,1, GL_DEPTH_COMPONENT,GL_FLOAT, &gldepth);
	dist = mapdepth(gldepth);
	if (dist >= .99*FHUGE)
		return(FHUGE);
	return((dist-odev.v.vfore)/DOT(direc,odev.v.vdir));
}


static
moveview(dx, dy, mov, orb)	/* move our view */
int	dx, dy, mov, orb;
{
	VIEW	nv;
	FVECT	odir, v1, wip;
	double	d,d1;
	register int	li;
				/* start with old view */
	copystruct(&nv, &odev.v);
				/* orient our motion */
	if (viewray(v1, odir, &odev.v,
			(dx+.5)/odev.hres, (dy+.5)/odev.vres) < -FTINY)
		return(0);		/* outside view */
	if (mov | orb) {	/* moving relative to geometry */
	        d = getdistance(dx, dy, odir);/*distance from front plane */
#ifdef DOBJ
		d1 = dobj_trace(NULL, v1, odir);	/* check objects */
						/* check holodeck */
		if (d1 < d)
			d = d1;
#endif
		if (d >= .99*FHUGE)
			d = 0.5*(dev_zmax+dev_zmin);	/* just guess */
		VSUM(wip, v1, odir, d);
		VSUB(odir, wip, odev.v.vp);
#if 0
		fprintf(stderr, "moveview: hit %s at (%f,%f,%f) (t=%f)\n",
				li < 0 ? "object" : "mesh",
				wip[0], wip[1], wip[2],
				(wip[0]-odev.v.vp[0])*odir[0] +
				(wip[1]-odev.v.vp[1])*odir[1] +
				(wip[2]-odev.v.vp[2])*odir[2]);
#endif
	} else			/* panning with constant viewpoint */
		VCOPY(nv.vdir, odir);
	if (orb && mov) {		/* orbit left/right */
		spinvector(odir, odir, nv.vup, d=MOVDEG*PI/180.*mov);
		VSUM(nv.vp, wip, odir, -1.);
		spinvector(nv.vdir, nv.vdir, nv.vup, d);
	} else if (orb) {		/* orbit up/down */
		fcross(v1, odir, nv.vup);
		if (normalize(v1) == 0.)
			return(0);
		spinvector(odir, odir, v1, d=MOVDEG*PI/180.*orb);
		VSUM(nv.vp, wip, odir, -1.);
		spinvector(nv.vdir, nv.vdir, v1, d);
	} else if (mov) {		/* move forward/backward */
		d = MOVPCT/100. * mov;
		VSUM(nv.vp, nv.vp, odir, d);
	}
	if (!mov ^ !orb && headlocked) {	/* restore head height */
		VSUM(v1, odev.v.vp, nv.vp, -1.);
		d = DOT(v1, nv.vup);
		VSUM(nv.vp, nv.vp, odev.v.vup, d);
	}
	if (setview(&nv) != NULL)
		return(0);	/* illegal view */
	dev_view(&nv);
	inpresflags |= DFL(DC_SETVIEW);
	return(1);
}


static
getframe(ebut)				/* get focus frame */
XButtonPressedEvent	*ebut;
{
	int	startx = ebut->x, starty = ebut->y;
	int	endx, endy;

	XMaskEvent(ourdisplay, ButtonReleaseMask, levptr(XEvent));
	endx = levptr(XButtonReleasedEvent)->x;
	endy = levptr(XButtonReleasedEvent)->y;
	if (endx == startx | endy == starty) {
		XBell(ourdisplay, 0);
		return;
	}
	if (endx < startx) {register int c = endx; endx = startx; startx = c;}
	if (endy < starty) {register int c = endy; endy = starty; starty = c;}
	sprintf(odev_args, "%.3f %.3f %.3f %.3f",
			(startx+.5)/odev.hres, 1.-(endy+.5)/odev.vres,
			(endx+.5)/odev.hres, 1.-(starty+.5)/odev.vres);
	inpresflags |= DFL(DC_FOCUS);
}


static
getmove(ebut)				/* get view change */
XButtonPressedEvent	*ebut;
{
	int	movdir = MOVDIR(ebut->button);
	int	movorb = MOVORB(ebut->state);
	int	qlevel = MAXQUALITY-1;
	time_t	lasttime, thistime;
	int	nframes;
	Window	rootw, childw;
	int	rootx, rooty, wx, wy;
	unsigned int	statemask;

	XNoOp(ourdisplay);		/* makes sure we're not idle */

	nframes = 0;
	while (!XCheckMaskEvent(ourdisplay,
			ButtonReleaseMask, levptr(XEvent))) {
					/* get cursor position */
		if (!XQueryPointer(ourdisplay, gwind, &rootw, &childw,
				&rootx, &rooty, &wx, &wy, &statemask))
			break;		/* on another screen */

		draw_grids(0);		/* clear old grid lines */
#ifdef STEREO
		pushright(); draw_grids(0); popright();
#endif
					/* compute view motion */
		if (!moveview(wx, odev.vres-1-wy, movdir, movorb)) {
			sleep(1);
			lasttime++;
			continue;	/* cursor in bad place */
		}
		draw_grids(1);		/* redraw grid */
#ifdef STEREO
		pushright();
		draw_grids(1);
		smUpdate(&vwright, qlevel);
#ifdef DOBJ
		dobj_render();
#endif
		popright();
#endif
					/* redraw mesh */
		smUpdate(&odev.v, qlevel);
#ifdef DOBJ
		dobj_render();		/* redraw object */
#endif
		glFlush();
					/* figure out good quality level */
		if (!nframes++) {		/* ignore first frame */
			lasttime = time(0);
			continue;
		}
		thistime = time(0);
		if (thistime - lasttime >= 6 ||
				nframes > (int)(3*3*TARGETFPS)) {
			qlevel = thistime<=lasttime ? 3*MAXQUALITY :
				(int)((double)nframes/(thistime-lasttime)
					/ TARGETFPS * qlevel + 0.5);
			nframes = 0;
			if (qlevel >= MAXQUALITY) {
				if (qlevel >= 3*MAXQUALITY) {	/* brakes!! */
					sleep(1);
					lasttime++;
				}
				qlevel = MAXQUALITY-1;
			} else if (qlevel < 1)
				qlevel = 1;
		}
	}
	if (!(inpresflags & DFL(DC_SETVIEW))) {	/* do final motion */
		movdir = MOVDIR(levptr(XButtonReleasedEvent)->button);
		wx = levptr(XButtonReleasedEvent)->x;
		wy = levptr(XButtonReleasedEvent)->y;
		moveview(wx, odev.vres-1-wy, movdir, movorb);
	}
}


static
setglpersp(vp)			/* set perspective view in GL */
register VIEW	*vp;
{
	double	d, xmin, xmax, ymin, ymax;

	if (mindpth >= maxdpth) {
		dev_zmin = 1.;
		dev_zmax = 100.;
	} else {
		dev_zmin = 0.1*mindpth;
		dev_zmax = 5.0*maxdpth;
	}
	if (odev.v.vfore > FTINY)
		dev_zmin = odev.v.vfore;
	if (odev.v.vaft > FTINY)
		dev_zmax = odev.v.vaft;
	if (dev_zmin*500. < dev_zmax)
		dev_zmax = dev_zmin*500.;
	xmax = dev_zmin * tan(PI/180./2. * odev.v.horiz);
	xmin = -xmax;
	d = odev.v.hoff * (xmax - xmin);
	xmin += d; xmax += d;
	ymax = dev_zmin * tan(PI/180./2. * odev.v.vert);
	ymin = -ymax;
	d = odev.v.voff * (ymax - ymin);
	ymin += d; ymax += d;
					/* set view matrix */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(xmin, xmax, ymin, ymax, dev_zmin, dev_zmax);
	gluLookAt(odev.v.vp[0], odev.v.vp[1], odev.v.vp[2],
		odev.v.vp[0] + odev.v.vdir[0],
		odev.v.vp[1] + odev.v.vdir[1],
		odev.v.vp[2] + odev.v.vdir[2],
		odev.v.vup[0], odev.v.vup[1], odev.v.vup[2]);
	checkglerr("setting perspective view");
}


static
wipeclean(tmflag)			/* prepare for redraw */
int	tmflag;
{
					/* clear depth buffer */
#ifdef STEREO
	setstereobuf(STEREO_BUFFER_RIGHT);
	glClear(GL_DEPTH_BUFFER_BIT);
	setstereobuf(STEREO_BUFFER_LEFT);
#endif
	glClear(GL_DEPTH_BUFFER_BIT);
	smClean(tmflag);		/* reset drawing routines */
	setglpersp(&odev.v);		/* reset view & clipping planes */
}


static
getkey(ekey)				/* get input key */
register XKeyPressedEvent  *ekey;
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
	case 'p':			/* pause computation */
		inpresflags |= DFL(DC_PAUSE);
		return;
	case 'v':			/* spit out view */
		inpresflags |= DFL(DC_GETVIEW);
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
	case '\n':
	case '\r':			/* resume computation */
		inpresflags |= DFL(DC_RESUME);
		return;
	case CTRL('R'):			/* redraw screen */
		wipeclean(1);
		return;
	case CTRL('L'):			/* refresh from server */
		if (inpresflags & DFL(DC_REDRAW))
			return;
		setglpersp(&odev.v);		/* reset clipping planes */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDisable(GL_DEPTH_TEST);	/* so grids will clear */
		draw_grids(1);
#ifdef STEREO
		pushright();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		draw_grids(1);
		popright();
#endif
		glEnable(GL_DEPTH_TEST);
		glFlush();
		smInit(rsL.max_samp);		/* get rid of old values */
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


static
fixwindow(eexp)				/* repair damage to window */
register XExposeEvent  *eexp;
{
	if (odev.hres == 0 | odev.vres == 0) {	/* first exposure */
		resizewindow((XConfigureEvent *)eexp);
		return;
	}
	if (eexp->count)		/* wait for final exposure */
		return;
	wipeclean(0);			/* clear depth */
}


static
resizewindow(ersz)			/* resize window */
register XConfigureEvent  *ersz;
{
	glViewport(0, 0, ersz->width, ersz->height);

	if (ersz->width == odev.hres && ersz->height == odev.vres)
		return;

	odev.hres = ersz->width;
	odev.vres = ersz->height;

	odev.v.horiz = 2.*180./PI * atan(0.5/VIEWDIST*pwidth*odev.hres);
	odev.v.vert = 2.*180./PI * atan(0.5/VIEWDIST*pheight*odev.vres);

	inpresflags |= DFL(DC_SETVIEW);
}
