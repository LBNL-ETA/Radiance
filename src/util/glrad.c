#ifndef lint
static const char	RCSid[] = "$Id: glrad.c,v 3.19 2004/03/26 23:34:23 schorsch Exp $";
#endif
/*
 * Program to display Radiance scene using OpenGL.
 */

#include <sys/types.h>
#include <GL/glx.h>
#ifndef NOSTEREO
#include <X11/extensions/SGIStereo.h>
#endif
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "radogl.h"
#include "view.h"
#include "paths.h"
#include "glradicon.h"
#include "rtprocess.h"

#ifndef MAXVIEW
#define MAXVIEW		63		/* maximum number of standard views */
#endif
#ifndef MAXSCENE
#define MAXSCENE	127		/* maximum number of scene files */
#endif

#define ZOOMPCT		9		/* percent to zoom at a time */

#define MOVPCT		7		/* percent distance to move /frame */
#define MOVDIR(b)	((b)==Button1 ? 1 : (b)==Button2 ? 0 : -1)
#define MOVDEG		(-5)		/* degrees to orbit CW/down /frame */
#define MOVORB(s)	((s)&ShiftMask ? 1 : (s)&ControlMask ? -1 : 0)

#define BORWIDTH	5		/* border width */

#define  ourscreen	DefaultScreen(ourdisplay)
#define  ourroot	RootWindow(ourdisplay,ourscreen)
#define  ourmask	(StructureNotifyMask|ExposureMask|KeyPressMask|\
			ButtonPressMask|ButtonReleaseMask)

#define  levptr(etype)	((etype *)&currentevent)

XEvent  currentevent;			/* current event */

int  mapped = 0;			/* window is mapped? */
unsigned long  ourblack=0, ourwhite=~0;

Display  *ourdisplay = NULL;		/* our display */
XVisualInfo  *ourvinf;			/* our visual information */
Window  gwind = 0;			/* our graphics window */
int	hres, vres;			/* rendering window dimensions */
int	maxhres, maxvres;		/* maximum given dimensions */
GLXContext	gctx;			/* our GLX context */

double	pwidth, pheight;		/* pixel dimensions (mm) */

int	headlocked = 0;			/* lock vertical motion */

struct {
	char	*nam;				/* view name (NULL if none) */
	VIEW	*v;				/* parameters (NULL term.) */
} vwl[MAXVIEW+1];			/* our list of views */

int	currentview = 0;		/* current view number */
VIEW	thisview = STDVIEW;		/* displayed view */
double	eyedist = 1;			/* interocular distance */
VIEW	lastview;			/* last recorded view */

char	*progname;			/* global argv[0] */
char	*radfile;			/* rad input file */
char	*scene[MAXSCENE+1];		/* material and scene file list */
int	nscenef = 0;			/* number of scene files */
char	*octree;			/* octree name (NULL if unnec.) */

SUBPROC	rtpd;			/* rtrace process descriptors */

int	silent = 0;			/* run rad silently? */
int	backvis = 1;			/* back faces visible? */
int	stereo = 0;			/* do stereo? */

#ifdef NOSTEREO
#define setstereobuf(bid)	
#else
#define setstereobuf(bid)	(glXWaitGL(), \
				XSGISetStereoBuffer(ourdisplay, gwind, bid), \
				glXWaitX())
#endif

int	displist;			/* our scene display list */

int	no_render = 0;			/* don't rerender */

extern int	nowarn;			/* turn warnings off? */

static void startrtrace(char	*octname);
static void runrad(int	ac, char	**av);
static int findvw(register char	*nm);
static int varmatch(register char	*s, register char	*vn);
static char * scan4var(char	*buf, int	buflen, char	*vname, FILE	*fp);
static void dev_open(char  *id);
static void dev_close(void);
static int dev_view(register VIEW	*nv);
static int dev_input(int	nsecs);
static void render(void);
static int moveview(int	dx, int	dy, int	mov, int	orb);
static void waitabit(void);
static void getmove(XButtonPressedEvent	*ebut);
static int getintersect(FVECT	wp, FVECT	org, FVECT	dir, double	md);
static void setglpersp(register VIEW	*vp);
static int getkey(register XKeyPressedEvent  *ekey);
static void zoomview(int	pct, int	dx, int	dy);
static void gotoview(int	vwnum);
static void appendview(char	*nm, VIEW	*vp);
static void copylastv(char	*cause);
static void fixwindow(register XExposeEvent  *eexp);
static void resizewindow(register XConfigureEvent  *ersz);


int
main(
	int	argc,
	char	*argv[]
)
{
	char	*viewsel = NULL;
	long	vwintvl = 0;
	int	i;

	progname = argv[0];
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'v':
			viewsel = argv[++i];
			break;
		case 'w':
			nowarn = !nowarn;
			break;
		case 's':
			silent = !silent;
			break;
		case 'S':
			stereo = !stereo;
			break;
		case 'c':
			vwintvl = atoi(argv[++i]);
			break;
		case 'b':
			backvis = !backvis;
			break;
		default:
			goto userr;
		}
	if (i >= argc)
		goto userr;
#ifdef NOSTEREO
	if (stereo)
		error(INTERNAL, "stereo not supported in this version");
#endif
					/* run rad and get views */
	runrad(argc-i, argv+i);
					/* check view */
	if (viewsel != NULL) {
		if (viewsel[0] == '-') {
			char	*ve = viewsel;
			if (!sscanview(&thisview, viewsel) ||
					(ve = setview(&thisview)) != NULL) {
				fprintf(stderr, "%s: bad view: %s\n",
						progname, ve);
				quit(1);
			}
			currentview = -1;
		} else if ((currentview = findvw(viewsel)) < 0) {
			fprintf(stderr, "%s: no such view: %s\n",
						progname, viewsel);
			quit(1);
		}
	}
					/* open GL */
	dev_open(radfile = argv[i]);
					/* load octree or scene files */
	if (octree != NULL) {
		displist = rgl_octlist(octree, NULL, NULL, NULL);
		startrtrace(octree);
	} else
		displist = rgl_filelist(nscenef, scene, NULL);
					/* set initial view */
	dev_view(currentview < 0 ? &thisview : vwl[currentview].v);
					/* input/render loop */
	while (dev_input(vwintvl))
		;
					/* all done */
	quit(0);
userr:
	fprintf(stderr,
		"Usage: %s [-w][-s][-b][-S][-v view] rfile [VAR=value]..\n",
			argv[0]);
	quit(1);
	return 1; /* pro forma return */
}


void
quit(				/* exit gracefully */
	int	code
)
{
	if (ourdisplay != NULL)
		dev_close();
	/* if (rtpd.pid > 0) { */
	if (rtpd.running) {
		if (close_process(&rtpd) > 0)
			wputs("bad exit status from rtrace\n");
		/* rtpd.pid = 0; */
	}
	exit(code);
}


static void
startrtrace(			/* start rtrace on octname */
	char	*octname
)
{
	static char	*av[12] = {"rtrace", "-h", "-fff", "-ld+",
					"-opL", "-x", "1"};
	int	ac = 7;

	if (nowarn) av[ac++] = "-w-";
	av[ac++] = octname;
	av[ac] = NULL;
	if (open_process(&rtpd, av) <= 0)
		error(SYSTEM, "cannot start rtrace process");
}


static void
runrad(				/* run rad and load variables */
	int	ac,
	char	**av
)
{
	static char	optfile[] = TEMPLATE;
	int	nvn = 0, nvv = 0;
	FILE	*fp;
	register char	*cp;
	char	radcomm[256], buf[128], nam[32];
					/* set rad commmand */
	strcpy(radcomm, "rad -w -v 0        ");	/* look out below! */
	cp = radcomm + 19;
	if (silent) {
		strcpy(cp, "-s ");
		cp += 3;
	}
	while (ac--) {
		strcpy(cp, *av++);
		while (*cp) cp++;
		*cp++ = ' ';
	}
	strcpy(cp, "OPTFILE=");		/* create temporary options file */
	strcpy(cp+8, mktemp(optfile));
	if (system(radcomm))		/* update octree */
		error(USER, "error executing rad command");
					/* replace "-v 0" with "-n -e -s -V" */
	strcpy(radcomm+7, "-n -e -s -V");
	radcomm[18] = ' ';
	if ((fp = popen(radcomm, "r")) == NULL)
		error(SYSTEM, "cannot start rad command");
	buf[0] = '\0';			/* read variables alphabetically */
						/* get exposure */
	if ((cp = scan4var(buf, sizeof(buf), "EXPOSURE", fp)) != NULL) {
		expval = atof(cp);
		if ((*cp == '-') | (*cp == '+'))
			expval = pow(2., expval);
		expval *= 0.5;		/* compensate for local shading */
	}
						/* look for eye separation */
	if ((cp = scan4var(buf, sizeof(buf), "EYESEP", fp)) != NULL)
		eyedist = atof(cp);
						/* look for materials */
	while ((cp = scan4var(buf, sizeof(buf), "materials", fp)) != NULL) {
		nscenef += wordstring(scene+nscenef, cp);
		buf[0] = '\0';
	}
						/* look for octree */
	if ((cp = scan4var(buf, sizeof(buf), "OCTREE", fp)) != NULL)
		octree = savqstr(cp);
						/* look for scene files */
	while ((cp = scan4var(buf, sizeof(buf), "scene", fp)) != NULL) {
		nscenef += wordstring(scene+nscenef, cp);
		buf[0] = '\0';
	}
						/* load view names */
	while ((cp = scan4var(buf, sizeof(buf), "view", fp)) != NULL) {
		if (nvn >= MAXVIEW)
			error(INTERNAL, "too many views in rad file");
		vwl[nvn++].nam = *cp == '-' ? (char *)NULL :
				savqstr(atos(nam, sizeof(nam), cp));
		buf[0] = '\0';
	}
						/* load actual views */
	do
		if (isview(buf)) {
			vwl[nvv].v = (VIEW *)bmalloc(sizeof(VIEW));
			*(vwl[nvv].v) = stdview;
			sscanview(vwl[nvv].v, buf);
			if ((cp = setview(vwl[nvv++].v)) != NULL) {
				fprintf(stderr, "%s: bad view %d - %s\n",
						progname, nvv, cp);
				quit(1);
			}
		}
	while (fgets(buf, sizeof(buf), fp) != NULL);
	if (nvv != nvn)
		error(INTERNAL, "view miscount in runrad");
	pclose(fp);
						/* open options file */
	if ((fp = fopen(optfile, "r")) == NULL)
		error(SYSTEM, "cannot open options file");
						/* get relevant options */
	while (fgets(buf, sizeof(buf), fp) != NULL)
		if (!strncmp(buf, "-av ", 4))
			setcolor(ambval, atof(buf+4),
					atof(sskip2(buf+4,1)),
					atof(sskip2(buf+4,2)));
		else if (backvis && !strncmp(buf, "-bv", 3) &&
				(!buf[3] || strchr("0-FfNn \n",buf[3])!=NULL))
			backvis = 0;
	fclose(fp);
	unlink(optfile);			/* delete options file */
}


static int
findvw(			/* find named view */
	register char	*nm
)
{
	register int	n;

	if ((*nm >= '1') & (*nm <= '9') &&
			(n = atoi(nm)-1) <= MAXVIEW && vwl[n].v != NULL)
		return(n);
	for (n = 0; vwl[n].v != NULL; n++)
		if (vwl[n].nam != NULL && !strcmp(nm, vwl[n].nam))
			return(n);
	return(-1);
}


static int
varmatch(				/* match line to variable */
	register char	*s,
	register char	*vn
)
{
	register int	c;

	for ( ; *vn && *s == *vn; s++, vn++)
		;
	while (isspace(*s))
		s++;
	if (*s == '=')
		return(*vn);
	while (!(c = toupper(*s++) - toupper(*vn)) && *vn++)
		;
	return(c);
}


static char *
scan4var(	/* scan for variable from fp */
	char	*buf,
	int	buflen,
	char	*vname,
	FILE	*fp
)
{
	int	cval;
	register char	*cp;
					/* search out matching line */
	while ((cval = varmatch(buf, vname))) {
		if (cval > 0)			/* gone too far? */
			return(NULL);
		buf[0] = '\0';			/* else get next line */
		if (fgetline(buf, buflen, fp) == NULL)
			return(NULL);
	}
					/* skip variable name and '=' */
	for (cp = buf; *cp++ != '='; )
		;
	while (isspace(*cp)) cp++;
	return(cp);
}


static void
dev_open(			/* initialize GLX driver */
	char  *id
)
{
	static int	atlBest[] = {GLX_RGBA, GLX_RED_SIZE,4,
				GLX_GREEN_SIZE,4, GLX_BLUE_SIZE,4,
				GLX_DOUBLEBUFFER, GLX_DEPTH_SIZE,15, None};
	XSetWindowAttributes	ourwinattr;
	XWMHints	ourxwmhints;
					/* open display server */
	ourdisplay = XOpenDisplay(NULL);
	if (ourdisplay == NULL)
		error(USER, "cannot open X-windows; DISPLAY variable set?\n");
					/* find a usable visual */
	ourvinf = glXChooseVisual(ourdisplay, ourscreen, atlBest);
	if (ourvinf == NULL)
		error(USER, "no suitable visuals available");
					/* get a context */
	gctx = glXCreateContext(ourdisplay, ourvinf, NULL, GL_TRUE);
					/* open window */
	ourwinattr.background_pixel = ourblack;
	ourwinattr.border_pixel = ourblack;
	ourwinattr.event_mask = ourmask;
					/* this is stupid */
	ourwinattr.colormap = XCreateColormap(ourdisplay, ourroot,
				ourvinf->visual, AllocNone);
	gwind = XCreateWindow(ourdisplay, ourroot, 0, 0,
		DisplayWidth(ourdisplay,ourscreen)-2*BORWIDTH,
		DisplayHeight(ourdisplay,ourscreen)-2*BORWIDTH,
		BORWIDTH, ourvinf->depth, InputOutput, ourvinf->visual,
		CWBackPixel|CWBorderPixel|CWColormap|CWEventMask, &ourwinattr);
	if (gwind == 0)
		error(SYSTEM, "cannot create window\n");
   	XStoreName(ourdisplay, gwind, id);
#ifndef NOSTEREO
	if (stereo)			/* check if stereo working */
		switch (XSGIQueryStereoMode(ourdisplay, gwind)) {
		case STEREO_TOP:
		case STEREO_BOTTOM:
			break;
		case STEREO_OFF:
			error(USER,
		"wrong video mode: run \"/usr/gfx/setmon -n STR_TOP\" first");
		case X_STEREO_UNSUPPORTED:
			error(USER, "stereo not supported on this screen");
		default:
			error(INTERNAL, "unknown stereo mode");
		}
#endif
					/* set window manager hints */
	ourxwmhints.flags = InputHint|IconPixmapHint;
	ourxwmhints.input = True;
	ourxwmhints.icon_pixmap = XCreateBitmapFromData(ourdisplay,
		gwind, glradicon_bits, glradicon_width, glradicon_height);
	XSetWMHints(ourdisplay, gwind, &ourxwmhints);
					/* set GLX context */
	glXMakeCurrent(ourdisplay, gwind, gctx);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_LIGHTING);
	glFrontFace(GL_CCW);
	glCullFace(GL_BACK);
	if (backvis)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	glDrawBuffer(GL_BACK);
					/* figure out sensible view */
	pwidth = (double)DisplayWidthMM(ourdisplay, ourscreen) /
			DisplayWidth(ourdisplay, ourscreen);
	pheight = (double)DisplayHeightMM(ourdisplay, ourscreen) /
			DisplayHeight(ourdisplay, ourscreen);
	if (stereo) {			/* set stereo mode */
		setstereobuf(STEREO_BUFFER_LEFT);
		pheight *= 2.;
	}
					/* map the window */
	XMapWindow(ourdisplay, gwind);
	no_render++;
	do
		dev_input(0);		/* get resize event */
	while ((hres == 0) & (vres == 0));
	no_render--;
	rgl_checkerr("initializing GLX");
}


static void
dev_close(void)			/* close our display and free resources */
{
	glXMakeCurrent(ourdisplay, None, NULL);
	glXDestroyContext(ourdisplay, gctx);
	XDestroyWindow(ourdisplay, gwind);
	gwind = 0;
	XCloseDisplay(ourdisplay);
	ourdisplay = NULL;
}


static int
dev_view(			/* assign new driver view */
	register VIEW	*nv
)
{
	int	newhres = hres, newvres = vres;
	double	wa, va;
					/* check view legality */
	if (nv->type != VT_PER) {
		error(COMMAND, "illegal view type");
		nv->type = VT_PER;
	}
	if ((nv->horiz > 160.) | (nv->vert > 160.)) {
		error(COMMAND, "illegal view angle");
		if (nv->horiz > 160.)
			nv->horiz = 160.;
		if (nv->vert > 160.)
			nv->vert = 160.;
	}
	if ((hres != 0) & (vres != 0)) {
		wa = (vres*pheight)/(hres*pwidth);
		va = viewaspect(nv);
		if (va > wa+.05) {
			newvres = (pwidth/pheight)*va*newhres + .5;
			if (newvres > maxvres) {
				newvres = maxvres;
				newhres = (pheight/pwidth)/va*newvres + .5;
			}
		} else if (va < wa-.05) {
			newhres = (pheight/pwidth)/va*newvres + .5;
			if (newhres > maxhres) {
				newhres = maxhres;
				newvres = (pwidth/pheight)*va*newhres + .5;
			}
		}
		if ((newhres != hres) | (newvres != vres)) {
			no_render++;
			XResizeWindow(ourdisplay, gwind, newhres, newvres);
			do
				dev_input(0);		/* get resize event */
			while ((newhres != hres) | (newvres != vres));
			no_render--;
		}
	}
	thisview = *nv;
	setglpersp(&thisview);
	render();
	return(1);
}


static int
dev_input(		/* get next input event */
	int	nsecs
)
{
#if 0
	static time_t	lasttime = 0;
	time_t	thistime;

	if (nsecs > 0) {
		thistime = time(0);
		nsecs -= (long)(thistime - lasttime);
		lasttime = thistime;
	}
	if (nsecs > 0)
		alarm(nsecs);
#endif
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
		return(getkey(levptr(XKeyPressedEvent)));
	case ButtonPress:
		getmove(levptr(XButtonPressedEvent));
		break;
	}
	return(1);
}


static void
render(void)			/* render our display list and swap buffers */
{
	double	d;

	if (!mapped | no_render)
		return;
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glCallList(displist);
	if (stereo) {				/* do right eye for stereo */
		setstereobuf(STEREO_BUFFER_RIGHT);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		d = -eyedist / sqrt(thisview.hn2);
		glTranslated(d*thisview.hvec[0], d*thisview.hvec[1],
				d*thisview.hvec[2]);
		glCallList(displist);
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
		setstereobuf(STEREO_BUFFER_LEFT);
	}
	glXSwapBuffers(ourdisplay, gwind);	/* calls glFlush() */
	rgl_checkerr("rendering display list");
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
	FVECT	odir, v1, wp;
	double	d;
				/* start with old view */
	nv = thisview;
				/* change view direction */
	if ((d = viewray(v1, odir, &thisview,
			(dx+.5)/hres, (dy+.5)/vres)) < -FTINY)
		return(0);		/* outside view */
	if (mov | orb) {
		if (!getintersect(wp, v1, odir, d))
			return(0);
		VSUM(odir, wp, nv.vp, -1.);
	} else
		VCOPY(nv.vdir, odir);
	if (orb && mov) {		/* orbit left/right */
		spinvector(odir, odir, nv.vup, d=MOVDEG*PI/180.*mov);
		VSUM(nv.vp, wp, odir, -1.);
		spinvector(nv.vdir, nv.vdir, nv.vup, d);
	} else if (orb) {		/* orbit up/down */
		fcross(v1, odir, nv.vup);
		if (normalize(v1) == 0.)
			return(0);
		spinvector(odir, odir, v1, d=MOVDEG*PI/180.*orb);
		VSUM(nv.vp, wp, odir, -1.);
		spinvector(nv.vdir, nv.vdir, v1, d);
	} else if (mov) {		/* move forward/backward */
		d = MOVPCT/100. * mov;
		VSUM(nv.vp, nv.vp, odir, d);
	}
	if (!mov ^ !orb && headlocked) {	/* restore head height */
		VSUM(v1, thisview.vp, nv.vp, -1.);
		d = DOT(v1, thisview.vup);
		VSUM(nv.vp, nv.vp, thisview.vup, d);
	}
	if (setview(&nv) != NULL)
		return(0);	/* illegal view */
	dev_view(&nv);
	return(1);
}


static void
waitabit(void)				/* pause a moment */
{
	struct timespec	ts;
	ts.tv_sec = 0;
	ts.tv_nsec = 5000000;
	nanosleep(&ts, NULL);
}


static void
getmove(				/* get view change */
	XButtonPressedEvent	*ebut
)
{
	int	movdir = MOVDIR(ebut->button);
	int	movorb = MOVORB(ebut->state);
	int	moved = 0;
	Window	rootw, childw;
	int	rootx, rooty, wx, wy;
	unsigned int	statemask;

	copylastv( movorb ? (movdir ? "left/right" : "up/down") :
			(movdir ? "fore/back" : "rotate") );
	XNoOp(ourdisplay);

	while (!XCheckMaskEvent(ourdisplay,
			ButtonReleaseMask, levptr(XEvent))) {
					/* pause so as not to move too fast */
		waitabit();

		if (!XQueryPointer(ourdisplay, gwind, &rootw, &childw,
				&rootx, &rooty, &wx, &wy, &statemask))
			break;		/* on another screen */

		if (!moveview(wx, vres-1-wy, movdir, movorb)) {
			sleep(1);
			continue;
		} else
			moved++;
	}
	if (!moved) {					/* do final motion */
		movdir = MOVDIR(levptr(XButtonReleasedEvent)->button);
		wx = levptr(XButtonReleasedEvent)->x;
		wy = levptr(XButtonReleasedEvent)->y;
		moveview(wx, vres-1-wy, movdir, movorb);
	}
}


static int
getintersect(		/* intersect ray with scene geometry */
	FVECT	wp,		/* returned world intersection point */
	FVECT	org,
	FVECT	dir,
	double	md
)
{
	float	fbuf[6];
				/* check to see if rtrace is running */
	/* if (rtpd.pid <= 0) */
	if (!rtpd.running)
		return(0);
				/* assign origin */
	fbuf[0] = org[0]; fbuf[1] = org[1]; fbuf[2] = org[2];
				/* compute clipping distance */
	if (md <= FTINY) md = FHUGE;
	fbuf[3] = dir[0]*md; fbuf[4] = dir[1]*md; fbuf[5] = dir[2]*md;
				/* trace that ray */
	if (process(&rtpd, (char *)fbuf, (char *)fbuf,
			4*sizeof(float), 6*sizeof(float)) != 4*sizeof(float))
		error(INTERNAL, "error getting data back from rtrace process");
	if (fbuf[3] >= .99*FHUGE)
		return(0);	/* missed local objects */
	wp[0] = fbuf[0]; wp[1] = fbuf[1]; wp[2] = fbuf[2];
	return(1);		/* else return world intersection */
}


static void
setglpersp(			/* set perspective view in GL */
	register VIEW	*vp
)
{
	double	d, xmin, xmax, ymin, ymax, zmin, zmax;

	zmin = 0.1;
	zmax = 1000.;
	if (thisview.vfore > FTINY)
		zmin = thisview.vfore;
	if (thisview.vaft > FTINY)
		zmax = thisview.vaft;
	xmax = zmin * tan(PI/180./2. * thisview.horiz);
	xmin = -xmax;
	d = thisview.hoff * (xmax - xmin);
	xmin += d; xmax += d;
	ymax = zmin * tan(PI/180./2. * thisview.vert);
	ymin = -ymax;
	d = thisview.voff * (ymax - ymin);
	ymin += d; ymax += d;
					/* set view matrix */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(xmin, xmax, ymin, ymax, zmin, zmax);
	gluLookAt(thisview.vp[0], thisview.vp[1], thisview.vp[2],
		thisview.vp[0] + thisview.vdir[0],
		thisview.vp[1] + thisview.vdir[1],
		thisview.vp[2] + thisview.vdir[2],
		thisview.vup[0], thisview.vup[1], thisview.vup[2]);
	rgl_checkerr("setting perspective view");
}


static int
getkey(				/* get input key */
	register XKeyPressedEvent  *ekey
)
{
	int  n;
	char	buf[8];

	n = XLookupString(ekey, buf, sizeof(buf), NULL, NULL);
	if (n != 1)
		return(1);
	switch (buf[0]) {
	case 'h':			/* turn on height motion lock */
		headlocked = 1;
		break;
	case 'H':			/* turn off height motion lock */
		headlocked = 0;
		break;
	case 'l':			/* retrieve last (premouse) view */
		if (lastview.type) {
			VIEW	vtmp;
			vtmp = thisview;
			dev_view(&lastview);
			lastview = vtmp;
		} else
			XBell(ourdisplay, 0);
		break;
	case 'n':			/* move to next standard view */
		gotoview(currentview+1);
		break;
	case 'p':			/* move to last standard view */
		gotoview(currentview-1);
		break;
	case '+':			/* zoom in */
		zoomview(100+ZOOMPCT, ekey->x, vres-1-ekey->y);
		break;
	case '-':			/* zoom out */
		zoomview(100-ZOOMPCT, ekey->x, vres-1-ekey->y);
		break;
	case 'v':			/* spit current view to stdout */
		fputs(VIEWSTR, stdout);
		fprintview(&thisview, stdout);
		fputc('\n', stdout);
		break;
	case 'V':			/* append view to rad file */
		appendview(NULL, &thisview);
		break;
	case 'q':			/* quit the program */
		return(0);
	default:
		XBell(ourdisplay, 0);
		break;
	}
	return(1);
}


static void
zoomview(			/* zoom in or out around (dx,dy) */
	int	pct,
	int	dx,
	int	dy
)
{
	double	h, v;

	if ((pct == 100) | (pct <= 0))
		return;
	copylastv("zooming");
	h = (dx+.5)/hres - 0.5;
	v = (dy+.5)/vres - 0.5;
	h *= (1. - 100./pct);
	v *= (1. - 100./pct);
	thisview.vdir[0] += h*thisview.hvec[0] + v*thisview.vvec[0];
	thisview.vdir[1] += h*thisview.hvec[1] + v*thisview.vvec[1];
	thisview.vdir[2] += h*thisview.hvec[2] + v*thisview.vvec[2];
	thisview.horiz = 2.*180./PI * atan( 100./pct *
					tan(PI/180./2.*thisview.horiz) );
	thisview.vert = 2.*180./PI * atan( 100./pct *
					tan(PI/180./2.*thisview.vert) );
	setview(&thisview);
	dev_view(&thisview);
}


static void
gotoview(				/* go to specified view number */
	int	vwnum
)
{
	if (vwnum < 0)
		for (vwnum = currentview; vwl[vwnum+1].v != NULL; vwnum++)
			;
	else if (vwnum >= MAXVIEW || vwl[vwnum].v == NULL)
		vwnum = 0;
	copylastv("standard view");
	dev_view(vwl[currentview=vwnum].v);
}


static void
appendview(			/* append standard view */
	char	*nm,
	VIEW	*vp
)
{
	FILE	*fp;
					/* check if already in there */
	if (!memcmp(&thisview, vwl[currentview].v, sizeof(VIEW))) {
		error(COMMAND, "view already in standard list");
		return;
	}
					/* append to file */
	if ((fp = fopen(radfile, "a")) == NULL) {
		error(COMMAND, "cannot append rad input file");
		return;
	}
	fputs("view=", fp);
	if (nm != NULL) {
		fputc(' ', fp); fputs(nm, fp);
	}
	fprintview(vp, fp); fputc('\n', fp);
	fclose(fp);
					/* append to our list */
	while (vwl[currentview].v != NULL)
		currentview++;
	if (currentview >= MAXVIEW)
		error(INTERNAL, "too many views in appendview");
	vwl[currentview].v = (VIEW *)bmalloc(sizeof(VIEW));
	*(vwl[currentview].v) = thisview;
	if (nm != NULL)
		vwl[currentview].nam = savqstr(nm);
}


static void
copylastv(			/* copy last view position */
	char	*cause
)
{
	static char	*lastvc;

	if (cause == lastvc)
		return;			/* only record one view per cause */
	lastvc = cause;
	lastview = thisview;
}


static void
fixwindow(				/* repair damage to window */
	register XExposeEvent  *eexp
)
{
	if ((hres == 0) | (vres == 0)) {	/* first exposure */
		resizewindow((XConfigureEvent *)eexp);
		return;
	}
	if (eexp->count)		/* wait for final exposure */
		return;
					/* rerender everything */
	render();
}


static void
resizewindow(			/* resize window */
	register XConfigureEvent  *ersz
)
{
	static char	resizing[] = "resizing window";
	double	wa, va;

	glViewport(0, 0, hres=ersz->width, vres=ersz->height);
	if (hres > maxhres) maxhres = hres;
	if (vres > maxvres) maxvres = vres;
	if (no_render)
		return;
	wa = (vres*pheight)/(hres*pwidth);
	va = viewaspect(&thisview);
	if (va > wa+.05) {
		copylastv(resizing);
		thisview.vert = 2.*180./PI *
				atan( tan(PI/180./2. * thisview.horiz) * wa );
	} else if (va < wa-.05) {
		copylastv(resizing);
		thisview.horiz = 2.*180./PI *
				atan( tan(PI/180./2. * thisview.vert) / wa );
	} else
		return;
	setview(&thisview);
	dev_view(&thisview);
}
