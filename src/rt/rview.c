#ifndef lint
static const char	RCSid[] = "$Id: rview.c,v 2.24 2003/09/24 14:55:54 greg Exp $";
#endif
/*
 *  rview.c - routines and variables for interactive view generation.
 *
 *  External symbols declared in rpaint.h
 */

#include "copyright.h"

#include  <signal.h>
#include  <ctype.h>

#include  "ray.h"
#include  "rpaint.h"


CUBE  thescene;				/* our scene */
OBJECT	nsceneobjs;			/* number of objects in our scene */

int  dimlist[MAXDIM];			/* sampling dimensions */
int  ndims = 0;				/* number of sampling dimensions */
int  samplendx = 0;			/* index for this sample */

extern void  ambnotify();
void  (*addobjnotify[])() = {ambnotify, NULL};

VIEW  ourview = STDVIEW;		/* viewing parameters */
int  hresolu, vresolu;			/* image resolution */

void  (*trace)() = NULL;		/* trace call */

int  do_irrad = 0;			/* compute irradiance? */

int  psample = 8;			/* pixel sample size */
double	maxdiff = .15;			/* max. sample difference */

double	exposure = 1.0;			/* exposure for scene */

double	dstrsrc = 0.0;			/* square source distribution */
double	shadthresh = .1;		/* shadow threshold */
double	shadcert = .25;			/* shadow certainty */
int  directrelay = 0;			/* number of source relays */
int  vspretest = 128;			/* virtual source pretest density */
int  directvis = 1;			/* sources visible? */
double	srcsizerat = 0.;		/* maximum ratio source size/dist. */

COLOR  cextinction = BLKCOLOR;		/* global extinction coefficient */
COLOR  salbedo = BLKCOLOR;		/* global scattering albedo */
double  seccg = 0.;			/* global scattering eccentricity */
double  ssampdist = 0.;			/* scatter sampling distance */

double	specthresh = .3;		/* specular sampling threshold */
double	specjitter = 1.;		/* specular sampling jitter */

int  backvis = 1;			/* back face visibility */

int  maxdepth = 6;			/* maximum recursion depth */
double	minweight = 1e-2;		/* minimum ray weight */

char  *ambfile = NULL;			/* ambient file name */
COLOR  ambval = BLKCOLOR;		/* ambient value */
int  ambvwt = 0;			/* initial weight for ambient value */
double	ambacc = 0.2;			/* ambient accuracy */
int  ambres = 32;			/* ambient resolution */
int  ambdiv = 128;			/* ambient divisions */
int  ambssamp = 32;			/* ambient super-samples */
int  ambounce = 0;			/* ambient bounces */
char  *amblist[128];			/* ambient include/exclude list */
int  ambincl = -1;			/* include == 1, exclude == 0 */

int  greyscale = 0;			/* map colors to brightness? */
char  *dvcname = dev_default;		/* output device name */

struct driver  *dev = NULL;		/* driver functions */

char  rifname[128];			/* rad input file name */

VIEW  oldview;				/* previous view parameters */

PNODE  ptrunk;				/* the base of our image */
RECT  pframe;				/* current frame boundaries */
int  pdepth;				/* image depth in current frame */

static char  *reserve_mem = NULL;	/* pre-allocated reserve memory */

#define RESERVE_AMT	32768		/* amount of memory to reserve */

#define	 CTRL(c)	((c)-'@')


void
quit(code)			/* quit program */
int  code;
{
#ifdef MSTATS
	if (code == 2 && errno == ENOMEM)
		printmemstats(stderr);
#endif
	devclose();
	exit(code);
}


void
devopen(dname)				/* open device driver */
char  *dname;
{
	extern char  *progname, *octname;
	char  *id;
	register int  i;

	id = octname!=NULL ? octname : progname;
						/* check device table */
	for (i = 0; devtable[i].name; i++)
		if (!strcmp(dname, devtable[i].name)) {
			if ((dev = (*devtable[i].init)(dname, id)) == NULL) {
				sprintf(errmsg, "cannot initialize %s", dname);
				error(USER, errmsg);
			} else
				return;
		}
						/* not there, try exec */
	if ((dev = comm_init(dname, id)) == NULL) {
		sprintf(errmsg, "cannot start device \"%s\"", dname);
		error(USER, errmsg);
	}
}


void
devclose()				/* close our device */
{
	if (dev != NULL)
		(*dev->close)();
	dev = NULL;
}


void
printdevices()				/* print list of output devices */
{
	register int  i;

	for (i = 0; devtable[i].name; i++)
		printf("%-16s # %s\n", devtable[i].name, devtable[i].descrip);
}


void
rview()				/* do a view */
{
	char  buf[32];

	devopen(dvcname);		/* open device */
	newimage();			/* start image (calls fillreserves) */

	for ( ; ; ) {			/* quit in command() */
		while (hresolu <= 1<<pdepth && vresolu <= 1<<pdepth)
			command("done: ");
		while (reserve_mem == NULL)
			command("out of memory: ");
		errno = 0;
		if (hresolu <= psample<<pdepth && vresolu <= psample<<pdepth) {
			sprintf(buf, "%d sampling...\n", 1<<pdepth);
			(*dev->comout)(buf);
			rsample();
		} else {
			sprintf(buf, "%d refining...\n", 1<<pdepth);
			(*dev->comout)(buf);
			refine(&ptrunk, 0, 0, hresolu, vresolu, pdepth+1);
		}
		if (errno == ENOMEM)		/* ran out of memory */
			freereserves();
		else if (dev->inpready)		/* noticed some input */
			command(": ");
		else				/* finished this depth */
			pdepth++;
	}
}


void
fillreserves()			/* fill memory reserves */
{
	if (reserve_mem != NULL)
		return;
	reserve_mem = (char *)malloc(RESERVE_AMT);
}


void
freereserves()			/* free memory reserves */
{
	if (reserve_mem == NULL)
		return;
	free(reserve_mem);
	reserve_mem = NULL;
}


void
command(prompt)			/* get/execute command */
char  *prompt;
{
#define	 badcom(s)	strncmp(s, inpbuf, args-inpbuf-1)
	char  inpbuf[256];
	char  *args;
again:
	(*dev->comin)(inpbuf, prompt);		/* get command + arguments */
	for (args = inpbuf; *args && *args != ' '; args++)
		;
	if (*args) *args++ = '\0';
	else *++args = '\0';
	
	switch (inpbuf[0]) {
	case 'f':				/* new frame (or free mem.) */
		if (badcom("frame")) {
			if (badcom("free"))
				goto commerr;
			free_objmem();
			break;
		}
		getframe(args);
		break;
	case 'v':				/* view */
		if (badcom("view"))
			goto commerr;
		getview(args);
		break;
	case 'l':				/* last view */
		if (badcom("last"))
			goto commerr;
		lastview(args);
		break;
	case 'V':				/* save view */
		if (badcom("V"))
			goto commerr;
		saveview(args);
		break;
	case 'L':				/* load view */
		if (badcom("L"))
			goto commerr;
		loadview(args);
		break;
	case 'e':				/* exposure */
		if (badcom("exposure"))
			goto commerr;
		getexposure(args);
		break;
	case 's':				/* set a parameter */
		if (badcom("set")) {
#ifdef	SIGTSTP
			if (!badcom("stop"))
				goto dostop;
#endif
			goto commerr;
		}
		setparam(args);
		break;
	case 'n':				/* new picture */
		if (badcom("new"))
			goto commerr;
		newimage();
		break;
	case 't':				/* trace a ray */
		if (badcom("trace"))
			goto commerr;
		traceray(args);
		break;
	case 'a':				/* aim camera */
		if (badcom("aim"))
			goto commerr;
		getaim(args);
		break;
	case 'm':				/* move camera (or memstats) */
		if (badcom("move"))
#ifdef	MSTATS
		{
			if (badcom("memory"))
				goto commerr;
			printmemstats(stderr);
			break;
		}
#else
			goto commerr;
#endif
		getmove(args);
		break;
	case 'r':				/* rotate/repaint */
		if (badcom("rotate")) {
			if (badcom("repaint")) {
				if (badcom("redraw"))
					goto commerr;
				redraw();
				break;
			}
			getrepaint(args);
			break;
		}
		getrotate(args);
		break;
	case 'p':				/* pivot view */
		if (badcom("pivot")) {
			if (badcom("pause"))
				goto commerr;
			goto again;
		}
		getpivot(args);
		break;
	case CTRL('R'):				/* redraw */
		redraw();
		break;
	case 'w':				/* write */
		if (badcom("write"))
			goto commerr;
		writepict(args);
		break;
	case 'q':				/* quit */
		if (badcom("quit"))
			goto commerr;
		quit(0);
	case CTRL('C'):				/* interrupt */
		goto again;
#ifdef	SIGTSTP
	case CTRL('Z'):;			/* stop */
dostop:
		devclose();
		kill(0, SIGTSTP);
		/* pc stops here */
		devopen(dvcname);
		redraw();
		break;
#endif
	case '\0':				/* continue */
		break;
	default:;
commerr:
		if (iscntrl(inpbuf[0]))
			sprintf(errmsg, "^%c: unknown control",
					inpbuf[0]|0100);
		else
			sprintf(errmsg, "%s: unknown command", inpbuf);
		error(COMMAND, errmsg);
		break;
	}
#undef	badcom
}


void
rsample()			/* sample the image */
{
	int  xsiz, ysiz, y;
	RECT  r;
	PNODE  *p;
	register RECT  *rl;
	register PNODE	**pl;
	register int  x;
	/*
	 *     We initialize the bottom row in the image at our current
	 * resolution.	During sampling, we check super-pixels to the
	 * right and above by calling bigdiff().  If there is a significant
	 * difference, we subsample the super-pixels.  The testing process
	 * includes initialization of the next row.
	 */
	xsiz = (((long)(pframe.r-pframe.l)<<pdepth)+hresolu-1) / hresolu;
	ysiz = (((long)(pframe.u-pframe.d)<<pdepth)+vresolu-1) / vresolu;
	rl = (RECT *)malloc(xsiz*sizeof(RECT));
	if (rl == NULL)
		return;
	pl = (PNODE **)malloc(xsiz*sizeof(PNODE *));
	if (pl == NULL) {
		free((void *)rl);
		return;
	}
	/*
	 * Initialize the bottom row.
	 */
	rl[0].l = rl[0].d = 0;
	rl[0].r = hresolu; rl[0].u = vresolu;
	pl[0] = findrect(pframe.l, pframe.d, &ptrunk, rl, pdepth);
	for (x = 1; x < xsiz; x++) {
		rl[x].l = rl[x].d = 0;
		rl[x].r = hresolu; rl[x].u = vresolu;
		pl[x] = findrect(pframe.l+((x*hresolu)>>pdepth),
				pframe.d, &ptrunk, rl+x, pdepth);
	}
						/* sample the image */
	for (y = 0; /* y < ysiz */ ; y++) {
		for (x = 0; x < xsiz-1; x++) {
			if (dev->inpready || errno == ENOMEM)
				goto escape;
			/*
			 * Test super-pixel to the right.
			 */
			if (pl[x] != pl[x+1] && bigdiff(pl[x]->v,
					pl[x+1]->v, maxdiff)) {
				refine(pl[x], rl[x].l, rl[x].d,
						rl[x].r, rl[x].u, 1);
				refine(pl[x+1], rl[x+1].l, rl[x+1].d,
						rl[x+1].r, rl[x+1].u, 1);
			}
		}
		if (y >= ysiz-1)
			break;
		for (x = 0; x < xsiz; x++) {
			if (dev->inpready || errno == ENOMEM)
				goto escape;
			/*
			 * Find super-pixel at this position in next row.
			 */
			r.l = r.d = 0;
			r.r = hresolu; r.u = vresolu;
			p = findrect(pframe.l+((x*hresolu)>>pdepth),
				pframe.d+(((y+1)*vresolu)>>pdepth),
					&ptrunk, &r, pdepth);
			/*
			 * Test super-pixel in next row.
			 */
			if (pl[x] != p && bigdiff(pl[x]->v, p->v, maxdiff)) {
				refine(pl[x], rl[x].l, rl[x].d,
						rl[x].r, rl[x].u, 1);
				refine(p, r.l, r.d, r.r, r.u, 1);
			}
			/*
			 * Copy into super-pixel array.
			 */
			rl[x].l = r.l; rl[x].d = r.d;
			rl[x].r = r.r; rl[x].u = r.u;
			pl[x] = p;
		}
	}
escape:
	free((void *)rl);
	free((void *)pl);
}


int
refine(p, xmin, ymin, xmax, ymax, pd)		/* refine a node */
register PNODE	*p;
int  xmin, ymin, xmax, ymax;
int  pd;
{
	int  growth;
	int  mx, my;
	int  i;

	if (dev->inpready)			/* quit for input */
		return(0);

	if (pd <= 0)				/* depth limit */
		return(0);

	mx = (xmin + xmax) >> 1;
	my = (ymin + ymax) >> 1;
	growth = 0;

	if (p->kid == NULL) {			/* subdivide */

		if ((p->kid = newptree()) == NULL)
			return(0);
		/*
		 *  The following paint order can leave a black pixel
		 *  if redraw() is called in (*dev->paintr)().
		 */
		if (p->x >= mx && p->y >= my)
			pcopy(p, p->kid+UR);
		else
			paint(p->kid+UR, mx, my, xmax, ymax);
		if (p->x < mx && p->y >= my)
			pcopy(p, p->kid+UL);
		else
			paint(p->kid+UL, xmin, my, mx, ymax);
		if (p->x >= mx && p->y < my)
			pcopy(p, p->kid+DR);
		else
			paint(p->kid+DR, mx, ymin, xmax, my);
		if (p->x < mx && p->y < my)
			pcopy(p, p->kid+DL);
		else
			paint(p->kid+DL, xmin, ymin, mx, my);

		growth++;
	}
						/* do children */
	if (mx > pframe.l) {
		if (my > pframe.d)
			growth += refine(p->kid+DL, xmin, ymin, mx, my, pd-1);
		if (my < pframe.u)
			growth += refine(p->kid+UL, xmin, my, mx, ymax, pd-1);
	}
	if (mx < pframe.r) {
		if (my > pframe.d)
			growth += refine(p->kid+DR, mx, ymin, xmax, my, pd-1);
		if (my < pframe.u)
			growth += refine(p->kid+UR, mx, my, xmax, ymax, pd-1);
	}
						/* recompute sum */
	if (growth) {
		setcolor(p->v, 0.0, 0.0, 0.0);
		for (i = 0; i < 4; i++)
			addcolor(p->v, p->kid[i].v);
		scalecolor(p->v, 0.25);
	}
	return(growth);
}
