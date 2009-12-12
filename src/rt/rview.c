#ifndef lint
static const char	RCSid[] = "$Id: rview.c,v 2.33 2009/12/12 23:08:13 greg Exp $";
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

#define	 CTRL(c)	((c)-'@')


void
quit(code)			/* quit program */
int  code;
{
#ifdef MSTATS
	if (code == 2 && errno == ENOMEM)
		printmemstats(stderr);
#endif
	if (ray_pnprocs > 0)	/* close children if any */
		ray_pclose(0);
	else if (!ray_pnprocs)	/* in parent */
		devclose();
	exit(code);
}


void
devopen(				/* open device driver */
	char  *dname
)
{
	extern char  *progname, *octname;
	char  *id;
	int  i;

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
devclose(void)				/* close our device */
{
	if (dev != NULL)
		(*dev->close)();
	dev = NULL;
}


void
printdevices(void)			/* print list of output devices */
{
	int  i;

	for (i = 0; devtable[i].name; i++)
		printf("%-16s # %s\n", devtable[i].name, devtable[i].descrip);
}


void
rview(void)				/* do a view */
{
	char  buf[32];

	devopen(dvcname);		/* open device */
	newimage(NULL);			/* start image */

	for ( ; ; ) {			/* quit in command() */
		while (hresolu <= 1<<pdepth && vresolu <= 1<<pdepth)
			command("done: ");
		errno = 0;
		if (hresolu <= psample<<pdepth && vresolu <= psample<<pdepth) {
			sprintf(buf, "%d sampling...\n", 1<<pdepth);
			(*dev->comout)(buf);
			rsample();
		} else {
			sprintf(buf, "%d refining...\n", 1<<pdepth);
			(*dev->comout)(buf);
			refine(&ptrunk, pdepth+1);
		}
		if (dev->inpready)		/* noticed some input */
			command(": ");
		else				/* finished this depth */
			pdepth++;
	}
}


void
command(			/* get/execute command */
	char  *prompt
)
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

	if (waitrays() < 0)			/* clear ray queue */
		quit(1);
	
	switch (inpbuf[0]) {
	case 'f':				/* new frame (|focus|free) */
		if (badcom("frame")) {
			if (badcom("focus")) {
				if (badcom("free"))
					goto commerr;
				free_objmem();
				break;
			}
			getfocus(args);
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
		newimage(args);
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
rsample(void)			/* sample the image */
{
	int  xsiz, ysiz, y;
	PNODE  *p;
	PNODE	**pl;
	int  x;
	/*
	 *     We initialize the bottom row in the image at our current
	 * resolution.	During sampling, we check super-pixels to the
	 * right and above by calling bigdiff().  If there is a significant
	 * difference, we subsample the super-pixels.  The testing process
	 * includes initialization of the next row.
	 */
	xsiz = (((long)(pframe.r-pframe.l)<<pdepth)+hresolu-1) / hresolu;
	ysiz = (((long)(pframe.u-pframe.d)<<pdepth)+vresolu-1) / vresolu;
	pl = (PNODE **)malloc(xsiz*sizeof(PNODE *));
	if (pl == NULL)
		return;
	/*
	 * Initialize the bottom row.
	 */
	pl[0] = findrect(pframe.l, pframe.d, &ptrunk, pdepth);
	for (x = 1; x < xsiz; x++) {
		pl[x] = findrect(pframe.l+((x*hresolu)>>pdepth),
				pframe.d, &ptrunk, pdepth);
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
				refine(pl[x], 1);
				refine(pl[x+1], 1);
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
			p = findrect(pframe.l+((x*hresolu)>>pdepth),
				pframe.d+(((y+1)*vresolu)>>pdepth),
					&ptrunk, pdepth);
			/*
			 * Test super-pixel in next row.
			 */
			if (pl[x] != p && bigdiff(pl[x]->v, p->v, maxdiff)) {
				refine(pl[x], 1);
				refine(p, 1);
			}
			/*
			 * Copy into super-pixel array.
			 */
			pl[x] = p;
		}
	}
escape:
	free((void *)pl);
}


int
refine(				/* refine a node */
	PNODE	*p,
	int  pd
)
{
	int  growth;
	int  mx, my;
	int  i;

	if (dev->inpready)			/* quit for input */
		return(0);

	if (pd <= 0)				/* depth limit */
		return(0);

	mx = (p->xmin + p->xmax) >> 1;
	my = (p->ymin + p->ymax) >> 1;
	growth = 0;

	if (p->kid == NULL) {			/* subdivide */

		if ((p->kid = newptree()) == NULL)
			return(0);

		p->kid[UR].xmin = mx;
		p->kid[UR].ymin = my;
		p->kid[UR].xmax = p->xmax;
		p->kid[UR].ymax = p->ymax;
		p->kid[UL].xmin = p->xmin;
		p->kid[UL].ymin = my;
		p->kid[UL].xmax = mx;
		p->kid[UL].ymax = p->ymax;
		p->kid[DR].xmin = mx;
		p->kid[DR].ymin = p->ymin;
		p->kid[DR].xmax = p->xmax;
		p->kid[DR].ymax = my;
		p->kid[DL].xmin = p->xmin;
		p->kid[DL].ymin = p->ymin;
		p->kid[DL].xmax = mx;
		p->kid[DL].ymax = my;
		/*
		 *  The following paint order can leave a black pixel
		 *  if redraw() is called in (*dev->paintr)().
		 */
		if (p->x >= mx && p->y >= my)
			pcopy(p, p->kid+UR);
		else if (paint(p->kid+UR) < 0)
			quit(1);
		if (p->x < mx && p->y >= my)
			pcopy(p, p->kid+UL);
		else if (paint(p->kid+UL) < 0)
			quit(1);
		if (p->x >= mx && p->y < my)
			pcopy(p, p->kid+DR);
		else if (paint(p->kid+DR) < 0)
			quit(1);
		if (p->x < mx && p->y < my)
			pcopy(p, p->kid+DL);
		else if (paint(p->kid+DL) < 0)
			quit(1);

		growth++;
	}
						/* do children */
	if (mx > pframe.l) {
		if (my > pframe.d)
			growth += refine(p->kid+DL, pd-1);
		if (my < pframe.u)
			growth += refine(p->kid+UL, pd-1);
	}
	if (mx < pframe.r) {
		if (my > pframe.d)
			growth += refine(p->kid+DR, pd-1);
		if (my < pframe.u)
			growth += refine(p->kid+UR, pd-1);
	}
	return(growth);
}
