#ifndef lint
static const char	RCSid[] = "$Id: raycalls.c,v 2.2 2003/02/25 02:47:23 greg Exp $";
#endif
/*
 *  raycalls.c - interface for running Radiance rendering as a library
 *
 *  External symbols declared in ray.h
 */

#include "copyright.h"

/*
 *  These routines are designed to aid the programmer who wishes
 *  to call Radiance as a library.  Unfortunately, the system was
 *  not originally intended to be run this way, and there are some
 *  awkward limitations to contend with.  The most irritating
 *  perhaps is that the global variables and functions do not have
 *  a prefix, and the symbols are a bit generic.  This results in a
 *  serious invasion of the calling application's name-space, and
 *  you may need to rename either some Radiance routines or some
 *  of your routines to avoid conflicts.  Another limitation is
 *  that the global variables are not gathered together into any
 *  sort of context, so it is impossible to simultaneously run
 *  this library on multiple scenes or in multiple threads.
 *  You get one scene and one thread, and if you want more, you
 *  will have to go with the process model used by the programs
 *  gen/mkillum, hd/rholo, and px/pinterp.  Finally, unrecoverable
 *  errors result in a call to the application-defined function
 *  quit().  The usual thing to do is to call exit().
 *  You might want to do something else instead, like
 *  call setjmp()/longjmp() to bring you back to the calling
 *  function for recovery.  You may also wish to define your own
 *  wputs(s) and eputs(s) functions to output warning and error
 *  messages, respectively.
 *
 *  With those caveats, we have attempted to make the interface
 *  as simple as we can.  Global variables and their defaults
 *  are defined below, and including "ray.h" declares these
 *  along with all the routines you are likely to need.  First,
 *  assign the global variable progname to your argv[0], then
 *  change the rendering parameters as you like.  If you have a set
 *  of option arguments you are working from, the getrenderopt(ac,av)
 *  call should be very useful.  Before tracing any rays, you
 *  must read in the octree with a call to ray_init(oct).
 *  Passing NULL for the file name causes ray_init() to read
 *  the octree from the standard input -- rarely a good idea.
 *  However, one may read an octree from a program (such as
 *  oconv) by preceding a shell command by a '!' character.
 *
 *  To trace a ray, define a RAY object myRay and assign:
 *
 *	myRay.rorg = ( ray origin point )
 *	myRay.rdir = ( normalized ray direction )
 *	myRay.rmax = ( maximum length, or zero for no limit )
 *
 *  If you are rendering from a VIEW structure, this can be
 *  accomplished with a single call for the ray at (x,y):
 *
 *	myRay.rmax = viewray(myRay.rorg, myRay.rdir, &myView, x, y);
 *
 *  Then, trace the primary ray with:
 *
 *	ray_trace(&myRay);
 *
 *  The resulting contents of myRay should provide you with
 *  more than enough information about what the ray hit,
 *  the computed value, etc.  For further clues of how to
 *  compute irradiance, how to get callbacks on the evaluated
 *  ray tree, etc., see the contents of rtrace.c.  See
 *  also the rpmain.c, rtmain.c, and rvmain.c modules
 *  to learn more how rendering options are processed.
 *
 *  When you are done, you may call ray_done(1) to clean
 *  up memory used by Radiance.  It doesn't free everything,
 *  but it makes a valiant effort.  If you call ray_done(0),
 *  it leaves data that is likely to be reused, including
 *  loaded data files and fonts.  The library may be
 *  restarted at any point by calling ray_init() on a new
 *  octree.
 *
 *  The call ray_save(rp) allocates and returns a buffer
 *  with the current global parameter settings, which may be
 *  restored at any time with a call to ray_restore(rp).
 *  This buffer contains no linked information, and thus
 *  may be passed between processes using write() and
 *  read() calls, so long as byte order is maintained.
 *  Calling ray_restore(NULL) restores the original
 *  default parameters, which is also retrievable with
 *  the call ray_defaults(rp).  (These  should be the
 *  same as the defaults for rtrace.)
 */

#include  "ray.h"

#include  "source.h"

#include  "ambient.h"

#include  "otypes.h"

#include  "random.h"

#include  "data.h"

#include  "font.h"

char	*progname = "unknown_app";	/* caller sets to argv[0] */

char	*octname;			/* octree name we are given */

char	*shm_boundary = NULL;		/* boundary of shared memory */

CUBE	thescene;			/* our scene */
OBJECT	nsceneobjs;			/* number of objects in our scene */

int	dimlist[MAXDIM];		/* sampling dimensions */
int	ndims = 0;			/* number of sampling dimensions */
int	samplendx = 0;			/* index for this sample */

void	(*trace)() = NULL;		/* trace call */

extern void	ambnotify();
void	(*addobjnotify[])() = {ambnotify, NULL};

int	do_irrad = 0;			/* compute irradiance? */

double	dstrsrc = 0.0;			/* square source distribution */
double	shadthresh = .05;		/* shadow threshold */
double	shadcert = .5;			/* shadow certainty */
int	directrelay = 2;		/* number of source relays */
int	vspretest = 512;		/* virtual source pretest density */
int	directvis = 1;			/* sources visible? */
double	srcsizerat = .2;		/* maximum ratio source size/dist. */

COLOR	cextinction = BLKCOLOR;		/* global extinction coefficient */
COLOR	salbedo = BLKCOLOR;		/* global scattering albedo */
double	seccg = 0.;			/* global scattering eccentricity */
double	ssampdist = 0.;			/* scatter sampling distance */

double	specthresh = .15;		/* specular sampling threshold */
double	specjitter = 1.;		/* specular sampling jitter */

int	backvis = 1;			/* back face visibility */

int	maxdepth = 6;			/* maximum recursion depth */
double	minweight = 4e-3;		/* minimum ray weight */

char	*ambfile = NULL;		/* ambient file name */
COLOR	ambval = BLKCOLOR;		/* ambient value */
int	ambvwt = 0;			/* initial weight for ambient value */
double	ambacc = 0.2;			/* ambient accuracy */
int	ambres = 128;			/* ambient resolution */
int	ambdiv = 512;			/* ambient divisions */
int	ambssamp = 0;			/* ambient super-samples */
int	ambounce = 0;			/* ambient bounces */
char	*amblist[AMBLLEN+1];		/* ambient include/exclude list */
int	ambincl = -1;			/* include == 1, exclude == 0 */


void
ray_init(otnm)			/* initialize ray-tracing calculation */
char	*otnm;
{
	if (nobjects > 0)		/* free old scene data */
		ray_done(0);
					/* initialize object types */
	if (ofun[OBJ_SPHERE].funp == o_default)
		initotypes();
					/* initialize urand */
	if (urperm == NULL)
		initurand(2048);
					/* read scene octree */
	readoct(octname = otnm, ~(IO_FILES|IO_INFO), &thescene, NULL);
	nsceneobjs = nobjects;
					/* find and mark sources */
	marksources();
					/* initialize ambient calculation */
	setambient();
					/* ready to go... */
}


void
ray_trace(RAY *r)		/* trace a primary ray */
{
	rayorigin(r, NULL, PRIMARY, 1.0);
	samplendx++;
	rayvalue(r);		/* assumes origin and direction are set */
}


void
ray_done(freall)		/* free ray-tracing data */
int	freall;
{
	retainfonts = 1;
	ambdone();
	ambnotify(OVOID);
	freesources();
	freeobjects(0, nobjects);
	donesets();
	octdone();
	thescene.cutree = EMPTY;
	octname = NULL;
	if (freall) {
		retainfonts = 0;
		freefont(NULL);
		freedata(NULL);
		initurand(0);
	}
	if (nobjects > 0) {
		sprintf(errmsg, "%d objects left after call to ray_done()",
				nobjects);
		error(WARNING, errmsg);
	}
}


void
ray_save(rp)			/* save current parameter settings */
RAYPARAMS	*rp;
{
	int	i, ndx;

	if (rp == NULL)
		return;
	rp->do_irrad = do_irrad;
	rp->dstrsrc = dstrsrc;
	rp->shadthresh = shadthresh;
	rp->shadcert = shadcert;
	rp->directrelay = directrelay;
	rp->vspretest = vspretest;
	rp->directvis = directvis;
	rp->srcsizerat = srcsizerat;
	copycolor(rp->cextinction, cextinction);
	copycolor(rp->salbedo, salbedo);
	rp->seccg = seccg;
	rp->ssampdist = ssampdist;
	rp->specthresh = specthresh;
	rp->specjitter = specjitter;
	rp->backvis = backvis;
	rp->maxdepth = maxdepth;
	rp->minweight = minweight;
	copycolor(rp->ambval, ambval);
	bzero(rp->ambfile, sizeof(rp->ambfile));
	if (ambfile != NULL)
		strncpy(rp->ambfile, ambfile, sizeof(rp->ambfile)-1);
	rp->ambvwt = ambvwt;
	rp->ambacc = ambacc;
	rp->ambres = ambres;
	rp->ambdiv = ambdiv;
	rp->ambssamp = ambssamp;
	rp->ambounce = ambounce;
	rp->ambincl = ambincl;
	bzero(rp->amblval, sizeof(rp->amblval));
	ndx = 0;
	for (i = 0; i < AMBLLEN && amblist[i] != NULL; i++) {
		int	len = strlen(amblist[i]);
		if (ndx+len >= sizeof(rp->amblval))
			break;
		strcpy(rp->amblval+ndx, amblist[i]);
		ndx += len+1;
	}
	while (i <= AMBLLEN)
		rp->amblndx[i++] = -1;
}


void
ray_restore(rp)			/* restore parameter settings */
RAYPARAMS	*rp;
{
	register int	i;

	if (rp == NULL) {		/* restore defaults */
		RAYPARAMS	dflt;
		ray_defaults(&dflt);
		ray_restore(&dflt);
		return;
	}
					/* restore saved settings */
	do_irrad = rp->do_irrad;
	dstrsrc = rp->dstrsrc;
	shadthresh = rp->shadthresh;
	shadcert = rp->shadcert;
	directrelay = rp->directrelay;
	vspretest = rp->vspretest;
	directvis = rp->directvis;
	srcsizerat = rp->srcsizerat;
	copycolor(cextinction, rp->cextinction);
	copycolor(salbedo, rp->salbedo);
	seccg = rp->seccg;
	ssampdist = rp->ssampdist;
	specthresh = rp->specthresh;
	specjitter = rp->specjitter;
	backvis = rp->backvis;
	maxdepth = rp->maxdepth;
	minweight = rp->minweight;
	copycolor(ambval, rp->ambval);
	ambvwt = rp->ambvwt;
	ambdiv = rp->ambdiv;
	ambssamp = rp->ambssamp;
	ambounce = rp->ambounce;
	for (i = 0; rp->amblndx[i] >= 0; i++)
		amblist[i] = rp->amblval + rp->amblndx[i];
	while (i <= AMBLLEN)
		amblist[i++] = NULL;
	ambincl = rp->ambincl;
					/* update ambient calculation */
	ambnotify(OVOID);
	if (thescene.cutree != EMPTY) {
		int	newamb = (ambfile == NULL) ?  rp->ambfile[0] :
					strcmp(ambfile, rp->ambfile) ;

		if (amblist[0] != NULL)
			for (i = 0; i < nobjects; i++)
				ambnotify(i);

		ambfile = (rp->ambfile[0]) ? rp->ambfile : (char *)NULL;
		if (newamb) {
			ambres = rp->ambres;
			ambacc = rp->ambacc;
			setambient();
		} else {
			setambres(rp->ambres);
			setambacc(rp->ambacc);
		}
	} else {
		ambfile = (rp->ambfile[0]) ? rp->ambfile : (char *)NULL;
		ambres = rp->ambres;
		ambacc = rp->ambacc;
	}
}


void
ray_defaults(rp)		/* get default parameter values */
RAYPARAMS	*rp;
{
	int	i;

	if (rp == NULL)
		return;

	rp->do_irrad = 0;
	rp->dstrsrc = 0.0;
	rp->shadthresh = .05;
	rp->shadcert = .5;
	rp->directrelay = 2;
	rp->vspretest = 512;
	rp->directvis = 1;
	rp->srcsizerat = .2;
	setcolor(rp->cextinction, 0., 0., 0.);
	setcolor(rp->salbedo, 0., 0., 0.);
	rp->seccg = 0.;
	rp->ssampdist = 0.;
	rp->specthresh = .15;
	rp->specjitter = 1.;
	rp->backvis = 1;
	rp->maxdepth = 6;
	rp->minweight = 4e-3;
	setcolor(rp->ambval, 0., 0., 0.);
	bzero(rp->ambfile, sizeof(rp->ambfile));
	rp->ambvwt = 0;
	rp->ambres = 128;
	rp->ambacc = 0.2;
	rp->ambdiv = 512;
	rp->ambssamp = 0;
	rp->ambounce = 0;
	rp->ambincl = -1;
	bzero(rp->amblval, sizeof(rp->amblval));
	for (i = AMBLLEN+1; i--; )
		rp->amblndx[i] = -1;
}
