#ifndef lint
static const char	RCSid[] = "$Id: glareval.c,v 2.12 2006/05/29 16:47:54 greg Exp $";
#endif
/*
 * Compute pixels for glare calculation
 */

#include "copyright.h"

#include <stdlib.h>
#include <string.h>

#include "rtprocess.h" /* Windows: must come first because of conflicts */
#include "glare.h"

					/* maximum rtrace buffer size */
#define MAXPIX		(4096/(6*sizeof(float)))

#define MAXSBUF		786432	/* maximum total size of scanline buffer */
#define HSIZE		317	/* size of scanline hash table */
#define NRETIRE		16	/* number of scanlines to retire at once */

static SUBPROC	rt_pd = SP_INACTIVE; /* process id & descriptors for rtrace */

FILE	*pictfp = NULL;		/* picture file pointer */
double	exposure;		/* picture exposure */
int	pxsiz, pysiz;		/* picture dimensions */

static int	curpos;		/* current scanline */
static long	*scanpos;	/* scanline positions */

typedef struct scan {
	int	y;		/* scanline position */
	long	lused;		/* for LRU replacement */
	struct scan	*next;	/* next in this hash or free list */
	/* followed by the scanline data */
} SCAN;			/* buffered scanline */

#define	scandata(sl)	((COLR *)((sl)+1))
#define shash(y)	((y)%HSIZE)

static int	maxpix;		/* maximum number of pixels to buffer */

static SCAN	*freelist;		/* scanline free list */
static SCAN	*hashtab[HSIZE];	/* scanline hash table */

static long	scanbufsiz;		/* size of allocated scanline buffer */

static long	ncall = 0L;	/* number of calls to getpictscan */
static long	nread = 0L;	/* number of scanlines read */
static long	nrecl = 0L;	/* number of scanlines reclaimed */

static int	wrongformat = 0;

static SCAN * claimscan(int	y); 
static COLR * getpictscan(int	y);
static double pict_val(FVECT	vd);
static void rt_compute(float	*pb, int	np);
static gethfunc getexpos;
static SCAN * scanretire(void);
static void initscans(void);
static void donescans(void);


static SCAN *
claimscan(			/* claim scanline from buffers */
	int	y
)
{
	int	hi = shash(y);
	SCAN	*slast;
	register SCAN	*sl;

	for (sl = hashtab[hi]; sl != NULL; sl = sl->next)
		if (sl->y == y)				/* active scanline */
			return(sl);
	for (slast = NULL, sl = freelist; sl != NULL; slast = sl, sl = sl->next)
		if (sl->y == -1 || sl->y == y || sl->next == NULL) {
			if (slast == NULL)		/* remove from free */
				freelist = sl->next;
			else
				slast->next = sl->next;
			if (sl->y == y) {		/* reclaim */
				sl->next = hashtab[hi];
				hashtab[hi] = sl;
				nrecl++;
			}
			return(sl);
		}
	return(scanretire());		/* need more free scanlines */
}


static COLR *
getpictscan(			/* get picture scanline */
	int	y
)
{
	register SCAN	*sl;
	register int	i;
					/* first check our buffers */
	sl = claimscan(y);
	if (sl == NULL)
		memerr("claimscan()");
	sl->lused = ncall++;
	if (sl->y == y)			/* scan hit */
		return(scandata(sl));
					/* else read in replacement */
	if (scanpos[y] < 0) {			/* need to search */
		for (i = y+1; i < curpos; i++)
			if (scanpos[i] >= 0) {
				if (fseek(pictfp, scanpos[i], 0) < 0)
					goto seekerr;
				curpos = i;
				break;
			}
		while (curpos >= y) {
			scanpos[curpos] = ftell(pictfp);
			if (freadcolrs(scandata(sl), pxsiz, pictfp) < 0)
				goto readerr;
			nread++;
			curpos--;
		}
	} else {
		if (curpos != y && fseek(pictfp, scanpos[y], 0) < 0)
			goto seekerr;
		if (freadcolrs(scandata(sl), pxsiz, pictfp) < 0)
			goto readerr;
		nread++;
		curpos = y-1;
	}
	sl->y = y;
	i = shash(y);			/* add to hash list */
	sl->next = hashtab[i];
	hashtab[i] = sl;
	return(scandata(sl));
readerr:
	fprintf(stderr, "%s: picture read error\n", progname);
	exit(1);
seekerr:
	fprintf(stderr, "%s: picture seek error\n", progname);
	exit(1);
}


#ifdef DEBUG
void
pict_stats(void)			/* print out picture read statistics */
{
	static long	lastcall = 0L;	/* ncall at last report */
	static long	lastread = 0L;	/* nread at last report */
	static long	lastrecl = 0L;	/* nrecl at last report */

	if (ncall == lastcall)
		return;
	fprintf(stderr, "%s: %ld scanlines read (%ld reclaimed) in %ld calls\n",
		progname, nread-lastread, nrecl-lastrecl, ncall-lastcall);
	lastcall = ncall;
	lastread = nread;
	lastrecl = nrecl;
}
#endif


static double
pict_val(			/* find picture value for view direction */
	FVECT	vd
)
{
	FVECT	pp;
	FVECT	ip;
	COLOR	res;

	if (pictfp == NULL)
		return(-1.0);
	pp[0] = pictview.vp[0] + vd[0];
	pp[1] = pictview.vp[1] + vd[1];
	pp[2] = pictview.vp[2] + vd[2];
	viewloc(ip, &pictview, pp);
	if (ip[2] <= FTINY || ip[0] < 0. || ip[0] >= 1. ||
			ip[1] < 0. || ip[1] >= 1.)
		return(-1.0);
	colr_color(res, getpictscan((int)(ip[1]*pysiz))[(int)(ip[0]*pxsiz)]);
	return(luminance(res)/exposure);
}


extern double
getviewpix(		/* compute single view pixel */
	int	vh,
	int	vv
)
{
	FVECT	dir;
	float	rt_buf[12];
	double	res;

	if (compdir(dir, vh, vv) < 0)
		return(-1.0);
	npixinvw++;
	if ((res = pict_val(dir)) >= 0.0)
		return(res);
	if (rt_pd.r == -1) {
		npixmiss++;
		return(-1.0);
	}
	rt_buf[0] = ourview.vp[0];
	rt_buf[1] = ourview.vp[1];
	rt_buf[2] = ourview.vp[2];
	rt_buf[3] = dir[0];
	rt_buf[4] = dir[1];
	rt_buf[5] = dir[2];
	rt_compute(rt_buf, 1);
	return(luminance(rt_buf));
}


extern void
getviewspan(		/* compute a span of view pixels */
	int	vv,
	float	*vb
)
{
	float	rt_buf[6*MAXPIX];	/* rtrace send/receive buffer */
	register int	n;		/* number of pixels in buffer */
	short	buf_vh[MAXPIX];		/* pixel positions */
	FVECT	dir;
	register int	vh;

#ifdef DEBUG
	if (verbose)
		fprintf(stderr, "%s: computing view span at %d...\n",
				progname, vv);
#endif
	n = 0;
	for (vh = -hsize; vh <= hsize; vh++) {
		if (compdir(dir, vh, vv) < 0) {		/* not in view */
			vb[vh+hsize] = -1.0;
			continue;
		}
		npixinvw++;
		if ((vb[vh+hsize] = pict_val(dir)) >= 0.0)
			continue;
		if (rt_pd.r == -1) {		/* missing information */
			npixmiss++;
			continue;
		}
						/* send to rtrace */
		if (n >= maxpix) {			/* flush */
			rt_compute(rt_buf, n);
			while (n-- > 0)
				vb[buf_vh[n]+hsize] = luminance(rt_buf+3*n);
		}
		rt_buf[6*n] = ourview.vp[0];
		rt_buf[6*n+1] = ourview.vp[1];
		rt_buf[6*n+2] = ourview.vp[2];
		rt_buf[6*n+3] = dir[0];
		rt_buf[6*n+4] = dir[1];
		rt_buf[6*n+5] = dir[2];
		buf_vh[n++] = vh;
	}
#ifdef DEBUG
	if (verbose)
		pict_stats();
#endif
	if (n > 0) {				/* process pending buffer */
		rt_compute(rt_buf, n);
		while (n-- > 0)
			vb[buf_vh[n]+hsize] = luminance(rt_buf+3*n);
	}
}


static void
rt_compute(		/* process buffer through rtrace */
	float	*pb,
	int	np
)
{
#ifdef DEBUG
	if (verbose && np > 1)
		fprintf(stderr, "%s: sending %d samples to rtrace...\n",
				progname, np);
#endif
	memset(pb+6*np, '\0', 6*sizeof(float));
	if (process(&rt_pd, (char *)pb, (char *)pb, 3*sizeof(float)*(np+1),
			6*sizeof(float)*(np+1)) < 3*sizeof(float)*(np+1)) {
		fprintf(stderr, "%s: rtrace communication error\n",
				progname);
		exit(1);
	}
}


static int
getexpos(			/* get exposure from header line */
	char	*s,
	void	*p
)
{
	char	fmt[32];

	if (isexpos(s))
		exposure *= exposval(s);
	else if (isformat(s)) {
		formatval(fmt, s);
		wrongformat = strcmp(fmt, COLRFMT);
	}
	return(0);
}


extern void
open_pict(			/* open picture file */
	char	*fn
)
{
	if ((pictfp = fopen(fn, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fn);
		exit(1);
	}
	exposure = 1.0;
	getheader(pictfp, getexpos, NULL);
	if (wrongformat || !fscnresolu(&pxsiz, &pysiz, pictfp)) {
		fprintf(stderr, "%s: incompatible picture format\n", fn);
		exit(1);
	}
	initscans();
}


extern void
close_pict(void)			/* done with picture */
{
	if (pictfp == NULL)
		return;
	fclose(pictfp);
	donescans();
	pictfp = NULL;
}


extern void
fork_rtrace(			/* open pipe and start rtrace */
	char	*av[]
)
{
	int	rval;

	rval = open_process(&rt_pd, av);
	if (rval < 0) {
		perror(progname);
		exit(1);
	}
	if (rval == 0) {
		fprintf(stderr, "%s: command not found\n", av[0]);
		exit(1);
	}
	maxpix = rval/(6*sizeof(float));
	if (maxpix > MAXPIX)
		maxpix = MAXPIX;
	maxpix--;
}


extern void
done_rtrace(void)			/* wait for rtrace to finish */
{
	int	status;

	status = close_process(&rt_pd);
	if (status > 0) {
		fprintf(stderr, "%s: bad status (%d) from rtrace\n",
				progname, status);
		exit(1);
	}
	rt_pd.r = -1;
}


static SCAN *
scanretire(void)			/* retire old scanlines to free list */
{
	SCAN	*sold[NRETIRE];
	int	n;
	int	h;
	register SCAN	*sl;
	register int	i;
					/* grab the NRETIRE oldest scanlines */
	sold[n = 0] = NULL;
	for (h = 0; h < HSIZE; h++)
		for (sl = hashtab[h]; sl != NULL; sl = sl->next) {
			for (i = n; i && sold[i-1]->lused > sl->lused; i--)
				if (i < NRETIRE)
					sold[i] = sold[i-1];
			if (i < NRETIRE) {
				sold[i] = sl;
				if (n < NRETIRE)	/* grow list */
					n++;
			}
		}
					/* put scanlines into free list */
	for (i = 0; i < n; i++) {
		h = shash(sold[i]->y);
		sl = hashtab[h];
		if (sl == sold[i])
			hashtab[h] = sl->next;
		else {
			while (sl->next != sold[i])	/* IS in list */
				sl = sl->next;
			sl->next = sold[i]->next;
		}
		if (i > 0) {		/* save oldest as return value */
			sold[i]->next = freelist;
			freelist = sold[i];
		}
	}
	return(sold[0]);
}


static char	*scan_buf;


static void
initscans(void)				/* initialize scanline buffers */
{
	int	scansize;
	register SCAN	*ptr;
	register int	i;
					/* initialize positions */
	scanpos = (long *)bmalloc(pysiz*sizeof(long));
	if (scanpos == NULL)
		memerr("scanline positions");
	for (i = pysiz-1; i >= 0; i--)
		scanpos[i] = -1L;
	curpos = pysiz-1;
					/* clear hash table */
	for (i = 0; i < HSIZE; i++)
		hashtab[i] = NULL;
					/* allocate scanline buffers */
	scansize = sizeof(SCAN) + pxsiz*sizeof(COLR);
#ifdef ALIGNT
	scansize = scansize+(sizeof(ALIGNT)-1) & ~(sizeof(ALIGNT)-1);
#endif
	i = MAXSBUF / scansize;		/* compute number to allocate */
	if (i > HSIZE)
		i = HSIZE;
	scanbufsiz = i*scansize;
	scan_buf = bmalloc(scanbufsiz);	/* get in one big chunk */
	if (scan_buf == NULL)
		memerr("scanline buffers");
	ptr = (SCAN *)scan_buf;
	freelist = NULL;		/* build our free list */
	while (i-- > 0) {
		ptr->y = -1;
		ptr->lused = -1;
		ptr->next = freelist;
		freelist = ptr;
		ptr = (SCAN *)((char *)ptr + scansize);	/* beware of C bugs */
	}
}


static void
donescans(void)				/* free up scanlines */
{
	bfree(scan_buf, scanbufsiz);
	bfree((char *)scanpos, pysiz*sizeof(long));
}
