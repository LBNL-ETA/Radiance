#ifndef lint
static const char	RCSid[] = "$Id: rhd_odraw.c,v 3.15 2003/06/20 00:25:49 greg Exp $";
#endif
/*
 * Routines for drawing samples using depth buffer checks.
 */

#include "standard.h"

#include <GL/glx.h>
#include <GL/glu.h>

#include "rhd_odraw.h"

#ifndef DEPTHEPS
#define DEPTHEPS	0.02		/* depth epsilon */
#endif
#ifndef SAMPSPERBLOCK
#define SAMPSPERBLOCK	1024		/* target samples per image block */
#endif
#ifndef SFREEFRAC
#define SFREEFRAC	0.2		/* fraction to free at a time */
#endif
#ifndef REDRAWTHRESH
#define REDRAWTHRESH	10240		/* number of samples for dissolve */
#endif
#ifndef MAXFAN
#define MAXFAN		32		/* maximum arms in a triangle fan */
#endif
#ifndef MINFAN
#define MINFAN		4		/* minimum arms in a triangle fan */
#endif
#ifndef FANSIZE
#define FANSIZE		3.5		/* fan sizing factor */
#endif

#define NEWMAP		01		/* need to recompute mapping */
#define NEWRGB		02		/* need to remap RGB values */
#define NEWHIST		04		/* clear histogram as well */

struct ODview	*odView;	/* our view list */
int	odNViews;		/* number of views in our list */

struct ODsamp	odS;		/* sample values */

static int	needmapping;	/* what needs doing with tone map */


#define SAMP32	(32*(2*sizeof(short)+sizeof(union ODfunion)+sizeof(TMbright)+\
			6*sizeof(BYTE))+sizeof(int32))

int
odInit(n)				/* initialize drawing routines */
int	n;
{
	int	nbytes, i, j, k, nextsamp, count, blockdiv;
	int	res[2];

	if (odNViews > 0) {		/* deallocate view structures */
		for (i = odNViews; i--; ) {
			free((void *)odView[i].bmap);
			free((void *)odView[i].pmap);
			if (odView[i].emap != NULL)
				free((void *)odView[i].emap);
		}
		free((void *)odView);
		odView = NULL;
		odNViews = 0;
	}
	if (n && n != odS.nsamp) {
				/* round space up to nearest power of 2 */
		nbytes = (n+31)/32 * SAMP32;
		for (i = 1024; nbytes > i-8; i <<= 1)
			;
		n = (i-8)/SAMP32 * 32;
		needmapping = NEWHIST;
	}
	if (n != odS.nsamp) {	/* (re)allocate sample array */
		if (odS.nsamp)
			free(odS.base);
		odS.nsamp = 0;
		if (!n)
			return(0);
		nbytes = (n+31)/32 * SAMP32;
		odS.base = (char *)malloc(nbytes);
		if (odS.base == NULL)
			return(0);
				/* assign larger alignment types earlier */
		odS.f = (union ODfunion *)odS.base;
		odS.redraw = (int32 *)(odS.f + n);
		odS.ip = (short (*)[2])(odS.redraw + n/32);
		odS.brt = (TMbright *)(odS.ip + n);
		odS.chr = (BYTE (*)[3])(odS.brt + n);
		odS.rgb = (BYTE (*)[3])(odS.chr + n);
		odS.nsamp = n;
	}
	if (!n)
		return(0);
				/* allocate view information */
	count = 0;			/* count pixels */
	for (i = 0; dev_auxview(i, res) != NULL; i++)
		count += res[0]*res[1];
	odView = (struct ODview *)malloc(i*sizeof(struct ODview));
	if (odView == NULL)
		return(0);
	odNViews = i;
	blockdiv = sqrt(count/(n/SAMPSPERBLOCK)) + 0.5;
	if (blockdiv < 8) blockdiv = 8;
	nextsamp = 0; count /= blockdiv*blockdiv;	/* # blocks */
	while (i--) {			/* initialize each view */
		dev_auxview(i, res);
		odView[i].hhi = res[0];
		odView[i].hlow = (res[0] + blockdiv/2) / blockdiv;
		if (odView[i].hlow < 1) odView[i].hlow = 1;
		odView[i].vhi = res[1];
		odView[i].vlow = (res[1] + blockdiv/2) / blockdiv;
		if (odView[i].vlow < 1) odView[i].vlow = 1;
		odView[i].emap = NULL;
		odView[i].dmap = NULL;
		odView[i].pmap = (int32 *)calloc(FL4NELS(res[0]*res[1]),
				sizeof(int32));
		if (odView[i].pmap == NULL)
			return(0);
		j = odView[i].hlow*odView[i].vlow;
		odView[i].bmap = (struct ODblock *)malloc(
				j * sizeof(struct ODblock));
		if (odView[i].bmap == NULL)
			return(0);
		DCHECK(count<=0 | nextsamp>=n,
				CONSISTENCY, "counter botch in odInit");
		if (!i) count = j;
		odView[i].sfirst = nextsamp;
		while (j--) {		/* initialize blocks & free lists */
			odView[i].bmap[j].pthresh = FHUGE;
			odView[i].bmap[j].first = k = nextsamp;
			nextsamp += odView[i].bmap[j].nsamp = 
				(n - nextsamp)/count--;
			odView[i].bmap[j].free = k;
			while (++k < nextsamp)
				odS.nextfree(k-1) = k;
			odS.nextfree(k-1) = ENDFREE;
			odView[i].bmap[j].nused = 0;
		}
		odView[i].snext = nextsamp;
		odView[i].n2redraw = 0;
	}
	CLR4ALL(odS.redraw, odS.nsamp);		/* clear redraw flags */
	for (i = odS.nsamp; i--; ) {		/* clear values */
		odS.ip[i][0] = odS.ip[i][1] = -1;
		odS.brt[i] = TM_NOBRT;
	}
	needmapping |= NEWMAP;			/* compute new map on update */
	return(odS.nsamp);			/* return number of samples */
}

#undef SAMP32


int
sampcmp(s0, s1)			/* sample order, descending proximity */
int	*s0, *s1;
{
	register double	diff = odS.closeness(*s1) - odS.closeness(*s0);

	return (diff > FTINY ? 1 : diff < -FTINY ? -1 : 0);
}


int
odAllocBlockSamp(vn, hh, vh, prox)	/* allocate sample from block */
int	vn, hh, vh;
double	prox;
{
	int	si[SAMPSPERBLOCK+SAMPSPERBLOCK/4];
	int	hl, vl;
	VIEW	*vw;
	FVECT	ro, rd;
	int	res[2];
	register struct ODblock	*bp;
	register int	i, j;
					/* get block */
	hl = hh*odView[vn].hlow/odView[vn].hhi;
	vl = vh*odView[vn].vlow/odView[vn].vhi;
	bp = odView[vn].bmap + vl*odView[vn].hlow + hl;
	if (prox > bp->pthresh)
		return(-1);		/* worse than free list occupants */
					/* check for duplicate pixel */
	if (CHK4(odView[vn].pmap, vh*odView[vn].hhi + hh))
		i = bp->first + bp->nsamp;
	else
		i = 0;
	while (i-- > bp->first)
		if (hh == odS.ip[i][0] && vh == odS.ip[i][1]) {	/* found it! */
						/* search free list for it */
			if (i == bp->free)
				break;		/* special case */
			if (bp->free != ENDFREE)
				for (j = bp->free; odS.nextfree(j) != ENDFREE;
						j = odS.nextfree(j))
					if (odS.nextfree(j) == i) {
						odS.nextfree(j) =
							odS.nextfree(i);
						bp->nused++;
						goto gotit;
					}
			if (prox >= 0.99*odS.closeness(i))
				return(-1);	/* previous sample is fine */
			goto gotit;
		}
	if (bp->free != ENDFREE) {	/* allocate from free list */
		i = bp->free;
		if (odS.ip[i][0] >= 0 & odS.ip[i][1] >= 0)
			CLR4(odView[vn].pmap, odS.ip[i][1]*odView[vn].hhi +
							odS.ip[i][0]);
		bp->free = odS.nextfree(i);
		bp->nused++;
		goto gotit;
	}
	DCHECK(bp->nsamp<=0, CONSISTENCY,
			"no available samples in odAllocBlockSamp");
	DCHECK(bp->nsamp > sizeof(si)/sizeof(si[0]), CONSISTENCY,
			"too many samples in odAllocBlockSamp");
					/* free some samples */
	if ((vw = dev_auxview(vn, res)) == NULL)
		error(CONSISTENCY, "bad view number in odAllocBlockSamp");
	for (i = bp->nsamp; i--; )	/* figure out which are worse */
		si[i] = bp->first + i;
	qsort((char *)si, bp->nsamp, sizeof(int), sampcmp);
	i = bp->nsamp*SFREEFRAC + .5;	/* put them into free list */
	if (i >= bp->nsamp) i = bp->nsamp-1;	/* paranoia */
	bp->pthresh = odS.closeness(si[i]);	/* new proximity threshold */
	while (--i > 0) {
		odS.nextfree(si[i]) = bp->free;
		bp->free = si[i];
		bp->nused--;
	}
	i = si[0];			/* use worst sample */
	CLR4(odView[vn].pmap, odS.ip[i][1]*odView[vn].hhi + odS.ip[i][0]);
gotit:
	odS.ip[i][0] = hh;
	odS.ip[i][1] = vh;
	odS.closeness(i) = prox;
	SET4(odView[vn].pmap, vh*odView[vn].hhi + hh);
	return(i);
}


odSample(c, d, p)			/* add a sample value */
COLR	c;
FVECT	d, p;
{
	FVECT	disp;
	double	d0, d1, h, v, prox;
	register VIEW	*vw;
	int	hh, vh;
	int	res[2];
	register int	i, id;

	DCHECK(odS.nsamp<=0, CONSISTENCY, "no samples allocated in odSample");
						/* add value to each view */
	for (i = 0; (vw = dev_auxview(i, res)) != NULL; i++) {
		DCHECK(i>=odNViews, CONSISTENCY, "too many views in odSample");
		CHECK(vw->type!=VT_PER, INTERNAL,
				"cannot handle non-perspective views");
		if (p != NULL) {		/* compute view position */
			VSUB(disp, p, vw->vp);
			d0 = DOT(disp, vw->vdir);
			if (d0 <= vw->vfore+FTINY)
				continue;		/* too close */
		} else {
			VCOPY(disp, d);
			d0 = DOT(disp, vw->vdir);
			if (d0 <= FTINY)		/* behind view */
				continue;
		}
		h = DOT(disp,vw->hvec)/(d0*vw->hn2) + 0.5 - vw->hoff;
		if (h < 0. || h >= 1.)
			continue;			/* left or right */
		v = DOT(disp,vw->vvec)/(d0*vw->vn2) + 0.5 - vw->voff;
		if (v < 0. || v >= 1.)
			continue;			/* above or below */
		hh = h * res[0];
		vh = v * res[1];
		if (odView[i].dmap != NULL) {		/* check depth */
			d1 = odView[i].dmap[vh*res[0] + hh];
			if (d1 < 0.99*FHUGE && (d0 > (1.+DEPTHEPS)*d1 ||
						(1.+DEPTHEPS)*d0 < d1))
			continue;			/* occlusion error */
		}
		if (p != NULL) {		/* compute closeness (sin^2) */
			d1 = DOT(disp, d);
			prox = 1. - d1*d1/DOT(disp,disp);
		} else
			prox = 0.;
						/* allocate sample */
		id = odAllocBlockSamp(i, hh, vh, prox);
		if (id < 0)
			continue;		/* not good enough */
							/* convert color */
		tmCvColrs(&odS.brt[id], odS.chr[id], (COLR *)c, 1);
		if (imm_mode | needmapping)		/* if immediate mode */
			needmapping |= NEWRGB;		/* map it later */
		else					/* else map it now */
			tmMapPixels(odS.rgb[id], &odS.brt[id], odS.chr[id], 1);
		SET4(odS.redraw, id);			/* mark for redraw */
		odView[i].n2redraw++;
	}
}


odRemap(newhist)			/* recompute tone mapping */
int	newhist;
{
	needmapping |= NEWMAP|NEWRGB;
	if (newhist)
		needmapping |= NEWHIST;
}


odRedrawAll()				/* mark all samples for redraw */
{
	register int	i;

	if ((needmapping&(NEWMAP|NEWRGB)) == (NEWMAP|NEWRGB))
		return;			/* will be called later, anyway */
	for (i = odS.nsamp; i--; )
		if (odS.ip[i][0] >= 0)
			SET4(odS.redraw, i);
					/* not right, but not important */
	for (i = 0; i < odNViews; i++)
		odView[i].n2redraw = odView[i].snext - odView[i].sfirst;
}


odRedraw(vn, hmin, vmin, hmax, vmax)	/* redraw view region */
int	vn, hmin, vmin, hmax, vmax;
{
	int	i, j;
	register struct ODblock	*bp;
	register int	k;

	if (vn<0 | vn>=odNViews)
		return;
				/* check view limits */
	if (hmin < 0) hmin = 0;
	if (hmax >= odView[vn].hhi) hmax = odView[vn].hhi-1;
	if (vmin < 0) vmin = 0;
	if (vmax >= odView[vn].vhi) vmax = odView[vn].vhi-1;
	if (hmax <= hmin | vmax <= vmin)
		return;
				/* convert to low resolution */
	hmin = hmin * odView[vn].hlow / odView[vn].hhi;
	hmax = hmax * odView[vn].hlow / odView[vn].hhi;
	vmin = vmin * odView[vn].vlow / odView[vn].vhi;
	vmax = vmax * odView[vn].vlow / odView[vn].vhi;
				/* mark block samples for redraw, inclusive */
	for (i = hmin; i <= hmax; i++)
		for (j = vmin; j <= vmax; j++) {
			bp = odView[vn].bmap + j*odView[vn].hlow + i;
			for (k = bp->nsamp; k--; )
				if (odS.ip[bp->first+k][0] >= 0) {
					SET4(odS.redraw, bp->first+k);
					odView[vn].n2redraw++;
				}
		}
}


odDepthMap(vn, dm)			/* assign depth map for view */
int	vn;
GLfloat	*dm;
{
	double	d0, d1;
	int	i, j, hmin, hmax, vmin, vmax;
	register int	k, l;

	if (dm == NULL) {			/* free edge map */
		if (vn<0 | vn>=odNViews)
			return;			/* too late -- they're gone! */
		if (odView[vn].emap != NULL)
			free((void *)odView[vn].emap);
		odView[vn].emap = NULL;
		odView[vn].dmap = NULL;
		return;
	}
	DCHECK(vn<0 | vn>=odNViews, CONSISTENCY,
			"bad view number in odDepthMap");
	odView[vn].dmap = dm;			/* initialize edge map */
	if (odView[vn].emap == NULL) {
		odView[vn].emap = (int32 *)malloc(
			FL4NELS(odView[vn].hlow*odView[vn].vlow)*sizeof(int32));
		if (odView[vn].emap == NULL)
			error(SYSTEM, "out of memory in odDepthMap");
	}
	CLR4ALL(odView[vn].emap, odView[vn].hlow*odView[vn].vlow);
						/* compute edge map */
	vmin = odView[vn].vhi;			/* enter loopsville */
	for (j = odView[vn].vlow; j--; ) {
		vmax = vmin;
		vmin = j*odView[vn].vhi/odView[vn].vlow;
		hmin = odView[vn].hhi;
		for (i = odView[vn].hlow; i--; ) {
			hmax = hmin;
			hmin = i*odView[vn].hhi/odView[vn].hlow;
			for (l = vmin; l < vmax; l++) {	/* vertical edges */
				d1 = dm[l*odView[vn].hhi+hmin];
				for (k = hmin+1; k < hmax; k++) {
					d0 = d1;
					d1 = dm[l*odView[vn].hhi+k];
					if (d0 > (1.+DEPTHEPS)*d1 ||
						(1.+DEPTHEPS)*d0 < d1) {
						SET4(odView[vn].emap,
							j*odView[vn].hlow + i);
						break;
					}
				}
				if (k < hmax)
					break;
			}
			if (l < vmax)
				continue;
			for (k = hmin; k < hmax; k++) {	/* horizontal edges */
				d1 = dm[vmin*odView[vn].hhi+k];
				for (l = vmin+1; l < vmax; l++) {
					d0 = d1;
					d1 = dm[l*odView[vn].hhi+k];
					if (d0 > (1.+DEPTHEPS)*d1 ||
						(1.+DEPTHEPS)*d0 < d1) {
						SET4(odView[vn].emap,
							j*odView[vn].hlow + i);
						break;
					}
				}
				if (l < vmax)
					break;
			}
		}
	}
}


odUpdate(vn)				/* update this view */
int	vn;
{
	static short	primes[] = {9431,6803,4177,2659,1609,887,587,251,47,1};
	int	myprime;
	register int	i, n;

	DCHECK(vn<0 | vn>=odNViews, CONSISTENCY,
			"bad view number in odUpdate");
					/* need to do some tone mapping? */
	if (needmapping & NEWRGB) {
		if (needmapping & NEWMAP) {
			if (needmapping & NEWHIST)
				tmClearHisto();
			needmapping &= ~NEWHIST;
			if (tmAddHisto(odS.brt,odS.nsamp,1) != TM_E_OK)
				return;
			if (tmComputeMapping(0.,0.,0.) != TM_E_OK)
				return;
			needmapping &= ~NEWMAP;
			odRedrawAll();			/* redraw everything */
		}
		if (tmMapPixels((BYTE *)(odS.rgb), odS.brt,
				(BYTE *)(odS.chr), odS.nsamp) != TM_E_OK)
			return;
		needmapping &= ~NEWRGB;
	}
	if (odView[vn].n2redraw <= 0)
		return;
#if REDRAWTHRESH
	if (odView[vn].n2redraw < REDRAWTHRESH)
		goto quickdraw;
					/* pick a good prime step size */
	n = odView[vn].snext - odView[vn].sfirst;
	for (i = 0; primes[i]<<5 >= n; i++)
		;
	while ((myprime = primes[i++]) > 1)
		if (n % myprime)
			break;
					/* dissolve in new samples */
	for (i = odView[vn].sfirst; n-- > 0; i += myprime) {
		if (i >= odView[vn].snext)
			i -= odView[vn].snext - odView[vn].sfirst;
		if (CHK4(odS.redraw, i)) {
			odDrawSamp(vn, i);
			CLR4(odS.redraw, i);
		}
	}
	odView[vn].n2redraw = 0;
	return;
quickdraw:				/* quicker sparse flag checking */
#endif
					/* redraw samples at end */
	for (i = odView[vn].snext-31; i < odView[vn].snext; i++)
		if (CHK4(odS.redraw, i)) {
			odDrawSamp(vn, i);
			CLR4(odS.redraw, i);
		}
					/* faster flag checks in middle */
	for (n = odView[vn].snext>>5; n-- > (odView[vn].sfirst+0x1f)>>5; )
		for (i = 0; odS.redraw[n]; i++)		/* skips faster */
			if (odS.redraw[n] & 1L<<i) {
				odDrawSamp(vn, (n<<5)+i);
				odS.redraw[n] &= ~(1L<<i);
			}
					/* redraw samples at beginning */
	for (i = odView[vn].sfirst; i < odView[vn].sfirst+31; i++)
		if (CHK4(odS.redraw, i)) {
			odDrawSamp(vn, i);
			CLR4(odS.redraw, i);
		}
	odView[vn].n2redraw = 0;
}


					/* this turned out to be unnecessary */
#if 0
static
clip_end(p, o, vp)			/* clip line segment to view */
GLshort	p[3];
short	o[2];
register struct ODview	*vp;
{
	if (p[0] < 0) {
		p[1] = -o[0]*(p[1]-o[1])/(p[0]-o[0]) + o[1];
		p[2] = -o[0]*p[2]/(p[0]-o[0]);
		p[0] = 0;
	} else if (p[0] >= vp->hhi) {
		p[1] = (vp->hhi-1-o[0])*(p[1]-o[1])/(p[0]-o[0]) + o[1];
		p[2] = (vp->hhi-1-o[0])*p[2]/(p[0]-o[0]);
		p[0] = vp->hhi-1;
	}
	if (p[1] < 0) {
		p[0] = -o[1]*(p[0]-o[0])/(p[1]-o[1]) + o[0];
		p[2] = -o[1]*p[2]/(p[1]-o[1]);
		p[1] = 0;
	} else if (p[1] >= vp->vhi) {
		p[0] = (vp->vhi-1-o[1])*(p[0]-o[0])/(p[1]-o[1]) + o[0];
		p[2] = (vp->vhi-1-o[1])*p[2]/(p[1]-o[1]);
		p[1] = vp->vhi-1;
	}
}
#endif


static int
make_arms(ar, cp, vp, sz)		/* make arms for triangle fan */
GLshort	ar[MAXFAN][3];
short	cp[2];
register struct ODview	*vp;
double	sz;
{
	int	na, dv;
	double	hrad, vrad, phi;
	register int	i;

	DCHECK(sz > 1, CONSISTENCY, "super-unary size in make_arms");
	na = MAXFAN*sz + 0.5;			/* keep arc length constant */
	if (na < MINFAN) na = MINFAN;
	hrad = FANSIZE*sz*vp->hhi/vp->hlow;
	vrad = FANSIZE*sz*vp->vhi/vp->vlow;
	if (hrad*vrad < 2.25)
		hrad = vrad = 1.5;
	dv = OMAXDEPTH*sz + 0.5;
	for (i = 0; i < na; i++) {
		phi = (2.*PI)*i/na;
		ar[i][0] = cp[0] + tcos(phi)*hrad + 0.5;
		ar[i][1] = cp[1] + tsin(phi)*vrad + 0.5;
		ar[i][2] = dv;
		/* clip_end(ar[i], cp, vp); */
	}
	return(na);
}


static int
depthchange(vp, x0, y0, x1, y1)		/* check depth discontinuity */
register struct ODview	*vp;
int	x0, y0, x1, y1;
{
	register double	d0, d1;

	DCHECK(x0<0 | x0>=vp->hhi | y0<0 | y0>=vp->vhi,
			CONSISTENCY, "coordinates off view in depthchange");

	if (x1<0 | x1>=vp->hhi | y1<0 | y1>=vp->vhi)
		return(1);

	d0 = vp->dmap[y0*vp->hhi + x0];
	d1 = vp->dmap[y1*vp->hhi + x1];

	return((1.+DEPTHEPS)*d0 < d1 || d0 > (1.+DEPTHEPS)*d1);
}


static
clip_edge(p, o, vp)			/* clip line segment to depth edge */
GLshort	p[3];
short	o[2];
register struct ODview	*vp;
{
	int	x, y, xstep, ystep, rise, rise2, run, run2, n;

	DCHECK(vp->dmap==NULL, CONSISTENCY,
			"clip_edge called with no depth map");
	x = o[0]; y = o[1];
	run = p[0] - x;
	xstep = run > 0 ? 1 : -1;
	run *= xstep;
	rise = p[1] - y;
	ystep = rise > 0 ? 1 : -1;
	rise *= ystep;
	rise2 = run2 = 0;
	if (rise > run) rise2 = 1;
	else run2 = 1;
	n = rise + run;
	while (n--)			/* run out arm, checking depth */
		if (run2 > rise2) {
			x += xstep;
			rise2 += rise;
			if (depthchange(vp, x-xstep, y, x, y))
				break;
		} else {
			y += ystep;
			run2 += run;
			if (depthchange(vp, x, y-ystep, x, y))
				break;
		}
	if (n < 0)			/* found something? */
		return;
	if (run > rise)
		p[2] = (x - o[0])*p[2]/(p[0] - o[0]);
	else
		p[2] = (y - o[1])*p[2]/(p[1] - o[1]);
	p[0] = x;
	p[1] = y;
}


static int
getblock(vp, h, v)			/* get block index */
register struct ODview	*vp;
register int	h, v;
{
	if (h<0 | h>=vp->hhi | v<0 | v>=vp->vhi)
		return(-1);
	return(h*vp->hlow/vp->hhi + v*vp->vlow/vp->vhi*vp->hlow);
}


static int
blockedge(vp, bi0, bi1)			/* check for edge between blocks? */
register struct ODview	*vp;
register int	bi0, bi1;
{
	if (bi1 == bi0)
		return(0);		/* same block */
	if (bi1 < 0)
		return(1);		/* end off view */
	if (CHK4(vp->emap, bi1))
		return(1);		/* end block has edges */
	if (bi1 == bi0+1 || bi1 == bi0-1 ||
			bi1 == bi0+vp->hlow || bi1 == bi0-vp->hlow)
		return(0);		/* end in adjacent block -- no edges */
	return(1);			/* conservative for rarer case */
}


odDrawSamp(vn, id)			/* draw view sample */
int	vn;
register int	id;
{
	GLshort	arm[MAXFAN][3];
	int	narms, blockindex;
	register struct ODview	*vp;
	double	size;
	int	home_edges;
	register int	i;

	vp = odView + vn;
	blockindex = getblock(vp, odS.ip[id][0], odS.ip[id][1]);
	DCHECK(blockindex<0, CONSISTENCY, "bad sample handed to odDrawSamp");
	DCHECK(vp->bmap[blockindex].nused <= 0,
			CONSISTENCY, "bad in-use count in odDrawSamp");
					/* create triangle fan */
	size = 1./sqrt((double)vp->bmap[blockindex].nused);
	narms = make_arms(arm, odS.ip[id], vp, size);
	if (vp->emap != NULL) {		/* check for edge collisions */
		home_edges = CHK4(vp->emap, blockindex);
		for (i = 0; i < narms; i++)
			if (home_edges || blockedge(vp, blockindex,
					getblock(vp, arm[i][0], arm[i][1])))
				clip_edge(arm[i], odS.ip[id], vp);
	}
					/* draw triangle fan */
	glColor3ub(odS.rgb[id][0], odS.rgb[id][1], odS.rgb[id][2]);
	glBegin(GL_TRIANGLE_FAN);
	glVertex3s((GLshort)odS.ip[id][0], (GLshort)odS.ip[id][1], (GLshort)0);
	for (i = 0; i < narms; i++)
		glVertex3sv(arm[i]);
	glVertex3sv(arm[0]);		/* connect last to first */
	glEnd();
}
