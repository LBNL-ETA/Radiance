#ifndef lint
static const char	RCSid[] = "$Id: ranimove2.c,v 3.2 2003/02/25 02:47:24 greg Exp $";
#endif
/*
 *  ranimove2.c
 *
 *  Frame refinement routines for ranimate(1).
 *
 *  Created by Gregory Ward on Wed Jan 08 2003.
 */

#include "copyright.h"

#include "ranimove.h"
#include "random.h"


#define	HL_ERR	0.32		/* highlight error threshold */

int	cerrzero;		/* is cerrmap all zeroes? */


int
refine_first()			/* initial refinement pass */
{
	int	*esamp = (int *)zprev;	/* OK to reuse */
	int	hl_erri = errori(HL_ERR);
	int	nextra = 0;
	int	x, y, xp, yp;
	int	neigh;
	register int	n, np;

	if (sizeof(int) < sizeof(*zprev))
		error(CONSISTENCY, "code error in refine_first");
	if (!silent) {
		printf("\tFirst refinement pass...");
		fflush(stdout);
	}
	bzero((void *)esamp, sizeof(int)*hres*vres);
	/*
	 * In our initial pass, we look for lower error pixels from
	 * the same objects in the previous frame, and copy them here.
	 */
	for (y = vres; y--; )
	    for (x = hres; x--; ) {
		n = fndx(x, y);
		if (obuffer[n] == OVOID)
			continue;
		if (xmbuffer[n] == MO_UNK)
			continue;
		xp = x + xmbuffer[n];
		if ((xp < 0 | xp >= hres))
			continue;
		yp = y + ymbuffer[n];
		if ((yp < 0 | yp >= vres))
			continue;
		np = fndx(xp, yp);
					/* make sure we hit same object */
		if (oprev[np] != obuffer[n])
			continue;
					/* is previous frame error lower? */
		if (aprev[np] < AMIN + ATIDIFF)
			continue;
		if (aprev[np] <= abuffer[n] + ATIDIFF)
			continue;
					/* shadow & highlight detection */
		if (abuffer[n] > hl_erri &&
				getclosest(&neigh, 1, x, y) &&
				bigdiff(cbuffer[neigh], cprev[np],
					HL_ERR*(.9+.2*frandom())))
			continue;
		abuffer[n] = aprev[np] - ATIDIFF;
		copycolor(cbuffer[n], cprev[np]);
		esamp[n] = 1;		/* record extrapolated sample */
		nextra++;
	    }
	for (n = hres*vres; n--; )	/* update sample counts */
		if (esamp[n])
			sbuffer[n] = 1;
	if (!silent)
		printf("extrapolated %d pixels\n", nextra);
	return(1);
}


/*
 * We use a recursive computation of the conspicuity
 * map to avoid associated memory costs and simplify
 * coding.  We create a virtual image pyramid, pooling
 * variance calculations, etc.  The top of the pyramid
 * corresponds to the foveal resolution, as there should
 * not be any interesting mechanisms above this level.
 */

#define CSF_C0	1.14
#define CSF_C1	0.67
#define CSF_C2	1.7
#define CSF_S1	6.1
#define CSF_S2	7.3
#define CSF_P1	45.9
#define CSF_PC	(30./45.9*CSF_P1)
#define CSF_VR0	0.15
#define CSF_VRC	80.

struct ConspSum {
	COLOR		vsum;		/* value sum */
	COLOR		v2sum;		/* value^2 sum */
	long		nsamp;		/* number of samples */
	long		xmsum;		/* x-motion sum */
	long		ymsum;		/* y-motion sum */
	int		npix;		/* number of pixels */
	double		hls;		/* high-level saliency */
};

static double	pixel_deg;	/* base pixel frequency */
static int	fhsiz, fvsiz;	/* foveal subimage size */

static void
clr_consp(cs)			/* initialize a conspicuity sum */
register struct ConspSum	*cs;
{
	if (cs == NULL)
		return;
	setcolor(cs->vsum, 0., 0., 0.);
	setcolor(cs->v2sum, 0., 0., 0.);
	cs->nsamp = 0;
	cs->xmsum = cs->ymsum = 0;
	cs->npix = 0;
	cs->hls = 0;
}

static void
sum_consp(cdest, cs)		/* sum in conspicuity result */
register struct ConspSum	*cdest, *cs;
{
	if ((cdest == NULL | cs == NULL))
		return;
	addcolor(cdest->vsum, cs->vsum);
	addcolor(cdest->v2sum, cs->v2sum);
	cdest->nsamp += cs->nsamp;
	cdest->xmsum += cs->xmsum;
	cdest->ymsum += cs->ymsum;
	cdest->npix += cs->npix;
	if (cs->hls > cdest->hls)
		cdest->hls = cs->hls;
}

static void
est_consp(x0,y0,x1,y1, cs)	/* estimate error conspicuity & update */
int	x0, y0, x1, y1;
register struct ConspSum	*cs;
{
	double	rad2, mtn2, cpd, vm, vr, csf, eest;
						/* do we care? */
	if (cs->hls <= FTINY)
		return;
						/* get relative error */
	if (cs->nsamp < NSAMPOK) {
		int	neigh[NSAMPOK];		/* gather neighbors */
		eest = comperr(neigh,
			getclosest(neigh, NSAMPOK, (x0+x1)>>1, (y0+y1)>>1),
				cs->nsamp);
	} else
		eest = estimaterr(cs->vsum, cs->v2sum, cs->nsamp, cs->nsamp);
	
	if ((x0 == x1-1 & y0 == y1-1)) {	/* update pixel error */
		int	n = fndx(x0, y0);
		int	ai;
		int	ne;
		if (sbuffer[n] >= 255) {
			abuffer[n] = ADISTANT;
		} else {
			ai = errori(eest);
			if (ai < AMIN) ai = AMIN;
			else if (ai >= ADISTANT) ai = ADISTANT-1;
			abuffer[n] = ai;
						/* can't improve on closest */
			if (!cs->nsamp && getclosest(&ne, 1, x0, y0) &&
					abuffer[ne] < ai &&
					abuffer[ne] >= AMIN)
				abuffer[n] = abuffer[ne];
		}
	}
						/* compute radius^2 */
	rad2 = 0.125*((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

						/* average motion^2 */
	mtn2 = (double)cs->xmsum*cs->xmsum + (double)cs->ymsum*cs->ymsum;
	mtn2 /= (double)(cs->npix*cs->npix);
						/* motion blur hides us? */
	if (mblur*mblur*mtn2 >= 4.*rad2)
		return;
						/* too small to see? */
	cpd = pixel_deg * pixel_deg / rad2;
	if (cpd > CSF_PC*CSF_PC)
		return;
	cpd = sqrt(cpd);
						/* compute CSF [Daley98] */
	vm = rate * sqrt(mtn2) / pixel_deg;
	vr = cs->hls/hlsmax*vm + CSF_VR0;	/* use hls tracking eff. */
	if (vr > CSF_VRC) vr = CSF_VRC;
	vr = vm - vr;
	if (vr < 0) vr = -vr;
	csf = log(CSF_C2*(1./3.)*vr);
	if (csf < 0) csf = -csf;
	csf = CSF_S1 + CSF_S2*csf*csf*csf;
	csf *= CSF_C0*CSF_C2*4.*PI*PI*CSF_C1*CSF_C1*cpd*cpd;
	csf *= exp(-CSF_C1*4.*PI/CSF_P1*(CSF_C2*vr + 2.)*cpd);
						/* compute visible error */
	eest = eest*csf/ndthresh - 1.;
	if (eest <= FTINY)
		return;
						/* scale by saleincy */
	eest *= cs->hls;
						/* worth the bother? */
	if (eest <= .01)
		return;
						/* sum into map */
	for ( ; y0 < y1; y0++) {
		float	*em0 = cerrmap + fndx(x0, y0);
		register float	*emp = em0 + (x1-x0);
		while (emp-- > em0)
			*emp += eest;
	}
	cerrzero = 0;
}

static void
subconspicuity(x0,y0,x1,y1, cs)	/* compute subportion of conspicuity */
int	x0, y0, x1, y1;
struct ConspSum	*cs;
{
	struct ConspSum	mysum;
	int	i;

	if ((x0 >= x1 | y0 >= y1))
		error(CONSISTENCY, "bad call to subconspicuity");

	clr_consp(&mysum);			/* prepare sum */

	if ((x0 == x1-1 & y0 == y1-1)) {	/* single pixel */
		double	hls;
		register int	n = fndx(x0, y0);
		if (sbuffer[n]) {
			copycolor(mysum.vsum, cbuffer[n]);
			copycolor(mysum.v2sum, val2map[n]);
			mysum.nsamp = sbuffer[n];
		}
		if ((mysum.xmsum = xmbuffer[n]) == MO_UNK)
			mysum.xmsum = 0;
		else
			mysum.ymsum = ymbuffer[n];
		mysum.npix = 1;
						/* max. hls in fovea */
		mysum.hls = obj_prio(obuffer[n]);
		if (x0 >= fhsiz) {
			hls = obj_prio(obuffer[fndx(x0-fhsiz,y0)]);
			if (hls > mysum.hls) mysum.hls = hls;
		}
		if (x0 < hres-fhsiz) {
			hls = obj_prio(obuffer[fndx(x0+fhsiz,y0)]);
			if (hls > mysum.hls) mysum.hls = hls;
		}
		if (y0 >= fvsiz) {
			hls = obj_prio(obuffer[fndx(x0,y0-fvsiz)]);
			if (hls > mysum.hls) mysum.hls = hls;
		}
		if (y0 < vres-fvsiz) {
			hls = obj_prio(obuffer[fndx(x0,y0+fvsiz)]);
			if (hls > mysum.hls) mysum.hls = hls;
		}
	} else if (x0 == x1-1) {		/* vertical pair */
		for (i = y0 ; i < y1; i++)
			subconspicuity(x0, i, x1, i+1, &mysum);
	} else if (y0 == y1-1) {		/* horizontal pair */
		for (i = x0 ; i < x1; i++)
			subconspicuity(i, y0, i+1, y1, &mysum);
	} else {				/* rectangle */
		subconspicuity(x0, y0, (x0+x1)>>1, (y0+y1)>>1, &mysum);
		subconspicuity((x0+x1)>>1, y0, x1, (y0+y1)>>1, &mysum);
		subconspicuity(x0, (y0+y1)>>1, (x0+x1)>>1, y1, &mysum);
		subconspicuity((x0+x1)>>1, (y0+y1)>>1, x1, y1, &mysum);
	}
						/* update conspicuity */
	est_consp(x0, y0, x1, y1, &mysum);
						/* sum into return value */
	sum_consp(cs, &mysum);
}

void
conspicuity()			/* compute conspicuous error map */
{
	int	fhres, fvres;
	int	fx, fy;
					/* reuse previous z-buffer */
	cerrmap = (float *)zprev;
	bzero((void *)cerrmap, sizeof(float)*hres*vres);
	cerrzero = 1;
					/* compute base pixel frequency */
	pixel_deg = .5*(hres/vw.horiz + vres/vw.vert);
					/* compute foveal resolution */
	fhres = vw.horiz/FOV_DEG + 0.5;
	if (fhres <= 0) fhres = 1;
	else if (fhres > hres) fhres = hres;
	fvres = vw.vert/FOV_DEG + 0.5;
	if (fvres <= 0) fvres = 1;
	else if (fvres > vres) fvres = vres;
	fhsiz = hres/fhres;
	fvsiz = vres/fvres;
					/* call our foveal subroutine */
	for (fy = fvres; fy--; )
		for (fx = fhres; fx--; )
			subconspicuity(hres*fx/fhres, vres*fy/fvres,
					hres*(fx+1)/fhres, vres*(fy+1)/fvres,
					NULL);
}


/*
 * The following structure is used to collect data on the
 * initial error in the ambient value estimate, in order
 * to correct for it in the subsequent frames.
 */
static struct AmbSum {
	double		diffsum[3];	/* sum of (correct - ambval) */
	long		nsamps;		/* number of values in sum */
}		*asump = NULL;


static int
ppri_cmp(pp1, pp2)		/* pixel priority comparison */
const void *pp1, *pp2;
{
	double	se1 = cerrmap[*(const int *)pp1];
	double	se2 = cerrmap[*(const int *)pp2];
	int	adiff;
					/* higher conspicuity to front */
	if (se1 < se2) return(1);
	if (se1 > se2) return(-1);
					/* else higher error to front */
	adiff = (int)abuffer[*(const int *)pp1] -
			(int)abuffer[*(const int *)pp2];
	if (adiff)
		return(adiff);
					/* else fewer samples to front */
	return((int)sbuffer[*(const int *)pp1] -
			(int)sbuffer[*(const int *)pp2]);
}


static int
ray_refine(n)			/* refine the given pixel by tracing a ray */
register int	n;
{
	RAY	ir;
	int	neigh[NSAMPOK];
	int	nc;
	COLOR	ctmp;
	int	i;

	if (n < 0) {			/* fetching stragglers */
		if (nprocs <= 1 || !ray_presult(&ir, 0))
			return(-1);
		n = ir.rno;
	} else {			/* else tracing a new ray */
		double	hv[2];
		if (sbuffer[n] >= 255)	/* reached limit? */
			return(-1);
		sample_pos(hv, n%hres, n/hres, sbuffer[n]);
		ir.rmax = viewray(ir.rorg, ir.rdir, &vw, hv[0], hv[1]);
		if (ir.rmax < -FTINY)
			return(-1);
		if (nprocs > 1) {
			int	rval;
			rayorigin(&ir, NULL, PRIMARY, 1.0);
			ir.rno = n;
			rval = ray_pqueue(&ir);
			if (!rval)
				return(-1);
			if (rval < 0)
				quit(1);
			n = ir.rno;
		} else
			ray_trace(&ir);
	}
	if (abuffer[n] == ALOWQ && asump != NULL) {
		if (sbuffer[n] != 1)
			error(CONSISTENCY, "bad code in ray_refine");
		if (getambcolor(ctmp, obuffer[n]) &&
				(colval(ctmp,RED) > 0.01 &
				 colval(ctmp,GRN) > 0.01 &
				 colval(ctmp,BLU) > 0.01)) {
			for (i = 0; i < 3; i++)
				asump->diffsum[i] +=
				    (colval(ir.rcol,i) - colval(cbuffer[n],i))
					/ colval(ctmp,i);
			asump->nsamps++;
		}
		sbuffer[n] = 0;
	}
	setcolor(ctmp,
		colval(ir.rcol,RED)*colval(ir.rcol,RED),
		colval(ir.rcol,GRN)*colval(ir.rcol,GRN),
		colval(ir.rcol,BLU)*colval(ir.rcol,BLU));
	if (!sbuffer[n]) {			/* first sample */
		copycolor(cbuffer[n], ir.rcol);
		copycolor(val2map[n], ctmp);
		abuffer[n] = AHIGHQ;
		sbuffer[n] = 1;
	} else {				/* else sum in sample */
		addcolor(cbuffer[n], ir.rcol);
		addcolor(val2map[n], ctmp);
		sbuffer[n]++;
	}
	return(n);
}


static long
refine_rays(nrays)		/* compute refinement rays */
long	nrays;
{
	int	*pord;
	int	ntodo;
	long	rdone;
	int	i;
					/* skip if nothing significant */
	if (ndtset && cerrzero)
		return;
					/* initialize priority list */
	pord = (int *)malloc(sizeof(int)*hres*vres);
	for (i = hres*vres; i--; )
		pord[i] = i;
					/* sort our priorities */
	ntodo = hres*vres;
	if (nrays < ntodo)
		qsort((void *)pord, hres*vres, sizeof(int), ppri_cmp);
	i = 0;
					/* trace rays in list */
	for (rdone = 0; rdone < nrays; rdone++) {
		if (ndtset && i >= 1000 && cerrmap[pord[i]] <= FTINY)
			ntodo = i;
		if (i >= ntodo) {	/* redo conspicuity & priority */
			while (ray_refine(-1) >= 0)
				;
			conspicuity();
			if (ndtset && cerrzero)
				break;
			qsort((void *)pord, hres*vres, sizeof(int), ppri_cmp);
			ntodo = hres*vres/8;
			i = 0;
		}
					/* sample next pixel */
		ray_refine(pord[i++]);
	}
					/* clean up and return */
	while (ray_refine(-1) >= 0)
		;
	free((void *)pord);
	return(rdone);
}


int
refine_frame(pass)		/* refine current frame */
int	pass;
{
	static double	rtime_used = 0;
	static long	ray_cnt = 0;
	static double	ctime_used = 0;
	static int	csp_cnt = 0;
	int	timed = (fcur > fbeg | pass > 0 | quickstart);
	double	time_start, rtime_start, time_done;
	struct AmbSum	myAmbSum;
	long	rays_todo, nr;
	register int	n;
					/* IBR refinement? */
	if ((pass == 0 & fcur > fbeg))
		return(refine_first());
					/* any time left? */
	time_start = getTime();
	if (timed) {
		if (time_start >= frm_stop)
			goto nomore;
		if (csp_cnt > 0 && time_start + ctime_used/csp_cnt >= frm_stop)
			goto nomore;
	}
	asump = NULL;			/* use resampling to update ambval? */
	if (!curparams->ambounce && hirendparams.ambounce) {
		myAmbSum.diffsum[RED] =
		myAmbSum.diffsum[GRN] =
		myAmbSum.diffsum[BLU] = 0;
		myAmbSum.nsamps = 0;
		asump = &myAmbSum;
	}
					/* initialize value-squared map */
	if (val2map == NULL) {
		val2map = cprev;	/* OK to reuse at this point */
		n = (asump == NULL) ? hres*vres : 0;
		while (n--)
			if (sbuffer[n])
				setcolor(val2map[n],
				colval(cbuffer[n],RED)*colval(cbuffer[n],RED),
				colval(cbuffer[n],GRN)*colval(cbuffer[n],GRN),
				colval(cbuffer[n],BLU)*colval(cbuffer[n],BLU));
			else
				setcolor(val2map[n], 0., 0., 0.);
	}
					/* compute conspicuity */
	if (!silent) {
		printf("\tComputing conspicuity map\n");
		fflush(stdout);
	}
	conspicuity();
	csp_cnt++;
#if 0
if (pass == 1) {
	char	fnm[256];
	sprintf(fnm, vval(BASENAME), fcur);
	strcat(fnm, "_incmap.pic");
	write_map(cerrmap, fnm);
}
#endif
					/* get ray start time */
	rtime_start = getTime();
	ctime_used += rtime_start - time_start;
	if (timed && rtime_start >= frm_stop)
		return(0);		/* error done but out of time */
	if (rtime_used <= FTINY) {
		if (quickstart)
			rays_todo = 1000;
		else
			rays_todo = hres*vres;
	} else {
		rays_todo = (long)((frm_stop - rtime_start) *
					ray_cnt / rtime_used);
		if (rays_todo < 1000)
			return(0);	/* let's call it a frame */
	}
					/* set higher rendering quality */
	if (twolevels && curparams != &hirendparams) {
		ray_restore(curparams = &hirendparams);
		if (nprocs > 1) {	/* need to update children */
			if (!silent) {
				printf("\tRestarting %d processes\n", nprocs);
				fflush(stdout);
			}
			ray_pclose(0);
			ray_popen(nprocs);
		}
	}
					/* compute refinement rays */
	if (!silent) {
		printf("\tRefinement pass %d...",
				pass+1, rays_todo);
		fflush(stdout);
	}
	if (asump != NULL)		/* flag low-quality samples */
		for (n = hres*vres; n--; )
			if (sbuffer[n])
				abuffer[n] = ALOWQ;
					/* trace those rays */
	nr = refine_rays(rays_todo);
	if (!silent)
		printf("traced %d HQ rays\n", nr);
	if (nr <= 0)
		return(0);
					/* update timing stats */
	while (ray_cnt >= 1L<<20) {
		ray_cnt >>= 1;
		rtime_used *= .5;
	}
	ray_cnt += nr;
	time_done = getTime();
	rtime_used += time_done - rtime_start;
	if (!timed && time_done > frm_stop)
		frm_stop = time_done;
					/* update ambient value */
	if (asump != NULL && asump->nsamps >= 1000) {
		double	sf = 1./(double)asump->nsamps;
		for (n = 3; n--; ) {
			asump->diffsum[n] *= sf;
			asump->diffsum[n] += colval(lorendparams.ambval,n);
			if (asump->diffsum[n] < 0) asump->diffsum[n] = 0;
		}
		setcolor(lorendparams.ambval,
				asump->diffsum[RED],
				asump->diffsum[GRN],
				asump->diffsum[BLU]);
		if (!silent)
			printf("\tUpdated parameter: -av %f %f %f\n",
					asump->diffsum[RED],
					asump->diffsum[GRN],
					asump->diffsum[BLU]);
		asump = NULL;
	}
	return(1);
nomore:
					/* make sure error map is updated */
	if ((fcur == fbeg | pass > 1))
		comp_frame_error();
	return(0);
}
