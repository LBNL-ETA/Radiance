#ifndef lint
static const char RCSid[] = "$Id: interp2d.c,v 2.11 2013/02/15 01:26:47 greg Exp $";
#endif
/*
 * General interpolation method for unstructured values on 2-D plane.
 *
 *	G.Ward Feb 2013
 */

#include "copyright.h"

/***************************************************************
 * This is a general method for 2-D interpolation similar to
 * radial basis functions but allowing for a good deal of local
 * anisotropy in the point distribution.  Each sample point
 * is examined to determine the closest neighboring samples in
 * each of NI2DIR surrounding directions.  To speed this
 * calculation, we sort the data into half-planes and apply
 * simple tests to see which neighbor is closest in each
 * angular slice.  Once we have our approximate neighborhood
 * for a sample, we can use it in a modified Gaussian weighting
 * with allowing local anisotropy.  Harmonic weighting is added
 * to reduce the influence of distant neighbors.  This yields a
 * smooth interpolation regardless of how the sample points are
 * initially distributed.  Evaluation is accelerated by use of
 * a fast approximation to the atan2(y,x) function and an array
 * of flags indicating where weights are (nearly) zero.
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "rtmath.h"
#include "interp2d.h"

#define DECODE_DIA(ip,ed)	((ip)->dmin*(1. + .5*(ed)))
#define ENCODE_DIA(ip,d)	((int)(2.*(d)/(ip)->dmin) - 2)

/* Sample order (private) */
typedef struct {
	int	si;		/* sample index */
	float	dm;		/* distance measure in this direction */
} SAMPORD;

/* Allocate a new set of interpolation samples (caller assigns spt[] array) */
INTERP2	*
interp2_alloc(int nsamps)
{
	INTERP2		*nip;

	if (nsamps <= 1)
		return(NULL);

	nip = (INTERP2 *)malloc(sizeof(INTERP2) + sizeof(float)*2*(nsamps-1));
	if (nip == NULL)
		return(NULL);

	nip->ns = nsamps;
	nip->dmin = 1;		/* default minimum diameter */
	nip->smf = NI2DSMF;	/* default smoothing factor */
	nip->da = NULL;
				/* caller must assign spt[] array */
	return(nip);
}

/* Resize interpolation array (caller must assign any new values) */
INTERP2	*
interp2_realloc(INTERP2 *ip, int nsamps)
{
	if (ip == NULL)
		return(interp2_alloc(nsamps));
	if (nsamps <= 1) {
		interp2_free(ip);
		return(NULL);
	}
	if (nsamps == ip->ns)
		return(ip);
	if (ip->da != NULL) {	/* will need to recompute distribution */
		free(ip->da);
		ip->da = NULL;
	}
	ip = (INTERP2 *)realloc(ip, sizeof(INTERP2)+sizeof(float)*2*(nsamps-1));
	if (ip == NULL)
		return(NULL);
	ip->ns = nsamps;
	return(ip);
}

/* Set minimum distance under which samples will start to merge */
void
interp2_spacing(INTERP2 *ip, double mind)
{
	if (mind <= 0)
		return;
	if ((.998*ip->dmin <= mind) & (mind <= 1.002*ip->dmin))
		return;
	if (ip->da != NULL) {	/* will need to recompute distribution */
		free(ip->da);
		ip->da = NULL;
	}
	ip->dmin = mind;
}

/* Modify smoothing parameter by the given factor */
void
interp2_smooth(INTERP2 *ip, double sf)
{
	if ((ip->smf *= sf) < NI2DSMF)
		ip->smf = NI2DSMF;
}

/* private call-back to sort position index */
static int
cmp_spos(const void *p1, const void *p2)
{
	const SAMPORD	*so1 = (const SAMPORD *)p1;
	const SAMPORD	*so2 = (const SAMPORD *)p2;

	if (so1->dm > so2->dm)
		return 1;
	if (so1->dm < so2->dm)
		return -1;
	return 0;
}

/* private routine to order samples in a particular direction */
static void
sort_samples(SAMPORD *sord, const INTERP2 *ip, double ang)
{
	const double	cosd = cos(ang);
	const double	sind = sin(ang);
	int		i;

	for (i = ip->ns; i--; ) {
		sord[i].si = i;
		sord[i].dm = cosd*ip->spt[i][0] + sind*ip->spt[i][1];
	}
	qsort(sord, ip->ns, sizeof(SAMPORD), &cmp_spos);
}

/* private routine to encode sample diameter with range checks */
static int
encode_diameter(const INTERP2 *ip, double d)
{
	const int	ed = ENCODE_DIA(ip, d);

	if (ed <= 0)
		return(0);
	if (ed >= 0xffff)
		return(0xffff);
	return(ed);
}

/* (Re)compute anisotropic basis function interpolant (normally automatic) */
int
interp2_analyze(INTERP2 *ip)
{
	SAMPORD	*sortord;
	int	*rightrndx, *leftrndx, *endrndx;
	int	i, bd;
					/* sanity checks */
	if (ip == NULL)
		return(0);
	if (ip->da != NULL) {		/* free previous data if any */
		free(ip->da);
		ip->da = NULL;
	}
	if ((ip->ns <= 1) | (ip->dmin <= 0))
		return(0);
					/* compute sample domain */
	ip->smin[0] = ip->smin[1] = FHUGE;
	ip->smax[0] = ip->smax[1] = -FHUGE;
	for (i = ip->ns; i--; ) {
		if (ip->spt[i][0] < ip->smin[0])
			ip->smin[0] = ip->spt[i][0];
		if (ip->spt[i][0] > ip->smax[0])
			ip->smax[0] = ip->spt[i][0];
		if (ip->spt[i][1] < ip->smin[1])
			ip->smin[1] = ip->spt[i][1];
		if (ip->spt[i][1] > ip->smax[1])
			ip->smax[1] = ip->spt[i][1];
	}
	ip->grid2 = ((ip->smax[0]-ip->smin[0])*(ip->smax[0]-ip->smin[0]) +
			(ip->smax[1]-ip->smin[1])*(ip->smax[1]-ip->smin[1])) *
				(1./NI2DIM/NI2DIM);
	if (ip->grid2 <= FTINY*ip->dmin*ip->dmin)
		return(0);
					/* allocate analysis data */
	ip->da = (struct interp2_samp *)calloc( ip->ns,
					sizeof(struct interp2_samp) );
	if (ip->da == NULL)
		return(0);
					/* get temporary arrays */
	sortord = (SAMPORD *)malloc(sizeof(SAMPORD)*ip->ns);
	rightrndx = (int *)malloc(sizeof(int)*ip->ns);
	leftrndx = (int *)malloc(sizeof(int)*ip->ns);
	endrndx = (int *)malloc(sizeof(int)*ip->ns);
	if ((sortord == NULL) | (rightrndx == NULL) |
			(leftrndx == NULL) | (endrndx == NULL))
		return(0);
					/* run through bidirections */
	for (bd = 0; bd < NI2DIR/2; bd++) {
	    const double	ang = 2.*PI/NI2DIR*bd;
	    int			*sptr;
					/* create right reverse index */
	    if (bd) {			/* re-use from previous iteration? */
		sptr = rightrndx;
		rightrndx = leftrndx;
		leftrndx = sptr;
	    } else {			/* else sort first half-plane */
		sort_samples(sortord, ip, PI/2. - PI/NI2DIR);
		for (i = ip->ns; i--; )
		    rightrndx[sortord[i].si] = i;
					/* & store reverse order for later */
		for (i = ip->ns; i--; )
		    endrndx[sortord[i].si] = ip->ns-1 - i;
	    }
					/* create new left reverse index */
	    if (bd == NI2DIR/2 - 1) {	/* use order from first iteration? */
		sptr = leftrndx;
		leftrndx = endrndx;
		endrndx = sptr;
	    } else {			/* else compute new half-plane */
		sort_samples(sortord, ip, ang + (PI/2. + PI/NI2DIR));
		for (i = ip->ns; i--; )
		    leftrndx[sortord[i].si] = i;
	    }
					/* sort grid values in this direction */
	    sort_samples(sortord, ip, ang);
					/* find nearest neighbors each side */
	    for (i = ip->ns; i--; ) {
		const int	ii = sortord[i].si;
		int		j;
					/* preload with large radii */
		ip->da[ii].dia[bd] =
		ip->da[ii].dia[bd+NI2DIR/2] = encode_diameter(ip,
				.5*(sortord[ip->ns-1].dm - sortord[0].dm));
		for (j = i; ++j < ip->ns; )	/* nearest above */
		    if (rightrndx[sortord[j].si] > rightrndx[ii] &&
				    leftrndx[sortord[j].si] < leftrndx[ii]) {
			ip->da[ii].dia[bd] = encode_diameter(ip,
						sortord[j].dm - sortord[i].dm);
			break;
		    }
		for (j = i; j-- > 0; )		/* nearest below */
		    if (rightrndx[sortord[j].si] < rightrndx[ii] &&
				    leftrndx[sortord[j].si] > leftrndx[ii]) {
			ip->da[ii].dia[bd+NI2DIR/2] = encode_diameter(ip,
						sortord[i].dm - sortord[j].dm);
			break;
		    }
	    }
	}
	free(sortord);			/* clean up */
	free(rightrndx);
	free(leftrndx);
	free(endrndx);
	return(1);
}

/* Compute unnormalized weight for a position relative to a sample */
double
interp2_wti(INTERP2 *ip, const int i, double x, double y)
{
	double		dir, rd, r2, d2;
	int		ri;
				/* get relative direction */
	x -= ip->spt[i][0];
	y -= ip->spt[i][1];
	dir = atan2a(y, x);
	dir += 2.*PI*(dir < 0);
				/* linear radius interpolation */
	rd = dir * (NI2DIR/2./PI);
	ri = (int)rd;
	rd -= (double)ri;
	rd = (1.-rd)*ip->da[i].dia[ri] + rd*ip->da[i].dia[(ri+1)%NI2DIR];
	rd = ip->smf * DECODE_DIA(ip, rd);
	r2 = 2.*rd*rd;
	d2 = x*x + y*y;
	if (d2 > 21.*r2)	/* result would be < 1e-9 */
		return(.0);
				/* Gaussian times harmonic weighting */
	return( exp(-d2/r2) * ip->dmin/(ip->dmin + sqrt(d2)) );
}

/* private call to get grid flag index */
static int
interp2_flagpos(int fgi[2], INTERP2 *ip, double x, double y)
{
	int	ingrid = 1;

	if (ip == NULL)		/* paranoia */
		return(-1);
				/* need to compute interpolant? */
	if (ip->da == NULL && !interp2_analyze(ip))
		return(-1);
				/* get grid position */
	fgi[0] = (x - ip->smin[0]) * NI2DIM / (ip->smax[0] - ip->smin[0]);
	if (fgi[0] >= NI2DIM) {
		fgi[0] = NI2DIM-1;
		ingrid = 0;
	} else if (fgi[0] < 0) {
		fgi[0] = 0;
		ingrid = 0;
	}
	fgi[1] = (y - ip->smin[1]) * NI2DIM / (ip->smax[1] - ip->smin[1]);
	if (fgi[1] >= NI2DIM) {
		fgi[1] = NI2DIM-1;
		ingrid = 0;
	} else if (fgi[1] < 0) {
		fgi[1] = 0;
		ingrid = 0;
	}
	return(ingrid);
}

/* private call to set black flag if not too close to the given sample */
static void
setblk(INTERP2 *ip, const int i, const int gi[2])
{
	double	dx = (gi[0]+.5)*(1./NI2DIM)*(ip->smax[0] - ip->smin[0]) +
					ip->smin[0] - ip->spt[i][0];
	double	dy = (gi[1]+.5)*(1./NI2DIM)*(ip->smax[1] - ip->smin[1]) +
					ip->smin[1] - ip->spt[i][1];

	if (dx*dx + dy*dy > 2.*ip->grid2)
		ip->da[i].blkflg[gi[1]] |= 1<<gi[0];
}

#define chkblk(ip,i,gi)	((ip)->da[i].blkflg[(gi)[1]]>>(gi)[0] & 1)

/* Assign full set of normalized weights to interpolate the given position */
int
interp2_weights(float wtv[], INTERP2 *ip, double x, double y)
{
	double	wnorm;
	int	fgi[2];
	int	ingrid;
	int	i;

	if (wtv == NULL)
		return(0);
					/* get flag position */
	if ((ingrid = interp2_flagpos(fgi, ip, x, y)) < 0)
		return(0);

	wnorm = 0;			/* compute raw weights */
	for (i = ip->ns; i--; )
	    if (chkblk(ip, i, fgi)) {
		wtv[i] = 0;
	    } else {
		double	wt = interp2_wti(ip, i, x, y);
		wtv[i] = wt;
		wnorm += wt;
		if (wt <= 1e-9 && ingrid)
			setblk(ip, i, fgi);
	    }
	if (wnorm <= 0)			/* too far from all our samples! */
		return(0);
	wnorm = 1./wnorm;		/* normalize weights */
	for (i = ip->ns; i--; )
		wtv[i] *= wnorm;
	return(ip->ns);			/* all done */
}


/* Get normalized weights and indexes for n best samples in descending order */
int
interp2_topsamp(float wt[], int si[], const int n, INTERP2 *ip, double x, double y)
{
	int	nn = 0;
	int	fgi[2];
	int	ingrid;
	double	wnorm;
	int	i, j;

	if ((n <= 0) | (wt == NULL) | (si == NULL))
		return(0);
					/* get flag position */
	if ((ingrid = interp2_flagpos(fgi, ip, x, y)) < 0)
		return(0);
					/* identify top n weights */
	for (i = ip->ns; i--; )
	    if (!chkblk(ip, i, fgi)) {
		const double	wti = interp2_wti(ip, i, x, y);
		if (wti <= 1e-9) {
			if (ingrid)
				setblk(ip, i, fgi);
			continue;
		}
		for (j = nn; j > 0; j--) {
			if (wt[j-1] >= wti)
				break;
			if (j < n) {
				wt[j] = wt[j-1];
				si[j] = si[j-1];
			}
		}
		if (j < n) {		/* add/insert sample */
			wt[j] = wti;
			si[j] = i;
			nn += (nn < n);
		}
	    }
	wnorm = 0;			/* normalize sample weights */
	for (j = nn; j--; )
		wnorm += wt[j];
	if (wnorm <= 0)
		return(0);
	wnorm = 1./wnorm;
	for (j = nn; j--; )
		wt[j] *= wnorm;
	return(nn);			/* return actual # samples */
}

/* Free interpolant */
void
interp2_free(INTERP2 *ip)
{
	if (ip == NULL)
		return;
	if (ip->da != NULL)
		free(ip->da);
	free(ip);
}
