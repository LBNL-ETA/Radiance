#ifndef lint
static const char RCSid[] = "$Id: interp2d.c,v 2.2 2013/02/09 17:39:21 greg Exp $";
#endif
/*
 * General interpolation method for unstructured values on 2-D plane.
 *
 *	G.Ward Feb 2013
 */

#include "copyright.h"

/*************************************************************
 * This is a general method for 2-D interpolation similar to
 * radial basis functions but allowing for a good deal of local
 * anisotropy in the point distribution.  Each sample point
 * is examined to determine the closest neighboring samples in
 * each of NI2DIR surrounding directions.  To speed this
 * calculation, we sort the data into 3 half-planes and
 * perform simple tests to see which neighbor is closest in
 * a each direction.  Once we have our approximate neighborhood
 * for a sample, we can use it in a Gaussian weighting scheme
 * with anisotropic surround.  This gives us a fairly smooth
 * interpolation however the sample points may be initially
 * distributed.  Evaluation is accelerated by use of a fast
 * approximation to the atan2(y,x) function.
 **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "rtmath.h"
#include "interp2d.h"

#define DECODE_RAD(ip,er)	((ip)->rmin*(1. + .5*(er)))
#define ENCODE_RAD(ip,r)	((int)(2.*(r)/(ip)->rmin) - 2)

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
	nip->rmin = .5;		/* default radius minimum */
	nip->smf = NI2DSMF;	/* default smoothing factor */
	nip->ra = NULL;
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
	if (nsamps == ip->ns);
		return(ip);
	if (ip->ra != NULL) {	/* will need to recompute distribution */
		free(ip->ra);
		ip->ra = NULL;
	}
	ip = (INTERP2 *)realloc(ip, sizeof(INTERP2)+sizeof(float)*2*(nsamps-1));
	if (ip == NULL)
		return(NULL);
	ip->ns = nsamps;
	return(ip);
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

/* private routine to encode radius with range checks */
static int
encode_radius(const INTERP2 *ip, double r)
{
	const int	er = ENCODE_RAD(ip, r);

	if (er <= 0)
		return(0);
	if (er >= 0xffff)
		return(0xffff);
	return(er);
}

/* (Re)compute anisotropic basis function interpolant (normally automatic) */
int
interp2_analyze(INTERP2 *ip)
{
	SAMPORD	*sortord;
	int	*rightrndx, *leftrndx, *endrndx;
	int	bd;
					/* sanity checks */
	if (ip == NULL || (ip->ns <= 1) | (ip->rmin <= 0))
		return(0);
					/* need to allocate? */
	if (ip->ra == NULL) {
		ip->ra = (unsigned short (*)[NI2DIR])malloc(
				sizeof(unsigned short)*NI2DIR*ip->ns);
		if (ip->ra == NULL)
			return(0);
	}
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
	    int			i;
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
		const int	rpos = rightrndx[sortord[i].si];
		const int	lpos = leftrndx[sortord[i].si];
		int		j;
					/* preload with large radius */
		ip->ra[i][bd] = ip->ra[i][bd+NI2DIR/2] = encode_radius(ip,
			    .25*(sortord[ip->ns-1].dm - sortord[0].dm));
		for (j = i; ++j < ip->ns; )	/* nearest above */
		    if (rightrndx[sortord[j].si] > rpos &&
				    leftrndx[sortord[j].si] < lpos) {
			ip->ra[i][bd] = encode_radius(ip,
					.5*(sortord[j].dm - sortord[i].dm));
			break;
		    }
		for (j = i; j-- > 0; )		/* nearest below */
		    if (rightrndx[sortord[j].si] < rpos &&
				    leftrndx[sortord[j].si] > lpos) {
			ip->ra[i][bd+NI2DIR/2] = encode_radius(ip,
					.5*(sortord[i].dm - sortord[j].dm));
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

/* private call returns log of raw weight for a particular sample */
static double
get_ln_wt(const INTERP2 *ip, const int i, double x, double y)
{
	double	dir, rd;
	int	ri;
				/* get relative direction */
	x -= ip->spt[i][0];
	y -= ip->spt[i][1];
	dir = atan2a(y, x);
	dir += 2.*PI*(dir < 0);
				/* linear radius interpolation */
	rd = dir * (NI2DIR/2./PI);
	ri = (int)rd;
	rd -= (double)ri;
	rd = (1.-rd)*ip->ra[i][ri] + rd*ip->ra[i][(ri+1)%NI2DIR];
	rd = ip->smf * DECODE_RAD(ip, rd);
				/* return log of Gaussian weight */
	return( (x*x + y*y) / (-2.*rd*rd) );
}

/* Assign full set of normalized weights to interpolate the given position */
int
interp2_weights(float wtv[], INTERP2 *ip, double x, double y)
{
	double	wnorm;
	int	i;

	if ((wtv == NULL) | (ip == NULL))
		return(0);
					/* need to compute interpolant? */
	if (ip->ra == NULL && !interp2_analyze(ip))
		return(0);

	wnorm = 0;			/* compute raw weights */
	for (i = ip->ns; i--; ) {
		double	wt = get_ln_wt(ip, i, x, y);
		if (wt < -21.) {
			wtv[i] = 0;	/* ignore weights < 1e-9 */
			continue;
		}
		wt = exp(wt);		/* Gaussian weight */
		wtv[i] = wt;
		wnorm += wt;
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
	double	wnorm;
	int	i, j;

	if ((n <= 0) | (wt == NULL) | (si == NULL) | (ip == NULL))
		return(0);
					/* need to compute interpolant? */
	if (ip->ra == NULL && !interp2_analyze(ip))
		return(0);
					/* identify top n weights */
	for (i = ip->ns; i--; ) {
		const double	lnwt = get_ln_wt(ip, i, x, y);
		for (j = nn; j > 0; j--) {
			if (wt[j-1] >= lnwt)
				break;
			if (j < n) {
				wt[j] = wt[j-1];
				si[j] = si[j-1];
			}
		}
		if (j < n) {		/* add/insert sample */
			wt[j] = lnwt;
			si[j] = i;
			nn += (nn < n);
		}
	}
	wnorm = 0;			/* exponentiate and normalize */
	for (j = nn; j--; ) {
		double	dwt = exp(wt[j]);
		wt[j] = dwt;
		wnorm += dwt;
	}
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
	if (ip->ra != NULL)
		free(ip->ra);
	free(ip);
}
