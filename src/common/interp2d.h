/* RCSid $Id: interp2d.h,v 2.9 2013/02/15 19:15:16 greg Exp $ */
/*
 * Header for interpolation of anisotropic samples on 2-D plane.
 *
 *	G.Ward Feb 2013
 */

#ifndef _RAD_INTERP2D_H_
#define _RAD_INTERP2D_H_

#ifdef __cplusplus
extern "C" {
#endif

#define	NI2DSMF	0.42f			/* minimal smoothing factor */

#define NI2DIR	(2*4)			/* # interpolation directions */
#define NI2DIM	16			/* size of influence map */

/* Data structure for interpolant */
typedef struct {
	int		ns;		/* number of sample positions */
	float		dmin;		/* minimum diameter (default=1) */
	float		smf;		/* smoothing factor (def=NI2DSMF) */
	float		smin[2];	/* sample minima */
	float		smax[2];	/* sample maxima */
	float		grid2;		/* grid diameter squared */
	struct interp2_samp {
		unsigned short	dia[NI2DIR];
		unsigned short	infl[NI2DIM];
	}		*da;		/* direction array (private) */
	float		spt[1][2];	/* sample positions (extends struct) */
} INTERP2;

/* Allocate a new set of interpolation samples (caller assigns spt[] array) */
extern INTERP2	*interp2_alloc(int nsamps);

/* Resize interpolation array (caller must assign any new values) */
extern INTERP2	*interp2_realloc(INTERP2 *ip, int nsamps);

/* Set minimum distance under which samples will start to merge */
extern void	interp2_spacing(INTERP2 *ip, double mind);

/* Modify smoothing parameter by the given factor */
extern void	interp2_smooth(INTERP2 *ip, double sf);

/* Assign full set of normalized weights to interpolate the given location */
extern int	interp2_weights(float wtv[], INTERP2 *ip, double x, double y);

/* Get normalized weights and indexes for n best samples in descending order */
extern int	interp2_topsamp(float wt[], int si[], const int n,
					INTERP2 *ip, double x, double y);
/* Free interpolant */
extern void	interp2_free(INTERP2 *ip);

/* (Re)compute anisotropic basis function interpolant (normally automatic) */
extern int	interp2_analyze(INTERP2 *ip);

/* Compute unnormalized weight for a position relative to a sample */
double		interp2_wti(INTERP2 *ip, const int i, double x, double y);

/***************************************************************
 * Typical use is to allocate an INTERP2 struct and assign the
 * spt[] array with the ordered sample locations in x & y.  (The
 * actual orientation of the axes is not important so long as
 * the application is consistent.)  During interpolation, either
 * interp2_weights() is called to obtain a normalized weighting
 * vector for all the samples, or interp2_topsamp() is called to
 * get the most important N samples for the specified location.
 * The weights (and indexes in the case of interp2_topsamp)
 * are then used as coefficients for corresponding sample
 * values in a vector sum that interpolates the function at
 * that location.
 * The minimum distance between sample positions defaults to 1.0.
 * Values spaced closer than this will be merged/averaged.
 * The interp2_spacing() call may be used to alter this distance,
 * causing the interpolant to be recalculated during the
 * next call to the sampling functions.
 * The default smoothing factor NI2DSMF provides near-optimal
 * interpolation when well-separated values are known
 * precisely.  Increase this setting by a factor > 1
 * with the interp2_smooth() call if greater mixing is desired.
 * A call of interp2_smooth(ip,0) resets to the minimum
 * default.  It is not possible to "sharpen" the data.
 **************************************************************/
 
#ifdef __cplusplus
}
#endif
#endif	/* !_RAD_INTERP2D_H_ */
