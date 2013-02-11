/* RCSid $Id: interp2d.h,v 2.4 2013/02/11 22:56:22 greg Exp $ */
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

/* Data structure for interpolant */
typedef struct {
	int		ns;		/* number of sample positions */
	float		dmin;		/* minimum diameter (default=1) */
	float		smf;		/* smoothing factor (def=NI2DSMF) */
	unsigned short	(*da)[NI2DIR];	/* anisotropic distances (private) */
	float		spt[1][2];	/* sample positions (extends struct) */
} INTERP2;

/* Allocate a new set of interpolation samples (caller assigns spt[] array) */
extern INTERP2	*interp2_alloc(int nsamps);

/* Resize interpolation array (caller must assign any new values) */
extern INTERP2	*interp2_realloc(INTERP2 *ip, int nsamps);

/* Assign full set of normalized weights to interpolate the given location */
extern int	interp2_weights(float wtv[], INTERP2 *ip, double x, double y);

/* Get normalized weights and indexes for n best samples in descending order */
extern int	interp2_topsamp(float wt[], int si[], const int n,
					INTERP2 *ip, double x, double y);
/* Free interpolant */
extern void	interp2_free(INTERP2 *ip);

/* (Re)compute anisotropic basis function interpolant (normally automatic) */
extern int	interp2_analyze(INTERP2 *ip);

/***************************************************************
 * Typical use is to allocate an INTERP2 struct and assign the
 * spt[] array with the ordered sample locations in x & y.  (The
 * actual orientation of the axes is not important so long as
 * the application is consistent.)  During interpolation, either
 * interp2_weights() is called to obtain a normalized weighting
 * vector for all the samples, or interp2_topsamp() is called to
 * get the most important N samples for the specified location.
 * The weights (and indexes in the case of interp2_topsamp)
 * are then used as coefficients for corresponding sample values
 * in a vector sum together that interpolates the function at that
 * location.  Between the initial allocation call and the first
 * weight evaluation, the dmin member may be adjusted to
 * the distance under which samples will start to merge.  If this
 * parameter is later changed, interp2_analyze() should be called
 * to recompute the interpolant.
 * The smf member sets the smoothing factor for interpolation.
 * The default setting of NI2DSMF produces near-optimal
 * interpolation when well-separated sample values are known
 * precisely.  If a greater degree of mixing is desired, this
 * paremter may be increased and it will affect the next call
 * to interp2_weights() or interp2_topsamp().
 **************************************************************/
 
#ifdef __cplusplus
}
#endif
#endif	/* !_RAD_INTERP2D_H_ */
