/* RCSid $Id: interp2d.h,v 2.2 2013/02/09 17:39:21 greg Exp $ */
/*
 * Header for interpolation of anisotropic values on 2-D plane.
 *
 *	G.Ward Feb 2013
 */

#define	NI2DSMF	0.42f			/* minimal smoothing radius */

#define NI2DIR	(2*4)			/* # interpolation directions */

/* Data structure for interpolant */
typedef struct {
	int		ns;		/* number of sample positions */
	float		rmin;		/* minimum valid radius (default=.5) */
	float		smf;		/* smoothing factor (def=NI2DSMF) */
	unsigned short	(*ra)[NI2DIR];	/* anisotropic radii (private) */
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
