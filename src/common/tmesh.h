/* RCSid $Id: tmesh.h,v 2.2 2003/03/11 17:08:55 greg Exp $ */
/*
 * Header file for triangle mesh routines using barycentric coordinates
 */

#define TCALNAME	"tmesh.cal"	/* the name of our auxiliary file */

typedef struct {
	int	ax;		/* major axis */
	FLOAT	tm[2][3];	/* transformation */
} BARYCCM;

#ifndef COSTOL
#define COSTOL		0.99985		/* cosine of tolerance for smoothing */
#endif

				/* flat_tri() return values */
#define ISBENT		0		/* is not flat */
#define ISFLAT		1		/* is flat */
#define RVBENT		2		/* reversed and not flat */
#define RVFLAT		3		/* reversed and flat */
#define DEGEN		-1		/* degenerate (zero area) */

#ifdef NOPROTO

int		flat_tri();
int		comp_baryc();
void		eval_baryc();
int		get_baryc();
void		put_baryc();

#else

int		flat_tri(FVECT v1, FVECT v2, FVECT v3,
				FVECT n1, FVECT n2, FVECT n3);
int		comp_baryc(BARYCCM *bcm,  FVECT v1, FVECT v2, FVECT v3);
void		eval_baryc(FLOAT wt[3], FVECT p, BARYCCM *bcm);
int		get_baryc(FLOAT wt[3], FVECT p, FVECT v1, FVECT v2, FVECT v3);
void		put_baryc(BARYCCM *bcm, FLOAT com[][3], int n);

#endif /* NOPROTO */
