/* RCSid: $Id: tmesh.h,v 2.3 2003/02/22 02:07:23 greg Exp $ */
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
