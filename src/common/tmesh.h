/* RCSid $Id: tmesh.h,v 2.5 2003/06/26 00:58:09 schorsch Exp $ */
/*
 * Header file for triangle mesh routines using barycentric coordinates
 */
#ifndef _RAD_TMESH_H_
#define _RAD_TMESH_H_
#ifdef __cplusplus
extern "C" {
#endif


#define TCALNAME	"tmesh.cal"	/* the name of our auxiliary file */

typedef struct {
	int	ax;		/* major axis */
	RREAL	tm[2][3];	/* transformation */
} BARYCCM;

#ifndef COSTOL
#define COSTOL		0.999995	/* cosine of tolerance for smoothing */
#endif

				/* flat_tri() return values */
#define ISBENT		0		/* is not flat */
#define ISFLAT		1		/* is flat */
#define RVBENT		2		/* reversed and not flat */
#define RVFLAT		3		/* reversed and flat */
#define DEGEN		-1		/* degenerate (zero area) */


int		flat_tri(FVECT v1, FVECT v2, FVECT v3,
				FVECT n1, FVECT n2, FVECT n3);
int		comp_baryc(BARYCCM *bcm,  FVECT v1, FVECT v2, FVECT v3);
void		eval_baryc(RREAL wt[3], FVECT p, BARYCCM *bcm);
int		get_baryc(RREAL wt[3], FVECT p, FVECT v1, FVECT v2, FVECT v3);
void		put_baryc(BARYCCM *bcm, RREAL com[][3], int n);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_TMESH_H_ */

