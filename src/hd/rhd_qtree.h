/* RCSid $Id: rhd_qtree.h,v 3.16 2003/07/14 22:24:00 schorsch Exp $ */
/*
 * Quadtree data structures for holodeck display drivers.
 */
#ifndef _RAD_RHD_QTREE_H_
#define _RAD_RHD_QTREE_H_

#include "color.h"
#include "tonemap.h"
#include "rhdriver.h"

#ifdef __cplusplus
extern "C" {
#endif

#define	DL		0		/* down left */
#define	DR		1		/* down right */
#define	UL		2		/* up left */
#define	UR		3		/* up right */

#define	BRF(i)		(0x1<<(i))	/* branch flag bit */
#define LFF(i)		(0x10<<(i))	/* leaf flag bit */
#define CHF(i)		(0x100<<(i))	/* change flag bit */
#define	CHBRF(i)	(0x101<<(i))	/* changed branch bit */
#define CHLFF(i)	(0x110<<(i))	/* changed leaf bit */
#define BR_ANY		0xf		/* flags for any branches */
#define LF_ANY		0xf0		/* flags for any leaves */
#define CH_ANY		0xf00		/* flags for any change */

typedef struct rtree {
	short	flgs;			/* content flags (defined above) */
	union {
		struct rtree	*b;		/* if branch */
		int	li;			/* if leaf */
	} k[4];				/* children */
}	RTREE;

extern struct rleaves {
	float		(*wp)[3];	/* world intersection point array */
	int32		*wd;		/* world direction array */
	TMbright	*brt;		/* encoded brightness array */
	BYTE		(*chr)[3];	/* encoded chrominance array */
	BYTE		(*rgb)[3];	/* tone-mapped color array */
	int		nl;		/* count of leaves in our pile */
	int		bl, tl;		/* bottom and top (next) leaf index */
	int		tml;		/* next leaf needing tone-mapping */
	char		*base;		/* base of allocated memory */
}	qtL;			/* our pile of leaves */

#define	is_stump(t)	(!((t)->flgs & (BR_ANY|LF_ANY)))

extern RTREE	qtrunk;		/* trunk of quadtree */
extern double	qtDepthEps;	/* epsilon to compare depths (z fraction) */
extern int	qtMinNodesiz;	/* minimum node dimension (pixels) */

extern int	rayqleft;	/* number of rays to queue before flush */

extern int32	encodedir();
extern double	fdir2diff(), dir2diff();

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RHD_QTREE_H_ */

