/* RCSid $Id: linregr.h,v 2.4 2003/02/25 02:47:21 greg Exp $ */
/*
 * Header file for linear regression calculation.
 */

#include "copyright.h"

typedef struct {
	double	xs, ys, xxs, yys, xys;
	int	n;
} LRSUM;

typedef struct {
	double	slope, intercept, correlation;
} LRLIN;

#define	lrpoint(x,y,l)	((l)->xs+=(x),(l)->ys+=(y),(l)->xxs+=(x)*(x), \
			(l)->yys+=(y)*(y),(l)->xys+=(x)*(y),++(l)->n)

#define	lrxavg(l)	((l)->xs/(l)->n)
#define	lryavg(l)	((l)->ys/(l)->n)
#define	lrxvar(l)	(((l)->xxs-(l)->xs*(l)->xs/(l)->n)/(l)->n)
#define	lryvar(l)	(((l)->yys-(l)->ys*(l)->ys/(l)->n)/(l)->n)
#define	lrxdev(l)	sqrt(((l)->xxs-(l)->xs*(l)->xs/(l)->n)/((l)->n-1))
#define	lrydev(l)	sqrt(((l)->yys-(l)->ys*(l)->ys/(l)->n)/((l)->n-1))

#ifdef NOPROTO

extern void	lrclear();
extern int	flrpoint();
extern int	lrfit();

#else

extern void	lrclear(LRSUM *l);
extern int	flrpoint(double x, double y, LRSUM *l);
extern int	lrfit(LRLIN *r, LRSUM *l);

#endif
