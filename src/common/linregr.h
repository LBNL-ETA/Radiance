/* RCSid $Id: linregr.h,v 2.6 2003/06/27 06:53:21 greg Exp $ */
/*
 * Header file for linear regression calculation.
 */
#ifndef _RAD_LINEGR_H_
#define _RAD_LINEGR_H_
#ifdef __cplusplus
extern "C" {
#endif

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


extern void	lrclear(LRSUM *l);
extern int	flrpoint(double x, double y, LRSUM *l);
extern int	lrfit(LRLIN *r, LRSUM *l);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_LINEGR_H_ */

