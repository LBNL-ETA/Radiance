/* RCSid $Id: fvect.h,v 2.15 2013/06/29 21:03:44 greg Exp $ */
/*
 * Declarations for floating-point vector operations.
 */
#ifndef _RAD_FVECT_H_
#define _RAD_FVECT_H_
#ifdef __cplusplus
extern "C" {
#endif

#ifdef  SMLFLT
#define  RREAL		float
#define  FTINY		(1e-3)
#else
#define  RREAL		double
#define  FTINY		(1e-6)
#endif
#define  FHUGE		(1e10)

typedef RREAL  FVECT[3];

#define  VCOPY(v1,v2)	((v1)[0]=(v2)[0],(v1)[1]=(v2)[1],(v1)[2]=(v2)[2])
#define  DOT(v1,v2)	((v1)[0]*(v2)[0]+(v1)[1]*(v2)[1]+(v1)[2]*(v2)[2])
#define  VLEN(v)	sqrt(DOT(v,v))
#define  VADD(vr,v1,v2)	((vr)[0]=(v1)[0]+(v2)[0], \
				(vr)[1]=(v1)[1]+(v2)[1], \
				(vr)[2]=(v1)[2]+(v2)[2])
#define  VSUB(vr,v1,v2)	((vr)[0]=(v1)[0]-(v2)[0], \
				(vr)[1]=(v1)[1]-(v2)[1], \
				(vr)[2]=(v1)[2]-(v2)[2])
#define  VSUM(vr,v1,v2,f)	((vr)[0]=(v1)[0]+(f)*(v2)[0], \
				(vr)[1]=(v1)[1]+(f)*(v2)[1], \
				(vr)[2]=(v1)[2]+(f)*(v2)[2])
#define  VCROSS(vr,v1,v2) \
			((vr)[0]=(v1)[1]*(v2)[2]-(v1)[2]*(v2)[1], \
			(vr)[1]=(v1)[2]*(v2)[0]-(v1)[0]*(v2)[2], \
			(vr)[2]=(v1)[0]*(v2)[1]-(v1)[1]*(v2)[0])

#define	GEOD_RAD	0	/* geodesic distance specified in radians */
#define	GEOD_ABS	1	/* absolute geodesic distance */
#define	GEOD_REL	2	/* relative geodesic distance */

extern double	Acos(double x);
extern double	Asin(double x);
extern double	fdot(const FVECT v1, const FVECT v2);
extern double	dist2(const FVECT v1, const FVECT v2);
extern double	dist2line(const FVECT p, const FVECT ep1, const FVECT ep2);
extern double	dist2lseg(const FVECT p, const FVECT ep1, const FVECT ep2);
extern void	fcross(FVECT vres, const FVECT v1, const FVECT v2);
extern void	fvsum(FVECT vres, const FVECT v0, const FVECT v1, double f);
extern double	normalize(FVECT v);
extern int	closestapproach(RREAL t[2],
			const FVECT rorg0, const FVECT rdir0,
			const FVECT rorg1, const FVECT rdir1);
extern void	spinvector(FVECT vres, const FVECT vorig,
			const FVECT vnorm, double theta);
extern double	geodesic(FVECT vres, const FVECT vorig,
			const FVECT vtarg, double t, int meas);
#ifdef __cplusplus
}
#endif
#endif /* _RAD_FVECT_H_ */

