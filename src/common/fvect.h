/* RCSid $Id: fvect.h,v 2.9 2003/06/06 16:38:47 schorsch Exp $ */
/*
 * Declarations for floating-point vector operations.
 */
#ifndef _RAD_FVECT_H_
#define _RAD_FVECT_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#ifdef  SMLFLT
#define  FLOAT		float
#define  FTINY		(1e-3)
#else
#define  FLOAT		double
#define  FTINY		(1e-6)
#endif
#define  FHUGE		(1e10)

typedef FLOAT  FVECT[3];

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


extern double	fdot(FVECT v1, FVECT v2);
extern double	dist2(FVECT v1, FVECT v2);
extern double	dist2line(FVECT p, FVECT ep1, FVECT ep2);
extern double	dist2lseg(FVECT p, FVECT ep1, FVECT ep2);
extern void	fcross(FVECT vres, FVECT v1, FVECT v2);
extern void	fvsum(FVECT vres, FVECT v0, FVECT v1, double f);
extern double	normalize(FVECT v);
extern void	spinvector(FVECT vres, FVECT vorig, FVECT vnorm, double theta);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_FVECT_H_ */

