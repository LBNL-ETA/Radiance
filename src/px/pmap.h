/* RCSid: $Id: pmap.h,v 2.4 2004/03/28 20:33:14 schorsch Exp $ */
/* Pmap return codes */
#ifndef _RAD_PMAP_H_
#define _RAD_PMAP_H_

#ifdef __cplusplus
extern "C" {
#endif

#define PMAP_BAD	-1
#define PMAP_LINEAR	0
#define PMAP_PERSP	1

/*  |a b|
 *  |c d|
 */
#define DET2(a,b, c,d) ((a)*(d) - (b)*(c))


	/* defined in pmapgen.c */
extern int pmap_quad_rect(double u0, double v0, double u1, double v1,
		double qdrl[4][2], double QR[3][3]);
extern int pmap_square_quad(double qdrl[4][2], double SQ[3][3]);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PMAP_H_ */

