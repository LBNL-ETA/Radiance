/* RCSid: $Id: pmap.h,v 2.3 2003/07/14 22:24:00 schorsch Exp $ */
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

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PMAP_H_ */

