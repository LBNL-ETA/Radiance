/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef __GBASIC_H
#define __GBASIC_H

#ifdef __cplusplus
extern "C" {
#endif


#include <math.h>

#ifndef g3Float
	#define g3Float		double
	#define	USE_DOUBLE	1
	#define	GB_EPSILON	1e-6
#else
	#define	GB_EPSILON	1e-4
#endif

typedef	unsigned char	uchar;
	
#define DEG2RAD(x) ((x)*M_PI/180.0)
#define RAD2DEG(x) ((x)*180.0/M_PI)

#define	GB_SWAP(a,b) {g3Float temp=(a);(a)=(b);(b)=temp;}
#define gb_swap(a,b) GB_SWAP(a,b)

#define	gb_random()	(((g3Float) rand())/(RAND_MAX-1))

int		gb_epsorder(g3Float a,g3Float b);

int		gb_epseq(g3Float a,g3Float b);

int		gb_inrange(g3Float a,g3Float rb,g3Float re);

int		gb_signum(g3Float a);

g3Float	gb_max(g3Float a,g3Float b);
g3Float	gb_min(g3Float a,g3Float b);

/* Roots of quadratic formula
 Returns zero if roots are not real, 1 for single root and 2 otherwise
*/
int		gb_getroots(g3Float* r1,g3Float* r2,g3Float a,g3Float b,g3Float c);

#ifdef __cplusplus
}
#endif


#endif
