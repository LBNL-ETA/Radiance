/* RCSid $Id: mat4.h,v 2.10 2003/06/27 06:53:21 greg Exp $ */
/*
 * Definitions for 4x4 matrix operations
 */
#ifndef _RAD_MAT4_H_
#define _RAD_MAT4_H_
#ifdef __cplusplus
extern "C" {
#endif

#include  <string.h>

#include  "fvect.h"

typedef RREAL  MAT4[4][4];

#ifdef  BSD
#define  copymat4(m4a,m4b)	bcopy((void *)m4b,(void *)m4a,sizeof(MAT4))
#else
#define  copymat4(m4a,m4b)	(void)memcpy((void *)m4a,(void *)m4b,sizeof(MAT4))
#endif

#define  MAT4IDENT		{ {1.,0.,0.,0.}, {0.,1.,0.,0.}, \
				{0.,0.,1.,0.}, {0.,0.,0.,1.} }

extern MAT4  m4ident;

#define  setident4(m4)		copymat4(m4, m4ident)


extern void	multmat4(MAT4 m4a, MAT4 m4b, MAT4 m4c);
extern void	multv3(FVECT v3a, FVECT v3b, MAT4 m4);
extern void	multp3(FVECT p3a, FVECT p3b, MAT4 m4);
extern int	invmat4(MAT4 inverse, MAT4 mat);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_MAT4_H_ */

