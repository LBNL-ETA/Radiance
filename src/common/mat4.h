/* RCSid $Id: mat4.h,v 2.5 2003/02/25 02:47:21 greg Exp $ */
/*
 * Definitions for 4x4 matrix operations
 */

#include "copyright.h"

#include  "fvect.h"

typedef FLOAT  MAT4[4][4];

#ifdef  BSD
#define  copymat4(m4a,m4b)	bcopy((char *)m4b,(char *)m4a,sizeof(MAT4))
#else
#define  copymat4(m4a,m4b)	(void)memcpy((char *)m4a,(char *)m4b,sizeof(MAT4))
#endif

#define  MAT4IDENT		{ {1.,0.,0.,0.}, {0.,1.,0.,0.}, \
				{0.,0.,1.,0.}, {0.,0.,0.,1.} }

extern MAT4  m4ident;

#define  setident4(m4)		copymat4(m4, m4ident)

#ifdef NOPROTO

extern void	multmat4();
extern void	multv3();
extern void	multp3();
extern int	invmat4();

#else

extern void	multmat4(MAT4 m4a, MAT4 m4b, MAT4 m4c);
extern void	multv3(FVECT v3a, FVECT v3b, MAT4 m4);
extern void	multp3(FVECT p3a, FVECT p3b, MAT4 m4);
extern int	invmat4(MAT4 inverse, MAT4 mat);

#endif
