/* RCSid $Id: plocate.h,v 2.5 2003/06/06 16:38:47 schorsch Exp $ */
/*
 *  plocate.h - header for 3D vector location.
 *
 *  Include after fvect.h
 */
#ifndef _RAD_PLOCATE_H_
#define _RAD_PLOCATE_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#define  EPSILON	FTINY		/* acceptable location error */

#define  XPOS		03		/* x position mask */
#define  YPOS		014		/* y position mask */
#define  ZPOS		060		/* z position mask */

#define  position(i)	(3<<((i)<<1))	/* macro version */

#define  BELOW		025		/* below bits */
#define  ABOVE		052		/* above bits */


extern int	clip(FLOAT *ep1, FLOAT *ep2, FVECT min, FVECT max);
extern int	plocate(FVECT p, FVECT min, FVECT max);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PLOCATE_H_ */

