/* RCSid $Id: plocate.h,v 2.4 2003/02/25 02:47:21 greg Exp $ */
/*
 *  plocate.h - header for 3D vector location.
 *
 *  Include after fvect.h
 */

#include "copyright.h"

#define  EPSILON	FTINY		/* acceptable location error */

#define  XPOS		03		/* x position mask */
#define  YPOS		014		/* y position mask */
#define  ZPOS		060		/* z position mask */

#define  position(i)	(3<<((i)<<1))	/* macro version */

#define  BELOW		025		/* below bits */
#define  ABOVE		052		/* above bits */

#ifdef NOPROTO

extern int	clip();

extern int	plocate();

#else

extern int	clip(FLOAT *ep1, FLOAT *ep2, FVECT min, FVECT max);

extern int	plocate(FVECT p, FVECT min, FVECT max);

#endif
