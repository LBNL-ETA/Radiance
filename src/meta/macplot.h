/* RCSid: $Id: macplot.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *   Definitions for MacIntosh plotting routines
 */

#undef  TRUE

#undef  FALSE

#include  <quickdraw.h>


#define  mapx(x)  CONV(x, dxsize)

#define  mapy(y)  (dysize - 1 - CONV(y, dysize))


extern int  pati[];

extern Pattern  macpat[];	/* fill patterns */

extern  dxsize, dysize;		/* plot dimensions */
