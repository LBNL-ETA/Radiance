/* RCSid: $Id: macplot.h,v 1.2 2003/07/14 22:24:00 schorsch Exp $ */
/*
 *   Definitions for MacIntosh plotting routines
 */
#ifndef _RAD_MACPLOT_H_
#define _RAD_MACPLOT_H_

#undef  TRUE
#undef  FALSE

#include  <quickdraw.h>

#ifdef __cplusplus
extern "C" {
#endif

#define  mapx(x)  CONV(x, dxsize)

#define  mapy(y)  (dysize - 1 - CONV(y, dysize))


extern int  pati[];

extern Pattern  macpat[];	/* fill patterns */

extern  dxsize, dysize;		/* plot dimensions */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_MACPLOT_H_ */

