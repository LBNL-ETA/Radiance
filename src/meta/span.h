/* RCSid: $Id: span.h,v 1.2 2003/07/14 22:24:00 schorsch Exp $ */
/*
 *   Structures for line segment output to dot matrix printers
 */
#ifndef _RAD_SPAN_H_
#define _RAD_SPAN_H_

#ifdef __cplusplus
extern "C" {
#endif


#define  MAXSPAN  5120			/* maximum span size in bytes */

struct span  {				/* line of printer output */
	      int  xleft, xright, ybot, ytop;
	      char  cols[MAXSPAN];
	      char  tcols[MAXSPAN];	
	      };




extern int  linwidt, linhite,		/* line size */
            nrows;			/* # of rows */

extern int  minwidth;			/* minimum line width */

extern struct span  outspan;		/* output span */

extern int  spanmin, spanmax;		/* current span dimensions */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_SPAN_H_ */

