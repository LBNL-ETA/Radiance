/* RCSid: $Id: span.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *   Structures for line segment output to dot matrix printers
 */


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
