/* RCSid: $Id: plot.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *   Definitions for plotting routines
 */


#define  NPATS  9		/* number of fill patterns */

#define  PATSIZE  16		/* pattern size (square) */


extern int  dxsize, dysize;	/* device size */

extern int  pati[];

extern unsigned char  pattern[][PATSIZE/8][PATSIZE];	/* fill patterns */
