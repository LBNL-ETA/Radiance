/* RCSid: $Id: plot.h,v 1.3 2003/07/14 22:24:00 schorsch Exp $ */
/*
 *   Definitions for plotting routines
 */
#ifndef _RAD_PLOT_H_
#define _RAD_PLOT_H_

#ifdef __cplusplus
extern "C" {
#endif

#define  NPATS  9		/* number of fill patterns */

#define  PATSIZE  16		/* pattern size (square) */


extern int  dxsize, dysize;	/* device size */

extern int  pati[];

extern unsigned char  pattern[][PATSIZE/8][PATSIZE];	/* fill patterns */

extern void set(int  attrib, char  *value);
extern void reset(int  attrib);
extern void unset(int  attrib);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PLOT_H_ */

