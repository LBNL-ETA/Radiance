/* RCSid: $Id: pfilt.h,v 2.4 2025/06/07 05:09:46 greg Exp $ */
/*
 * Header file for picture filtering
 */
#ifndef _RAD_PFILT_H_
#define _RAD_PFILT_H_

#include  "standard.h"
#include  "paths.h"
#include  "color.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double  CHECKRAD;	/* radius over which gaussian is summed */
extern double  rad;		/* output pixel radius for filtering */
extern double  thresh;		/* maximum contribution for subpixel */
extern int  nrows, ncols;	/* number of rows and columns for output */
extern int  xres, yres;		/* x and y resolution */
extern int  avghot;		/* true means average in avgbrt spots */
extern double  hotlvl;		/* brightness considered "hot" */
extern int  npts;		/* # of points for stars */
extern double  spread;		/* spread for star points */
extern char  *progname;
extern COLOR  exposure;		/* exposure for frame */

extern double  x_c, y_r;	/* conversion factors */

extern int  xrad;		/* x search radius */
extern int  yrad;		/* y search radius */
extern int  xbrad;		/* x box size */
extern int  ybrad;		/* y box size */

extern int  barsize;		/* size of input scan bar */
extern COLORV  **scanin;	/* input scan bar */
extern COLORV  *scanout;	/* output scan line */
extern COLORV  **scoutbar;	/* output scan bar (if thresh > 0) */
extern float  **greybar;	/* grey-averaged input values */
extern int  obarsize;		/* size of output scan bar */
extern int  orad;		/* output window radius */

extern int  wrapfilt;		/* wrap filter horizontally? */

typedef double brightfunc_t(SCOLOR c);
extern brightfunc_t *ourbright;	/* brightness calculation function */

	/* defined in pf2.c */
extern void pass1init(void);		/* prepare for first pass */
extern void pass1default(void);		/* for single pass */
extern void pass1scan(COLORV *scan, int  y);	/* process first pass scanline */
extern void pass2init(void);		/* prepare for final pass */
extern void pass2scan(COLORV *scan, int  y);	/* process final pass scanline */

	/* defined in pf3.c */
extern void initmask(void);			/* initialize gaussian lookup table */
extern void dobox(SCOLOR  csum, int  xcent, int  ycent,
		int  c, int  r);		/* simple box filter */
extern void dogauss(SCOLOR  csum, int  xcent, int  ycent,
		int  c, int  r);	/* gaussian filter */
extern void dothresh(int  xcent, int  ycent,
		int  ccent, int  rcent);	/* gaussian threshold filter */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PFILT_H_ */
