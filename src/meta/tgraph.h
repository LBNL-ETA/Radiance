/* RCSid: $Id: tgraph.h,v 1.4 2003/11/15 02:13:37 schorsch Exp $ */

#ifndef _RAD_TGRAPH_H_
#define _RAD_TGRAPH_H_

#include  <stdlib.h>
#include  <math.h>

#include  "meta.h"

#ifdef __cplusplus
extern "C" {
#endif

#define  MINDIVS  4		/* Minimum number of divisions for axis */

#define  FTINY  1e-10		/* tiny float */

#define  FHUGE  1e10		/* large float */

#define  PI  3.14159265358979323846

#define  XBEG  2048		/* starting x coordinate */

#define  XSIZ  10240		/* x axis coordinate size */

#define  YBEG  2048		/* starting y coordinate */

#define  YSIZ  10240		/* y axis coordinate size */

#define  XCONV(x)  ((int)((x-xmin)*XSIZ/xsize+XBEG))

#define  YCONV(y)  ((int)((y-ymin)*YSIZ/ysize+YBEG))

#define  TSIZ  200		/* tick size */

#define  SYMRAD  100		/* default symbol radius */

#define  XTICS  1		/* flags for making axes */
#define  YTICS  2
#define  XNUMS  4
#define  YNUMS  8
#define  XGRID  16
#define  YGRID  32
#define  ORIGIN  64
#define  BOX  128

#define  RADIANS  1		/* flags for polar coordinates */
#define  DEGREES  2

#define  NCUR  16		/* number of supported curve types */
 


extern short  usecurve[];		/* booleans for curve usage */

extern int  symrad;			/* current symbol radius */

extern double  xmin, xmax,		/* extrema */
               ymin, ymax;

extern double  xsize, ysize;		/* coordinate dimensions */

extern double  xmnset, xmxset,		/* domain settings */
               ymnset, ymxset;
       
extern short  logx, logy;		/* flags for logarithmic axes */

extern short  grid;			/* flag for grid */

extern short  polar;			/* flag for polar coordiates */

extern int    ncurves;			/* number of curves in file */

extern char *instr(char  *s, char  *t);
extern char *snagquo(register char  *s);
extern void normalize(FILE  *fp, FILE *fout);
extern void initialize(void);
extern void option(char *s);
extern int istitle(char *s);
extern int isdata(char *s);
extern int islabel(char *s);
extern int isdivlab(char *s);
extern int isxlabel(char *s);
extern int isylabel(char *s);
extern void makeaxis(int flag);
extern int getdata(char  *s, double  *xp, double  *yp);
extern void boxstring(int a0, int xmn, int ymn, int xmx, int ymx, char  *s);
extern void symout(int  a0, int  x, int  y, char  *sname);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_TGRAPH_H_ */

