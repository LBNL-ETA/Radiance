/* RCSid: $Id: tgraph.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
#include  "meta.h"

#include  <math.h>


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

extern char  *snagquo(), *instr();

extern double  atof(), floor(), ceil(), log();
