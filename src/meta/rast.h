/* RCSid: $Id: rast.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *   Structures for line segment output to raster files
 */

#define  NUMSCANS	16		/* number of scanlines per block */

typedef struct {			/* raster scanline block */
	      unsigned char  *cols[NUMSCANS];
	      int  ybot, ytop;			/* ybot is scan[0] */
	      int  xleft, xright;
	      int  width;
	      } SCANBLOCK;

extern int  ydown;			/* y going down? */

extern int  minwidth;			/* minimum line width */

extern SCANBLOCK  outblock;		/* output span */

#define  IBLK	0		/* index for black */
#define  IRED	1
#define  IGRN	2
#define  IBLU	3
#define  IYEL	4
#define  IMAG	5
#define  ICYN	6
#define  IWHT	7

#define  pixtog(x,y,c)	{ register unsigned char \
		*cp = outblock.cols[(y)-outblock.ybot]+(x); \
		*cp = ((~*cp ^ (c)<<3)&070) | (*cp&07); }

#define  pixmix(x,y,c)	(outblock.cols[(y)-outblock.ybot][x] &= 070|(c))

#define  someabove(p,y)	(CONV((p)->xy[YMX],dysize) > (y))
#define  somebelow(p,y)	(CONV((p)->xy[YMN],dysize) < (y))

#define  inthis(p)	(ydown ? someabove(p,outblock.ybot-1) : \
				somebelow(p,outblock.ytop+1))
#define  innext(p)	(ydown ? somebelow(p,outblock.ybot) : \
				someabove(p,outblock.ytop))
