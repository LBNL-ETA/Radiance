/* RCSid: $Id: pic.h,v 2.2 2003/02/22 02:07:27 greg Exp $ */
/* the following three structures are used by ciq */

typedef int colormap[3][256];

typedef unsigned char pixel;

typedef struct {
    pixel r,g,b;
} rgbpixel;

#define  rgb_bright(p)		(int)((77L*(p)->r+151L*(p)->g+28L*(p)->b)/256)

/* image resolution */
extern int	xmax,ymax;
