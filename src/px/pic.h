/* RCSid: $Id: pic.h,v 2.3 2003/07/14 22:24:00 schorsch Exp $ */
/* the following three structures are used by ciq */

#ifndef _RAD_PIC_H_
#define _RAD_PIC_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef int colormap[3][256];

typedef unsigned char pixel;

typedef struct {
    pixel r,g,b;
} rgbpixel;

#define  rgb_bright(p)		(int)((77L*(p)->r+151L*(p)->g+28L*(p)->b)/256)

/* image resolution */
extern int	xmax,ymax;

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PIC_H_ */

