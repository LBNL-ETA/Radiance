/* RCSid: $Id: ciq.h,v 2.4 2003/07/27 22:12:03 schorsch Exp $ */

#ifndef _RAD_CIQ_H_
#define _RAD_CIQ_H_

#include <stdio.h>

#include "pic.h"

#ifdef __cplusplus
extern "C" {
#endif

#define red(i) (((i)>>7&0xf8)|4)	/* 5 bits red, 5 bits green, 5 bits blue */
#define gre(i) (((i)>>2&0xf8)|4)
#define blu(i) (((i)<<3&0xf8)|4)
#define len 32768

extern int hist[len];	/* list of frequencies or pixelvalues for coded color */

extern colormap color;	/* quantization colormap */
extern int n;		/* number of colors in it */

#define linealloc(xdim)		(pixel *)emalloc(sizeof(pixel)*xdim)
#define line3alloc(xdim)	(rgbpixel *)emalloc(sizeof(rgbpixel)*xdim)

extern char	*emalloc();


#ifdef __cplusplus
}
#endif
#endif /* _RAD_CIQ_H_ */

