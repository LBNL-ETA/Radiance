/* RCSid: $Id: pic.h,v 2.4 2004/03/28 20:33:14 schorsch Exp $ */
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

	/* defined in closest.c */
extern void initializeclosest(void);
extern int closest(int r, int g, int b);

	/* defined in cut.c */
extern int makecm(int nw, int *na);

	/* defined in ciq.c */
extern void ciq(int dith, int nw, int synth, colormap cm);

	/* defined in biq.c */
extern void biq(int dith, int nw, int synth, colormap cm);

	/* defined in in the calling program */
extern void picreadcm(colormap  map);
extern void picwritecm(colormap  cm);
extern void picwriteline(int  y, pixel  *l);
extern void picreadline3(int  y, rgbpixel  *l3);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PIC_H_ */

