/* RCSid: $Id: x11raster.h,v 2.5 2011/05/20 02:06:39 greg Exp $ */
/*
 * x11raster.h - header file for X routines using images.
 *
 *	3/1/90
 */
#ifndef _RAD_X11RASTER_H_
#define _RAD_X11RASTER_H_

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
	Display	*disp;				/* the display */
	int	screen;				/* the screen */
	Visual	*visual;			/* pointer to visual used */
	XImage	*image;				/* the X image */
	GC	gc;				/* private graphics context */
	int	ncolors;			/* number of colors */
	XColor	*cdefs;				/* color definitions */
	short	*pmap;				/* inverse pixel mapping */
	unsigned long	*pixels;		/* allocated table entries */
	Colormap	cmap;			/* installed color map */
	Pixmap	pm;				/* storage on server side */
}	XRASTER;

extern Colormap newcmap(Display	*disp, int	scrn, Visual	*vis);
extern int init_rcolors(XRASTER	*xr, uby8	cmap[][3]);
extern unsigned long * map_rcolors(XRASTER	*xr, Window	w);
extern Pixmap make_rpixmap(XRASTER	*xr, Window	w);
extern XRASTER * make_raster( Display	*disp, XVisualInfo	*vis,
	int	npixbits, char	*data, int	width, int height, int	bm_pad);
extern void patch_raster(Drawable	d, int	xsrc, int	ysrc,
	int	xdst, int	ydst, int	width, int	height, register XRASTER	*xr);
extern void unmap_rcolors(XRASTER *xr);
extern void free_rpixmap(XRASTER *xr);
extern void free_raster(XRASTER *xr);

#define put_raster(d,xdst,ydst,xr) patch_raster(d,0,0,xdst,ydst, \
				(xr)->image->width,(xr)->image->height,xr)

#ifdef __cplusplus
}
#endif
#endif /* _RAD_X11RASTER_H_ */

