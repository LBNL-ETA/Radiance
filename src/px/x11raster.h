/* RCSid: $Id: x11raster.h,v 2.2 2003/02/22 02:07:28 greg Exp $ */
/*
 * x11raster.h - header file for X routines using images.
 *
 *	3/1/90
 */

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

extern Colormap	newcmap();

extern unsigned long	*map_rcolors();

extern Pixmap	make_rpixmap();

extern XRASTER	*make_raster();

#define put_raster(d,xdst,ydst,xr) patch_raster(d,0,0,xdst,ydst, \
				(xr)->image->width,(xr)->image->height,xr)
