#ifndef lint
static const char	RCSid[] = "$Id: x11raster.c,v 2.12 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
 * x11raster.c - routines to handle images for X windows.
 *
 *	2/18/88
 */

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "color.h"
#include "x11raster.h"


XRASTER *
make_raster(disp, vis, npixbits, data, width, height, bm_pad)
Display	*disp;
XVisualInfo	*vis;
int	npixbits;
char	*data;
int	width, height;
int	bm_pad;
{
	static long	swaptest = 1;
	register XRASTER	*xr;

	if ((xr = (XRASTER *)calloc(1, sizeof(XRASTER))) == NULL)
		return(NULL);
	xr->disp = disp;
	xr->screen = vis->screen;
	xr->visual = vis->visual;
	if (npixbits == 1)
		xr->image = XCreateImage(disp,vis->visual,1,
				XYBitmap,0,data,width,height,bm_pad,0);
	else
		xr->image = XCreateImage(disp,vis->visual,vis->depth,
				ZPixmap,0,data,width,height,bm_pad,0);
	xr->image->bitmap_bit_order = MSBFirst;
	xr->image->bitmap_unit = bm_pad;
	xr->image->byte_order = *(char *)&swaptest ? LSBFirst : MSBFirst;
	if (vis->depth >= 24 && (xr->image->red_mask != 0xff ||
			xr->image->green_mask != 0xff00 ||
			xr->image->blue_mask != 0xff0000) &&
			(xr->image->red_mask != 0xff0000 ||
			xr->image->green_mask != 0xff00 ||
			xr->image->blue_mask != 0xff)) {
		xr->image->red_mask = 0xff;
		xr->image->green_mask = 0xff00;
		xr->image->blue_mask = 0xff0000;
	}
	xr->gc = 0;
	return(xr);
}


int
init_rcolors(xr, cmap)			/* initialize colors */
register XRASTER	*xr;
BYTE	cmap[][3];
{
	register unsigned char	*p;
	register int	i;

	if ((xr->image->depth > 8) | (xr->ncolors != 0))
		return(xr->ncolors);
	xr->pmap = (short *)malloc(256*sizeof(short));
	if (xr->pmap == NULL)
		return(0);
	xr->cdefs = (XColor *)malloc(256*sizeof(XColor));
	if (xr->cdefs == NULL)
		return(0);
	for (i = 0; i < 256; i++)
		xr->pmap[i] = -1;
	for (p = (unsigned char *)xr->image->data,
			i = xr->image->width*xr->image->height;
			i--; p++)
		if (xr->pmap[*p] == -1) {
			xr->cdefs[xr->ncolors].red = cmap[*p][RED] << 8;
			xr->cdefs[xr->ncolors].green = cmap[*p][GRN] << 8;
			xr->cdefs[xr->ncolors].blue = cmap[*p][BLU] << 8;
			xr->cdefs[xr->ncolors].pixel = *p;
			xr->cdefs[xr->ncolors].flags = DoRed|DoGreen|DoBlue;
			xr->pmap[*p] = xr->ncolors++;
		}
	xr->cdefs = (XColor *)realloc((void *)xr->cdefs,
			xr->ncolors*sizeof(XColor));
	if (xr->cdefs == NULL)
		return(0);
	return(xr->ncolors);
}


Colormap
newcmap(disp, scrn, vis)		/* get colormap and fix b & w */
Display	*disp;
int	scrn;
Visual	*vis;
{
	XColor	thiscolor;
	unsigned long	*pixels;
	Colormap	cmap;
	int	n;
	register int	i, j;

	cmap = XCreateColormap(disp, RootWindow(disp,scrn), vis, AllocNone);
	if (cmap == 0 || vis->class != PseudoColor)
		return(cmap);
	pixels=(unsigned long *)malloc(vis->map_entries*sizeof(unsigned long));
	if (pixels == NULL)
		return(0);
	for (n = vis->map_entries; n > 0; n--)
		if (XAllocColorCells(disp, cmap, 0, NULL, 0, pixels, n) != 0)
			break;
	if (n == 0)
		return(0);
					/* reset black and white */
	for (i = 0; i < n; i++) {
		if (pixels[i] != BlackPixel(disp,scrn)
				&& pixels[i] != WhitePixel(disp,scrn))
			continue;
		thiscolor.pixel = pixels[i];
		thiscolor.flags = DoRed|DoGreen|DoBlue;
		XQueryColor(disp, DefaultColormap(disp,scrn), &thiscolor);
		XStoreColor(disp, cmap, &thiscolor);
		for (j = i; j+1 < n; j++)
			pixels[j] = pixels[j+1];
		n--;
		i--;
	}
	XFreeColors(disp, cmap, pixels, n, 0);
	free((void *)pixels);
	return(cmap);
}


unsigned long *
map_rcolors(xr, w)			/* get and assign pixels */
register XRASTER	*xr;
Window	w;
{
	register int	i;
	register unsigned char	*p;

	if (xr->ncolors == 0 || xr->image->depth > 8)
		return(NULL);
	if (xr->pixels != NULL)
		return(xr->pixels);
	xr->pixels = (unsigned long *)malloc(xr->ncolors*sizeof(unsigned long));
	if (xr->pixels == NULL)
		return(NULL);
	if (xr->visual == DefaultVisual(xr->disp, xr->screen))
		xr->cmap = DefaultColormap(xr->disp, xr->screen);
	else
		xr->cmap = newcmap(xr->disp, xr->screen, xr->visual);
	while (XAllocColorCells(xr->disp, xr->cmap, 0,
			NULL, 0, xr->pixels, xr->ncolors) == 0)
		if (xr->cmap == DefaultColormap(xr->disp, xr->screen))
			xr->cmap = newcmap(xr->disp, xr->screen, xr->visual);
		else {
			free((void *)xr->pixels);
			xr->pixels = NULL;
			return(NULL);
		}
	for (i = 0; i < xr->ncolors; i++)
		if (xr->pmap[xr->pixels[i]] == -1)
			break;
	if (i < xr->ncolors) {			/* different pixels */
		for (p = (unsigned char *)xr->image->data,
				i = xr->image->width*xr->image->height;
				i--; p++)
			*p = xr->pixels[xr->pmap[*p]];
		for (i = 0; i < 256; i++)
			xr->pmap[i] = -1;
		for (i = 0; i < xr->ncolors; i++) {
			xr->cdefs[i].pixel = xr->pixels[i];
			xr->pmap[xr->pixels[i]] = i;
		}
		free_rpixmap(xr);		/* Pixmap invalid */
	}
	XStoreColors(xr->disp, xr->cmap, xr->cdefs, xr->ncolors);
	XSetWindowColormap(xr->disp, w, xr->cmap);
	return(xr->pixels);
}


Pixmap
make_rpixmap(xr, w)			/* make pixmap for raster */
register XRASTER	*xr;
Window	w;
{
	XWindowAttributes	xwattr;
	Pixmap	pm;

	if (xr->pm != 0)
		return(xr->pm);
	XGetWindowAttributes(xr->disp, w, &xwattr);
	pm = XCreatePixmap(xr->disp, w,
			xr->image->width, xr->image->height,
			xwattr.depth);
	if (pm == 0)
		return(0);
	put_raster(pm, 0, 0, xr);
	return(xr->pm = pm);
}


patch_raster(d, xsrc, ysrc, xdst, ydst, width, height, xr)	/* redraw */
Drawable	d;
int	xsrc, ysrc, xdst, ydst;
int	width, height;
register XRASTER	*xr;
{
	if (xsrc >= xr->image->width || ysrc >= xr->image->height)
		return;
	if (xsrc < 0) {
		xdst -= xsrc; width += xsrc;
		xsrc = 0;
	}
	if (ysrc < 0) {
		ydst -= ysrc; height += ysrc;
		ysrc = 0;
	}
	if (width <= 0 || height <= 0)
		return;
	if (xsrc + width > xr->image->width)
		width = xr->image->width - xsrc;
	if (ysrc + height > xr->image->height)
		height = xr->image->height - ysrc;

	if (xr->gc == 0) {
		xr->gc = XCreateGC(xr->disp, d, 0, 0);
		XSetState(xr->disp, xr->gc, BlackPixel(xr->disp,xr->screen),
			WhitePixel(xr->disp,xr->screen), GXcopy, AllPlanes);
	}
	if (xr->pm == 0)
		XPutImage(xr->disp, d, xr->gc, xr->image, xsrc, ysrc,
				xdst, ydst, width, height);
	else
		XCopyArea(xr->disp, xr->pm, d, xr->gc, xsrc, ysrc,
				width, height, xdst, ydst);
}


unmap_rcolors(xr)			/* free colors */
register XRASTER	*xr;
{
	if (xr->pixels == NULL)
		return;
	XFreeColors(xr->disp, xr->cmap, xr->pixels, xr->ncolors, 0);
	if (xr->cmap != DefaultColormap(xr->disp, xr->screen))
		XFreeColormap(xr->disp, xr->cmap);
	free((void *)xr->pixels);
	xr->pixels = NULL;
}


free_rpixmap(xr)			/* release Pixmap */
register XRASTER	*xr;
{
	if (xr->pm == 0)
		return;
	XFreePixmap(xr->disp, xr->pm);
	xr->pm = 0;
}


free_raster(xr)				/* free raster data */
register XRASTER	*xr;
{
	free_rpixmap(xr);
	if (xr->ncolors > 0) {
		unmap_rcolors(xr);
		free((void *)xr->pmap);
		free((void *)xr->cdefs);
	}
	XDestroyImage(xr->image);
	if (xr->gc != 0)
		XFreeGC(xr->disp, xr->gc);
	free((void *)xr);
}
