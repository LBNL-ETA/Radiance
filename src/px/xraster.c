#ifndef lint
static const char	RCSid[] = "$Id: xraster.c,v 3.2 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * xraster.c - routines to handle images for X windows.
 */

#include "copyright.h"

#include <stdio.h>

#include <X/Xlib.h>

#include "xraster.h"


int *
map_rcolors(xr)					/* get and assign pixels */
register XRASTER	*xr;
{
	register int	i;
	register unsigned char	*p;
	int	j;

	if (xr->ncolors == 0)
		return(NULL);
	if (xr->pixels != NULL)
		return(xr->pixels);
	xr->pixels = (int *)malloc(xr->ncolors*sizeof(int));
	if (xr->pixels == NULL)
		return(NULL);
	if (XGetColorCells(0, xr->ncolors, 0, &j, xr->pixels) == 0) {
		free((void *)xr->pixels);
		xr->pixels = NULL;
		return(NULL);
	}
	for (i = 0; i < xr->ncolors; i++)
		if (xr->pmap[xr->pixels[i]] == -1)
			break;
	if (i < xr->ncolors) {			/* different pixels */
		for (p = xr->data.bz, i = BZPixmapSize(xr->width,xr->height);
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
	XStoreColors(xr->ncolors, xr->cdefs);
	return(xr->pixels);
}


Pixmap
make_rpixmap(xr)			/* make pixmap for raster */
register XRASTER	*xr;
{
	Bitmap	bm;

	if (xr->pm != 0)
		return(xr->pm);
	if (xr->ncolors > 0) {
		if (DisplayPlanes() < 2 || DisplayPlanes() > 8)
			return(0);
		xr->pm = XStorePixmapZ(xr->width, xr->height, xr->data.bz);
	} else {
		bm = XStoreBitmap(xr->width, xr->height, xr->data.m);
		if (bm == 0)
			return(0);
		xr->pm = XMakePixmap(bm, BlackPixel, WhitePixel);
		XFreeBitmap(bm);
	}
	return(xr->pm);
}


unmap_rcolors(xr)			/* free colors */
register XRASTER	*xr;
{
	if (xr->pixels == NULL)
		return;
	XFreeColors(xr->pixels, xr->ncolors, 0);
	free((void *)xr->pixels);
	xr->pixels = NULL;
}


free_rpixmap(xr)			/* release Pixmap */
register XRASTER	*xr;
{
	if (xr->pm == 0)
		return;
	XFreePixmap(xr->pm);
	xr->pm = 0;
}


free_raster(xr)				/* free raster data */
register XRASTER	*xr;
{
	free_rpixmap(xr);
	if (xr->ncolors > 0) {
		unmap_rcolors(xr);
		free((void *)xr->data.bz);
		free((void *)xr->pmap);
		free((void *)xr->cdefs);
	} else
		free((void *)xr->data.m);
	free((void *)xr);
}


put_raster(w, x, y, xr)			/* put raster into window */
Window	w;
int	x, y;
register XRASTER	*xr;
{
	if (xr->pm != 0) {
		XPixmapPut(w, 0, 0, x, y, xr->width, xr->height,
				xr->pm, GXcopy, AllPlanes);
	} else if (xr->ncolors > 0) {
		if (DisplayPlanes() < 2 || DisplayPlanes() > 8)
			return(0);
		XPixmapBitsPutZ(w, x, y, xr->width, xr->height,
				xr->data.bz, 0, GXcopy, AllPlanes);
	} else {
		XBitmapBitsPut(w, x, y, xr->width, xr->height,
				xr->data.m, BlackPixel, WhitePixel,
				0, GXcopy, AllPlanes);
	}
	return(1);
}


patch_raster(w, xsrc, ysrc, xdst, ydst, width, height, xr)	/* piece */
Window	w;
int	xsrc, ysrc, xdst, ydst;
int	width, height;
register XRASTER	*xr;
{
	register char	*p;

	if (xsrc >= xr->width || ysrc >= xr->height)
		return(1);
	if (xsrc < 0) {
		xdst -= xsrc; width += xsrc;
		xsrc = 0;
	}
	if (ysrc < 0) {
		ydst -= ysrc; height += ysrc;
		ysrc = 0;
	}
	if (width <= 0 || height <= 0)
		return(1);
	if (xsrc + width > xr->width)
		width = xr->width - xsrc;
	if (ysrc + height > xr->height)
		height = xr->height - ysrc;
	if (xr->pm != 0) {
		XPixmapPut(w, xsrc, ysrc, xdst, ydst, width, height,
				xr->pm, GXcopy, AllPlanes);
	} else if (xr->ncolors > 0) {
		if (DisplayPlanes() < 2 || DisplayPlanes() > 8)
			return(0);
		p = (char *)xr->data.bz + BZPixmapSize(xr->width,ysrc)
				+ BZPixmapSize(xsrc,1);
		while (height--) {
			XPixmapBitsPutZ(w, xdst, ydst, width, 1, p,
					0, GXcopy, AllPlanes);
			p += BZPixmapSize(xr->width,1);
			ydst++;
		}
	} else {
		xdst -= xsrc&15; width += xsrc&15;
		xsrc &= ~15;
		p = (char *)xr->data.m + BitmapSize(xr->width,ysrc)
				+ BitmapSize(xsrc,1);
		while (height--) {
			XBitmapBitsPut(w, xdst, ydst, width, 1,
					p, BlackPixel, WhitePixel,
					0, GXcopy, AllPlanes);
			p += BitmapSize(xr->width,1);
			ydst++;
		}
	}
	return(1);
}
