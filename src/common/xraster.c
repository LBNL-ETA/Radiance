/*
 * xraster.c - routines to handle images for X windows.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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
