/* RCSid $Id: xraster.h,v 3.4 2003/07/14 22:24:00 schorsch Exp $ */
/*
 * xraster.h - header file for X routines using images.
 */
#ifndef _RAD_XRASTER_H_
#define _RAD_XRASTER_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int	width, height;			/* image size */
	int	ncolors;			/* number of colors */
	union {
		unsigned short	*m;			/* monochrome */
		unsigned char	*bz;			/* color */
	}	data;				/* storage on our side */
	Pixmap	pm;				/* storage on server side */
	Color	*cdefs;				/* color definitions */
	int	*pmap;				/* inverse pixel mapping */
	int	*pixels;			/* allocated table entries */
}	XRASTER;

extern int	*map_rcolors();

extern Pixmap	make_rpixmap();

#ifdef __cplusplus
}
#endif
#endif /* _RAD_XRASTER_H_ */

