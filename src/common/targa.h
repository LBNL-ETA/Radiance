/* RCSid $Id: targa.h,v 2.5 2003/06/27 06:53:22 greg Exp $ */
/*
 *  tardev.h - header file for reading and writing Targa format files.
 */
#ifndef _RAD_TARGA_H_
#define _RAD_TARGA_H_
#ifdef __cplusplus
extern "C" {
#endif
 
			/* header structure adapted from tardev.h */
struct hdStruct {
	char textSize;			/* size of info. line ( < 256) */
	char mapType;			/* color map type */
	char dataType;			/* data type */
	int  mapOrig;			/* first color index */
	int  mapLength;			/* length of file map */
	char CMapBits;			/* bits per map entry */
	int  XOffset;			/* picture offset */
	int  YOffset;
	int  x;				/* picture size */
	int  y;
	int  dataBits;			/* bits per pixel */
	int  imType;			/* image descriptor byte */
};

#define  IM_NODATA	0		/* no data included */
#define  IM_CMAP	1		/* color-mapped */
#define  IM_RGB		2		/* straight RGB */
#define  IM_MONO	3		/* straight monochrome */
#define  IM_CCMAP	9		/* compressed color-mapped */
#define  IM_CRGB	10		/* compressed RGB */
#define  IM_CMONO	11		/* compressed monochrome */

				/* color map types */
#define  CM_NOMAP	0		/* no color map */
#define  CM_HASMAP	1		/* has color map */

#define  bits_bytes(n)	(((n)+7)>>3)	/* number of bits to number of bytes */


#ifdef __cplusplus
}
#endif
#endif /* _RAD_TARGA_H_ */
