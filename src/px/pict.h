/* RCSid: $Id: pict.h,v 2.4 2003/07/14 22:24:00 schorsch Exp $ */
/*
 Header file for ra_pict.h
*/

#ifndef _RAD_PICT_H_
#define _RAD_PICT_H_

#ifdef __cplusplus
extern "C" {
#endif

#define HEADER_SIZE	    512
#define PICT_picVersion	    0x0011
#define PICT_reservedHeader 0x0C00
#define PICT_clipRgn	    0x0001
#define PICT_EndOfPicture   0x00FF
#define PICT_PackBitsRect   0x0098
#define PICT_Pack32BitsRect 0x009A

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PICT_H_ */
