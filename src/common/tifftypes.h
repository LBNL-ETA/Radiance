/* RCSid $Id: tifftypes.h,v 3.1 2003/06/20 00:25:49 greg Exp $ */
/* The following is lifted straight out of tiff.h */
#ifndef _TIFF_DATA_TYPEDEFS_
#define _TIFF_DATA_TYPEDEFS_
/*
 * Intrinsic data types required by the file format:
 *
 * 8-bit quantities	int8/uint8
 * 16-bit quantities	int16/uint16
 * 32-bit quantities	int32/uint32
 */
#ifdef __STDC__
typedef	signed char int8;	/* NB: non-ANSI compilers may not grok */
#else
typedef	char int8;
#endif
typedef	unsigned char uint8;
typedef	short int16;
typedef	unsigned short uint16;	/* sizeof (uint16) must == 2 */
#if defined(__alpha) || (defined(_MIPS_SZLONG) && _MIPS_SZLONG == 64)
typedef	int int32;
typedef	unsigned int uint32;	/* sizeof (uint32) must == 4 */
#else
typedef	long int32;
typedef	unsigned long uint32;	/* sizeof (uint32) must == 4 */
#endif
#endif /* _TIFF_DATA_TYPEDEFS_ */
