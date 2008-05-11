/* RCSid $Id: tifftypes.h,v 3.3 2008/05/11 04:55:01 greg Exp $ */
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
#ifndef HAVE_INT8
#ifdef __STDC__
typedef	signed char int8;	/* NB: non-ANSI compilers may not grok */
#else
typedef	char int8;
#endif
typedef	unsigned char uint8;
#define HAVE_INT8
#endif
#ifndef HAVE_INT16
typedef	short int16;
typedef	unsigned short uint16;	/* sizeof (uint16) must == 2 */
#define HAVE_INT16
#endif
#ifndef HAVE_INT32
#if defined(__alpha) || (defined(_MIPS_SZLONG) && _MIPS_SZLONG == 64) || defined(__LP64__) || defined(__x86_64)
typedef	int int32;
typedef	unsigned int uint32;	/* sizeof (uint32) must == 4 */
#else
typedef	long int32;
typedef	unsigned long uint32;	/* sizeof (uint32) must == 4 */
#endif
#define HAVE_INT32
#endif
#endif /* _TIFF_DATA_TYPEDEFS_ */
