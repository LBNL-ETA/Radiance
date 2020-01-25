/* RCSid $Id: depthcodec.h,v 2.7 2020/01/25 05:35:34 greg Exp $ */
/*
 * Definitions and declarations for 16-bit depth encode/decode
 *
 *  Include after rtio.h and fvect.h
 *  Includes view.h
 */

#ifndef _RAD_DEPTHCODEC_H_
#define _RAD_DEPTHCODEC_H_

#include "view.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEPTHSTR	"REFDEPTH="
#define LDEPTHSTR	9
#define DEPTH16FMT	"16-bit_encoded_depth"

#define	HF_HEADIN	0x1			/* expect input header */
#define HF_HEADOUT	0x2			/* write header to stdout */
#define HF_RESIN	0x4			/* expect resolution string */
#define HF_RESOUT	0x8			/* write resolution to stdout */
#define HF_STDERR	0x10			/* report errors to stderr */
#define HF_ALL		0x1f			/* all flags above */
#define HF_ENCODE	0x20			/* we are encoding */

/* Structure for encoding/decoding depths and world points */
typedef struct {
	FILE		*finp;			/* input stream */
	const char	*inpname;		/* input name */
	short		format;			/* decoded format */
	short		swapped;		/* byte-swapped input */
	short		last_dc;		/* last depth-code read */
	short		use_last;		/* use last depth-code next */
	long		dstart;			/* start of data */
	long		curpos;			/* current input position */
	double		refdepth;		/* reference depth */
	char		depth_unit[32];		/* string including units */
	short		hdrflags;		/* header i/o flags */
	char		inpfmt[MAXFMTLEN];	/* format from header */
	short		gotview;		/* got input view? */
	VIEW		vw;			/* input view parameters */
	RESOLU		res;			/* input resolution */
} DEPTHCODEC;

/* Encode depth as 16-bit signed integer */
#if 1
#define	depth2code(d, dref) \
		( (d) > (dref) ? (int)(32768.001 - 32768.*(dref)/(d))-1 : \
		  (d) > .0 ? (int)(32767.*(d)/(dref) - 32768.999) : -32768 )
#else
extern int	depth2code(double d, double dref);
#endif

/* Decode depth from 16-bit signed integer */
#if 1
#define code2depth(c, dref) \
		( (c) <= -32768 ? .0 : (c) >= 32767 ? FHUGE : \
		  (c) < -1 ? (dref)*(32768.5 + (c))*(1./32767.) : \
				(dref)*32768./(32766.5 - (c)) )
#else
extern double	code2depth(int c, double dref);
#endif

/* Set codec defaults */
extern void	set_dc_defaults(DEPTHCODEC *dcp);

/* Load/copy header */
extern int	process_dc_header(DEPTHCODEC *dcp, int ac, char *av[]);

/* Check that we have what we need to decode depths */
extern int	check_decode_depths(DEPTHCODEC *dcp);

/* Decode next depth pixel from input */
extern double	decode_depth_next(DEPTHCODEC *dcp);

/* Seek to the indicated pixel position */
extern int	seek_dc_pix(DEPTHCODEC *dcp, int x, int y);

/* Read and decode depth for the given pixel */
extern double	decode_depth_pix(DEPTHCODEC *dcp, int x, int y);

/* Check that we have what we need to decode world positions */
extern int	check_decode_worldpos(DEPTHCODEC *dcp);

/* Compute world position from pixel position and depth */
extern int	compute_worldpos(FVECT wpos, DEPTHCODEC *dcp,
					int x, int y, double d);

/* Decode the next world position from input */
int		decode_worldpos_next(FVECT wpos, DEPTHCODEC *dcp);

/* Decode depth and compute world position for the given pixel */
extern int	get_worldpos_pix(FVECT wpos, DEPTHCODEC *dcp, int x, int y);

extern char	*progname;	/* global argv[0] (set by main) */

#ifdef __cplusplus
}
#endif
#endif		/* _RAD_DEPTHCODEC_H_ */
