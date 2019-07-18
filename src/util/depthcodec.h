/* RCSid $Id: depthcodec.h,v 2.2 2019/07/18 22:33:34 greg Exp $ */
/*
 * Definitions and declarations for 16-bit depth encode/decode
 *
 *  Include after stdio.h and fvect.h
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

#define	HF_HEADIN	0x1
#define HF_HEADOUT	0x2
#define HF_RESIN	0x4
#define HF_RESOUT	0x8
#define HF_ALL		0xf
#define HF_ENCODE	0x10

/* Structure for encoding/decoding depths and world points */
typedef struct {
	FILE		*finp;			/* input stream */
	const char	*inpname;		/* input name */
	int		format;			/* decoded format */
	long		dstart;			/* start of data */
	long		curpos;			/* current input position */
	double		refdepth;		/* reference depth */
	char		depth_unit[32];		/* string including units */
	int		hdrflags;		/* header i/o flags */
	char		inpfmt[MAXFMTLEN];	/* format from header */
	VIEW		vw;			/* input view parameters */
	int		gotview;		/* got input view? */
	RESOLU		res;			/* input resolution */
} DEPTHCODEC;

/* Encode depth as 16-bit signed integer */
#define	depth2code(d, dref) \
		( (d) > (dref) ? (int)(32768 - 32768*(dref)/(d))-1 : \
		  (d) > .0 ? (int)(32767*(d)/(dref) - 32768) : -32768 )

/* Decode depth from 16-bit signed integer */
extern double	code2depth(int c, double dref);

/* Set codec defaults */
extern void	set_dc_defaults(DEPTHCODEC *dcp);

/* Load/copy header */
extern int	process_dc_header(DEPTHCODEC *dcp, int ac, char *av[]);

/* Check that we have what we need to decode depths */
extern int	check_decode_depths(DEPTHCODEC *dcp);

/* Decode next depth pixel */
extern double	decode_depth_next(DEPTHCODEC *dcp);

/* Seek to the indicated pixel position */
extern int	seek_dc_pix(DEPTHCODEC *dcp, int x, int y);

/* Read and decode depth for the given pixel */
extern double	decode_depth_pix(DEPTHCODEC *dcp, int x, int y);

/* Check that we have what we need to decode world positions */
extern int	check_decode_worldpos(DEPTHCODEC *dcp);

/* Compute world position from depth */
extern int	compute_worldpos(FVECT wpos, DEPTHCODEC *dcp,
					int x, int y, double d);

/* Decode the next world position */
int		decode_worldpos_next(FVECT wpos, DEPTHCODEC *dcp);

/* Read and decode the world position at the given pixel */
extern int	get_worldpos_pix(FVECT wpos, DEPTHCODEC *dcp, int x, int y);

extern char	*progname;	/* global argv[0] (set by main) */

#ifdef __cplusplus
}
#endif
#endif		/* _RAD_DEPTHCODEC_H_ */
