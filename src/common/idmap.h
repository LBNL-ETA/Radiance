/* RCSid $Id: idmap.h,v 2.2 2019/07/26 18:37:21 greg Exp $ */
/*
 * Definitions and delcarations for loading identifier maps
 *
 *  Include after stdio.h
 *  Includes resolu.h
 */

#ifndef _RAD_IDMAP_H_
#define _RAD_IDMAP_H_

#include "resolu.h"

#ifdef __cplusplus
extern "C" {
#endif

#define	IDMAP8FMT	"8-bit_indexed_name"
#define	IDMAP16FMT	"16-bit_indexed_name"
#define IDMAP24FMT	"24-bit_indexed_name"
#define	IDMAPFMT	"*-bit_indexed_name"

#define HF_HEADOUT	0x2			/* copy header to stdout */
#define HF_RESOUT	0x8			/* copy resolution */
#define HF_STDERR	0x10			/* report errors to stderr */

/* Structure for reading identifier maps */
typedef struct {
	FILE		*finp;			/* input file pointer */
	long		dstart;			/* start of data */
	long		curpos;			/* current input position */
	RESOLU		res;			/* input resolution */
	int		bytespi;		/* 1, 2, or 3 bytes per index */
	int		nids;			/* ID count */
	int		*idoffset;		/* offsets to identifiers */
	char		idtable[1];		/* nul-terminated ID strings */
} IDMAP;

/* Get indexed ID */
#define mapID(mp,i)	( (i) >= (mp)->nids ? (char *)NULL : \
				(mp)->idtable+(mp)->idoffset[i] )

/* Open ID map file for reading, copying info to stdout based on hflags */
extern IDMAP		*idmap_ropen(const char *fname, int hflags);

/* Read the next ID index from input */
extern int		idmap_next_i(IDMAP *idmp);

/* Read the next ID from input */
extern const char	*idmap_next(IDMAP *idmp);

/* Seek to a specific pixel position in ID map */
extern int		idmap_seek(IDMAP *idmp, int x, int y);

/* Read ID at a pixel position */
extern const char	*idmap_pix(IDMAP *idmp, int x, int y);

/* Close ID map and free resources */
extern void		idmap_close(IDMAP *idmp);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_IDMAP_H_ */
