/* RCSid: $Id: holo.h,v 3.25 2004/09/09 01:06:19 greg Exp $ */
/*
 * Header file for holodeck programs
 *
 *	9/29/97	GWLarson
 */
#ifndef _RAD_HOLO_H_
#define _RAD_HOLO_H_

#include "standard.h"
#include "color.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HDMAX
#define HDMAX		128	/* maximum active holodeck sections */
#endif

#ifndef MAXDIRSE
#define MAXDIRSE	32	/* maximum seeks per directory write */
#endif

#define DCINF	(unsigned)((1L<<16)-1)	/* special value for infinity */
#define DCLIN	(unsigned)(1L<<11)	/* linear depth limit */

typedef struct {
	BYTE 	r[2][2];	/* ray direction index */
	COLR 	v;		/* color value */
	uint16	d;		/* depth code (from entry wall) */
} RAYVAL;		/* ray value */

/*
 * Walls are ordered:		X0	X1	X2	WN
 *				0	?	?	0
 *				1	?	?	1
 *				?	0	?	2
 *				?	1	?	3
 *				?	?	0	4
 *				?	?	1	5
 *
 * Grid on wall WN corresponds to X(WN/2+1)%3 and X(WN/2+2)%3, resp.
 */

typedef struct {
	short	w;		/* wall number */
	short	i[2];		/* index on wall grid */
} GCOORD;		/* grid coordinates (two for beam) */

typedef struct {
	uint32	nrd;		/* number of beam rays bundled on disk */
	off_t	fo;		/* position in file */
} BEAMI;		/* beam index */

typedef struct {
	uint32	nrm;		/* number of beam rays bundled in memory */
	unsigned long	tick;	/* clock tick for LRU replacement */
} BEAM;			/* followed by nrm RAYVAL's */

#define hdbray(bp)	((RAYVAL *)((bp)+1))
#define hdbsiz(nr)	(sizeof(BEAM)+(nr)*sizeof(RAYVAL))

typedef struct {
	FVECT	orig;		/* prism origin (first) */
	FVECT	xv[3];		/* side vectors (second) */
	int16	grid[3];	/* grid resolution (third) */
} HDGRID;		/* holodeck section grid (must match HOLO struct) */

typedef struct holo {
	FVECT	orig;		/* prism origin (first) */
	FVECT	xv[3];		/* side vectors (second) */
	int16	grid[3];	/* grid resolution (third) */
	int	fd;		/* file descriptor */
	struct {
		int	s, n;		/* dirty section start, length */
	} dirseg[MAXDIRSE+1];	/* dirty beam index segments */
	short	dirty;		/* number of dirty segments */
	double	tlin;		/* linear range for depth encoding */
	FVECT	wg[3];		/* wall grid vectors (derived) */
	int	wi[6];		/* wall super-indices (derived) */
	char	*priv;		/* pointer to private client data */
	BEAM	**bl;		/* beam pointers (memory cache) */
	BEAMI	bi[1];		/* beam index (extends struct) */
} HOLO;			/* holodeck section */

typedef struct {
	HOLO	*h;		/* pointer to holodeck */
	int	b;		/* beam index */
} HDBEAMI;		/* holodeck beam index */

#define nbeams(hp)	(((hp)->wi[5]-1)<<1)
#define biglob(hp)	((hp)->bi)
#define blglob(hp)	(*(hp)->bl)

#define bnrays(hp,i)	((hp)->bl[i]!=NULL ? (hp)->bl[i]->nrm : (hp)->bi[i].nrd)

#define hdflush(hp)	(hdfreebeam(hp,0), hdsync(hp,0))
#define hdclobber(hp)	(hdkillbeam(hp,0), hdsync(hp,0))

/*
extern HOLO	*hdinit(), *hdalloc();
extern BEAM	*hdgetbeam();
extern RAYVAL	*hdnewrays();
extern unsigned	hdmemuse();
extern off_t	hdfiluse(), hdfilen(), hdallocfrag();
extern double	hdray(), hdinter();
extern unsigned	hdcode();
extern int	hdfilord();
*/

#define FF_NEVER	0		/* never free fragments */
#define FF_WRITE	01		/* free fragment on write */
#define FF_ALLOC	02		/* free fragment on ray alloc */
#define FF_READ		04		/* free fragment on read */
#define FF_KILL		010		/* free fragment on beam kill */

extern int	hdfragflags;		/* tells when to free fragments */
extern unsigned	hdcachesize;		/* target cache size (bytes) */
extern unsigned long	hdclock;	/* holodeck system clock */
extern HOLO	*hdlist[HDMAX+1];	/* holodeck pointers (NULL term.) */

extern float	hd_depthmap[];		/* depth conversion map */

extern int	hdwg0[6];		/* wall grid 0 index */
extern int	hdwg1[6];		/* wall grid 1 index */

#define hddepth(hp,dc)	( (dc) >= DCINF ? FHUGE : \
				(hp)->tlin * ( (dc) >= DCLIN ? \
					hd_depthmap[(dc)-DCLIN] : \
					((dc)+.5)/DCLIN ) )

#define HOLOFMT		"Holodeck"	/* file format identifier */
#define HOLOVERS	0		/* file format version number */
#define HOLOMAGIC	(323+sizeof(off_t)+8*HOLOVERS)	/* file magic number */

/*
 * A holodeck file consists of an information header terminated by a
 * blank line, with "FORMAT=Holodeck" somewhere in it.
 * The first integer after the information header is the
 * above magic number, which includes the file format version number.
 * The first longword after the magic number is a pointer to the pointer
 * just before the SECOND holodeck section, or 0 if there is only one.
 * This longword is immediately followed by the first holodeck
 * section header and directory.
 * Similarly, every holodeck section in the file is preceeded by
 * a pointer to the following section, or 0 for the final section.
 * Since holodeck files consist of directly written C data structures, 
 * they are not generally portable between different machine architectures.
 * In particular, different floating point formats or bit/byte ordering
 * will make the data unusable.  This is unfortunate, and may be changed
 * in future versions, but we thought this would be best for paging speed
 * in our initial implementation.
 */

	/* clumpbeams.c */
extern void clumpbeams(HOLO *hp, int maxcnt, int maxsiz,
int (*cf)(HOLO *hp, int *bqueue, int bqlen));
	/* holo.c */
extern void hdcompgrid(HOLO *hp);
extern int hdbcoord(GCOORD gc[2], HOLO *hp, int i);
extern int hdbindex(HOLO *hp, GCOORD gc[2]);
extern void hdcell(FVECT cp[4], HOLO *hp, GCOORD *gc);
extern int hdlseg(int	lseg[2][3], HOLO	*hp, GCOORD	gc[2]);
extern unsigned int hdcode(HOLO *hp, double d);
extern void hdgrid( FVECT gp, HOLO *hp, FVECT wp);
extern void hdworld(FVECT wp, HOLO *hp, FVECT gp);
extern double hdray(FVECT ro, FVECT rd, HOLO *hp, GCOORD gc[2], BYTE r[2][2]);
extern double hdinter(GCOORD gc[2], BYTE r[2][2], double *ed, HOLO *hp,
		FVECT ro, FVECT rd);
	/* holofile.c */
extern HOLO * hdinit(int fd, HDGRID *hproto);
extern void hddone(HOLO *hp);
extern int hdsync(HOLO *hp, int all);
extern off_t hdfilen(int fd);
extern off_t hdfiluse(int fd);
extern RAYVAL * hdnewrays(HOLO *hp, int i, int nr);
extern BEAM * hdgetbeam(HOLO *hp, int i);
extern void hdloadbeams(HDBEAMI *hb, int n, void (*bf)(BEAM *bp, HDBEAMI *hb));
extern int hdfreebeam(HOLO *hp, int i);
extern int hdfreefrag(HOLO *hp, int i);
extern int hdfragOK(int fd, int *listlen, int32 *listsiz);
extern int hdkillbeam(HOLO *hp, int i);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_HOLO_H_ */

