/* RCSid: $Id: holo.h,v 3.21 2003/05/29 16:26:21 greg Exp $ */
/*
 * Header file for holodeck programs
 *
 *	9/29/97	GWLarson
 */

#include "standard.h"
#include "color.h"

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
	unsigned int2	d;	/* depth code (from entry wall) */
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
	unsigned int4	nrd;	/* number of beam rays bundled on disk */
	off_t	fo;		/* position in file */
} BEAMI;		/* beam index */

typedef struct {
	unsigned int4	nrm;	/* number of beam rays bundled in memory */
	unsigned long	tick;	/* clock tick for LRU replacement */
} BEAM;			/* followed by nrm RAYVAL's */

#define hdbray(bp)	((RAYVAL *)((bp)+1))
#define hdbsiz(nr)	(sizeof(BEAM)+(nr)*sizeof(RAYVAL))

typedef struct {
	FVECT	orig;		/* prism origin (first) */
	FVECT	xv[3];		/* side vectors (second) */
	int2	grid[3];	/* grid resolution (third) */
} HDGRID;		/* holodeck section grid (must match HOLO struct) */

typedef struct holo {
	FVECT	orig;		/* prism origin (first) */
	FVECT	xv[3];		/* side vectors (second) */
	int2	grid[3];	/* grid resolution (third) */
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

extern HOLO	*hdinit(), *hdalloc();
extern BEAM	*hdgetbeam();
extern RAYVAL	*hdnewrays();
extern unsigned	hdmemuse();
extern off_t	hdfiluse(), hdfilen(), hdallocfrag();
extern double	hdray(), hdinter();
extern unsigned	hdcode();
extern int	hdfilord();

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
