/* RCSid $Id: source.h,v 2.11 2003/07/14 22:24:00 schorsch Exp $ */
/*
 *  source.h - header file for ray tracing sources.
 *
 *  Include after ray.h
 */
#ifndef _RAD_SOURCE_H_
#define _RAD_SOURCE_H_

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define  AIMREQT	100		/* required aim success/failure */

#define  SDISTANT	01		/* source distant flag */
#define  SSKIP		02		/* source skip flag */
#define  SPROX		04		/* source proximity flag */
#define  SSPOT		010		/* source spotlight flag */
#define  SVIRTUAL	020		/* source virtual flag */
#define  SFLAT		040		/* source flat flag */
#define  SCYL		0100		/* source cylindrical flag */
#define  SFOLLOW	0200		/* source follow path flag */

typedef struct {
	FVECT  aim;		/* aim direction or center */
	float  siz;		/* output solid angle or area */
	float  flen;		/* focal length (negative if distant source) */
} SPOT;			/* spotlight */

typedef struct {
	int  sflags;		/* source flags */
	FVECT  sloc;		/* direction or position of source */
	FVECT  ss[3];		/* source dimension vectors, U, V, and W */
	float  srad;		/* maximum source radius */
	float  ss2;		/* solid angle or projected area */
	struct {
		float  prox;		/* proximity */
		SPOT  *s;		/* spot */
	} sl;			/* localized source information */
	union {
		long  success;		/* successes - AIMREQT*failures */
		struct {
			short  pn;		/* projection number */
			int  sn;		/* next source to aim for */
		}  sv;			/* virtual source */
	} sa;			/* source aiming information */
	unsigned long
		ntests, nhits;	/* shadow tests and hits */
	OBJREC  *so;		/* source destination object */
}  SRCREC;		/* light source */

#define MAXSPART	64		/* maximum partitions per source */

#define SU		0		/* U vector or partition */
#define SV		1		/* V vector or partition */
#define SW		2		/* W vector or partition */
#define S0		3		/* leaf partition */

#define snorm		ss[SW]		/* normal vector for flat source */

typedef struct {
	int  sn;				/* source number */
	short  np;				/* number of partitions */
	short  sp;				/* this partition number */
	double  dom;				/* solid angle of partition */
	unsigned char  spt[MAXSPART/2];		/* source partitioning */
}  SRCINDEX;		/* source index structure */

#define initsrcindex(s)	((s)->sn = (s)->sp = -1, (s)->np = 0)

#define clrpart(pt)	memset((char *)(pt), '\0', MAXSPART/2)
#define setpart(pt,i,v)	((pt)[(i)>>2] |= (v)<<(((i)&3)<<1))
#define spart(pt,pi)	((pt)[(pi)>>2] >> (((pi)&3)<<1) & 3)

/*
 * Special support functions for sources
 */

/*
 * Virtual source materials must define the following.
 *
 *	vproj(pm, op, sp, i)	Compute i'th virtual projection
 *				of source sp in object op and assign
 *				the 4x4 transformation matrix pm.
 *				Return 1 on success, 0 if no i'th projection.
 *
 *	nproj			The number of projections.  The value of
 *				i passed to vproj runs from 0 to nproj-1.
 */

typedef struct {
	int  (*vproj)();	/* project virtual sources */
	int  nproj;		/* number of possible projections */
} VSMATERIAL;		/* virtual source material functions */

typedef struct {
	void	(*setsrc)();	/* set light source for object */
	void	(*partit)();	/* partition light source object */
	double  (*getpleq)();	/* plane equation for surface */
	double  (*getdisk)();	/* maximum disk for surface */
} SOBJECT;		/* source object functions */

typedef union {
	VSMATERIAL  *mf;	/* material functions */
	SOBJECT  *of;		/* object functions */
} SRCFUNC;		/* source functions */

extern SRCFUNC  sfun[];			/* source dispatch table */

extern SRCREC  *source;			/* our source list */
extern int  nsources;			/* the number of sources */

#define  sflatform(sn,dir)	-DOT(source[sn].snorm, dir)

#define  getplaneq(c,o)		(*sfun[(o)->otype].of->getpleq)(c,o)
#define  getmaxdisk(c,o)	(*sfun[(o)->otype].of->getdisk)(c,o)
#define  setsource(s,o)		(*sfun[(o)->otype].of->setsrc)(s,o)

					/* defined in source.c */
extern void	marksources(void);
extern void	freesources(void);
extern int	srcray(RAY *sr, RAY *r, SRCINDEX *si);
extern void	srcvalue(RAY *r);
extern int	sourcehit(RAY *r);
extern void	direct(RAY *r, void (*f)(), char *p);
extern void	srcscatter(RAY *r);
extern int	m_light(OBJREC *m, RAY *r);
					/* defined in srcsamp.c */
extern double	nextssamp(RAY *r, SRCINDEX *si);
extern int	skipparts(int ct[3], int sz[3], int pp[2], unsigned char *pt);
extern void	nopart(SRCINDEX *si, RAY *r);
extern void	cylpart(SRCINDEX *si, RAY *r);
extern void	flatpart(SRCINDEX *si, RAY *r);
extern double	scylform(int sn, FVECT dir);
					/* defined in srcsupp.c */
extern void	initstypes(void);
extern int	newsource(void);
extern void	setflatss(SRCREC *src);
extern void	fsetsrc(SRCREC *src, OBJREC *so);
extern void	ssetsrc(SRCREC *src, OBJREC *so);
extern void	sphsetsrc(SRCREC *src, OBJREC *so);
extern void	rsetsrc(SRCREC *src, OBJREC *so);
extern void	cylsetsrc(SRCREC *src, OBJREC *so);
extern SPOT	*makespot(OBJREC *m);
extern int	spotout(RAY *r, SPOT *s);
extern double	fgetmaxdisk(FVECT ocent, OBJREC *op);
extern double	rgetmaxdisk(FVECT ocent, OBJREC *op);
extern double	fgetplaneq(FVECT nvec, OBJREC *op);
extern double	rgetplaneq(FVECT nvec, OBJREC *op);
extern int	commonspot(SPOT *sp1, SPOT *sp2, FVECT org);
extern int	commonbeam(SPOT *sp1, SPOT *sp2, FVECT org);
extern int	checkspot(SPOT *sp, FVECT nrm);
extern double	spotdisk(FVECT oc, OBJREC *op, SPOT *sp, FVECT pos);
extern double	beamdisk(FVECT oc, OBJREC *op, SPOT *sp, FVECT dir);
extern double	intercircle(FVECT cc, FVECT c1, FVECT c2,
			double r1s, double r2s);
					/* defined in virtuals.c */
extern void	markvirtuals(void);
extern void	addvirtuals(int sn, int nr);
extern void	vproject(OBJREC *o, int sn, int n);
extern OBJREC	*vsmaterial(OBJREC *o);
extern int	makevsrc(OBJREC *op, int sn, MAT4 pm);
extern double	getdisk(FVECT oc, OBJREC *op, int sn);
extern int	vstestvis(int f, OBJREC *o, FVECT oc, double or2, int sn);
extern void	virtverb(int sn, FILE *fp);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_SOURCE_H_ */

