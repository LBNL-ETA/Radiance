/* RCSid: $Id: mkillum.h,v 2.7 2004/09/17 21:43:50 greg Exp $ */
/*
 * Common definitions for mkillum
 */
#ifndef _RAD_MKILLUM_H_
#define _RAD_MKILLUM_H_

#include  "standard.h"
#include  "object.h"
#include  "rtprocess.h"

#ifdef __cplusplus
extern "C" {
#endif
				/* illum flags */
#define  IL_LIGHT	0x1		/* light rather than illum */
#define  IL_COLDST	0x2		/* use color distribution */
#define  IL_COLAVG	0x4		/* use average color */
#define  IL_DATCLB	0x8		/* OK to clobber data file */

struct illum_args {
	int	flags;			/* flags from list above */
	char	matname[MAXSTR];	/* illum material name */
	char	datafile[MAXSTR];	/* distribution data file name */
	int	dfnum;			/* data file number */
	char	altmat[MAXSTR];		/* alternate material name */
	int	sampdens;		/* point sample density */
	int	nsamps;			/* # of samples in each direction */
	float	minbrt;			/* minimum average brightness */
	float	col[3];			/* computed average color */
};				/* illum options */

struct rtproc {
	struct rtproc	*next;		/* next in list of processes */
	SUBPROC pd;			/* rtrace pipe descriptors */
	float	*buf;			/* rtrace i/o buffer */
	int	bsiz;			/* maximum rays for rtrace buffer */
	float	**dest;			/* destination for each ray result */
	int	nrays;			/* current length of rtrace buffer */
};				/* rtrace process */

extern void printobj(char  *mod, register OBJREC  *obj);
extern int average(register struct illum_args  *il, float  *da, int  n);
extern void flatout(struct illum_args  *il, float  *da, int  n, int  m,
	FVECT  u, FVECT  v, FVECT  w);
extern void illumout(register struct illum_args  *il, OBJREC  *ob);
extern void roundout(struct illum_args  *il, float  *da, int  n, int  m);

/* The header file otypes.h has to follow definition of our struct's */
#define FUN_ARGLIST	OBJREC *, struct illum_args *, struct rtproc *, char *
#ifdef __cplusplus
}
#include  "otypes.h"
extern "C" {
#else
#include  "otypes.h"
#endif

extern int o_default(FUN_ARGLIST);
extern int o_face(FUN_ARGLIST);
extern int o_sphere(FUN_ARGLIST);
extern int o_ring(FUN_ARGLIST);
extern void raysamp(float res[3], FVECT org, FVECT dir, struct rtproc *rt0);
extern void rayflush(struct rtproc *rt, int doall);
extern struct rtproc *raywait(struct rtproc *rt0);


#ifdef __cplusplus
}
#endif

#endif /* _RAD_MKILLUM_H_ */

