/* RCSid: $Id: mkillum.h,v 2.6 2003/11/21 07:15:30 greg Exp $ */
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

#ifdef __cplusplus
}
#endif

#define FUN_ARGLIST	OBJREC *, struct illum_args *, struct rtproc *, char *
#include  "otypes.h"

#endif /* _RAD_MKILLUM_H_ */

