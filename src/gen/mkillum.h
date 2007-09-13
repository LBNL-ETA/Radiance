/* RCSid: $Id: mkillum.h,v 2.8 2007/09/13 06:31:21 greg Exp $ */
/*
 * Common definitions for mkillum
 */
#ifndef _RAD_MKILLUM_H_
#define _RAD_MKILLUM_H_

#include  "ray.h"
#include  "otypes.h"

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

extern void printobj(char  *mod, register OBJREC  *obj);
extern int average(register struct illum_args  *il, float  *da, int  n);
extern void flatout(struct illum_args  *il, float  *da, int  n, int  m,
	FVECT  u, FVECT  v, FVECT  w);
extern void illumout(register struct illum_args  *il, OBJREC  *ob);
extern void roundout(struct illum_args  *il, float  *da, int  n, int  m);

extern int my_default(OBJREC *, struct illum_args *, char *);
extern int my_face(OBJREC *, struct illum_args *, char *);
extern int my_sphere(OBJREC *, struct illum_args *, char *);
extern int my_ring(OBJREC *, struct illum_args *, char *);

extern char	*progname;

#ifdef __cplusplus
}
#endif

#endif /* _RAD_MKILLUM_H_ */

