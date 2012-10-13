/* RCSid: $Id: mkillum.h,v 2.19 2012/10/13 20:15:43 greg Exp $ */
/*
 * Common definitions for mkillum
 */
#ifndef _RAD_MKILLUM_H_
#define _RAD_MKILLUM_H_

#include  "ray.h"
#include  "otypes.h"
#include  "random.h"

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
	double	thick;			/* object thickness */
	char	matname[MAXSTR];	/* illum material name */
	char	datafile[MAXSTR];	/* distribution data file name */
	int	dfnum;			/* data file number */
	char	altmat[MAXSTR];		/* alternate material name */
	int	sampdens;		/* point sample density */
	int	nsamps;			/* # of samples in each direction */
	float	minbrt;			/* minimum average brightness */
	COLOR	col;			/* computed average color */
};				/* illum options */

extern void printobj(char *mod, OBJREC *obj);
extern int average(struct illum_args *il, COLORV *da, int n);
extern void flatout(struct illum_args *il, COLORV *da, int n, int m,
	FVECT u, FVECT v, FVECT w);
extern void illumout(struct illum_args *il, OBJREC *ob);
extern void roundout(struct illum_args *il, COLORV *da, int n, int m);

extern void newdist(int siz);
extern int process_ray(RAY *r, int rv);
extern void raysamp(int ndx, FVECT org, FVECT dir);
extern void rayclean(void);

extern void flatdir(FVECT  dv, double  alt, double  azi);
extern int flatindex(FVECT dv, int nalt, int nazi);

extern int my_default(OBJREC *, struct illum_args *, char *);
extern int my_face(OBJREC *, struct illum_args *, char *);
extern int my_sphere(OBJREC *, struct illum_args *, char *);
extern int my_ring(OBJREC *, struct illum_args *, char *);

extern COLORV *	distarr;		/* distribution array */
extern int	distsiz;

extern char	*progname;

#ifdef __cplusplus
}
#endif

#endif /* _RAD_MKILLUM_H_ */

