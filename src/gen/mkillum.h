/* RCSid: $Id: mkillum.h,v 2.13 2007/12/05 20:07:34 greg Exp $ */
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

				/* up directions */
typedef enum {
	UDzneg=-3,
	UDyneg=-2,
	UDxneg=-1,
	UDunknown=0,
	UDxpos=1,
	UDypos=2,
	UDzpos=3
} UpDir;
				/* BSDF coordinate calculation routines */
				/* vectors always point away from surface */
typedef int	b_vecf(FVECT, int, char *);
typedef int	b_ndxf(FVECT, char *);
typedef double	b_radf(int, char *);

/* Bidirectional Scattering Distribution Function */
struct BSDF_data {
	int	ninc;			/* number of incoming directions */
	int	nout;			/* number of outgoing directions */
	char	ib_priv[32];		/* input basis private data */
	b_vecf	*ib_vec;		/* get input vector from index */
	b_ndxf	*ib_ndx;		/* get input index from vector */
	b_radf	*ib_rad;		/* get input radius for index */
	char	ob_priv[32];		/* output basis private data */
	b_vecf	*ob_vec;		/* get output vector from index */
	b_ndxf	*ob_ndx;		/* get output index from vector */
	b_radf	*ob_rad;		/* get output radius for index */
	float	*bsdf;			/* scattering distribution data */
};				/* bidirectional scattering distrib. func. */

struct illum_args {
	int	flags;			/* flags from list above */
	UpDir	udir;			/* up direction */
	double	thick;			/* object thickness */
	char	matname[MAXSTR];	/* illum material name */
	char	datafile[MAXSTR];	/* distribution data file name */
	int	dfnum;			/* data file number */
	char	altmat[MAXSTR];		/* alternate material name */
	int	sampdens;		/* point sample density */
	int	nsamps;			/* # of samples in each direction */
	struct BSDF_data
		*sd;			/* scattering data (if set) */
	float	minbrt;			/* minimum average brightness */
	COLOR	col;			/* computed average color */
};				/* illum options */

#define getBSDF_incvec(v,b,i)	(*(b)->ib_vec)(v,i,(b)->ib_priv)
#define getBSDF_incndx(b,v)	(*(b)->ib_ndx)(v,(b)->ib_priv)
#define getBSDF_incrad(b,i)	(*(b)->ib_rad)(i,(b)->ib_priv)
#define getBSDF_outvec(v,b,o)	(*(b)->ob_vec)(v,o,(b)->ob_priv)
#define getBSDF_outndx(b,v)	(*(b)->ob_ndx)(v,(b)->ob_priv)
#define getBSDF_outrad(b,o)	(*(b)->ob_rad)(o,(b)->ob_priv)
#define BSDF_visible(b,i,o)	(b)->bsdf[(o)*(b)->ninc + (i)]

extern struct BSDF_data *load_BSDF(char *fname);
extern void free_BSDF(struct BSDF_data *b);
extern void r_BSDF_incvec(FVECT v, struct BSDF_data *b, int i,
				double rv, MAT4 xm);
extern void r_BSDF_outvec(FVECT v, struct BSDF_data *b, int o,
				double rv, MAT4 xm);
extern int getBSDF_xfm(MAT4 xm, FVECT nrm, UpDir ud);
extern void redistribute(struct BSDF_data *b, int nalt, int nazi,
				FVECT u, FVECT v, FVECT w, MAT4 xm);

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

