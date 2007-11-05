/* RCSid: $Id: mkillum.h,v 2.12 2007/11/05 23:40:26 greg Exp $ */
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

/* XXX need to add per-band data */
struct BSDF_data {
	double	om_scale;		/* maximum solid angle (sr/256) */
	int	ninc;			/* number of incoming directions */
	int32	*inc_dir;		/* incoming direction codes */
	BYTE	*inc_rad;		/* incoming radians to neighbor */
	int	nout;			/* number of outgoing directions */
	int32	*out_dir;		/* outgoing direction codes */
	BYTE	*out_rad;		/* outgoing radians to neighbor */
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

#define getBSDF_incvec(v,b,i)	decodedir(v, (b)->inc_dir[i])
#define getBSDF_outvec(v,b,o)	decodedir(v, (b)->out_dir[o])
#define getBSDF_incrad(b,i)	((b)->om_scale*((b)->inc_rad[i] + .5));
#define getBSDF_outrad(b,o)	((b)->om_scale*((b)->out_rad[o] + .5));
#define BSDF_data(b,i,o)	(b)->bsdf[(o)*(b)->ninc + (i)]

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

