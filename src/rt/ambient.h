/* RCSid $Id: ambient.h,v 2.14 2005/05/28 22:27:54 greg Exp $ */
/*
 * Common definitions for interreflection routines.
 *
 * Include after ray.h
 */
#ifndef _RAD_AMBIENT_H_
#define _RAD_AMBIENT_H_
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Since we've defined our vectors as float below to save space,
 * watch out for changes in the definitions of VCOPY() and DOT()
 * and don't pass these vectors to fvect routines.
 */
typedef struct ambrec {
	unsigned long  latick;	/* last accessed tick */
	float  pos[3];		/* position in space */
	float  dir[3];		/* normal direction */
	int    lvl;		/* recursion level of parent ray */
	float  weight;		/* weight of parent ray */
	float  rad;		/* validity radius */
	COLOR  val;		/* computed ambient value */
	float  gpos[3];		/* gradient wrt. position */
	float  gdir[3];		/* gradient wrt. direction */
	struct ambrec  *next;	/* next in list */
}  AMBVAL;			/* ambient value */

typedef struct ambtree {
	AMBVAL	*alist;		/* ambient value list */
	struct ambtree	*kid;	/* 8 child nodes */
}  AMBTREE;			/* ambient octree */

typedef struct {
	short  t, p;		/* theta, phi indices */
	COLOR  v;		/* value sum */
	float  r;		/* 1/distance sum */
	float  k;		/* variance for this division */
	int  n;			/* number of subsamples */
}  AMBSAMP;		/* ambient sample division */

typedef struct {
	COLOR  acoef;		/* ambient contribution coefficient */
	float  drc;		/* daughter ray coefficient */
	FVECT  ux, uy, uz;	/* x, y and z axis directions */
	short  nt, np;		/* number of theta and phi directions */
	int    ns;		/* number of super-samples */
}  AMBHEMI;		/* ambient sample hemisphere */

extern double  maxarad;		/* maximum ambient radius */
extern double  minarad;		/* minimum ambient radius */

#ifndef AVGREFL
#define  AVGREFL	0.5	/* assumed average reflectance */
#endif

#define  AMBVALSIZ	75	/* number of bytes in portable AMBVAL struct */
#define  AMBMAGIC	557	/* magic number for ambient value files */
#define  AMBFMT		"Radiance_ambval"	/* format id string */

					/* defined in ambient.c */
extern void	setambres(int ar);
extern void	setambacc(double newa);
extern void	setambient(void);
extern void	multambient(COLOR aval, RAY *r, FVECT nrm);
extern void	ambdone(void);
extern void	ambnotify(OBJECT obj);
extern double	sumambient(COLOR acol, RAY *r, double rw, FVECT rn, int al,
				AMBTREE *at, FVECT c0, double s);
extern double	makeambient(COLOR acol, RAY *r, COLOR ac, FVECT rn, int al);
extern void	extambient(COLOR cr, AMBVAL *ap, FVECT pv, FVECT nv);
extern int	ambsync(void);
					/* defined in ambcomp.c */
extern double	doambient(COLOR acol, RAY *r, COLOR ac, double wt,
					FVECT pg, FVECT dg);
extern void	inithemi(AMBHEMI *hp, RAY *r, COLOR ac, double wt);
extern int	divsample(AMBSAMP *dp, AMBHEMI *h, RAY *r);
extern void	comperrs(AMBSAMP *da, AMBHEMI *hp);
extern void	posgradient(FVECT gv, AMBSAMP *da, AMBHEMI *hp);
extern void	dirgradient(FVECT gv, AMBSAMP *da, AMBHEMI *hp);
					/* defined in ambio.c */
extern void	putambmagic(FILE *fp);
extern int	hasambmagic(FILE *fp);
extern int	writambval(AMBVAL *av, FILE *fp);
extern int	ambvalOK(AMBVAL *av);
extern int	readambval(AMBVAL *av, FILE *fp);
					/* defined in lookamb.c */
extern void	lookamb(FILE *fp);
extern void	writamb(FILE *fp);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_AMBIENT_H_ */

