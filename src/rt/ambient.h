/* RCSid $Id: ambient.h,v 2.8 2003/02/25 02:47:22 greg Exp $ */
/*
 * Common definitions for interreflection routines.
 *
 * Include after ray.h
 */

#include "copyright.h"

/*
 * Since we've defined our vectors as float below to save space,
 * watch out for changes in the definitions of VCOPY() and DOT()
 * and don't pass these vectors to fvect routines.
 */
typedef struct ambrec {
	unsigned long  latick;	/* last accessed tick */
	float  pos[3];		/* position in space */
	float  dir[3];		/* normal direction */
	int  lvl;		/* recursion level of parent ray */
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
	FVECT  ux, uy, uz;	/* x, y and z axis directions */
	short  nt, np;		/* number of theta and phi directions */
}  AMBHEMI;		/* ambient sample hemisphere */

extern double  maxarad;		/* maximum ambient radius */
extern double  minarad;		/* minimum ambient radius */

#define  AVGREFL	0.5	/* assumed average reflectance */

#define  AMBVALSIZ	75	/* number of bytes in portable AMBVAL struct */
#define  AMBMAGIC	557	/* magic number for ambient value files */
#define  AMBFMT		"Radiance_ambval"	/* format id string */

#ifdef NOPROTO

extern int	divsample();
extern double	doambient();
extern void	inithemi();
extern void	comperrs();
extern void	posgradient();
extern void	dirgradient();
extern void	setambres();
extern void	setambacc();
extern void	setambient();
extern void	ambdone();
extern void	ambnotify();
extern void	ambient();
extern double	sumambient();
extern double	makeambient();
extern void	extambient();
extern int	ambsync();
extern void	putambmagic();
extern int	hasambmagic();
extern int	writambval();
extern int	ambvalOK();
extern int	readambval();
extern void	lookamb();
extern void	writamb();

#else
					/* defined in ambcomp.c */
extern int	divsample(AMBSAMP *dp, AMBHEMI *h, RAY *r);
extern double	doambient(COLOR acol, RAY *r, double wt, FVECT pg, FVECT dg);
extern void	inithemi(AMBHEMI *hp, RAY *r, double wt);
extern void	comperrs(AMBSAMP *da, AMBHEMI *hp);
extern void	posgradient(FVECT gv, AMBSAMP *da, AMBHEMI *hp);
extern void	dirgradient(FVECT gv, AMBSAMP *da, AMBHEMI *hp);
					/* defined in ambient.c */
extern void	setambres(int ar);
extern void	setambacc(double newa);
extern void	setambient(void);
extern void	ambdone(void);
extern void	ambnotify(OBJECT obj);
extern void	ambient(COLOR acol, RAY *r, FVECT nrm);
extern double	sumambient(COLOR acol, RAY *r, FVECT rn, int al,
				AMBTREE *at, FVECT c0, double s);
extern double	makeambient(COLOR acol, RAY *r, FVECT rn, int al);
extern void	extambient(COLOR cr, AMBVAL *ap, FVECT pv, FVECT nv);
extern int	ambsync(void);
					/* defined in ambio.c */
extern void	putambmagic(FILE *fp);
extern int	hasambmagic(FILE *fp);
extern int	writambval(AMBVAL *av, FILE *fp);
extern int	ambvalOK(AMBVAL *av);
extern int	readambval(AMBVAL *av, FILE *fp);
					/* defined in lookamb.c */
extern void	lookamb(FILE *fp);
extern void	writamb(FILE *fp);

#endif
