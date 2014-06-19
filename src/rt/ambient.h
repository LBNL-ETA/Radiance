/* RCSid $Id: ambient.h,v 2.23 2014/06/19 16:26:55 greg Exp $ */
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

#ifndef OLDAMB

/*
 * Normal and u-vector directions encoded using dircode.c
 */
typedef struct ambrec {
	struct ambrec  *next;	/* next in list */
	unsigned long  latick;	/* last accessed tick */
	float  pos[3];		/* position in space */
	int32  ndir;		/* encoded surface normal */
	int32  udir;		/* u-vector direction */
	int  lvl;		/* recursion level of parent ray */
	float  weight;		/* weight of parent ray */
	float  rad[2];		/* anisotropic radii (rad[0] <= rad[1]) */
	COLOR  val;		/* computed ambient value */
	float  gpos[2];		/* (u,v) gradient wrt. position */
	float  gdir[2];		/* (u,v) gradient wrt. direction */
	uint32  corral;		/* potential light leak direction flags */
}  AMBVAL;			/* ambient value */

typedef struct ambtree {
	AMBVAL	*alist;		/* ambient value list */
	struct ambtree	*kid;	/* 8 child nodes */
}  AMBTREE;			/* ambient octree */

extern double  maxarad;		/* maximum ambient radius */
extern double  minarad;		/* minimum ambient radius */

#ifndef AVGREFL
#define  AVGREFL	0.5	/* assumed average reflectance */
#endif

#define  AMBVALSIZ	67	/* number of bytes in portable AMBVAL struct */
#define  AMBMAGIC	559	/* magic number for ambient value files */
#define  AMBFMT		"Radiance_ambval"	/* format id string */

					/* defined in ambient.c */
extern void	setambres(int ar);
extern void	setambacc(double newa);
extern void	setambient(void);
extern void	multambient(COLOR aval, RAY *r, FVECT nrm);
extern void	ambdone(void);
extern void	ambnotify(OBJECT obj);
extern int	ambsync(void);
					/* defined in ambcomp.c */
extern int	doambient(COLOR acol, RAY *r, double wt,
				FVECT uv[2], float rad[2],
				float gpos[2], float gdir[2], uint32 *crlp);
					/* defined in ambio.c */
extern void	putambmagic(FILE *fp);
extern int	hasambmagic(FILE *fp);
extern int	writambval(AMBVAL *av, FILE *fp);
extern int	readambval(AMBVAL *av, FILE *fp);
extern int	ambvalOK(AMBVAL *av);

#else /* ! NEWAMB */

/*
 * Since we've defined our vectors as float below to save space,
 * watch out for changes in the definitions of VCOPY() and DOT()
 * and don't pass these vectors to fvect routines.
 */
typedef struct ambrec {
	struct ambrec  *next;	/* next in list */
	unsigned long  latick;	/* last accessed tick */
	float  pos[3];		/* position in space */
	float  dir[3];		/* normal direction */
	int    lvl;		/* recursion level of parent ray */
	float  weight;		/* weight of parent ray */
	float  rad;		/* validity radius */
	COLOR  val;		/* computed ambient value */
	float  gpos[3];		/* gradient wrt. position */
	float  gdir[3];		/* gradient wrt. direction */
}  AMBVAL;			/* ambient value */

typedef struct ambtree {
	AMBVAL	*alist;		/* ambient value list */
	struct ambtree	*kid;	/* 8 child nodes */
}  AMBTREE;			/* ambient octree */

typedef struct {
	COLOR  v;		/* division sum (partial) */
	float  r;		/* 1/distance sum */
	float  k;		/* variance for this division */
	int    n;		/* number of subsamples */
	unsigned short  t, p;	/* theta, phi indices */
}  AMBSAMP;		/* ambient sample division */

typedef struct {
	FVECT  ux, uy, uz;	/* x, y and z axis directions */
	COLOR  acoef;		/* division contribution coefficient */
	int    ns;		/* number of super-samples */
	int    nt, np;		/* number of theta and phi directions */
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
extern int	ambsync(void);
					/* defined in ambcomp.c */
extern double	doambient(COLOR acol, RAY *r, double wt,
					FVECT pg, FVECT dg);
extern void	inithemi(AMBHEMI *hp, COLOR ac, RAY *r, double wt);
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

#endif	/* ! NEWAMB */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_AMBIENT_H_ */

