/* RCSid $Id: ambient.h,v 2.25 2019/05/14 17:39:10 greg Exp $ */
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
 * Normal and u-vector directions encoded using dircode.c
 */
typedef struct ambrec {
	struct ambrec  *next;	/* next in list */
	unsigned long  latick;	/* last accessed tick */
	float  pos[3];		/* position in space */
	int32  ndir;		/* encoded surface normal */
	int32  udir;		/* u-vector direction */
	short  lvl;		/* recursion level of parent ray */
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

#ifdef __cplusplus
}
#endif
#endif /* _RAD_AMBIENT_H_ */

