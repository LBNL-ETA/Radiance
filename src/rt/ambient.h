/* RCSid $Id: ambient.h,v 2.27 2023/11/15 18:02:52 greg Exp $ */
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
	float  pos[3];		/* position in space */
	int32  ndir;		/* encoded surface normal */
	int32  udir;		/* u-vector direction */
	short  lvl;		/* recursion level of parent ray */
	float  weight;		/* weight of parent ray */
	float  rad[2];		/* anisotropic radii (rad[0] <= rad[1]) */
	float  gpos[2];		/* (u,v) gradient wrt. position */
	float  gdir[2];		/* (u,v) gradient wrt. direction */
	uint32  corral;		/* potential light leak direction flags */
	SCOLOR  val;		/* computed indirect irradiance (last!) */
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

#define  AMBVALSIZ	(64+AMB_CNDX[3])	/* number of bytes in portable AMBVAL */
#define  AMBMAGIC	561			/* magic number for ambient value file */
#define  AMBFMT		"Radiance_ambval"	/* format id string */

					/* defined in ambient.c */
extern void	setambres(int ar);
extern void	setambacc(double newa);
extern void	setambient(void);
extern void	multambient(SCOLOR aval, RAY *r, FVECT nrm);
extern void	ambdone(void);
extern void	ambnotify(OBJECT obj);
extern int	ambsync(void);
					/* defined in ambcomp.c */
extern int	doambient(SCOLOR acol, RAY *r, double wt,
				FVECT uv[2], float rad[2],
				float gpos[2], float gdir[2], uint32 *crlp);
					/* defined in ambio.c */
extern int	*AMB_CNDX;	/* open ambient file RGBE indices */
extern float	*AMB_WLPART;	/* open ambient file spectral range */
extern int	amb_headline(char *hl, void *p);
extern void	putambmagic(FILE *fp);
extern int	hasambmagic(FILE *fp);
extern int	writambval(AMBVAL *av, FILE *fp);
extern int	readambval(AMBVAL *av, FILE *fp);
extern int	ambvalOK(AMBVAL *av);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_AMBIENT_H_ */

