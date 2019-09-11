/* RCSid $Id: bsdf_m.h,v 3.10 2019/09/11 00:24:03 greg Exp $ */
/*
 *  bsdf_m.h
 *  
 *  Support for BSDF matrices.
 *  Assumes "bsdf.h" already included.
 *  Include after "ezxml.h" for SDloadMtx() declaration.
 *
 *  Created by Greg Ward on 2/2/11.
 *
 */

#ifndef _BSDF_M_H_
#define	_BSDF_M_H_

#ifdef __cplusplus
extern "C" {
#endif
				/* Fixed-position coordinate functions */
typedef int	b_vecf(FVECT vec, double ndxr, void *c_data);
typedef int	b_ndxf(const FVECT vec, void *c_data);
typedef double	b_ohmf(int ndx, void *c_data);

/* Rectangular matrix format BSDF */
typedef struct {
	int		ninc;		/* number of incoming directions */
	int		nout;		/* number of outgoing directions */
	void		*ib_priv;	/* input basis private data */
	b_vecf		*ib_vec;	/* get input vector from index */
	b_ndxf		*ib_ndx;	/* get input index from vector */
	b_ohmf		*ib_ohm;	/* get input proj. SA for index */
	void		*ob_priv;	/* output basis private data */
	b_vecf		*ob_vec;	/* get output vector from index */
	b_ndxf		*ob_ndx;	/* get output index from vector */
	b_ohmf		*ob_ohm;	/* get output proj. SA for index */
	C_CHROMA	*chroma;	/* chromaticity data */
	float		bsdf[1];	/* scattering data (extends struct) */
} SDMat;

/* Matrix BSDF accessors */
#define mBSDF_incvec(v,b,ix)	(*(b)->ib_vec)(v,ix,(b)->ib_priv)
#define mBSDF_incndx(b,v)	(*(b)->ib_ndx)(v,(b)->ib_priv)
#define mBSDF_incohm(b,i)	(*(b)->ib_ohm)(i,(b)->ib_priv)
#define mBSDF_outvec(v,b,ox)	(*(b)->ob_vec)(v,ox,(b)->ob_priv)
#define mBSDF_outndx(b,v)	(*(b)->ob_ndx)(v,(b)->ob_priv)
#define mBSDF_outohm(b,o)	(*(b)->ob_ohm)(o,(b)->ob_priv)
#define mBSDF_value(b,o,i)	(b)->bsdf[(o)*(b)->ninc + (i)]
#define mBSDF_chroma(b,o,i)	(b)->chroma[(o)*(b)->ninc + (i)]

/* Holder for cumulative distribution (sum of BSDF * projSA) */
typedef struct SDMatCDst_s {
	SD_CDIST_BASE(SDMatCDst_s);	/* base fields; must come first */
	int		indx;		/* incident angle index */
	void		*ob_priv;	/* private data for generator */
	b_vecf		*ob_vec;	/* outbound vector generator */
	int		calen;		/* cumulative array length */
	unsigned	carr[1];	/* cumulative array (extends struct) */
} SDMatCDst;	

#ifdef _EZXML_H
/* Load a set of BSDF matrices from an open XML file */
extern SDError		SDloadMtx(SDData *sd, ezxml_t wtl);
#endif

/* Our matrix handling routines */
extern const SDFunc	SDhandleMtx;

/******** Klems basis declarations for more intimate access ********/

#define MAXLATS		46		/* maximum number of latitudes */

/* BSDF angle specification */
typedef struct {
	char	name[64];		/* basis name */
	int	nangles;		/* total number of directions */
	struct {
		float	tmin;			/* starting theta */
		int	nphis;			/* number of phis (0 term) */
	} lat[MAXLATS+1];		/* latitudes */
} ANGLE_BASIS;

#define	MAXABASES	7		/* limit on defined bases */

extern ANGLE_BASIS	abase_list[MAXABASES];

extern int		nabases;	/* current number of defined bases */

extern C_COLOR	mtx_RGB_prim[3];	/* matrix RGB primaries  */
extern float	mtx_RGB_coef[3];	/* corresponding Y coefficients */

/* Get color or grayscale value for BSDF in the given directions */
extern int		mBSDF_color(float coef[], const SDMat *b, int i, int o);

/* Get vector for this angle basis index (front exiting) */
extern b_vecf		fo_getvec;

/* Get index corresponding to the given vector (front exiting) */
extern b_ndxf		fo_getndx;

/* Get projected solid angle for this angle basis index (universal) */
extern b_ohmf		io_getohm;

/* Get vector for this angle basis index (back incident) */
extern b_vecf		bi_getvec;

/* Get index corresponding to the vector (back incident) */
extern b_ndxf		bi_getndx;

/* Get vector for this angle basis index (back exiting) */
extern b_vecf		bo_getvec;

/* Get index corresponding to the vector (back exiting) */
extern b_ndxf		bo_getndx;

/* Get vector for this angle basis index (front incident) */
extern b_vecf		fi_getvec;

/* Get index corresponding to the vector (front incident) */
extern b_ndxf		fi_getndx;

#ifdef __cplusplus
}
#endif
#endif	/* ! _BSDF_M_H_ */
