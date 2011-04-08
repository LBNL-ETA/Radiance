/* RCSid $Id: bsdf_m.h,v 3.4 2011/04/08 18:13:48 greg Exp $ */
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
	float		bsdf[1];	/* scattering data (extends struct) */
} SDMat;

/* Matrix BSDF accessors */
#define mBSDF_incvec(v,b,ix)	(*(b)->ib_vec)(v,ix,(b)->ib_priv)
#define mBSDF_incndx(b,v)	(*(b)->ib_ndx)(v,(b)->ib_priv)
#define mBSDF_incohm(b,i)	(*(b)->ib_ohm)(i,(b)->ib_priv)
#define mBSDF_outvec(v,b,ox)	(*(b)->ob_vec)(v,ox,(b)->ob_priv)
#define mBSDF_outndx(b,v)	(*(b)->ob_ndx)(v,(b)->ob_priv)
#define mBSDF_outohm(b,o)	(*(b)->ob_ohm)(o,(b)->ob_priv)
#define mBSDF_value(b,i,o)	(b)->bsdf[(o)*(b)->ninc + (i)]

/* Holder for cumulative distribution (sum of BSDF * projSA) */
typedef struct {
	SD_CDIST_BASE;			/* base fields; must come first */
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
extern SDFunc		SDhandleMtx;

#ifdef __cplusplus
}
#endif
#endif	/* ! _BSDF_M_H_ */
