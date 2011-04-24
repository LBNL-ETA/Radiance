/* RCSid $Id: bsdf_t.h,v 3.6 2011/04/24 19:39:21 greg Exp $ */
/*
 *  bsdf_t.h
 *  
 *  Support for variable-resolution BSDF trees.
 *  Assumes "bsdf.h" already included.
 *  Include after "ezxml.h" for SDloadTre() declaration.
 *
 *  Created by Greg Ward on 2/2/11.
 *
 */

#ifndef _BSDF_T_H_
#define	_BSDF_T_H_

#ifdef __cplusplus
extern "C" {
#endif

#define SD_MAXDIM	4	/* maximum expected # dimensions */

/* Basic node structure for variable-resolution BSDF data */
typedef struct SDNode_s {
	short	ndim;		/* number of dimensions */
	short	log2GR;		/* log(2) of grid resolution (< 0 for tree) */
	union {
		struct SDNode_s	*t[1];		/* subtree pointers */
		float		v[1];		/* scattering value(s) */
	} u;			/* subtrees or values (extends struct) */
} SDNode;

/* Variable-resolution BSDF holder */
typedef struct {
	int	isxmit;		/* transmitted component? */
	SDNode	*st;		/* BSDF tree */
} SDTre;

/* Holder for cumulative distribution (sum of BSDF * projSA) */
typedef struct {
	SD_CDIST_BASE;		/* base fields; must come first */
	double	clim[2][2];	/* input coordinate limits */
	double	max_psa;	/* maximum projected solid angle */
	int	isxmit;		/* transmitted component? */
	int	calen;		/* cumulative array length */
	struct {
		unsigned	hndx;	/* hilbert index */
		unsigned	cuml;	/* cumulative value */
	}		carr[1];	/* cumulative array (extends struct) */
} SDTreCDst;	

#ifdef _EZXML_H
/* Load a variable-resolution BSDF tree from an open XML file */
extern SDError		SDloadTre(SDData *sd, ezxml_t wtl);
#endif

/* Our matrix handling routines */
extern SDFunc	SDhandleTre;

#ifdef __cplusplus
}
#endif
#endif	/* ! _BSDF_T_H_ */
