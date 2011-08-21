/* RCSid $Id: bsdf_t.h,v 3.9 2011/08/21 22:38:12 greg Exp $ */
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

#define SD_UFRONT	0x1	/* flag for BSDF uses front side */
#define SD_UBACK	0x2	/* flag for BSDF uses back side */
#define SD_XMIT		0x3	/* combine the two for transmission */

/* Variable-resolution BSDF holder */
typedef struct {
	int	sidef;		/* transmitted component? */
	SDNode	*st;		/* BSDF tree */
} SDTre;

/* Holder for cumulative distribution (sum of BSDF * projSA) */
typedef struct SDTreCDst_s {
				/* base fields; must come first */
	SD_CDIST_BASE(SDTreCDst_s);
	double	clim[2][2];	/* input coordinate limits */
	double	max_psa;	/* maximum projected solid angle */
	short	sidef;		/* which side(s) to use */
	short	isodist;	/* isotropic distribution? */
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
