/* RCSid $Id: bsdf_t.h,v 3.2 2011/02/18 00:41:44 greg Exp $ */
/*
 *  bsdf_t.h
 *  
 *  Support for variable-resolution BSDF trees
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

/* Load a variable-resolution BSDF tree from an open XML file */
extern SDError		SDloadTre(SDData *sd, ezxml_t fl);

/* Our matrix handling routines */
extern const SDFunc	SDhandleTre;

#ifdef __cplusplus
}
#endif
#endif	/* ! _BSDF_T_H_ */
