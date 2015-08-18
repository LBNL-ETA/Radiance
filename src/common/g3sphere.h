/* RCSid $Id: g3sphere.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef __G3SPHERE_H
#define __G3SPHERE_H
#ifndef M_PI
    #define M_PI       3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "gbasic.h"
#include "g3vector.h"


/**
***	Functions for doing some geometry on spheres.
***	Subscripts in names indicates: "cc" cartesian coords,"sph" spherical
*** coords, normalization on cart. coords ist done, "chk" checks are made.
***	If the function name contains no subsript 
***	no checks are made and normalized cartesian coordinates are expected.
*** All functions allow for "res" to be equivalent to some other param.
*/

#define	G3S_RAD		0
#define	G3S_THETA	1
#define	G3S_PHI		2
#define	G3S_MY		1
#define	G3S_MZ		2

#define	G3S_RESCOPY(res,co) int copy = 0; \
							if (res == co) { \
								res = g3v_create(); \
								copy = 1; \
							}		
#define	G3S_RESFREE(res,co)	if (copy) { \
								g3v_copy(co,res); \
								g3v_free(res); \
								res = co; \
							}

/*
** Some transformation between coord. systems:
** cc: Cartesian 
** sph: spherical
** mtr: metrical spherical
** tr: transversal spherical
*/
g3Vec	g3s_sphtocc(g3Vec res,g3Vec sph);
g3Vec	g3s_cctosph(g3Vec res,g3Vec cc);
g3Vec	g3s_sphwrap(g3Vec sph);
g3Vec	g3s_trwrap(g3Vec tr);
g3Vec	g3s_mtrtocc(g3Vec res,g3Vec mtr);
g3Vec	g3s_cctomtr(g3Vec res,g3Vec cc);
g3Vec	g3s_trtocc(g3Vec res,g3Vec mtr);
g3Vec	g3s_cctotr(g3Vec res,g3Vec cc);
g3Vec	g3s_sphtotr(g3Vec res,g3Vec sph);
g3Vec	g3s_trtosph(g3Vec res,g3Vec tr);
g3Float	g3s_dist(const g3Vec cc1,const g3Vec cc2);
g3Float	g3s_dist_norm(const g3Vec cc1,const g3Vec cc2);


#ifdef __cplusplus
}
#endif

#endif
