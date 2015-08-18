/* RCSid $Id: g3affine.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef	__G3AFFINE_H
#define __G3AFFINE_H



#ifdef __cplusplus
extern "C" {
#endif

#include "g3vector.h"

#ifdef GL_CONVERSIONS
#include "GL/gl.h"
#endif


typedef	g3Float		g3AVec[4];
typedef g3Float		g3ATrans[4][4];


void	g3av_copy(g3AVec res,g3AVec v);
void	g3at_tounit(g3ATrans t);
void	g3av_vectohom(g3AVec res,g3Vec v);
void	g3av_homtovec(g3Vec res,g3AVec h);
int	 	g3av_tonorm(g3AVec v);

#ifdef GL_CONVERSIONS
void	g3at_get_openGL(g3ATrans t,GLdouble* glm);
void	g3at_from_openGL(g3ATrans t,GLdouble* glm);
#endif


void	g3at_copy(g3ATrans t,g3ATrans tsrc);
void	g3at_print(g3ATrans t,FILE* outf);
void	g3at_comb(g3ATrans t,g3ATrans tf);
void	g3at_comb_to(g3ATrans res,g3ATrans t,g3ATrans tf);
void	g3at_translate(g3ATrans t,g3Vec tv);
void	g3at_add_trans(g3ATrans t,g3Vec tv);
void	g3at_prepend_trans(g3ATrans t,g3Vec tv);
void	g3at_rotate(g3ATrans t,g3Vec rnorm,g3Float ang);
void	g3at_add_rot(g3ATrans t,g3Vec rnorm,g3Float ang);
void	g3at_prepend_rot(g3ATrans t,g3Vec rnorm,g3Float ang);

void	g3at_btrans(g3ATrans t,g3Vec xv,g3Vec yv,g3Vec zv);
int		g3at_sph_ctrans(g3ATrans t,g3Vec vdir,g3Vec vup);
int		g3at_ctrans(g3ATrans t,g3Vec vdir,g3Vec vup);
int		g3at_inverse(g3ATrans t);
void	g3at_apply_h(g3ATrans t,g3AVec v);
void	g3at_iapply_h(g3Vec v,g3ATrans t);
void	g3at_apply(g3ATrans t,g3Vec v);
void	g3at_iapply(g3Vec v,g3ATrans t);

#ifdef __cplusplus
}
#endif


#endif
