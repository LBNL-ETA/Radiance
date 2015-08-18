/* RCSid $Id: g3vector.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef __G3VECTOR_H
#define	__G3VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>

#include "gbasic.h"
#include "fvect.h"

typedef g3Float*	g3Vec;

#define	g3v_angle(v1,v2)	(acos(g3v_dot((v1),(v2))))


g3Vec	g3v_create();
void	g3v_free(g3Vec v);
g3Vec	g3v_copy(g3Vec dest,g3Vec src);
g3Vec	g3v_zero(g3Vec v);

g3Vec	g3v_set(g3Vec v,g3Float x,g3Float y,g3Float z);
g3Vec	g3v_sane(g3Vec v);

g3Vec	g3v_fromrad(g3Vec v,FVECT rv);
void	g3v_torad(FVECT rv,g3Vec v);

/*
**	Orderering is done component wise, with decreasing significance
**	Returns -1, 0 or 1 if v1 is less, equal or greater than v2
*/
int		g3v_epsorder(g3Vec v1,g3Vec v2);

/* 
** 	Returns true if v1 and v2 are epsilon equal 
*/
int		g3v_epseq(g3Vec v1,g3Vec v2);

int		g3v_epszero(g3Vec v);

/* Euclidean length of vector */
g3Float	g3v_length(g3Vec v);

/*
**	Common operations on vectors. All routines returns zero on errors or
**	the result as a g3Float
*/
g3Float	g3v_dot(g3Vec v1,g3Vec v2);
g3Vec	g3v_scale(g3Vec res,g3Vec v,g3Float sc);
g3Vec	g3v_normalize(g3Vec v);
g3Vec	g3v_cross(g3Vec res,g3Vec v1,g3Vec v2);
g3Vec	g3v_sub(g3Vec res,g3Vec v1,g3Vec v2);
g3Vec	g3v_add(g3Vec res,g3Vec v1,g3Vec v2);

g3Vec	g3v_rotate(g3Vec res,g3Vec v,g3Vec rnorm,g3Float ang);

g3Vec	g3v_print(g3Vec v,FILE* fp);
int		g3v_read(g3Vec v,FILE* fp);

#ifdef __cplusplus
}
#endif

#endif
