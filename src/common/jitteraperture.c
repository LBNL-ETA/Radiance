#ifndef lint
static const char	RCSid[] = "$Id: jitteraperture.c,v 2.3 2022/01/14 22:57:40 greg Exp $";
#endif
/*
 *  jitteraperture.c - routine to sample depth-of-field
 *
 *  External symbols declared in view.h
 */

#include "copyright.h"

#include  <stdio.h>
#include  "rtmath.h"
#include  "random.h"
#include  "view.h"


int
jitteraperture(			/* random aperture shift for depth-of-field */
FVECT  orig,		/* assigned previously by... */
FVECT  direc,		/* ...viewray() call */
VIEW  *v,
double  dia
)
{
	RREAL  df[2];
	double  adj;
	int  i;
					/* are we even needed? */
	if (dia <= FTINY)
		return(1);
					/* get random point on disk */
	square2disk(df, frandom(), frandom());
	df[0] *= .5*dia;
	df[1] *= .5*dia;

	if ((v->type == VT_PER) | (v->type == VT_PAR)) {
		adj = 1.0;		/* basic view cases */
		if (v->type == VT_PER)
			adj /= DOT(direc, v->vdir);
		df[0] /= sqrt(v->hn2);
		df[1] /= sqrt(v->vn2);
		for (i = 3; i--; ) {
			orig[i] = v->vp[i] + df[0]*v->hvec[i] +
						df[1]*v->vvec[i] ;
			direc[i] = v->vp[i] + adj*v->vdist*direc[i]
					- orig[i];
		}
	} else {			/* difficult view cases */
		adj = PI/4.*dia*(.5 - frandom());
		if ((v->type != VT_ANG) & (v->type != VT_PLS)) {
			if (v->type != VT_CYL)
				df[0] /= sqrt(v->hn2);
			df[1] /= sqrt(v->vn2);
		}
		for (i = 3; i--; ) {
			orig[i] = v->vp[i] + df[0]*v->hvec[i] +
						df[1]*v->vvec[i] +
						adj*v->vdir[i] ;
			direc[i] = v->vp[i] + v->vdist*direc[i]
					- orig[i];
		}
	}
	if (normalize(direc) == 0.0)
		return(0);
	if ((adj = v->vfore) <= FTINY)
		return(1);
	if (v->type == VT_PER)
		adj /= DOT(direc, v->vdir);
	VSUM(orig, orig, direc, adj);
	return(1);
}
