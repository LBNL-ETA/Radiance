/* RCSid $Id: view.h,v 2.15 2003/10/27 10:19:31 schorsch Exp $ */
/*
 *  view.h - header file for image generation.
 *
 *  Include after stdio.h and fvect.h
 *  Includes resolu.h
 */
#ifndef _RAD_VIEW_H_
#define _RAD_VIEW_H_

#include  <time.h>

#include  "rtmath.h"
#include  "resolu.h"

#ifdef __cplusplus
extern "C" {
#endif

				/* view types */
#define  VT_PER		'v'		/* perspective */
#define  VT_PAR		'l'		/* parallel */
#define  VT_ANG		'a'		/* angular fisheye */
#define  VT_HEM		'h'		/* hemispherical fisheye */
#define  VT_CYL		'c'		/* cylindrical panorama */

typedef struct {
	int  type;		/* view type */
	FVECT  vp;		/* view origin */
	FVECT  vdir;		/* view direction */
	FVECT  vup;		/* view up */
	double  horiz;		/* horizontal view size */
	double  vert;		/* vertical view size */
	double  hoff;		/* horizontal image offset */
	double  voff;		/* vertical image offset */
	double  vfore;		/* fore clipping plane */
	double  vaft;		/* aft clipping plane (<=0 for inf) */
	FVECT  hvec;		/* computed horizontal image vector */
	FVECT  vvec;		/* computed vertical image vector */
	double  hn2;		/* DOT(hvec,hvec) */
	double  vn2;		/* DOT(vvec,vvec) */
} VIEW;			/* view parameters */

extern VIEW  stdview;

#define  viewaspect(v)	sqrt((v)->vn2/(v)->hn2)

#define  STDVIEW	{VT_PER,{0.,0.,0.},{0.,1.,0.},{0.,0.,1.}, \
				45.,45.,0.,0.,0.,0., \
				{0.,0.,0.},{0.,0.,0.},0.,0.}

#define  VIEWSTR	"VIEW="
#define  VIEWSTRL	5


extern char	*setview(VIEW *v);
extern void	normaspect(double va, double *ap, int *xp, int *yp);
extern double	viewray(FVECT orig, FVECT direc, VIEW *v, double x, double y);
extern void	viewloc(FVECT ip, VIEW *v, FVECT p);
extern void	pix2loc(RREAL loc[2], RESOLU *rp, int px, int py);
extern void	loc2pix(int pp[2], RESOLU *rp, double lx, double ly);
extern int	getviewopt(VIEW *v, int ac, char *av[]);
extern int	sscanview(VIEW *vp, char *s);
extern void	fprintview(VIEW *vp, FILE *fp);
extern char	*viewopt(VIEW *vp);
extern int	isview(char *s);
extern int	viewfile(char *fname, VIEW *vp, RESOLU *rp);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_VIEW_H_ */

