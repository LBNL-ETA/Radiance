/* RCSid $Id: radogl.h,v 3.8 2003/06/26 00:58:09 schorsch Exp $ */
/*
 * Header file for Radiance - OpenGL routines.
 */
#ifndef _RAD_RADOGL_H_
#define _RAD_RADOGL_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "copyright.h"

#include "standard.h"
#include <GL/glu.h>
#include "color.h"
#include "object.h"
#include "otypes.h"
#include "lookup.h"

#define MAXLIGHTS	8	/* number of OGL light sources */

#define MAXSPECEXP	128.	/* maximum allowed specular exponent */
#define UNKSPECEXP	25.	/* value to use when exponent unknown */

typedef struct {
	short	type;		/* material type (from otypes.h) */
	short	nlinks;		/* number of links to this material */
	union {
		struct {
			COLOR	ambdiff;	/* ambient and diffuse color */
			COLOR	specular;	/* specular color */
			GLfloat	specexp;	/* specular exponent */
		}	m;		/* regular material */
		struct {
			COLOR	emission;	/* emitting component */
			GLfloat	spotdir[3];	/* spot direction */
			GLfloat	spotang;	/* spot cutoff angle */
		}	l;		/* light source */
	}	u;		/* union of types */
}	MATREC;		/* OGL material properties */

extern double	expval;		/* global exposure value */
extern COLOR	ambval;		/* global ambient value */
extern int	glightid[MAXLIGHTS];	/* OpenGL GL_LIGHTi values */
extern int	dolights;	/* are we outputting light sources? */
extern int	domats;		/* are we doing materials? */

extern LUTAB	mtab;		/* material/modifier lookup table */

#define issrcmat(m)	((m) != NULL && islight((m)->type) && \
				(m)->type != MAT_GLOW)

				/* defined in rgldomat.c */
extern void	domatobj(MATREC *mp, FVECT cent);
extern void	domatvert(MATREC *mp, FVECT v, FVECT n);
				/* defined in rglfile.c */
extern int	newglist(void);
extern void	rgl_checkerr(char *where);
extern int	rgl_filelist(int ic, char **inp, int *nl);
extern int	rgl_octlist(char *fname, FVECT cent, RREAL *radp, int *nl);
extern void	rgl_load(char *inpspec);
extern void	rgl_object(char *name, FILE *fp);
				/* defined in rglinst.c */
extern int	o_instance(OBJREC *o);
extern int	loadoctrees(void);
extern double	checkoct(char *fname, FVECT cent);
extern int	loadoct(char *fname);
				/* defined in rglmat.c */
extern void	rgl_matclear(void);
extern MATREC	*getmatp(char *nam);
extern int	o_default(OBJREC *o);
extern int	o_unsupported(OBJREC *o);
extern MATREC	*newmaterial(char *nam);
extern void	freemtl(MATREC *mp);
extern int	m_normal(OBJREC *o);
extern int	m_aniso(OBJREC *o);
extern int	m_glass(OBJREC *o);
extern int	m_brdf(OBJREC *o);
extern int	m_brdf2(OBJREC *o);
extern int	m_light(OBJREC *o);
extern int	m_mirror(OBJREC *o);
extern int	m_prism(OBJREC *o);
				/* defined in rglsrc.c */
extern void	lightinit(void);
extern void	lightclean(void);
extern void	lightdefs(void);
extern int	o_source(OBJREC *o);
extern int	doflatsrc(MATREC *m, FVECT pos, FVECT norm, double area);
extern int	dosphsrc(MATREC *m, FVECT pos, double area);
				/* defined in rglsurf.c */
extern void	setmaterial(MATREC *mp, FVECT cent, int ispoly);
extern double	polyarea(FVECT cent, FVECT norm, int n, FVECT v[]);
extern int	o_face(OBJREC *o);
extern void	surfclean(void);
extern int	o_sphere(OBJREC *o);
extern int	o_cone(OBJREC *o);
extern int	o_ring(OBJREC *o);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_RADOGL_H_ */

