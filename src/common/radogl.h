/* RCSid: $Id: radogl.h,v 3.4 2003/02/22 02:07:22 greg Exp $ */
/*
 * Header file for Radiance - OpenGL routines.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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

#ifdef NOPROTO

extern void	domatobj();
extern void	domatvert();
extern int	newglist();
extern void	rgl_checkerr();
extern int	rgl_filelist();
extern int	rgl_octlist();
extern void	rgl_load();
extern void	rgl_object();
extern int	o_instance();
extern int	loadoctrees();
extern double	checkoct();
extern int	loadoct();
extern void	rgl_matclear();
extern MATREC	*getmatp();
extern int	o_default();
extern MATREC	*newmaterial();
extern void	freemtl();
extern int	m_normal();
extern int	m_aniso();
extern int	m_glass();
extern int	m_brdf();
extern int	m_brdf2();
extern int	m_light();
extern int	m_mirror();
extern int	m_prism();
extern void	lightinit();
extern void	lightclean();
extern void	lightdefs();
extern int	o_source();
extern int	doflatsrc();
extern int	dosphsrc();
extern void	setmaterial();
extern double	polyarea();
extern int	o_face();
extern void	surfclean();
extern int	o_sphere();
extern int	o_cone();
extern int	o_ring();

#else
				/* defined in rgldomat.c */
extern void	domatobj(MATREC *mp, FVECT cent);
extern void	domatvert(MATREC *mp, FVECT v, FVECT n);
				/* defined in rglfile.c */
extern int	newglist(void);
extern void	rgl_checkerr(char *where);
extern int	rgl_filelist(int ic, char **inp, int *nl);
extern int	rgl_octlist(char *fname, FVECT cent, FLOAT *radp, int *nl);
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

#endif
