#ifndef lint
static const char	RCSid[] = "$Id: rglmat.c,v 3.2 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Routines for Radiance -> OpenGL materials.
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

#include "radogl.h"

int	domats = 1;			/* are we doing materials? */

LUTAB	mtab = LU_SINIT(free,freemtl);


void
rgl_matclear()			/* clean up materials */
{
	lu_done(&mtab);
	domats = 1;
}


MATREC *
getmatp(nam)			/* find material record for modifier name */
char	*nam;
{
	register LUENT	*lup;

	if (nam == NULL)
		return(NULL);
	if ((lup = lu_find(&mtab, nam)) == NULL)
		return(NULL);
	return((MATREC *)lup->data);
}


int
o_default(o)			/* default object is non-material modifier */
register OBJREC	*o;
{
	register LUENT	*lup;
#ifdef DEBUG
	if (o->otype >= 0 && !ismodifier(o->otype))
		error(CONSISTENCY, "o_default handed non-modifier");
#endif
					/* find name in lookup table */
	if ((lup = lu_find(&mtab, o->oname)) == NULL)
		goto memerr;
	if (lup->key == NULL) {		/* new entry? */
		lup->key = (char *)malloc(strlen(o->oname)+1);
		if (lup->key == NULL)
			goto memerr;
		strcpy(lup->key, o->oname);
	} else if (lup->data != NULL)
		freemtl((MATREC *)lup->data);
	if ((lup->data = o->os) != NULL)	/* make material reference */
		((MATREC *)lup->data)->nlinks++;
	return;
memerr:
	error(SYSTEM, "out of memory in o_default");
}


MATREC *
newmaterial(nam)		/* get an entry for a new material */
char	*nam;
{
	register LUENT	*lup;
					/* look it up (assign entry) */
	if ((lup = lu_find(&mtab, nam)) == NULL)
		goto memerr;
	if (lup->key == NULL) {		/* new entry? */
		lup->key = (char *)malloc(strlen(nam)+1);
		if (lup->key == NULL)
			goto memerr;
		strcpy(lup->key, nam);
	} else if (lup->data != NULL)
		freemtl((MATREC *)lup->data);
	lup->data = (char *)malloc(sizeof(MATREC));
	if (lup->data == NULL)
		goto memerr;
	((MATREC *)lup->data)->nlinks = 1;
	return((MATREC *)lup->data);
memerr:
	error(SYSTEM, "out of memory in newmaterial");
}


void
freemtl(mp)			/* free a material */
register MATREC	*mp;
{
	if (!--mp->nlinks)
		free((void *)mp);
}


int
m_normal(o)			/* compute normal material parameters */
register OBJREC	*o;
{
	register MATREC	*m;
					/* check arguments */
	if (o->oargs.nfargs != (o->otype == MAT_TRANS ? 7 : 5))
		objerror(o, USER, "bad # of real arguments");
					/* allocate/insert material */
	m = newmaterial(o->oname);
					/* assign parameters */
	setcolor(m->u.m.ambdiff, o->oargs.farg[0],
			o->oargs.farg[1], o->oargs.farg[2]);
	if ((m->type = o->otype) == MAT_METAL)
		copycolor(m->u.m.specular, m->u.m.ambdiff);
	else
		setcolor(m->u.m.specular, 1., 1., 1.);
	scalecolor(m->u.m.specular, o->oargs.farg[3]);
	scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[3]);
	if (m->type == MAT_TRANS) {
		scalecolor(m->u.m.specular, 1.-o->oargs.farg[5]);
		scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[5]);
	}
	if (o->oargs.farg[4] <= FTINY)
		m->u.m.specexp = MAXSPECEXP;
	else
		m->u.m.specexp = 2./(o->oargs.farg[4]*o->oargs.farg[4]);
	if (m->u.m.specexp > MAXSPECEXP)
		m->u.m.specexp = MAXSPECEXP;
}


int
m_aniso(o)			/* anisotropic material */
register OBJREC	*o;
{
	register MATREC	*m;
					/* check arguments */
	if (o->oargs.nfargs < (o->otype == MAT_TRANS2 ? 8 : 6))
		objerror(o, USER, "bad # of real arguments");
					/* allocate/insert material */
	m = newmaterial(o->oname);
					/* assign parameters */
	setcolor(m->u.m.ambdiff, o->oargs.farg[0],
			o->oargs.farg[1], o->oargs.farg[2]);
	if ((m->type = o->otype) == MAT_METAL2)
		copycolor(m->u.m.specular, m->u.m.ambdiff);
	else
		setcolor(m->u.m.specular, 1., 1., 1.);
	scalecolor(m->u.m.specular, o->oargs.farg[3]);
	scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[3]);
	if (m->type == MAT_TRANS2) {
		scalecolor(m->u.m.specular, 1.-o->oargs.farg[6]);
		scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[6]);
	}
	if (o->oargs.farg[4]*o->oargs.farg[5] <= FTINY*FTINY)
		m->u.m.specexp = MAXSPECEXP;
	else
		m->u.m.specexp = 2./(o->oargs.farg[4]*o->oargs.farg[5]);
	if (m->u.m.specexp > MAXSPECEXP)
		m->u.m.specexp = MAXSPECEXP;
}


int
m_glass(o)			/* glass material (hopeless) */
OBJREC	*o;
{
	register MATREC	*m;

	m = newmaterial(o->oname);
	m->type = o->otype;
	setcolor(m->u.m.ambdiff, 0., 0., 0.);
	setcolor(m->u.m.specular, .08, .08, .08);
	m->u.m.specexp = MAXSPECEXP;
}


int
m_brdf(o)			/* convert functional material */
register OBJREC	*o;
{
	register MATREC	*m;
					/* check arguments */
	if (o->oargs.nfargs < (o->otype == MAT_TFUNC ? 6 : 4))
		objerror(o, USER, "bad # of real arguments");
					/* allocate/insert material */
	m = newmaterial(o->oname);
					/* assign parameters */
	setcolor(m->u.m.ambdiff, o->oargs.farg[0],
			o->oargs.farg[1], o->oargs.farg[2]);
	if ((m->type = o->otype) == MAT_MFUNC)
		copycolor(m->u.m.specular, m->u.m.ambdiff);
	else
		setcolor(m->u.m.specular, 1., 1., 1.);
	scalecolor(m->u.m.specular, o->oargs.farg[3]);
	scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[3]);
	if (m->type == MAT_TFUNC) {
		scalecolor(m->u.m.specular, 1.-o->oargs.farg[4]);
		scalecolor(m->u.m.ambdiff, 1.-o->oargs.farg[4]);
	}
	m->u.m.specexp = UNKSPECEXP;
}


int
m_brdf2(o)			/* convert advanced functional material */
register OBJREC	*o;
{
	register MATREC	*m;

	if (o->oargs.nfargs < 9)
		objerror(o, USER, "bad # of real arguments");
	m = newmaterial(o->oname);
	m->type = o->otype;
					/* assign average diffuse front+back */
	setcolor(m->u.m.ambdiff, (o->oargs.farg[0]+o->oargs.farg[3])*.5,
				 (o->oargs.farg[1]+o->oargs.farg[4])*.5,
				 (o->oargs.farg[2]+o->oargs.farg[5])*.5);
					/* guess the rest */
	setcolor(m->u.m.specular, .1, .1, .1);
	m->u.m.specexp = UNKSPECEXP;
}


int
m_light(o)			/* convert light type */
register OBJREC	*o;
{
	FVECT	v;
	register MATREC	*m;

	if (o->oargs.nfargs < (o->otype == MAT_SPOT ? 7 : 3))
		objerror(o, USER, "bad # of real arguments");
	m = newmaterial(o->oname);
	setcolor(m->u.l.emission, o->oargs.farg[0],
			o->oargs.farg[1], o->oargs.farg[2]);
	if ((m->type = o->otype) == MAT_SPOT) {
		if ((m->u.l.spotang = o->oargs.farg[3]/2.) > 90.)
			m->u.l.spotang = 180.;
		v[0] = o->oargs.farg[4];
		v[1] = o->oargs.farg[5];
		v[2] = o->oargs.farg[6];
		if (normalize(v) == 0.)
			objerror(o, USER, "illegal direction");
		VCOPY(m->u.l.spotdir, v);
	} else {
		m->u.l.spotang = 180.;
		m->u.l.spotdir[0] = m->u.l.spotdir[1] = 0.;
		m->u.l.spotdir[2] = -1.;
	}
}


int
m_mirror(o)			/* convert mirror type */
register OBJREC	*o;
{
	register MATREC	*m;

	if (o->oargs.nfargs != 3)
		objerror(o, USER, "bad # real arguments");
	m = newmaterial(o->oname);
	m->type = o->otype;
	setcolor(m->u.m.ambdiff, 0., 0., 0.);
	setcolor(m->u.m.specular, o->oargs.farg[0],
			o->oargs.farg[1], o->oargs.farg[2]);
	m->u.m.specexp = MAXSPECEXP;
}


int
m_prism(o)			/* convert prism type */
register OBJREC	*o;
{
	register MATREC	*m;
					/* can't really deal with this type */
	m = newmaterial(o->oname);
	m->type = o->otype;
	setcolor(m->u.m.ambdiff, 0.2, 0.2, 0.2);
	setcolor(m->u.m.specular, 0.1, 0.1, 0.1);
	m->u.m.specexp = UNKSPECEXP;
}
