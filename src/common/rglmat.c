#ifndef lint
static const char	RCSid[] = "$Id: rglmat.c,v 3.5 2003/09/16 06:31:48 greg Exp $";
#endif
/*
 * Routines for Radiance -> OpenGL materials.
 */

#include "copyright.h"

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
		freemtl(lup->data);
	if ((lup->data = o->os) != NULL)	/* make material reference */
		((MATREC *)lup->data)->nlinks++;
	return(0);
memerr:
	error(SYSTEM, "out of memory in o_default");
	return(0);
}


int
o_unsupported(o)		/* unsupported object primitive */
OBJREC	*o;
{
	objerror(o, WARNING, "unsupported type");
	return(0);
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
		freemtl(lup->data);
	lup->data = (char *)malloc(sizeof(MATREC));
	if (lup->data == NULL)
		goto memerr;
	((MATREC *)lup->data)->nlinks = 1;
	return((MATREC *)lup->data);
memerr:
	error(SYSTEM, "out of memory in newmaterial");
}


void
freemtl(p)			/* free a material */
void	*p;
{
	register MATREC	*mp = (MATREC *)p;

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
	return(0);
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
	return(0);
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
	return(0);
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
	return(0);
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
	return(0);
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
	return(0);
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
	return(0);
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
	return(0);
}
