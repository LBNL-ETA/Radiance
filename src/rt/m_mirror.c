#ifndef lint
static const char	RCSid[] = "$Id: m_mirror.c,v 2.22 2020/01/21 21:31:58 greg Exp $";
#endif
/*
 * Routines for mirror material supporting virtual light sources
 */

#include "copyright.h"

#include  "ray.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "source.h"

/*
 * The real arguments for MAT_MIRROR are simply:
 *
 *	3 rrefl grefl brefl
 *
 * Additionally, the user may specify a single string argument
 * which is interpreted as the name of the material to use
 * instead of the mirror if the ray being considered is not
 * part of the direct calculation.
 */

static int mir_proj(MAT4  pm, OBJREC  *o, SRCREC  *s, int  n);
static void mirrorproj(MAT4  m, FVECT  nv, double  offs);

VSMATERIAL  mirror_vs = {mir_proj, 1};


int
m_mirror(			/* shade mirrored ray */
	OBJREC  *m,
	RAY  *r
)
{
	COLOR  mcolor;
	RAY  nr;
	int  rpure = 1;
	int  i;
					/* check arguments */
	if (m->oargs.nfargs != 3 || m->oargs.nsargs > 1)
		objerror(m, USER, "bad number of arguments");
					/* check for substitute material */
					/* but avoid double-counting */
	if ( m->oargs.nsargs > 0 &&
			(r->rsrc < 0 || source[r->rsrc].so != r->ro) &&
			!(r->crtype & (AMBIENT|SPECULAR) && r->rod > 0.) ) {
		if (!strcmp(m->oargs.sarg[0], VOIDID)) {
			raytrans(r);
			return(1);
		}
		return(rayshade(r, lastmod(objndx(m), m->oargs.sarg[0])));
	}
					/* check for bad source ray */
	if (r->rsrc >= 0 && source[r->rsrc].so != r->ro)
		return(1);

	if (r->rod < 0.) {		/* back is black */
		if (!backvis)
			raytrans(r);	/* unless back visibility is off */
		return(1);
	}
					/* get modifiers */
	raytexture(r, m->omod);
					/* assign material color */
	setcolor(mcolor, m->oargs.farg[0],
			m->oargs.farg[1],
			m->oargs.farg[2]);
	multcolor(mcolor, r->pcol);
					/* compute reflected ray */
	if (r->rsrc >= 0) {			/* relayed light source */
		rayorigin(&nr, REFLECTED, r, mcolor);
					/* ignore textures */
		for (i = 0; i < 3; i++)
			nr.rdir[i] = r->rdir[i] + 2.*r->rod*r->ron[i];
					/* source we're aiming for next */
		nr.rsrc = source[r->rsrc].sa.sv.sn;
	} else {				/* ordinary reflection */
		FVECT  pnorm;
		double  pdot;

		if (rayorigin(&nr, REFLECTED, r, mcolor) < 0)
			return(1);
		if (!(r->crtype & AMBIENT) &&
				DOT(r->pert,r->pert) > FTINY*FTINY) {
			pdot = raynormal(pnorm, r);	/* use textures */
			for (i = 0; i < 3; i++)
				nr.rdir[i] = r->rdir[i] + 2.*pdot*pnorm[i];
			rpure = 0;
		}
						/* check for penetration */
		if (rpure || DOT(nr.rdir, r->ron) <= FTINY)
			for (i = 0; i < 3; i++)
				nr.rdir[i] = r->rdir[i] + 2.*r->rod*r->ron[i];
	}
	checknorm(nr.rdir);
	rayvalue(&nr);
	multcolor(nr.rcol, nr.rcoef);
	copycolor(r->mcol, nr.rcol);
	addcolor(r->rcol, nr.rcol);
	r->rmt = r->rot;
	if (rpure && r->ro != NULL && isflat(r->ro->otype))
		r->rmt += raydistance(&nr);
	return(1);
}


static int
mir_proj(		/* compute a mirror's projection */
	MAT4  pm,
	OBJREC  *o,
	SRCREC  *s,
	int  n
)
{
	double	corr = 1.;
	FVECT  nv, sc;
	double  od, offs;
	int  i;
				/* get surface normal and offset */
	offs = od = getplaneq(nv, o);
	if (s->sflags & SDISTANT)
		offs = 0.;
				/* check for extreme point behind */
	if (s->sflags & SCIR) {
		if (s->sflags & (SFLAT|SDISTANT))
			corr = 1.12837917;	/* correct setflatss() */
		else
			corr = 1.0/0.7236;	/* correct sphsetsrc() */
	}
	VCOPY(sc, s->sloc);
	for (i = s->sflags & SFLAT ? SV : SW; i >= 0; i--)
		if (DOT(nv, s->ss[i]) > offs)
			VSUM(sc, sc, s->ss[i], corr);
		else
			VSUM(sc, sc, s->ss[i], -corr);
	if (DOT(sc, nv) <= offs+FTINY)
		return(0);
				/* everything OK -- compute projection */
	mirrorproj(pm, nv, od);
	return(1);
}


static void
mirrorproj(		/* get mirror projection for surface */
	MAT4  m,
	FVECT  nv,
	double  offs
)
{
	int  i, j;
					/* assign matrix */
	setident4(m);
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 3; i++)
			m[i][j] -= 2.*nv[i]*nv[j];
		m[3][j] = 2.*offs*nv[j];
	}
}
