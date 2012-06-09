#ifndef lint
static const char	RCSid[] = "$Id: m_direct.c,v 2.14 2012/06/09 07:16:47 greg Exp $";
#endif
/*
 * Routines for light-redirecting materials and
 *   their associated virtual light sources
 */

#include "copyright.h"

#include  "ray.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "source.h"
#include  "func.h"

/*
 * The arguments for MAT_DIRECT1 are:
 *
 *	5+ coef1 dx1 dy1 dz1 funcfile transform..
 *	0
 *	n A1 A2 .. An
 *
 * The arguments for MAT_DIRECT2 are:
 *
 *	9+ coef1 dx1 dy1 dz1 coef2 dx2 dy2 dz2 funcfile transform..
 *	0
 *	n A1 A2 .. An
 */

static int redirect(OBJREC  *m, RAY  *r, int  n);
static int dir_proj(MAT4  pm, OBJREC  *o, SRCREC  *s, int  n);

VSMATERIAL  direct1_vs = {dir_proj, 1};
VSMATERIAL  direct2_vs = {dir_proj, 2};

#define getdfunc(m)	( (m)->otype == MAT_DIRECT1 ? \
				getfunc(m, 4, 0xf, 1) : \
				getfunc(m, 8, 0xff, 1) )


extern int
m_direct(			/* shade redirected ray */
	register OBJREC  *m,
	register RAY  *r
)
{
					/* check if source ray */
	if (r->rsrc >= 0 && source[r->rsrc].so != r->ro)
		return(1);			/* got the wrong guy */
					/* compute first projection */
	if (m->otype == MAT_DIRECT1 ||
			(r->rsrc < 0 || source[r->rsrc].sa.sv.pn == 0))
		redirect(m, r, 0);
					/* compute second projection */
	if (m->otype == MAT_DIRECT2 &&
			(r->rsrc < 0 || source[r->rsrc].sa.sv.pn == 1))
		redirect(m, r, 1);
	return(1);
}


static int
redirect(		/* compute n'th ray redirection */
	OBJREC  *m,
	RAY  *r,
	int  n
)
{
	MFUNC  *mf;
	register EPNODE  **va;
	FVECT  nsdir;
	RAY  nr;
	double  coef;
	register int  j;
					/* set up function */
	mf = getdfunc(m);
	setfunc(m, r);
					/* assign direction variable */
	if (r->rsrc >= 0) {
		register SRCREC  *sp = source + source[r->rsrc].sa.sv.sn;

		if (sp->sflags & SDISTANT)
			VCOPY(nsdir, sp->sloc);
		else {
			for (j = 0; j < 3; j++)
				nsdir[j] = sp->sloc[j] - r->rop[j];
			normalize(nsdir);
		}
		multv3(nsdir, nsdir, funcxf.xfm);
		varset("DxA", '=', nsdir[0]/funcxf.sca);
		varset("DyA", '=', nsdir[1]/funcxf.sca);
		varset("DzA", '=', nsdir[2]/funcxf.sca);
	} else {
		varset("DxA", '=', 0.0);
		varset("DyA", '=', 0.0);
		varset("DzA", '=', 0.0);
	}
					/* compute coefficient */
	errno = 0;
	va = mf->ep + 4*n;
	coef = evalue(va[0]);
	if ((errno == EDOM) | (errno == ERANGE))
		goto computerr;
	setcolor(nr.rcoef, coef, coef, coef);
	if (rayorigin(&nr, TRANS, r, nr.rcoef) < 0)
		return(0);
	va++;				/* compute direction */
	for (j = 0; j < 3; j++) {
		nr.rdir[j] = evalue(va[j]);
		if (errno == EDOM || errno == ERANGE)
			goto computerr;
	}
	if (mf->fxp != &unitxf)
		multv3(nr.rdir, nr.rdir, mf->fxp->xfm);
	if (r->rox != NULL)
		multv3(nr.rdir, nr.rdir, r->rox->f.xfm);
	if (normalize(nr.rdir) == 0.0)
		goto computerr;
					/* compute value */
	if (r->rsrc >= 0)
		nr.rsrc = source[r->rsrc].sa.sv.sn;
	rayvalue(&nr);
	multcolor(nr.rcol, nr.rcoef);
	addcolor(r->rcol, nr.rcol);
	if (r->ro != NULL && isflat(r->ro->otype))
		r->rt = r->rot + nr.rt;
	return(1);
computerr:
	objerror(m, WARNING, "compute error");
	return(-1);
}


static int
dir_proj(		/* compute a director's projection */
	MAT4  pm,
	OBJREC  *o,
	SRCREC  *s,
	int  n
)
{
	RAY  tr;
	OBJREC  *m;
	MFUNC  *mf;
	EPNODE  **va;
	FVECT  cent, newdir, nv, h;
	double  coef, olddot, newdot, od;
	register int  i, j;
				/* initialize test ray */
	getmaxdisk(cent, o);
	if (s->sflags & SDISTANT)
		for (i = 0; i < 3; i++) {
			tr.rdir[i] = -s->sloc[i];
			tr.rorg[i] = cent[i] - tr.rdir[i];
		}
	else {
		for (i = 0; i < 3; i++) {
			tr.rdir[i] = cent[i] - s->sloc[i];
			tr.rorg[i] = s->sloc[i];
		}
		if (normalize(tr.rdir) == 0.0)
			return(0);		/* at source! */
	}
	od = getplaneq(nv, o);
	olddot = DOT(tr.rdir, nv);
	if (olddot <= FTINY && olddot >= -FTINY)
		return(0);		/* old dir parallels plane */
	tr.rmax = 0.0;
	rayorigin(&tr, PRIMARY, NULL, NULL);
	if (!(*ofun[o->otype].funp)(o, &tr))
		return(0);		/* no intersection! */
				/* compute redirection */
	m = vsmaterial(o);
	mf = getdfunc(m);
	setfunc(m, &tr);
	varset("DxA", '=', 0.0);
	varset("DyA", '=', 0.0);
	varset("DzA", '=', 0.0);
	errno = 0;
	va = mf->ep + 4*n;
	coef = evalue(va[0]);
	if (errno == EDOM || errno == ERANGE)
		goto computerr;
	if (coef <= FTINY)
		return(0);		/* insignificant */
	va++;
	for (i = 0; i < 3; i++) {
		newdir[i] = evalue(va[i]);
		if (errno == EDOM || errno == ERANGE)
			goto computerr;
	}
	if (mf->fxp != &unitxf)
		multv3(newdir, newdir, mf->fxp->xfm);
					/* normalization unnecessary */
	newdot = DOT(newdir, nv);
	if (newdot <= FTINY && newdot >= -FTINY)
		return(0);		/* new dir parallels plane */
				/* everything OK -- compute shear */
	for (i = 0; i < 3; i++)
		h[i] = newdir[i]/newdot - tr.rdir[i]/olddot;
	setident4(pm);
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 3; i++)
			pm[i][j] += nv[i]*h[j];
		pm[3][j] = -od*h[j];
	}
	if ((newdot > 0.0) ^ (olddot > 0.0))	/* add mirroring */
		for (j = 0; j < 3; j++) {
			for (i = 0; i < 3; i++)
				pm[i][j] -= 2.*nv[i]*nv[j];
			pm[3][j] += 2.*od*nv[j];
		}
	return(1);
computerr:
	objerror(m, WARNING, "projection compute error");
	return(0);
}
