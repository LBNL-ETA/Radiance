#ifndef lint
static const char	RCSid[] = "$Id: glass.c,v 2.11 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 *  glass.c - simpler shading function for thin glass surfaces.
 */

#include "copyright.h"

#include  "ray.h"

#include  "otypes.h"

/*
 *  This definition of glass provides for a quick calculation
 *  using a single surface where two closely spaced parallel
 *  dielectric surfaces would otherwise be used.  The chief
 *  advantage to using this material is speed, since internal
 *  reflections are avoided.
 *
 *  The specification for glass is as follows:
 *
 *	modifier glass id
 *	0
 *	0
 *	3+ red grn blu [refractive_index]
 *
 *  The color is used for the transmission at normal incidence.
 *  To compute transmissivity (tn) from transmittance (Tn) use:
 *
 *	tn = (sqrt(.8402528435+.0072522239*Tn*Tn)-.9166530661)/.0036261119/Tn
 *
 *  The transmissivity of standard 88% transmittance glass is 0.96.
 *  A refractive index other than the default can be used by giving
 *  it as the fourth real argument.  The above formula no longer applies.
 *
 *  If we appear to hit the back side of the surface, then we
 *  turn the normal around.
 */

#define  RINDEX		1.52		/* refractive index of glass */


m_glass(m, r)		/* color a ray which hit a thin glass surface */
OBJREC  *m;
register RAY  *r;
{
	COLOR  mcolor;
	double  pdot;
	FVECT  pnorm;
	double  rindex, cos2;
	COLOR  trans, refl;
	int  hastexture;
	double  d, r1e, r1m;
	double  transtest, transdist;
	double  mirtest, mirdist;
	RAY  p;
	register int  i;
						/* check arguments */
	if (m->oargs.nfargs == 3)
		rindex = RINDEX;		/* default value of n */
	else if (m->oargs.nfargs == 4)
		rindex = m->oargs.farg[3];	/* use their value */
	else
		objerror(m, USER, "bad arguments");

	setcolor(mcolor, m->oargs.farg[0], m->oargs.farg[1], m->oargs.farg[2]);

	if (r->rod < 0.0)			/* reorient if necessary */
		flipsurface(r);
	mirtest = transtest = 0;
	mirdist = transdist = r->rot;
						/* get modifiers */
	raytexture(r, m->omod);
	if (hastexture = DOT(r->pert,r->pert) > FTINY*FTINY)
		pdot = raynormal(pnorm, r);
	else {
		VCOPY(pnorm, r->ron);
		pdot = r->rod;
	}
						/* angular transmission */
	cos2 = sqrt( (1.0-1.0/(rindex*rindex)) +
		     pdot*pdot/(rindex*rindex) );
	setcolor(mcolor, pow(colval(mcolor,RED), 1.0/cos2),
			 pow(colval(mcolor,GRN), 1.0/cos2),
			 pow(colval(mcolor,BLU), 1.0/cos2));

						/* compute reflection */
	r1e = (pdot - rindex*cos2) / (pdot + rindex*cos2);
	r1e *= r1e;
	r1m = (1.0/pdot - rindex/cos2) / (1.0/pdot + rindex/cos2);
	r1m *= r1m;
						/* compute transmittance */
	for (i = 0; i < 3; i++) {
		d = colval(mcolor, i);
		colval(trans,i) = .5*(1.0-r1e)*(1.0-r1e)*d/(1.0-r1e*r1e*d*d);
		colval(trans,i) += .5*(1.0-r1m)*(1.0-r1m)*d/(1.0-r1m*r1m*d*d);
	}
						/* transmitted ray */
	if (rayorigin(&p, r, TRANS, bright(trans)) == 0) {
		if (!(r->crtype & SHADOW) && hastexture) {
			for (i = 0; i < 3; i++)		/* perturb direction */
				p.rdir[i] = r->rdir[i] +
						2.*(1.-rindex)*r->pert[i];
			if (normalize(p.rdir) == 0.0) {
				objerror(m, WARNING, "bad perturbation");
				VCOPY(p.rdir, r->rdir);
			}
		} else {
			VCOPY(p.rdir, r->rdir);
			transtest = 2;
		}
		rayvalue(&p);
		multcolor(p.rcol, r->pcol);	/* modify */
		multcolor(p.rcol, trans);
		addcolor(r->rcol, p.rcol);
		transtest *= bright(p.rcol);
		transdist = r->rot + p.rt;
	}

	if (r->crtype & SHADOW) {		/* skip reflected ray */
		r->rt = transdist;
		return(1);
	}
						/* compute reflectance */
	for (i = 0; i < 3; i++) {
		d = colval(mcolor, i);
		d *= d;
		colval(refl,i) = .5*r1e*(1.0+(1.0-2.0*r1e)*d)/(1.0-r1e*r1e*d);
		colval(refl,i) += .5*r1m*(1.0+(1.0-2.0*r1m)*d)/(1.0-r1m*r1m*d);
	}
						/* reflected ray */
	if (rayorigin(&p, r, REFLECTED, bright(refl)) == 0) {
		for (i = 0; i < 3; i++)
			p.rdir[i] = r->rdir[i] + 2.0*pdot*pnorm[i];
		rayvalue(&p);
		multcolor(p.rcol, refl);
		addcolor(r->rcol, p.rcol);
		if (!hastexture && r->ro != NULL && isflat(r->ro->otype)) {
			mirtest = 2.0*bright(p.rcol);
			mirdist = r->rot + p.rt;
		}
	}
                                        /* check distance */
	d = bright(r->rcol);
	if (transtest > d)
		r->rt = transdist;
	else if (mirtest > d)
		r->rt = mirdist;
	return(1);
}
