#ifndef lint
static const char	RCSid[] = "$Id: glass.c,v 2.10 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  glass.c - simpler shading function for thin glass surfaces.
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
