#ifndef lint
static const char	RCSid[] = "$Id: sphere.c,v 2.3 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 *  sphere.c - compute ray intersection with spheres.
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


o_sphere(so, r)			/* compute intersection with sphere */
OBJREC  *so;
register RAY  *r;
{
	double  a, b, c;	/* coefficients for quadratic equation */
	double  root[2];	/* quadratic roots */
	int  nroots;
	double  t;
	register FLOAT  *ap;
	register int  i;

	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad # arguments");
	ap = so->oargs.farg;
	if (ap[3] < -FTINY) {
		objerror(so, WARNING, "negative radius");
		so->otype = so->otype == OBJ_SPHERE ?
				OBJ_BUBBLE : OBJ_SPHERE;
		ap[3] = -ap[3];
	} else if (ap[3] <= FTINY)
		objerror(so, USER, "zero radius");

	/*
	 *	We compute the intersection by substituting into
	 *  the surface equation for the sphere.  The resulting
	 *  quadratic equation in t is then solved for the
	 *  smallest positive root, which is our point of
	 *  intersection.
	 *	Since the ray is normalized, a should always be
	 *  one.  We compute it here to prevent instability in the
	 *  intersection calculation.
	 */
				/* compute quadratic coefficients */
	a = b = c = 0.0;
	for (i = 0; i < 3; i++) {
		a += r->rdir[i]*r->rdir[i];
		t = r->rorg[i] - ap[i];
		b += 2.0*r->rdir[i]*t;
		c += t*t;
	}
	c -= ap[3] * ap[3];

	nroots = quadratic(root, a, b, c);	/* solve quadratic */
	
	for (i = 0; i < nroots; i++)		/* get smallest positive */
		if ((t = root[i]) > FTINY)
			break;
	if (i >= nroots)
		return(0);			/* no positive root */

	if (t >= r->rot)
		return(0);			/* other is closer */

	r->ro = so;
	r->rot = t;
					/* compute normal */
	a = ap[3];
	if (so->otype == OBJ_BUBBLE)
		a = -a;			/* reverse */
	for (i = 0; i < 3; i++) {
		r->rop[i] = r->rorg[i] + r->rdir[i]*t;
		r->ron[i] = (r->rop[i] - ap[i]) / a;
	}
	r->rod = -DOT(r->rdir, r->ron);
	r->rox = NULL;

	return(1);			/* hit */
}
