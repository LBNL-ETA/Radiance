#ifndef lint
static const char	RCSid[] = "$Id: o_cone.c,v 2.2 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 *  o_cone.c - routine to determine ray intersection with cones.
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

#include  "cone.h"


o_cone(o, r)			/* intersect ray with cone */
OBJREC  *o;
register RAY  *r;
{
	FVECT  rox, rdx;
	double  a, b, c;
	double  root[2];
	int  nroots, rn;
	register CONE  *co;
	register int  i;

						/* get cone structure */
	co = getcone(o, 1);

	/*
	 *     To intersect a ray with a cone, we transform the
	 *  ray into the cone's normalized space.  This greatly
	 *  simplifies the computation.
	 *     For a cone or cup, normalization results in the
	 *  equation:
	 *
	 *		x*x + y*y - z*z == 0
	 *
	 *     For a cylinder or tube, the normalized equation is:
	 *
	 *		x*x + y*y - r*r == 0
	 *
	 *     A normalized ring obeys the following set of equations:
	 *
	 *		z == 0			&&
	 *		x*x + y*y >= r0*r0	&&
	 *		x*x + y*y <= r1*r1
	 */

					/* transform ray */
	multp3(rox, r->rorg, co->tm);
	multv3(rdx, r->rdir, co->tm);
					/* compute intersection */

	if (o->otype == OBJ_CONE || o->otype == OBJ_CUP) {

		a = rdx[0]*rdx[0] + rdx[1]*rdx[1] - rdx[2]*rdx[2];
		b = 2.0*(rdx[0]*rox[0] + rdx[1]*rox[1] - rdx[2]*rox[2]);
		c = rox[0]*rox[0] + rox[1]*rox[1] - rox[2]*rox[2];

	} else if (o->otype == OBJ_CYLINDER || o->otype == OBJ_TUBE) {

		a = rdx[0]*rdx[0] + rdx[1]*rdx[1];
		b = 2.0*(rdx[0]*rox[0] + rdx[1]*rox[1]);
		c = rox[0]*rox[0] + rox[1]*rox[1] - CO_R0(co)*CO_R0(co);

	} else { /* OBJ_RING */

		if (rdx[2] <= FTINY && rdx[2] >= -FTINY)
			return(0);			/* parallel */
		root[0] = -rox[2]/rdx[2];
		if (root[0] <= FTINY || root[0] >= r->rot)
			return(0);			/* distance check */
		b = root[0]*rdx[0] + rox[0];
		c = root[0]*rdx[1] + rox[1];
		a = b*b + c*c;
		if (a < CO_R0(co)*CO_R0(co) || a > CO_R1(co)*CO_R1(co))
			return(0);			/* outside radii */
		r->ro = o;
		r->rot = root[0];
		for (i = 0; i < 3; i++)
			r->rop[i] = r->rorg[i] + r->rdir[i]*r->rot;
		VCOPY(r->ron, co->ad);
		r->rod = -rdx[2];
		r->rox = NULL;
		return(1);				/* good */
	}
					/* roots for cone, cup, cyl., tube */
	nroots = quadratic(root, a, b, c);

	for (rn = 0; rn < nroots; rn++) {	/* check real roots */
		if (root[rn] <= FTINY)
			continue;		/* too small */
		if (root[rn] >= r->rot)
			break;			/* too big */
						/* check endpoints */
		for (i = 0; i < 3; i++) {
			rox[i] = r->rorg[i] + root[rn]*r->rdir[i];
			rdx[i] = rox[i] - CO_P0(co)[i];
		}
		b = DOT(rdx, co->ad); 
		if (b < 0.0)
			continue;		/* before p0 */
		if (b > co->al)
			continue;		/* after p1 */
		r->ro = o;
		r->rot = root[rn];
		VCOPY(r->rop, rox);
						/* get normal */
		if (o->otype == OBJ_CYLINDER)
			a = CO_R0(co);
		else if (o->otype == OBJ_TUBE)
			a = -CO_R0(co);
		else { /* OBJ_CONE || OBJ_CUP */
			c = CO_R1(co) - CO_R0(co);
			a = CO_R0(co) + b*c/co->al;
			if (o->otype == OBJ_CUP) {
				c = -c;
				a = -a;
			}
		}
		for (i = 0; i < 3; i++)
			r->ron[i] = (rdx[i] - b*co->ad[i])/a;
		if (o->otype == OBJ_CONE || o->otype == OBJ_CUP)
			for (i = 0; i < 3; i++)
				r->ron[i] = (co->al*r->ron[i] - c*co->ad[i])
						/co->sl;
		r->rod = -DOT(r->rdir, r->ron);
		r->rox = NULL;
		return(1);			/* good */
	}
	return(0);
}
