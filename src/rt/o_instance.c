#ifndef lint
static const char	RCSid[] = "$Id: o_instance.c,v 2.4 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 *  o_instance.c - routines for computing ray intersections with octrees.
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

#include  "instance.h"


o_instance(o, r)		/* compute ray intersection with octree */
OBJREC  *o;
register RAY  *r;
{
	RAY  rcont;
	register INSTANCE  *in;
	register int  i;
					/* get the octree */
	in = getinstance(o, IO_ALL);
					/* copy and transform ray */
	copystruct(&rcont, r);
	multp3(rcont.rorg, r->rorg, in->x.b.xfm);
	multv3(rcont.rdir, r->rdir, in->x.b.xfm);
	for (i = 0; i < 3; i++)
		rcont.rdir[i] /= in->x.b.sca;
	rcont.rmax *= in->x.b.sca;
					/* clear and trace it */
	rayclear(&rcont);
	if (!localhit(&rcont, &in->obj->scube))
		return(0);			/* missed */
	if (rcont.rot * in->x.f.sca >= r->rot)
		return(0);			/* not close enough */

	if (o->omod != OVOID) {		/* if we have modifier, use it */
		r->ro = o;
		r->rox = NULL;
	} else {			/* else use theirs */
		r->ro = rcont.ro;
		if (rcont.rox != NULL) {
			newrayxf(r);		/* allocate transformation */
					/* NOTE: r->rox may equal rcont.rox! */
			multmat4(r->rox->f.xfm, rcont.rox->f.xfm, in->x.f.xfm);
			r->rox->f.sca = rcont.rox->f.sca * in->x.f.sca;
			multmat4(r->rox->b.xfm, in->x.b.xfm, rcont.rox->b.xfm);
			r->rox->b.sca = in->x.b.sca * rcont.rox->b.sca;
		} else
			r->rox = &in->x;
	}
					/* transform it back */
	r->rot = rcont.rot * in->x.f.sca;
	multp3(r->rop, rcont.rop, in->x.f.xfm);
	multv3(r->ron, rcont.ron, in->x.f.xfm);
	for (i = 0; i < 3; i++)
		r->ron[i] /= in->x.f.sca;
	r->rod = rcont.rod;
					/* return hit */
	return(1);
}
