#ifndef lint
static const char	RCSid[] = "$Id: o_face.c,v 2.2 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 *  o_face.c - compute ray intersection with faces.
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

#include  "face.h"


o_face(o, r)		/* compute intersection with polygonal face */
OBJREC  *o;
register RAY  *r;
{
	double  rdot;		/* direction . normal */
	double  t;		/* distance to intersection */
	FVECT  pisect;		/* intersection point */
	register FACE  *f;	/* face record */
	register int  i;

	f = getface(o);
		
	/*
	 *  First, we find the distance to the plane containing the
	 *  face.  If this distance is less than zero or greater
	 *  than a previous intersection, we return.  Otherwise,
	 *  we determine whether in fact the ray intersects the
	 *  face.  The ray intersects the face if the
	 *  point of intersection with the plane of the face
	 *  is inside the face.
	 */
						/* compute dist. to plane */
	rdot = -DOT(r->rdir, f->norm);
	if (rdot <= FTINY && rdot >= -FTINY)	/* ray parallels plane */
		t = FHUGE;
	else
		t = (DOT(r->rorg, f->norm) - f->offset) / rdot;
	
	if (t <= FTINY || t >= r->rot)		/* not good enough */
		return(0);
						/* compute intersection */
	for (i = 0; i < 3; i++)
		pisect[i] = r->rorg[i] + r->rdir[i]*t;

	if (!inface(pisect, f))			/* ray intersects face? */
		return(0);

	r->ro = o;
	r->rot = t;
	VCOPY(r->rop, pisect);
	VCOPY(r->ron, f->norm);
	r->rod = rdot;
	r->rox = NULL;

	return(1);				/* hit */
}
