#ifndef lint
static const char	RCSid[] = "$Id: m_clip.c,v 2.5 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  m_clip.c - routine for clipped (cut) objects.
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

/*
 *  Clipping objects permit holes and sections to be taken out
 *  of other objects.  The method is simple:  
 *
 *  The argument is the clipped materials;
 *  the first is used to shade upon exit.
 */


int
m_clip(m, r)			/* clip objects from ray */
register OBJREC  *m;
register RAY  *r;
{
	OBJECT  cset[MAXSET+1], *modset;
	OBJECT  obj, mod;
	int  entering;
	register int  i;

	obj = objndx(m);
	if ((modset = (OBJECT *)m->os) == NULL) {
		if (m->oargs.nsargs < 1 || m->oargs.nsargs > MAXSET)
			objerror(m, USER, "bad # arguments");
		modset = (OBJECT *)malloc((m->oargs.nsargs+1)*sizeof(OBJECT));
		if (modset == NULL)
			error(SYSTEM, "out of memory in m_clip");
		modset[0] = 0;
		for (i = 0; i < m->oargs.nsargs; i++) {
			if (!strcmp(m->oargs.sarg[i], VOIDID))
				continue;
			if ((mod = lastmod(obj, m->oargs.sarg[i])) == OVOID) {
				sprintf(errmsg, "unknown modifier \"%s\"",
						m->oargs.sarg[i]);
				objerror(m, WARNING, errmsg);
				continue;
			}
			if (inset(modset, mod)) {
				objerror(m, WARNING, "duplicate modifier");
				continue;
			}
			insertelem(modset, mod);
		}
		m->os = (char *)modset;
	}
	if (r->clipset != NULL)
		setcopy(cset, r->clipset);
	else
		cset[0] = 0;

	entering = r->rod > 0.0;		/* entering clipped region? */

	for (i = modset[0]; i > 0; i--) {
		if (entering) {
			if (!inset(cset, modset[i])) {
				if (cset[0] >= MAXSET)
					error(INTERNAL, "set overflow in m_clip");
				insertelem(cset, modset[i]);
			}
		} else {
			if (inset(cset, modset[i]))
				deletelem(cset, modset[i]);
		}
	}
					/* compute ray value */
	r->newcset = cset;
	if (strcmp(m->oargs.sarg[0], VOIDID)) {
		int  inside = 0;
		register RAY  *rp;
					/* check for penetration */
		for (rp = r; rp->parent != NULL; rp = rp->parent)
			if (!(rp->rtype & RAYREFL) && rp->parent->ro != NULL
					&& inset(modset, rp->parent->ro->omod))
				if (rp->parent->rod > 0.0)
					inside++;
				else
					inside--;
		if (inside > 0) {	/* we just hit the object */
			flipsurface(r);
			return(rayshade(r, lastmod(obj, m->oargs.sarg[0])));
		}
	}
	raytrans(r);			/* else transfer ray */
	return(1);
}
