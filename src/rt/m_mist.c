#ifndef lint
static const char	RCSid[] = "$Id: m_mist.c,v 2.10 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 * Mist volumetric material.
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

#include  "source.h"

/*
 *  A mist volume is used to specify a region in the scene where a certain
 *  light source (or sources) is going to contribute to scattering.  The
 *  material can add to the existing global medium, and override any ray
 *  settings for scattering albedo and eccentricity.  Overlapping mist
 *  regions should agree w.r.t. albedo and eccentricity, and
 *  should have disjoint source lists.
 *
 *  A pattern, if used, should compute the line integral of extinction,
 *  and will modify the first three arguments directly.  This will tend
 *  to invalidate results when there are other objects intersected within
 *  the mist region.
 *
 *  The string arguments for MAT_MIST are the identifiers for the important
 *  light sources, which will be looked up in the source array.  The last
 *  source found matching a name is the one used.  A relayed light source
 *  may be indicated by the relay surface name, followed by a '>' character,
 *  followed by the relayed source name (which may be another relay).
 *
 *  Up to five real arguments may be given for MAT_MIST:
 *
 *	[ext_r  ext_g  ext_b  [albedo_r albedo_g albedo_b [gecc]]]
 *
 *  The primaries indicate medium extinction per unit length (absorption
 *  plus scattering), which is added to the global extinction coefficient, set
 *  by the -me option.  The albedo is the ratio of scattering to extinction,
 *  and is set globally by the -ma option (salbedo) and overridden here.
 *  The Heyney-Greenstein eccentricity parameter (-mg seccg) indicates how much
 *  scattering favors the forward direction.  A value of 0 means isotropic
 *  scattering.  A value approaching 1 indicates strong forward scattering.
 */

#ifndef  MAXSLIST
#define  MAXSLIST	32	/* maximum sources to check */
#endif

#define RELAYDELIM	'>'		/* relay delimiter character */


static int
inslist(sl, n)		/* return index of source n if it's in list sl */
register int  *sl;
register int  n;
{
	register int  i;

	for (i = sl[0]; i > 0; i--)
		if (sl[i] == n)
			return(i);
	return(0);
}


static int
srcmatch(sp, id)	/* check for an id match on a light source */
register SRCREC  *sp;
register char  *id;
{
	extern char  *index();
	register char  *cp;
						/* check for relay sources */
	while ((cp = index(id, RELAYDELIM)) != NULL) {
		if (!(sp->sflags & SVIRTUAL) || sp->so == NULL)
			return(0);
		if (strncmp(id, sp->so->oname, cp-id) || sp->so->oname[cp-id])
			return(0);
		id = cp + 1;				/* relay to next */
		sp = source + sp->sa.sv.sn;
	}
	if (sp->sflags & SVIRTUAL || sp->so == NULL)
		return(0);
	return(!strcmp(id, sp->so->oname));
}


static void
add2slist(r, sl)	/* add source list to ray's */
register RAY  *r;
register int  *sl;
{
	static int  slspare[MAXSLIST+1];	/* in case of emergence */
	register int  i;

	if (sl == NULL || sl[0] == 0)		/* nothing to add */
		return;
	if (r->slights == NULL)
		(r->slights = slspare)[0] = 0;	/* just once per ray path */
	for (i = sl[0]; i > 0; i--)
		if (!inslist(r->slights, sl[i])) {
			if (r->slights[0] >= MAXSLIST)
				error(INTERNAL,
					"scattering source list overflow");
			r->slights[++r->slights[0]] = sl[i];
		}
}


int
m_mist(m, r)		/* process a ray entering or leaving some mist */
OBJREC  *m;
register RAY  *r;
{
	RAY  p;
	int  *myslist = NULL;
	int  newslist[MAXSLIST+1];
	COLOR  mext;
	double  re, ge, be;
	register int  i, j;
					/* check arguments */
	if (m->oargs.nfargs > 7)
		objerror(m, USER, "bad arguments");
					/* get source indices */
	if (m->oargs.nsargs > 0 && (myslist = (int *)m->os) == NULL) {
		if (m->oargs.nsargs > MAXSLIST)
			objerror(m, INTERNAL, "too many sources in list");
		myslist = (int *)malloc((m->oargs.nsargs+1)*sizeof(int));
		if (myslist == NULL)
			goto memerr;
		myslist[0] = 0;			/* size is first in list */
		for (j = 0; j < m->oargs.nsargs; j++) {
			i = nsources;		/* look up each source id */
			while (i--)
				if (srcmatch(source+i, m->oargs.sarg[j]))
					break;
			if (i < 0) {
				sprintf(errmsg, "unknown source \"%s\"",
						m->oargs.sarg[j]);
				objerror(m, WARNING, errmsg);
			} else if (inslist(myslist, i)) {
				sprintf(errmsg, "duplicate source \"%s\"",
						m->oargs.sarg[j]);
				objerror(m, WARNING, errmsg);
			} else
				myslist[++myslist[0]] = i;
		}
		m->os = (char *)myslist;
	}
	if (m->oargs.nfargs > 2) {		/* compute extinction */
		setcolor(mext, m->oargs.farg[0], m->oargs.farg[1],
				m->oargs.farg[2]);
		raytexture(r, m->omod);			/* get modifiers */
		multcolor(mext, r->pcol);
	} else
		setcolor(mext, 0., 0., 0.);
						/* start transmitted ray */
	if (rayorigin(&p, r, TRANS, 1.) < 0)
		return(1);
	VCOPY(p.rdir, r->rdir);
	p.slights = newslist;
	if (r->slights != NULL)			/* copy old list if one */
		for (j = r->slights[0]; j >= 0; j--)
			p.slights[j] = r->slights[j];
	else
		p.slights[0] = 0;
	if (r->rod > 0.) {			/* entering ray */
		addcolor(p.cext, mext);
		if (m->oargs.nfargs > 5)
			setcolor(p.albedo, m->oargs.farg[3],
					m->oargs.farg[4], m->oargs.farg[5]);
		if (m->oargs.nfargs > 6)
			p.gecc = m->oargs.farg[6];
		add2slist(&p, myslist);			/* add to list */
	} else {				/* leaving ray */
		if (myslist != NULL) {			/* delete from list */
			for (j = myslist[0]; j > 0; j--)
				if (i = inslist(p.slights, myslist[j]))
					p.slights[i] = -1;
			for (i = 0, j = 1; j <= p.slights[0]; j++)
				if (p.slights[j] != -1)
					p.slights[++i] = p.slights[j];
			if (p.slights[0] - i < myslist[0]) {	/* fix old */
				addcolor(r->cext, mext);
				if (m->oargs.nfargs > 5)
					setcolor(r->albedo, m->oargs.farg[3],
					m->oargs.farg[4], m->oargs.farg[5]);
				if (m->oargs.nfargs > 6)
					r->gecc = m->oargs.farg[6];
				add2slist(r, myslist);
			}
			p.slights[0] = i;
		}
		if ((re = colval(r->cext,RED) - colval(mext,RED)) <
				colval(cextinction,RED))
			re = colval(cextinction,RED);
		if ((ge = colval(r->cext,GRN) - colval(mext,GRN)) <
				colval(cextinction,GRN))
			ge = colval(cextinction,GRN);
		if ((be = colval(r->cext,BLU) - colval(mext,BLU)) <
				colval(cextinction,BLU))
			be = colval(cextinction,BLU);
		setcolor(p.cext, re, ge, be);
		if (m->oargs.nfargs > 5)
			copycolor(p.albedo, salbedo);
		if (m->oargs.nfargs > 6)
			p.gecc = seccg;
	}
	rayvalue(&p);				/* calls rayparticipate() */
	copycolor(r->rcol, p.rcol);		/* return value */
	r->rt = r->rot + p.rt;
	return(1);
memerr:
	error(SYSTEM, "out of memory in m_mist");
}
