#ifndef lint
static const char	RCSid[] = "$Id: srcsamp.c,v 2.7 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 * Source sampling routines
 *
 *  External symbols declared in source.h
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

#include  "random.h"


static int  cyl_partit(), flt_partit();


double
nextssamp(r, si)		/* compute sample for source, rtn. distance */
register RAY  *r;		/* origin is read, direction is set */
register SRCINDEX  *si;		/* source index (modified to current) */
{
	int  cent[3], size[3], parr[2];
	FVECT  vpos;
	double  d;
	register int  i;
nextsample:
	while (++si->sp >= si->np) {	/* get next sample */
		if (++si->sn >= nsources)
			return(0.0);	/* no more */
		if (source[si->sn].sflags & SSKIP)
			si->np = 0;
		else if (srcsizerat <= FTINY)
			nopart(si, r);
		else {
			for (i = si->sn; source[i].sflags & SVIRTUAL;
					i = source[i].sa.sv.sn)
				;		/* partition source */
			(*sfun[source[i].so->otype].of->partit)(si, r);
		}
		si->sp = -1;
	}
					/* get partition */
	cent[0] = cent[1] = cent[2] = 0;
	size[0] = size[1] = size[2] = MAXSPART;
	parr[0] = 0; parr[1] = si->sp;
	if (!skipparts(cent, size, parr, si->spt))
		error(CONSISTENCY, "bad source partition in nextssamp");
					/* compute sample */
	if (dstrsrc > FTINY) {			/* jitter sample */
		dimlist[ndims] = si->sn + 8831;
		dimlist[ndims+1] = si->sp + 3109;
		d = urand(ilhash(dimlist,ndims+2)+samplendx);
		if (source[si->sn].sflags & SFLAT) {
			multisamp(vpos, 2, d);
			vpos[2] = 0.5;
		} else
			multisamp(vpos, 3, d);
		for (i = 0; i < 3; i++)
			vpos[i] = dstrsrc * (1. - 2.*vpos[i]) *
					(double)size[i]/MAXSPART;
	} else
		vpos[0] = vpos[1] = vpos[2] = 0.0;

	for (i = 0; i < 3; i++)
		vpos[i] += (double)cent[i]/MAXSPART;
					/* compute direction */
	for (i = 0; i < 3; i++)
		r->rdir[i] = source[si->sn].sloc[i] +
				vpos[SU]*source[si->sn].ss[SU][i] +
				vpos[SV]*source[si->sn].ss[SV][i] +
				vpos[SW]*source[si->sn].ss[SW][i];

	if (!(source[si->sn].sflags & SDISTANT))
		for (i = 0; i < 3; i++)
			r->rdir[i] -= r->rorg[i];
					/* compute distance */
	if ((d = normalize(r->rdir)) == 0.0)
		goto nextsample;		/* at source! */

					/* compute sample size */
	if (source[si->sn].sflags & SFLAT) {
		si->dom = sflatform(si->sn, r->rdir);
		si->dom *= size[SU]*size[SV]/(MAXSPART*(double)MAXSPART);
	} else if (source[si->sn].sflags & SCYL) {
		si->dom = scylform(si->sn, r->rdir);
		si->dom *= size[SU]/(double)MAXSPART;
	} else {
		si->dom = size[SU]*size[SV]*(double)size[SW] /
				(MAXSPART*MAXSPART*(double)MAXSPART) ;
	}
	if (source[si->sn].sflags & SDISTANT) {
		si->dom *= source[si->sn].ss2;
		return(FHUGE);
	}
	if (si->dom <= 1e-4)
		goto nextsample;		/* behind source? */
	si->dom *= source[si->sn].ss2/(d*d);
	return(d);		/* sample OK, return distance */
}


int
skipparts(ct, sz, pp, pt)		/* skip to requested partition */
int  ct[3], sz[3];		/* center and size of partition (returned) */
register int  pp[2];		/* current index, number to skip (modified) */
unsigned char  *pt;		/* partition array */
{
	register int  p;
					/* check this partition */
	p = spart(pt, pp[0]);
	pp[0]++;
	if (p == S0)			/* leaf partition */
		if (pp[1]) {
			pp[1]--;
			return(0);	/* not there yet */
		} else
			return(1);	/* we've arrived */
				/* else check lower */
	sz[p] >>= 1;
	ct[p] -= sz[p];
	if (skipparts(ct, sz, pp, pt))
		return(1);	/* return hit */
				/* else check upper */
	ct[p] += sz[p] << 1;
	if (skipparts(ct, sz, pp, pt))
		return(1);	/* return hit */
				/* else return to starting position */
	ct[p] -= sz[p];
	sz[p] <<= 1;
	return(0);		/* return miss */
}


void
nopart(si, r)			/* single source partition */
register SRCINDEX  *si;
RAY  *r;
{
	clrpart(si->spt);
	setpart(si->spt, 0, S0);
	si->np = 1;
}


void
cylpart(si, r)			/* partition a cylinder */
SRCINDEX  *si;
register RAY  *r;
{
	double  dist2, safedist2, dist2cent, rad2;
	FVECT  v;
	register SRCREC  *sp;
	int  pi;
					/* first check point location */
	clrpart(si->spt);
	sp = source + si->sn;
	rad2 = 1.365 * DOT(sp->ss[SV],sp->ss[SV]);
	v[0] = r->rorg[0] - sp->sloc[0];
	v[1] = r->rorg[1] - sp->sloc[1];
	v[2] = r->rorg[2] - sp->sloc[2];
	dist2 = DOT(v,sp->ss[SU]);
	safedist2 = DOT(sp->ss[SU],sp->ss[SU]);
	dist2 *= dist2 / safedist2;
	dist2cent = DOT(v,v);
	dist2 = dist2cent - dist2;
	if (dist2 <= rad2) {		/* point inside extended cylinder */
		si->np = 0;
		return;
	}
	safedist2 *= 4.*r->rweight*r->rweight/(srcsizerat*srcsizerat);
	if (dist2 <= 4.*rad2 ||		/* point too close to subdivide */
			dist2cent >= safedist2) {	/* or too far */
		setpart(si->spt, 0, S0);
		si->np = 1;
		return;
	}
	pi = 0;
	si->np = cyl_partit(r->rorg, si->spt, &pi, MAXSPART,
			sp->sloc, sp->ss[SU], safedist2);
}


static int
cyl_partit(ro, pt, pi, mp, cent, axis, d2)	/* slice a cylinder */
FVECT  ro;
unsigned char  *pt;
register int  *pi;
int  mp;
FVECT  cent, axis;
double  d2;
{
	FVECT  newct, newax;
	int  npl, npu;

	if (mp < 2 || dist2(ro, cent) >= d2) {	/* hit limit? */
		setpart(pt, *pi, S0);
		(*pi)++;
		return(1);
	}
					/* subdivide */
	setpart(pt, *pi, SU);
	(*pi)++;
	newax[0] = .5*axis[0];
	newax[1] = .5*axis[1];
	newax[2] = .5*axis[2];
	d2 *= 0.25;
					/* lower half */
	newct[0] = cent[0] - newax[0];
	newct[1] = cent[1] - newax[1];
	newct[2] = cent[2] - newax[2];
	npl = cyl_partit(ro, pt, pi, mp/2, newct, newax, d2);
					/* upper half */
	newct[0] = cent[0] + newax[0];
	newct[1] = cent[1] + newax[1];
	newct[2] = cent[2] + newax[2];
	npu = cyl_partit(ro, pt, pi, mp/2, newct, newax, d2);
					/* return total */
	return(npl + npu);
}


void
flatpart(si, r)				/* partition a flat source */
register SRCINDEX  *si;
register RAY  *r;
{
	register FLOAT  *vp;
	FVECT  v;
	double  du2, dv2;
	int  pi;

	clrpart(si->spt);
	vp = source[si->sn].sloc;
	v[0] = r->rorg[0] - vp[0];
	v[1] = r->rorg[1] - vp[1];
	v[2] = r->rorg[2] - vp[2];
	vp = source[si->sn].snorm;
	if (DOT(v,vp) <= FTINY) {	/* behind source */
		si->np = 0;
		return;
	}
	dv2 = 2.*r->rweight/srcsizerat;
	dv2 *= dv2;
	vp = source[si->sn].ss[SU];
	du2 = dv2 * DOT(vp,vp);
	vp = source[si->sn].ss[SV];
	dv2 *= DOT(vp,vp);
	pi = 0;
	si->np = flt_partit(r->rorg, si->spt, &pi, MAXSPART,
		source[si->sn].sloc,
		source[si->sn].ss[SU], source[si->sn].ss[SV], du2, dv2);
}


static int
flt_partit(ro, pt, pi, mp, cent, u, v, du2, dv2)	/* partition flatty */
FVECT  ro;
unsigned char  *pt;
register int  *pi;
int  mp;
FVECT  cent, u, v;
double  du2, dv2;
{
	double  d2;
	FVECT  newct, newax;
	int  npl, npu;

	if (mp < 2 || ((d2 = dist2(ro, cent)) >= du2
			&& d2 >= dv2)) {	/* hit limit? */
		setpart(pt, *pi, S0);
		(*pi)++;
		return(1);
	}
	if (du2 > dv2) {			/* subdivide in U */
		setpart(pt, *pi, SU);
		(*pi)++;
		newax[0] = .5*u[0];
		newax[1] = .5*u[1];
		newax[2] = .5*u[2];
		u = newax;
		du2 *= 0.25;
	} else {				/* subdivide in V */
		setpart(pt, *pi, SV);
		(*pi)++;
		newax[0] = .5*v[0];
		newax[1] = .5*v[1];
		newax[2] = .5*v[2];
		v = newax;
		dv2 *= 0.25;
	}
					/* lower half */
	newct[0] = cent[0] - newax[0];
	newct[1] = cent[1] - newax[1];
	newct[2] = cent[2] - newax[2];
	npl = flt_partit(ro, pt, pi, mp/2, newct, u, v, du2, dv2);
					/* upper half */
	newct[0] = cent[0] + newax[0];
	newct[1] = cent[1] + newax[1];
	newct[2] = cent[2] + newax[2];
	npu = flt_partit(ro, pt, pi, mp/2, newct, u, v, du2, dv2);
				/* return total */
	return(npl + npu);
}


double
scylform(sn, dir)		/* compute cosine for cylinder's projection */
int  sn;
register FVECT  dir;		/* assume normalized */
{
	register FLOAT  *dv;
	double  d;

	dv = source[sn].ss[SU];
	d = DOT(dir, dv);
	d *= d / DOT(dv,dv);
	return(sqrt(1. - d));
}
