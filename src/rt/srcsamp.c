#ifndef lint
static const char	RCSid[] = "$Id: srcsamp.c,v 2.22 2024/11/15 20:47:42 greg Exp $";
#endif
/*
 * Source sampling routines
 *
 *  External symbols declared in source.h
 */

#include "copyright.h"

#include  "ray.h"

#include  "source.h"

#include  "random.h"

#ifdef SSKIPOPT
/* The following table is used for skipping sources */
static uby8	*srcskipflags = NULL;		/* source inclusion lookup */
static int	ssf_count = 0;			/* number of flag entries */
static int	ssf_max = 0;			/* current array size */
static uby8	*ssf_noskip = NULL;		/* set of zero flags */

uby8		*ssf_select = NULL;		/* sources we may skip */

/* Find/allocate source skip flag entry (free all if NULL) */
int
sskip_rsi(uby8 *flags)
{
	uby8	*flp;
	int	i;

	if (flags == NULL) {		/* means clear all */
		efree(srcskipflags); srcskipflags = NULL;
		ssf_count = ssf_max = 0;
		sskip_free(ssf_noskip);
		sskip_free(ssf_select);
		return(0);
	}
	if (ssf_noskip == NULL)		/* first call? */
		ssf_noskip = sskip_new();

	if (sskip_eq(flags, ssf_noskip))
		return(-1);		/* nothing to skip */
					/* search recent entries */
	flp = srcskipflags + ssf_count*SSKIPFLSIZ;
	for (i = ssf_count; i-- > 0; )
		if (sskip_eq(flp -= SSKIPFLSIZ, flags))
			return(-2-i);	/* found it! */
					/* else tack on new entry */
	if (ssf_count >= ssf_max) {	/* need more space? */
fprintf(stderr, "DEBUG: skip flag array > %d entries (%.2f MBytes)\n",
ssf_count, SSKIPFLSIZ/1024./1024.*ssf_count);
		ssf_max = ssf_count + (ssf_count>>2) + 64;
		if (ssf_max <= ssf_count &&
				(ssf_max = ssf_count+1024) <= ssf_count)
			error(SYSTEM, "out of space in sskip_rsi()");

		srcskipflags = (uby8 *)erealloc(srcskipflags,
						ssf_max*SSKIPFLSIZ);
	}
	sskip_cpy(srcskipflags + ssf_count*SSKIPFLSIZ, flags);

	return(-2 - ssf_count++);	/* return index (< -1) */
}

/* Get skip flags associated with RAY rsrc index (or NULL) */
uby8 *
sskip_flags(int rsi)
{
	if (rsi >= -1)
		return(ssf_noskip);

	if ((rsi = -2 - rsi) >= ssf_count)
		error(CONSISTENCY, "bad index to sskip_flags()");

	return(srcskipflags + rsi*SSKIPFLSIZ);
}

/* OR in a second set of flags into a first */
void
sskip_addflags(uby8 *dfl, const uby8 *sfl)
{
	int	nb = SSKIPFLSIZ;

	while (nb--)
		*dfl++ |= *sfl++;
}
#endif

int
srcskip(			/* pre-emptive test for source to skip */
	int  sn,
	RAY  *r
)
{
	SRCREC  *sp = source + sn;

	if (sp->sflags & SSKIP)
		return(1);
#ifdef SSKIPOPT
	if (r->rsrc < -1 &&	/* ray has custom skip flags? */
			sskip_chk(sskip_flags(r->rsrc), sn))
		return(1);
#endif
	if ((sp->sflags & (SPROX|SDISTANT)) != SPROX)
		return(0);

	return(dist2(r->rorg, sp->sloc) >
			(sp->sl.prox + sp->srad)*(sp->sl.prox + sp->srad));
}

double
nextssamp(			/* compute sample for source, rtn. distance */
	RAY  *r,		/* origin is read, direction is set */
	SRCINDEX  *si		/* source index (modified to current) */
)
{
	int  cent[3], size[3], parr[2];
	SRCREC  *srcp;
	double  vpos[3];
	double  d;
	int  i;
nextsample:
	while (++si->sp >= si->np) {	/* get next sample */
		if (++si->sn >= nsources)
			return(0.0);	/* no more */
		if (srcskip(si->sn, r))
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
	srcp = source + si->sn;
	if (dstrsrc > FTINY) {			/* jitter sample */
		dimlist[ndims] = si->sn + 8831;
		dimlist[ndims+1] = si->sp + 3109;
		d = urand(ilhash(dimlist,ndims+2)+samplendx);
		if (srcp->sflags & SFLAT) {
			multisamp(vpos, 2, d);
			vpos[SW] = 0.5;
		} else
			multisamp(vpos, 3, d);
		for (i = 0; i < 3; i++)
			vpos[i] = dstrsrc * (1. - 2.*vpos[i]) *
					(double)size[i]*(1.0/MAXSPART);
	} else
		vpos[0] = vpos[1] = vpos[2] = 0.0;

	VSUM(vpos, vpos, cent, 1.0/MAXSPART);
					/* avoid circular aiming failures */
	if ((srcp->sflags & SCIR) && (si->np > 1) | (dstrsrc > 0.7)) {
		FVECT	trim;
		if (srcp->sflags & (SFLAT|SDISTANT)) {
			d = 1.12837917;		/* correct setflatss() */
			trim[SU] = d*sqrt(1.0 - 0.5*vpos[SV]*vpos[SV]);
			trim[SV] = d*sqrt(1.0 - 0.5*vpos[SU]*vpos[SU]);
			trim[SW] = 0.0;
		} else {
			trim[SW] = trim[SU] = vpos[SU]*vpos[SU];
			d = vpos[SV]*vpos[SV];
			if (d > trim[SW]) trim[SW] = d;
			trim[SU] += d;
			d = vpos[SW]*vpos[SW];
			if (d > trim[SW]) trim[SW] = d;
			trim[SU] += d;
			if (trim[SU] > FTINY*FTINY) {
				d = 1.0/0.7236;	/* correct sphsetsrc() */
				trim[SW] = trim[SV] = trim[SU] =
						d*sqrt(trim[SW]/trim[SU]);
			} else
				trim[SW] = trim[SV] = trim[SU] = 0.0;
		}
		for (i = 0; i < 3; i++)
			vpos[i] *= trim[i];
	}
					/* compute direction */
	for (i = 0; i < 3; i++)
		r->rdir[i] = srcp->sloc[i] +
				vpos[SU]*srcp->ss[SU][i] +
				vpos[SV]*srcp->ss[SV][i] +
				vpos[SW]*srcp->ss[SW][i];

	if (!(srcp->sflags & SDISTANT))
		VSUB(r->rdir, r->rdir, r->rorg);
					/* compute distance */
	if ((d = normalize(r->rdir)) == 0.0)
		goto nextsample;		/* at source! */

					/* compute sample size */
	if (srcp->sflags & SFLAT) {
		si->dom = sflatform(si->sn, r->rdir);
		si->dom *= size[SU]*size[SV]*(1.0/MAXSPART/MAXSPART);
	} else if (srcp->sflags & SCYL) {
		si->dom = scylform(si->sn, r->rdir);
		si->dom *= size[SU]*(1.0/MAXSPART);
	} else {
		si->dom = size[SU]*size[SV]*(double)size[SW] *
				(1.0/MAXSPART/MAXSPART/MAXSPART) ;
	}
	if (srcp->sflags & SDISTANT) {
		si->dom *= srcp->ss2;
		return(FHUGE);
	}
	if (si->dom <= 1e-4)
		goto nextsample;		/* behind source? */
	si->dom *= srcp->ss2/(d*d);
	return(d);		/* sample OK, return distance */
}


int
skipparts(			/* skip to requested partition */
	int  ct[3],
	int  sz[3],		/* center and size of partition (returned) */
	int  pp[2],		/* current index, number to skip (modified) */
	unsigned char  *pt	/* partition array */
)
{
	int  p;
					/* check this partition */
	p = spart(pt, pp[0]);
	pp[0]++;
	if (p == S0) {			/* leaf partition */
		if (pp[1]) {
			pp[1]--;
			return(0);	/* not there yet */
		} else
			return(1);	/* we've arrived */
	}
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
nopart(				/* single source partition */
	SRCINDEX  *si,
	RAY  *r
)
{
	clrpart(si->spt);
	setpart(si->spt, 0, S0);
	si->np = 1;
}


static int
cyl_partit(				/* slice a cylinder */
	FVECT  ro,
	unsigned char  *pt,
	int  *pi,
	int  mp,
	FVECT  cent,
	FVECT  axis,
	double  d2
)
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
cylpart(			/* partition a cylinder */
	SRCINDEX  *si,
	RAY  *r
)
{
	double  dist2, safedist2, dist2cent, rad2;
	FVECT  v;
	SRCREC  *sp;
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
flt_partit(				/* partition flatty */
	FVECT  ro,
	unsigned char  *pt,
	int  *pi,
	int  mp,
	FVECT  cent,
	FVECT  u,
	FVECT  v,
	double  du2,
	double  dv2
)
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


void
flatpart(				/* partition a flat source */
	SRCINDEX  *si,
	RAY  *r
)
{
	RREAL  *vp;
	FVECT  v;
	double  du2, dv2;
	int  pi;

	clrpart(si->spt);
	vp = source[si->sn].sloc;
	v[0] = r->rorg[0] - vp[0];
	v[1] = r->rorg[1] - vp[1];
	v[2] = r->rorg[2] - vp[2];
	vp = source[si->sn].snorm;
	if (DOT(v,vp) <= 0.) {		/* behind source */
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


double
scylform(			/* compute cosine for cylinder's projection */
	int  sn,
	FVECT  dir		/* assume normalized */
)
{
	RREAL  *dv;
	double  d;

	dv = source[sn].ss[SU];
	d = DOT(dir, dv);
	d *= d / DOT(dv,dv);
	return(sqrt(1. - d));
}
