#ifndef lint
static const char RCSid[] = "$Id: srcobstr.c,v 2.12 2005/04/19 01:15:07 greg Exp $";
#endif
/*
 * Source occlusion caching routines
 */

#include  "ray.h"

#include  "otypes.h"

#include  "otspecial.h"

#include  "source.h"

#define ABS(x)  ((x)>0 ? (x) : -(x))


#if  SHADCACHE			/* preemptive shadow checking */


static int				/* cast source ray to first blocker */
castshadow(int sn, FVECT rorg, FVECT rdir)
{
	RAY     rt;
	
	VCOPY(rt.rorg, rorg);
	VCOPY(rt.rdir, rdir);
	rt.rmax = 0;
	rayorigin(&rt, PRIMARY, NULL, NULL);
					/* check for intersection */
	while (localhit(&rt, &thescene)) {
		RAY	rt1 = rt;	/* pretend we were aimed at source */
		rt1.crtype |= rt1.rtype = SHADOW;
		rt1.rdir[0] = -rt.rdir[0];
		rt1.rdir[1] = -rt.rdir[1];
		rt1.rdir[2] = -rt.rdir[2];
		rt1.rod = -rt.rod;
		VSUB(rt1.rorg, rt.rop, rt.rdir);
		rt1.rot = 1.;
		rt1.rsrc = sn;
					/* record blocker */
		if (srcblocker(&rt1))
			return(1);
					/* move past failed blocker */
		VSUM(rt.rorg, rt.rop, rt.rdir, FTINY);
		rayclear(&rt);		/* & try again... */
	}
	return(0);			/* found no blockers */
}


void					/* initialize occlusion cache */
initobscache(int sn)
{
	register SRCREC *srcp = &source[sn];
	int		cachelen;
	FVECT		rorg, rdir;
	RREAL		d;
	int		i, j, k;
	int		ax, ax1, ax2;

	if (srcp->sflags & (SSKIP|SPROX|SSPOT|SVIRTUAL))
		return;			/* don't cache these */
	if (srcp->sflags & SDISTANT)
		cachelen = 4*SHADCACHE*SHADCACHE;
	else if (srcp->sflags & SFLAT)
		cachelen = SHADCACHE*SHADCACHE*3 + (SHADCACHE&1)*SHADCACHE*4;
	else /* spherical distribution */
		cachelen = SHADCACHE*SHADCACHE*6;
					/* allocate cache */
	srcp->obscache = (OBSCACHE *)malloc(sizeof(OBSCACHE) +
						sizeof(OBJECT)*(cachelen-1));
	if (srcp->obscache == NULL)
		error(SYSTEM, "out of memory in initobscache()");
					/* set parameters */
	if (srcp->sflags & SDISTANT) {
		RREAL   amax = 0;
		for (ax1 = 3; ax1--; )
			if (ABS(srcp->sloc[ax1]) > amax) {
				amax = ABS(srcp->sloc[ax1]);
				ax = ax1;
			}
		srcp->obscache->p.d.ax = ax;
		ax1 = (ax+1)%3;
		ax2 = (ax+2)%3;
		VCOPY(srcp->obscache->p.d.o, thescene.cuorg);
		if (srcp->sloc[ax] > 0)
			srcp->obscache->p.d.o[ax] += thescene.cusize;
		if (srcp->sloc[ax1] < 0)
			srcp->obscache->p.d.o[ax1] += thescene.cusize *
					srcp->sloc[ax1] / amax;
		if (srcp->sloc[ax2] < 0)
			srcp->obscache->p.d.o[ax2] += thescene.cusize *
					srcp->sloc[ax2] / amax;
		srcp->obscache->p.d.e1 = 1. / (thescene.cusize*(1. +
				fabs(srcp->sloc[ax1])/amax));
		srcp->obscache->p.d.e2 = 1. / (thescene.cusize*(1. +
				fabs(srcp->sloc[ax2])/amax));
	} else if (srcp->sflags & SFLAT) {
		VCOPY(srcp->obscache->p.f.u, srcp->ss[SU]);
		normalize(srcp->obscache->p.f.u);
		fcross(srcp->obscache->p.f.v,
				srcp->snorm, srcp->obscache->p.f.u);
	}
					/* clear cache */
	for (i = cachelen; i--; )
		srcp->obscache->obs[i] = OVOID;
					/* cast shadow rays */
	if (srcp->sflags & SDISTANT) {
		for (k = 3; k--; )
			rdir[k] = -srcp->sloc[k];
		for (i = 2*SHADCACHE; i--; )
			for (j = 2*SHADCACHE; j--; ) {
				VCOPY(rorg, srcp->obscache->p.d.o);
				rorg[ax1] += (i+.5) /
					(2*SHADCACHE*srcp->obscache->p.d.e1);
				rorg[ax2] += (j+.5) /
					(2*SHADCACHE*srcp->obscache->p.d.e2);
				castshadow(sn, rorg, rdir);
			}
	} else if (srcp->sflags & SFLAT) {
		d = 0.01*srcp->srad;
		VSUM(rorg, srcp->sloc, srcp->snorm, d);
		for (i = SHADCACHE; i--; )
			for (j = SHADCACHE; j--; ) {
				d = 2./SHADCACHE*(i+.5) - 1.;
				VSUM(rdir, srcp->snorm,
						srcp->obscache->p.f.u, d);
				d = 2./SHADCACHE*(j+.5) - 1.;
				VSUM(rdir, rdir, srcp->obscache->p.f.v, d);
				normalize(rdir);
				castshadow(sn, rorg, rdir);
			}
		for (k = 2; k--; )
		    for (i = SHADCACHE; i--; )
			for (j = SHADCACHE>>1; j--; ) {
				d = 2./SHADCACHE*(i+.5) - 1.;
				if (k)
					VSUM(rdir, srcp->obscache->p.f.u,
						srcp->obscache->p.f.v, d);
				else
					VSUM(rdir, srcp->obscache->p.f.v,
						srcp->obscache->p.f.u, d);
				d = 1. - 2./SHADCACHE*(j+.5);
				VSUM(rdir, rdir, srcp->snorm, d);
				normalize(rdir);
				castshadow(sn, rorg, rdir);
				d = 2.*DOT(rdir, srcp->snorm);
				rdir[0] = d*srcp->snorm[0] - rdir[0];
				rdir[1] = d*srcp->snorm[1] - rdir[1];
				rdir[2] = d*srcp->snorm[2] - rdir[2];
				castshadow(sn, rorg, rdir);
			}
	} else /* spherical distribution */
	    for (k = 6; k--; ) {
		ax = k%3;
		ax1 = (k+1)%3;
		ax2 = (k+2)%3;
		for (i = SHADCACHE; i--; )
			for (j = SHADCACHE; j--; ) {
				rdir[0]=rdir[1]=rdir[2] = 0.;
				rdir[ax] = k<3 ? 1. : -1.;
				rdir[ax1] = 2./SHADCACHE*(i+.5) - 1.;
				rdir[ax2] = 2./SHADCACHE*(j+.5) - 1.;
				normalize(rdir);
				d = 1.05*srcp->srad;
				VSUM(rorg, srcp->sloc, rdir, d);
				castshadow(sn, rorg, rdir);
			}
	    }
}


static OBJECT *			/* return occluder cache entry */
srcobstructp(register RAY *r)
{
	static unsigned long	lastrno = ~0;
	static OBJECT   noobs;
	static OBJECT	*lastobjp;
	SRCREC		*srcp;
	int		ondx;

	noobs = OVOID;
	if (r->rno == lastrno)
		return lastobjp;	/* just recall last pointer */
	DCHECK(r->rsrc < 0, CONSISTENCY,
			"srcobstructp() called with unaimed ray");
	lastrno = r->rno;
	lastobjp = &noobs;
	srcp = &source[r->rsrc];
	if (srcp->sflags & (SSKIP|SPROX|SSPOT|SVIRTUAL))
		return(&noobs);		/* don't cache these */
	if (srcp->obscache == NULL)     /* initialize cache */
		initobscache(r->rsrc);
					/* compute cache index */
	if (srcp->sflags & SDISTANT) {
		int     ax, ax1, ax2;
		double  t;
		ax = srcp->obscache->p.d.ax;
		if ((ax1 = ax+1) >= 3) ax1 -= 3;
		if ((ax2 = ax+2) >= 3) ax2 -= 3;
		t = (srcp->obscache->p.d.o[ax] - r->rorg[ax]) / srcp->sloc[ax];
		if (t <= FTINY)
			return(&noobs); /* could happen if ray is outside */
		ondx = 2*SHADCACHE*(int)(2*SHADCACHE*srcp->obscache->p.d.e1 *
				(r->rorg[ax1] + t*srcp->sloc[ax1] -
					srcp->obscache->p.d.o[ax1]));
		ondx += (int)(2*SHADCACHE*srcp->obscache->p.d.e2 *
				(r->rorg[ax2] + t*srcp->sloc[ax2] -
					srcp->obscache->p.d.o[ax2]));
		if ((ondx < 0) | (ondx >= 4*SHADCACHE*SHADCACHE))
			return(&noobs); /* could happen if ray is outside */
	} else if (srcp->sflags & SFLAT) {
		FVECT   sd;
		RREAL   sd0m, sd1m;
		sd[0] = -DOT(r->rdir, srcp->obscache->p.f.u);
		sd[1] = -DOT(r->rdir, srcp->obscache->p.f.v);
		sd[2] = -DOT(r->rdir, srcp->snorm);
		if (sd[2] < 0)
			return(&noobs); /* shouldn't happen */
		sd0m = ABS(sd[0]);
		sd1m = ABS(sd[1]);
		if (sd[2] >= sd0m && sd[2] >= sd1m) {
			ondx = SHADCACHE*(int)(SHADCACHE*(.5-FTINY) *
					(1. + sd[0]/sd[2]));
			ondx += (int)(SHADCACHE*(.5-FTINY) *
					(1. + sd[1]/sd[2]));
		} else if (sd0m >= sd1m) {
			ondx = SHADCACHE*SHADCACHE;
			if (sd[0] < 0)
				ondx += ((SHADCACHE+1)>>1)*SHADCACHE;
			ondx += SHADCACHE*(int)(SHADCACHE*(.5-FTINY) *
					(1. - sd[2]/sd0m));
			ondx += (int)(SHADCACHE*(.5-FTINY) *
					(1. + sd[1]/sd0m));
		} else /* sd1m > sd0m */ {
			ondx = SHADCACHE*SHADCACHE +
					((SHADCACHE+1)>>1)*SHADCACHE*2;
			if (sd[1] < 0)
				ondx += ((SHADCACHE+1)>>1)*SHADCACHE;
			ondx += SHADCACHE*(int)(SHADCACHE*(.5-FTINY) *
					(1. - sd[2]/sd1m));
			ondx += (int)(SHADCACHE*(.5-FTINY) *
					(1. + sd[0]/sd1m));
		}
		DCHECK((ondx < 0) | (ondx >= SHADCACHE*SHADCACHE*3 +
				(SHADCACHE&1)*SHADCACHE*4), CONSISTENCY,
				"flat source cache index out of bounds");
	} else /* spherical distribution */ {
		int     ax, ax1, ax2;
		RREAL   amax = 0;
		for (ax1 = 3; ax1--; )
			if (ABS(r->rdir[ax1]) > amax) {
				amax = ABS(r->rdir[ax1]);
				ax = ax1;
			}
		if ((ax1 = ax+1) >= 3) ax1 -= 3;
		if ((ax2 = ax+2) >= 3) ax2 -= 3;
		ondx = 2*SHADCACHE*SHADCACHE * ax;
		if (r->rdir[ax] < 0)
			ondx += SHADCACHE*SHADCACHE;
		ondx += SHADCACHE*(int)(SHADCACHE*(.5-FTINY) *
					(1. + r->rdir[ax1]/amax));
		ondx += (int)(SHADCACHE*(.5-FTINY) *
				(1. + r->rdir[ax2]/amax));
		DCHECK((ondx < 0) | (ondx >= SHADCACHE*SHADCACHE*6), CONSISTENCY,
				"radial source cache index out of bounds");
	}
					/* return cache pointer */
	return(lastobjp = &srcp->obscache->obs[ondx]);
}


void				/* free obstruction cache */
freeobscache(SRCREC *srcp)
{
	if (srcp->obscache == NULL)
		return;
	free((void *)srcp->obscache);
	srcp->obscache = NULL;
}

	
int				/* record a source blocker */
srcblocker(register RAY *r)
{
	OBJREC  *m;

	if (r->robj == OVOID || objptr(r->robj) != r->ro ||
			isvolume(r->ro->otype))
		return(0);		/* don't record complex blockers */
	if (r->rsrc < 0 || source[r->rsrc].so == r->ro)
		return(0);		/* just a mistake, that's all */
	m = findmaterial(r->ro);
	if (m == NULL)
		return(0);		/* no material?! */
	if (!isopaque(m->otype))
		return(0);		/* material not a reliable blocker */
	*srcobstructp(r) = r->robj;     /* else record obstructor */
	return(1);
}


int				/* check ray against cached blocker */
srcblocked(RAY *r)
{
	OBJECT  obs = *srcobstructp(r);
	OBJREC  *op;

	if (obs == OVOID)
		return(0);
	op = objptr(obs);		/* check blocker intersection */
	if (!(*ofun[op->otype].funp)(op, r))
		return(0);
	if (source[r->rsrc].sflags & SDISTANT)
		return(1);
	op = source[r->rsrc].so;	/* check source intersection */
	if (!(*ofun[op->otype].funp)(op, r))
		return(1);
	rayclear(r);
	return(0);			/* source in front */
}


#endif  /* SHADCACHE */
