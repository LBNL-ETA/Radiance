#ifndef lint
static const char RCSid[] = "$Id: RtraceSimulManager.cpp,v 2.1 2023/02/08 17:41:48 greg Exp $";
#endif
/*
 *  RtraceSimulManager.cpp
 *
 *	Rtrace simulation manager class implementation
 *
 *  Created by Greg Ward on 2/2/2023.
 */

#include "RtraceSimulManager.h"
#include "source.h"

extern int	castonly;	// doing ray-casting only?

// Load octree and prepare renderer
bool
RadSimulManager::LoadOctree(const char *octn)
{
	if (octname) {		// already running?
		if (octn && !strcmp(octn, octname))
			return true;
		Cleanup();
	}
	if (!octn)
		return false;

	ray_init((char *)octn);
	return true;
}

// Set number of computation threads (0 => #cores)
int
RadSimulManager::SetThreadCount(int nt)
{
	return nThreads = 1;	// XXX temporary
}

// Close octree, free data, return status
int
RadSimulManager::Cleanup()
{
	ray_done(0);
	return 0;
}

// How many threads are currently unoccupied?
int
RadSimulManager::ThreadsAvailable() const
{
	return 1;	// XXX temporary
}

// Global pointer to simulation manager for trace call-back (only one)
static const RtraceSimulManager *	ourRTsimMan = NULL;

void	// static call-back
RtraceSimulManager::RTracer(RAY *r)
{
	(*ourRTsimMan->traceCall)(r, ourRTsimMan->tcData);
}

// Check for changes to render flags & adjust accordingly
bool
RtraceSimulManager::UpdateMode()
{
	rtFlags &= RTmask;
	if (!cookedCall)
		rtFlags &= ~RTdoFIFO;
	if (!traceCall)
		rtFlags &= ~RTtraceSources;
	if (rtFlags & RTimmIrrad)
		rtFlags &= ~RTlimDist;

	int	misMatch = rtFlags ^ curFlags;
				// updates based on toggled flags
	if (misMatch & RTtraceSources) {
		if (rtFlags & RTtraceSources) {
			for (int sn = 0; sn < nsources; sn++)
				source[sn].sflags |= SFOLLOW;
		} else		// cannot undo this...
			rtFlags |= RTtraceSources;
	}
	if (misMatch & RTdoFIFO) {
		if (!FlushQueue())
			return false;
	}
	curFlags = rtFlags;
				// update trace callback
	if (traceCall) {
		if (ourRTsimMan && ourRTsimMan != this)
			error(WARNING, "Competing top-level simulation managers?");
		ourRTsimMan = this;
		trace = RTracer;
	} else if (ourRTsimMan == this) {
		trace = NULL;
		ourRTsimMan = NULL;
	}
	return true;
}

extern "C" int	m_normal(OBJREC *m, RAY *r);

/* compute irradiance rather than radiance */
static void
rayirrad(RAY *r)
{
					/* pretend we hit surface */
	r->rxt = r->rot = 1e-5;
	VSUM(r->rop, r->rorg, r->rdir, r->rot);
	r->ron[0] = -r->rdir[0];
	r->ron[1] = -r->rdir[1];
	r->ron[2] = -r->rdir[2];
	r->rod = 1.0;
					/* compute result */
	r->revf = raytrace;
	m_normal(&Lamb, r);
	r->revf = rayirrad;
}

/* compute first ray intersection only */
static void
raycast(RAY *r)
{
	if (!localhit(r, &thescene)) {
		if (r->ro == &Aftplane) {	/* clipped */
			r->ro = NULL;
			r->rot = FHUGE;
		} else
			sourcehit(r);
	}
}

// Add ray bundle to queue w/ optional 1st ray ID
int
RtraceSimulManager::EnqueueBundle(const FVECT orig_direc[], int n, RNUMBER rID0)
{
	int	nqueued = 0;
	RAY	res;

	if (!Ready())
		return -1;

	if (castonly && !cookedCall)
		error(CONSISTENCY, "EnqueueBundle() called in castonly mode without cookedCall");

	if (!UpdateMode())		// update rendering mode if requested
		return -1;

	while (n-- > 0) {		// queue each ray
		double	d;
		VCOPY(res.rorg, orig_direc[0]);
		VCOPY(res.rdir, orig_direc[1]);
		orig_direc += 2;
		rayorigin(&res, PRIMARY, NULL, NULL);
		if (rID0) res.rno = rID0++;
		else res.rno = ++lastRayID;
		if (curFlags & RTimmIrrad)
			res.revf = rayirrad;
		else if (castonly)
			res.revf = raycast;
		d = normalize(res.rdir);
		if (d > 0) {		// direction vector is valid?
			if (curFlags & RTlimDist)
				res.rmax = d;
			samplendx++;
			rayvalue(&res);		// XXX single-threaded for now
			++nqueued;
		} else if (ThreadsAvailable() < NThreads() &&
				!FlushQueue())
			return -1;
		if (cookedCall)
			(*cookedCall)(&res, ccData);
	}
	return nqueued;
}

// Finish pending rays and complete callbacks
bool
RtraceSimulManager::FlushQueue()
{
	return true;		// XXX no-op for now
}
