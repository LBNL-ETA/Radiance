#ifndef lint
static const char RCSid[] = "$Id: RtraceSimulManager.cpp,v 2.8 2024/05/02 22:10:43 greg Exp $";
#endif
/*
 *  RtraceSimulManager.cpp
 *
 *	Rtrace simulation manager class implementation
 *
 *  Created by Greg Ward on 2/2/2023.
 */

#include <unistd.h>
#include "RtraceSimulManager.h"
#include "source.h"

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

// How many processors are available?
int
RadSimulManager::GetNCores()
{
	return sysconf(_SC_NPROCESSORS_ONLN);
}

// Set number of computation threads (0 => #cores)
int
RadSimulManager::SetThreadCount(int nt)
{
	if (!Ready())
		return 0;

	if (nt <= 0) nt = castonly ? 1 : GetNCores();

	if (nt == 1)
		ray_pclose(ray_pnprocs);
	else if (nt < ray_pnprocs)
		ray_pclose(ray_pnprocs - nt);
	else if (nt > ray_pnprocs)
		ray_popen(nt - ray_pnprocs);

	return NThreads();
}

// Assign ray to subthread (fails if NThreads()<2)
bool
RadSimulManager::SplitRay(RAY *r)
{
	if (!ray_pnprocs || ThreadsAvailable() < 1)
		return false;

	return (ray_psend(r) > 0);
}

// Process a ray (in subthread), optional result
bool
RadSimulManager::ProcessRay(RAY *r)
{
	if (!Ready()) return false;

	if (!ray_pnprocs) {	// single-threaded mode?
		samplendx++;
		rayvalue(r);
		return true;
	}
	int	rv = ray_pqueue(r);
	if (rv < 0) {
		error(WARNING, "ray tracing process(es) died");
		return false;
	}
	return (rv > 0);
}

// Wait for next result (or fail)
bool
RadSimulManager::WaitResult(RAY *r)
{
	if (!ray_pnprocs)
		return false;

	return (ray_presult(r, 0) > 0);
}

// Close octree, free data, return status
int
RadSimulManager::Cleanup(bool everything)
{
	if (!ray_pnprocs)
		ray_pdone(everything);
	else
		ray_done(everything);
	return 0;
}

// How many threads are currently unoccupied?
int
RadSimulManager::ThreadsAvailable() const
{
	if (!ray_pnprocs) return 1;

	return ray_pnidle;
}

// Global pointer to simulation manager for trace call-back (only one)
static const RtraceSimulManager *	ourRTsimMan = NULL;

// Call-back for trace output
void
RtraceSimulManager::RTracer(RAY *r)
{
	(*ourRTsimMan->traceCall)(r, ourRTsimMan->tcData);
}

// Call-back for FIFO output
int
RtraceSimulManager::Rfifout(RAY *r)
{
	return (*ourRTsimMan->cookedCall)(r, ourRTsimMan->ccData);
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
	if (misMatch & RTdoFIFO && FlushQueue() < 0)
		return false;
	curFlags = rtFlags;
				// update callbacks
	if (traceCall)
		trace = RTracer;
	else if (trace == RTracer)
		trace = NULL;
	if (rtFlags & RTdoFIFO)
		ray_fifo_out = Rfifout;
	else if (ray_fifo_out == Rfifout)
		ray_fifo_out = NULL;
	if ((trace != RTracer) & (ray_fifo_out != Rfifout)) {
		ourRTsimMan = NULL;
	} else if (ourRTsimMan != this) {
		if (ourRTsimMan)
			error(WARNING, "Competing top-level simulation managers?");
		ourRTsimMan = this;
	}
	return true;
}

extern "C" int	m_normal(OBJREC *m, RAY *r);

// compute irradiance rather than radiance
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

// compute first ray intersection only
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
		error(INTERNAL, "EnqueueBundle() called in castonly mode without cookedCall");

	if (!UpdateMode())		// update rendering mode if requested
		return -1;

	while (n-- > 0) {		// queue each ray
		VCOPY(res.rorg, orig_direc[0]);
		VCOPY(res.rdir, orig_direc[1]);
		orig_direc += 2;
		rayorigin(&res, PRIMARY, NULL, NULL);
		if (!rID0)
			res.rno = ++lastRayID;
		else if (curFlags & RTdoFIFO)
			error(INTERNAL, "Ray number assignment unsupported with FIFO");
		else
			res.rno = lastRayID = rID0++;
		if (curFlags & RTimmIrrad)
			res.revf = rayirrad;
		else if (castonly)
			res.revf = raycast;
		double	d = normalize(res.rdir);
		bool	sendRes = (cookedCall != NULL);
		if (d > 0) {		// direction vector is valid?
			if (curFlags & RTlimDist)
				res.rmax = d;
			if (((curFlags&RTdoFIFO) != 0) & (ray_pnprocs > 0)) {
				if (ray_fifo_in(&res) < 0)
					return -1;
				sendRes = false;
			} else
				sendRes &= ProcessRay(&res);
		} else if (ThreadsAvailable() < NThreads() &&
				FlushQueue() < 0)
			return -1;
					// may be dummy ray
		if (sendRes && (*cookedCall)(&res, ccData) < 0)
			return -1;
		nqueued++;
	}
	return nqueued;
}

// Finish pending rays and complete callbacks
int
RtraceSimulManager::FlushQueue()
{
	if (curFlags & RTdoFIFO) {
		if (ray_pnprocs)
			return ray_fifo_flush();
		return 0;
	}
	int	nsent = 0;
	RAY	res;

	while (WaitResult(&res)) {
		if (!cookedCall) continue;
		int	rv = (*cookedCall)(&res, ccData);
		if (rv < 0) return -1;
		nsent += rv;
	}
	return nsent;
}
