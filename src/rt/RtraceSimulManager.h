/* RCSid $Id: RtraceSimulManager.h,v 2.18 2024/11/19 20:39:40 greg Exp $ */
/*
 *  RtraceSimulManager.h
 *
 *	Rtrace simulation manager class declaration (along with base class)
 *	Enqueuing rays will block caller iff #rays >= ThreadsAvail()
 *	Reporting call-backs made from EnqueBundle() and FlushQueue()
 *
 *  Created by Greg Ward on 11/10/22.
 */

#ifndef RtraceSimulManager_h
#define RtraceSimulManager_h

#include "ray.h"
#include "abitmap.h"

extern char *	octname;	// global octree name

extern int	castonly;	// doing ray-casting only?

/// Ray reporting callback method -- returns # successfully reported, -1 to abort
typedef int	RayReportCall(RAY *r, void *cd);

/// Multi-threaded simulation manager base class
class RadSimulManager {
	char *			header;			// header (less intro and format)
	int			hlen;			// header string length
protected:
	bool			SplitRay(RAY *r) {
					return (ray_pnprocs && ray_psend(r) > 0);
				}
public:
				RadSimulManager(const char *octn = NULL) {
					header = NULL; hlen = 0;
					LoadOctree(octn);
				}
				~RadSimulManager() {
					Cleanup();
				}
				/// Load octree and prepare renderer
	bool			LoadOctree(const char *octn);
				/// Prepare header from previous input (or clear)
				/// Normally called during octree load
	bool			NewHeader(const char *inspec = NULL);
				/// Add a line to header (adds newline if none)
	bool			AddHeader(const char *str);
				/// Append program line to header
	bool			AddHeader(int ac, char *av[]);
				/// Get current header length in bytes
	int			GetHeadLen() const {
					return hlen;
				}
				/// Get header lines or empty string
	const char *		GetHeadStr() const {
					return hlen ? header : "";
				}
				/// Look for specific header keyword, return value
	const char *		GetHeadStr(const char *key, bool inOK = false) const;
				/// How many cores are available?
	static int		GetNCores();
				/// Set number of computation threads (0 => #cores)
	int			SetThreadCount(int nt = 0);
				/// Check thread count (1 means no multi-threading)
	int			NThreads() const {
					return ray_pnprocs + !ray_pnprocs;
				}
				/// How many threads are currently unoccupied?
	int			ThreadsAvailable() const {
					return ray_pnprocs ? ray_pnidle : 1;
				}
				/// Are we ready?
	bool			Ready() const {
					return (octname && nobjects > 0);
				}
				/// Process a ray (in subthread), optional result
	int			ProcessRay(RAY *r);
				/// Wait for next result (or fail)
	bool			WaitResult(RAY *r);
				/// Close octree, free data, return status
	int			Cleanup(bool everything = false);
};

/// Flags to control rendering operations
enum {RTdoFIFO=1, RTtraceSources=2, RTlimDist=4, RTimmIrrad=8, RTmask=15};

/// rtrace-like simulation manager (at most one such object)
class RtraceSimulManager : public RadSimulManager {
	RayReportCall *		cookedCall;	// callback for cooked primary rays
	void *			ccData;		// client data for cooked primary rays
	RayReportCall *		traceCall;	// call for every ray in tree
	void *			tcData;		// client data for traced rays
	int			curFlags;	// current operating flags
	ABitMap			srcFollowed;	// source flags changed
				// Call-back for global ray-tracing context
	static void		RTracer(RAY *r);
				// Call-back for FIFO
	static int		Rfifout(RAY *r);
				// Check for changes to render flags, etc.
	bool			UpdateMode();
protected:
	RNUMBER			lastRayID;	// last ray ID assigned
public:
	int			rtFlags;	// operation (RT*) flags
				RtraceSimulManager(RayReportCall *cb = NULL, void *cd = NULL,
						const char *octn = NULL) : RadSimulManager(octn) {
					lastRayID = 0;
					rtFlags = curFlags = 0;
					SetCookedCall(cb, cd);
					traceCall = NULL; tcData = NULL;
				}
				~RtraceSimulManager() {
					FlushQueue();
				}
				/// Set number of computation threads (0 => #cores)
	int			SetThreadCount(int nt = 0) {
					if (nt <= 0) nt = castonly ? 1 : GetNCores();
					if (nt == NThreads()) return nt;
					if (nt < NThreads() && FlushQueue() < 0) return 0;
					return RadSimulManager::SetThreadCount(nt);
				}
				/// Add ray bundle to queue w/ optional 1st ray ID
	int			EnqueueBundle(const FVECT orig_direc[], int n,
						RNUMBER rID0 = 0);
				/// Enqueue a single ray w/ optional ray ID
	int			EnqueueRay(const FVECT org, const FVECT dir,
						RNUMBER rID = 0) {
					if (dir == org+1)
						return(EnqueueBundle((const FVECT *)org, 1, rID) > 0);
					FVECT	orgdir[2];
					VCOPY(orgdir[0], org); VCOPY(orgdir[1], dir);
					return EnqueueBundle(orgdir, 1, rID);
				}
				/// Set/change cooked ray callback
	void			SetCookedCall(RayReportCall *cb, void *cd = NULL) {
					if (cookedCall && (cookedCall != cb) | (ccData != cd))
						FlushQueue();
					cookedCall = cb;
					ccData = cb ? cd : NULL;
				}
				/// Set/change trace callback
	void			SetTraceCall(RayReportCall *cb, void *cd = NULL) {
					if (cb == traceCall) {
						if (cb) tcData = cd;
						return;
					}
					int	nt = NThreads();
					if (nt > 1) SetThreadCount(1);
					traceCall = cb;
					tcData = cb ? cd : NULL;
					if (nt > 1) SetThreadCount(nt);
				}
				/// Are we ready?
	bool			Ready() const {
					return (cookedCall != NULL) | (traceCall != NULL) &&
						RadSimulManager::Ready();
				}
				/// Finish pending rays and complete callbacks (return #sent)
	int			FlushQueue();
				/// Close octree, free data, return status
	int			Cleanup(bool everything = false) {
					SetCookedCall(NULL);
					SetTraceCall(NULL);
					rtFlags = 0;
					UpdateMode();
					lastRayID = 0;
					return RadSimulManager::Cleanup(everything);
				}
};

/// Determine if vector is all zeroes
inline bool
IsZeroVec(const FVECT vec)
{
	return (vec[0] == 0.0) & (vec[1] == 0.0) & (vec[2] == 0.0);
}

#endif /* RtraceSimulManager_h */
