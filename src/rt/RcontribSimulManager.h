/* RCSid $Id: RcontribSimulManager.h,v 2.7 2024/11/08 15:49:34 greg Exp $ */
/*
 *  RcontribSimulManager.h
 *
 *	Rcontrib simulation manager class declaration
 *
 *  Created by Greg Ward on 10/11/2024.
 */

#ifndef RcontribSimulManager_h
#define RcontribSimulManager_h

#include "RtraceSimulManager.h"
#include "RdataShare.h"
#include "lookup.h"
#include "rtprocess.h"

/*
 * As with other single-object classes, many global variable affect
 * behavior.  Besides rendering parameters, there are spectral parameters
 * and output dimensions taken from the environment.
 *
 * NOTE:  "record" and "row" are used interchangeably throughout.
 */

extern char		RCCONTEXT[];		// global rcontrib context

class RcontribSimulManager;			// need forward decl

/// Shared data object for record output (includes header; may be write-only)
class RcontribOutput {
	RcontribOutput *	next;		// next in sorted list
	char *			ofname;		// output file name
	uint32			rowCountPos;	// row count position in header
	void *			rowp;		// current row memory
	bool			NewHeader(const RcontribSimulManager *rcp);
	int			CheckHeader(const RcontribSimulManager *rcp);
public:
	RdataShare *		rData;		// data sharing object
	size_t			rowBytes;	// byte count per row
	const char *		omod;		// single modifier (or NULL)
	int32			obin;		// single output bin (or -1)
	uint32			begData;	// start of data (type-aligned)
	int32			curRow;		// current output row
	uint32			nRows;		// total number of rows
				RcontribOutput(const char *fnm = NULL) {
					next = NULL;
					omod = NULL;
					obin = -1;
					rData = NULL;
					ofname = savqstr(fnm);
					rowBytes = 0;
					nRows = 0;
					rowp = NULL; curRow = -1;
				}
				~RcontribOutput() {
					DoneRow();
					delete rData;
					delete next;
					freeqstr(ofname);
				}
				/// Return output channel name
	const char *		GetName() const {
					if (rData) return rData->GetName();
					return ofname;
				}
				/// Update output row count
	bool			SetRowsDone(int r) {
					if (!rData | (0 >= r) | (r > nRows)) return false;
					char *	rbuf = (char *)rData->GetMemory(rowCountPos, 17, 0);
					sprintf(rbuf, "%-16d", r);
					rbuf[16] = '\n';	// replaces nul byte
					return rData->ReleaseMemory(rbuf, RDSwrite);
				}
				/// Get buffer for indicated row (contents may be random)
	void *			GetRow(int r) {
					if (!rData | (r < 0)) return NULL;
					if (r != curRow) {
						DoneRow();
						if (r < nRows)
							rowp = rData->GetMemory(begData + r*rowBytes,
										 rowBytes, 0);
						if (rowp) curRow = r;
					}
					return rowp;
				}
				/// Current row with byte offset
	void *			InsertionP(int coffset) const {
					if (!rowp | (coffset < 0) | (coffset >= rowBytes))
						return NULL;
					return (char *)rowp + coffset;
				}
				/// Release current row, writing contents
	void			DoneRow() {
					if (rowp) rData->ReleaseMemory(rowp, RDSwrite);
					rowp = NULL; curRow = -1;
				}
				/// Get next in list
	const RcontribOutput *	Next() const {
					return next;
				}
	RcontribOutput *	Next() {
					return next;
				}
				/// RcontribSimulManager gets full access
	friend class		RcontribSimulManager;
};

typedef double		DCOLORV;	// color accumulator type

/// Modifier channel for recording contributions (no constructor/destructor)
struct RcontribMod;

/// Allocate rcontrib accumulator
extern RcontribMod *	NewRcMod(const char *prms = NULL, const char *binexpr = NULL, int ncbins = 1);
/// Free an RcontribMod
extern lut_free_t	FreeRcMod;

/*
 * General RcontribSimulManager class operation:
 *
 *  1)  Call LoadOctree(), then alter the header as desired
 *  2)  Set number of spectral samples (NCSAMP) and call SetDataFormat()
 *  3)  Set xres and yres to desired dimensions (xres>0 for picture output)
 *  4)  Call AddModifier() and AddModFile() to indicate tracked modifiers
 *  5)  Set outOp and cdsF according to desired output/recovery
 *  6)  Set desired computation flags via SetFlag()
 *  7)  Call PrepOutput() to open output channels
 *  8)  Call SetThreadCount() to fork children if desired
 *  9)  Set accum to the number of ray samples per record
 * 10)  Call ComputeRecord() with accum ray samples
 * 11)  Continue until GetRowMax() records have been sent
 * 12)  Call Cleanup()
 *
 * The order of some of these calls may be changed.  Technically, the octree
 * may be loaded anytime before PrepOutput() is called.  Also, SetThreadCount()
 * may be called anytime *after* PrepOutput(), and may be interleaved with
 * calls to ComputeRecord().  The accum setting may be changed at any time.
 * Finally, it is possible to restart the output using ResetRow(), and
 * a zero argument will rewind to the beginning, whence all records
 * may be recalculated.  The previous output rows are not zeroed or deleted,
 * but are overwritten as the calculation proceeds from the new starting point.
 * However, the output file(s) will indicate in the NROWS= line in the header
 * that only the newly calculated rows are present.  If you wish to start over
 * with a different set of modifiers or outputs, call ClearModifiers() instead,
 * which keeps the current octree in memory.  This call also returns to single
 * process mode if any children were running.
 *
 * It is not possible to write to standard output, but the output
 * model is quite flexible thanks to the RdataShare polymorphic class.
 * The current default output class creates a shared, memory-mapped file,
 * which is the most efficient object on most systems.
 *
 * ASCII output is not supported, so full data recovery is.
 */

/// Output channel opening options: new/exclusive, overwrite if exists, or recover data
enum RCOutputOp {RCOnew=0, RCOforce, RCOrecover};

/// Converts above to RdataShare open flags (may be adjusted by calling program)
extern int	RSDOflags[];

/// Call-back function type to create named data channel (freed using "delete" operator)
typedef RdataShare *	RcreateDataShareF(const char *name, RCOutputOp op, size_t siz);

/// Our default data share function
extern RcreateDataShareF	defDataShare;

/// Modifiable ray-tracing flags for rcontrib
#define RCcontrib		(RTmask+1)	// compute contributions? (r.t. coefficients)
#define RCmask			(RTlimDist|RTimmIrrad|RCcontrib)

/// rcontrib-like simulation manager (at most one such object)
class RcontribSimulManager : protected RtraceSimulManager {
protected:
	static RayReportCall	RctCall;	// our callback for traced rays
	ABitMap			rowsDone;	// bit mask of completed rows
	uint32			rInPos;		// which row (record) is next on input?
	uby8			nChan;		// NCSAMP setting for this calculation
	char			dtyp;		// data type ('f', 'd', or 'c')
	uint16			dsiz;		// N-component element size in bytes
	RcontribOutput *	outList;	// ordered list of output channels
	LUTAB			modLUT;		// modifier lookup table
	SUBPROC *		kid;		// array of child processes
	int32 *			kidRow;		// row assigned to each child
	int			nkids;		// child process count (-1 in child)
	bool			UpdateRowsDone(int r);
	int			GetChild(bool forceWait = false);
	bool			StartKids(int n2go);
	int			StopKids(int n2end = 0);
	void			RunChild();
public:
	RCOutputOp		outOp;		// output operation
	RcreateDataShareF *	cdsF;		// data share creator
	int			xres, yres;	// output (picture) size
	uint32			accum;		// # rays to accumulate per record
				RcontribSimulManager(const char *octn = NULL)
						: RtraceSimulManager(NULL, NULL, octn) {
					rInPos = 0;
					nChan = 0;
					dtyp = 'f';
					dsiz = 0;
					outList = NULL;
					memset(&modLUT, 0, sizeof(modLUT));
					modLUT.hashf = lu_shash;
					modLUT.keycmp = strcmp;
					modLUT.freek = efree;
					modLUT.freed = FreeRcMod;
					kid = NULL; kidRow = NULL; nkids = 0;
					rtFlags = RTtraceSources;
					SetTraceCall(&RctCall, this);
					outOp = RCOnew;
					cdsF = &defDataShare;
					xres = yres = 0;
					accum = 1;
				}
				~RcontribSimulManager() {
					if (nkids >= 0) ClearModifiers();
				}
				/// Check modifiable ray-tracing computation flag(s)
	bool			HasFlag(int fl) const {
					return ((rtFlags & RCmask & fl) != 0);
				}
				/// Set/reset modifiable ray-tracing computation flag(s)
	bool			SetFlag(int fl, bool val = true) {
					if (!(fl &= RCmask)) return false;
					if (val) rtFlags |= fl;
					else rtFlags &= ~fl;
					return true;
				}
				/// Load octree and prepare renderer
	bool			LoadOctree(const char *octn) {
					return RtraceSimulManager::LoadOctree(octn);
				}
				/// Prepare header from previous input (or clear)
	bool			NewHeader(const char *inspec=NULL) {
					return RtraceSimulManager::NewHeader(inspec);
				}
				/// Add a string to header (adds newline if none)
	bool			AddHeader(const char *str) {
					return RtraceSimulManager::AddHeader(str);
				}
				/// Append program line to header
	bool			AddHeader(int ac, char *av[]) {
					return RtraceSimulManager::AddHeader(ac, av);
				}
				/// Get current header length in bytes
	int			GetHeadLen() const {
					return RtraceSimulManager::GetHeadLen();
				}
				/// Get header lines if any
	const char *		GetHeadStr() const {
					return RtraceSimulManager::GetHeadStr();
				}
				/// Look for specific header keyword, return value
	const char *		GetHeadStr(const char *key, bool inOK = false) const {
					return RtraceSimulManager::GetHeadStr(key, inOK);
				}
				/// Set output format ('f', 'd', or 'c'), call before mods
	bool			SetDataFormat(int ty);
				/// Get current format (and element size in bytes)
	int			GetFormat(int *siz = NULL) const {
					if (siz) *siz = dsiz;
					return dtyp;
				}
				/// Add a modifier and arguments, create output(s)
	bool			AddModifier(const char *modn, const char *outspec,
						const char *prms = NULL,
						const char *binval = NULL, int bincnt = 1);
				/// Add a file of modifiers with associated arguments
	bool			AddModFile(const char *modfn, const char *outspec,
						const char *prms = NULL,
						const char *binval = NULL, int bincnt = 1);
				/// Get named rcontrib output (or list)
	const RcontribOutput *	GetOutput(const char *nm = NULL) const {
					if (!nm) return outList;
					const RcontribOutput *	op = outList;
					while (op && strcmp(op->GetName(), nm))
						op = op->next;
					return op;
				}
				/// Open output channels and return # completed rows
	int			PrepOutput();
				/// Are we ready to compute some records?
	bool			Ready() const {
					return (rowsDone.Length() > 0) & (accum > 0);
				}
				/// Set number of computation threads (0 => #cores)
	int			SetThreadCount(int nt = 0);
				/// Check thread count (1 means no multi-processing)
	int			NThreads() const {
					return nkids + !nkids;
				}
				/// What is maximum row?
	int			GetRowMax() const {
					if (!outList) return yres * (xres + !xres);
					return outList->nRows;
				}
				/// Get current row count (# rows sent for computation)
	int			GetRowCount() const {
					return rInPos;
				}
				/// Get # rows completed
	int			GetRowFinished() const {
					if (!nkids) return rInPos;
					uint32	nDone = rowsDone.Find(0, false);
					if (nDone == ABMend)
						return rowsDone.Length();
					return nDone;
				}
				/// Add a ray/bundle to compute next record (n=accum)
	int			ComputeRecord(const FVECT orig_direc[]);
				/// Finish pending rays if multi-processing
	bool			FlushQueue() {
					if (nkids <= 0) return true;
					while (GetChild(true) >= 0)
						;
					return true;
				}
				/// Rewind calculation (previous results unchanged)
	bool			ResetRow(int r);
				/// Clear the modifiers and close all outputs
	void			ClearModifiers() {
					if (rowsDone.Length()) {
						SetThreadCount(1);
						cow_doneshare();
						rowsDone.NewBitMap(0);
					}
					lu_done(&modLUT);
					delete outList; outList = NULL;
					nChan = 0;
				}
				/// Close octree, free data, return status
	int			Cleanup(bool everything = false) {
					ClearModifiers();
					return RtraceSimulManager::Cleanup(everything);
				}
};

#endif /* RcontribSimulManager_h */
