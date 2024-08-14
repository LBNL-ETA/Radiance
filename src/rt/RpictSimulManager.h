/* RCSid $Id: RpictSimulManager.h,v 2.1 2024/08/14 20:05:23 greg Exp $ */
/*
 *  RpictSimulManager.h
 *
 *	Rpict simulation manager class declaration
 *
 *  Created by Greg Ward on 07/11/2024.
 */

#ifndef RpictSimulManager_h
#define RpictSimulManager_h

#include "RtraceSimulManager.h"
#include "view.h"
#include "depthcodec.h"
#include "abitmap.h"

/// Data type flags for pixel access and output
enum RenderDataType {
	RDTnone=0,
	RDTscolor=0x1, RDTrgb=0x2, RDTxyz=0x3, RDTscolr=0x4, RDTrgbe=0x5, RDTxyze=0x6,
	RDTcolorM=0x7,
	RDTdfloat=0x8, RDTdshort=0x10,
	RDTdepthM=0x18
};

#define	RDTcolorT(f)	RenderDataType((f) & RDTcolorM)
#define RDTdepthT(f)	RenderDataType((f) & RDTdepthM)
#define RDTcommonE(f)	(RDTcolorT(f) >= RDTscolr)
#define RDTnewCT(f,c)	RenderDataType((f) & ~RDTcolorM | (c))
#define RDTnewDT(f,d)	RenderDataType((f) & ~RDTdepthM | (d))

/// Pixel accessor (read/write to caller's buffer with possible conversion)
class PixelAccess {
	union {
		COLORV *	f;
		COLRV *		b;
	}		pbase;		// pixel base pointer
	union {
		float *		f;
		short *		s;
	}		dbase;		// depth base pointer
	long		rowStride;	// # values to next y position
	int		dtyp;		// data type flags
	RGBPRIMP	primp;		// color primaries if tristimulus
	COLORMAT	xyz2myrgbmat;	// custom color conversion matrix
	COLORV *	CF3(int x, int y) {
				return pbase.f + (rowStride*y + x)*3;
			}
	COLORV *	GetCF3(int x, int y) const {
				return const_cast<COLORV *>(pbase.f + (rowStride*y + x)*3);
			}
	COLORV *	SCF(int x, int y) {
				return pbase.f + (rowStride*y + x)*NCSAMP;
			}
	COLORV *	GetSCF(int x, int y) const {
				return const_cast<COLORV *>(pbase.f + (rowStride*y + x)*NCSAMP);
			}
	COLRV *		CB3(int x, int y) {
				return pbase.b + (rowStride*y + x)*4;
			}
	COLRV *		GetCB3(int x, int y) const {
				return const_cast<COLRV *>(pbase.b + (rowStride*y + x)*4);
			}
	COLRV *		SCB(int x, int y) {
				return pbase.b + (rowStride*y + x)*(NCSAMP+1);
			}
	COLRV *		GetSCB(int x, int y) const {
				return const_cast<COLRV *>(pbase.b + (rowStride*y + x)*(NCSAMP+1));
			}
public:
	double		refDepth;	// reference depth
			PixelAccess() {
				refDepth = 1.;
				Init();
			}
			PixelAccess(COLORV *rp, int ystride, float *zp=NULL) {
				refDepth = 1.;
				Init(rp, ystride, zp);
			}
			PixelAccess(COLRV *bp, int ystride, float *zp=NULL) {
				refDepth = 1.;
				Init(bp, ystride, zp);
			}
			PixelAccess(COLRV *bp, int ystride, short *dp) {
				refDepth = 1.;
				Init(bp, ystride, dp);
			}
	void		Init() {
				pbase.f = NULL; dbase.f = NULL;
				rowStride = 0;
				primp = NULL;
				dtyp = 0;
			}
			/// Initializers default to rendering color space
	void		Init(COLORV *rp, int ystride, float *zp=NULL) {
				pbase.f = rp; dbase.f = zp;
				rowStride = ystride;
				if (NCSAMP > 3) {
					dtyp = RDTscolor; primp = NULL;
				} else {
					dtyp = RDTrgb; primp = stdprims;
				}
				if (zp) dtyp |= RDTdfloat;
			}
	void		Init(COLRV *bp, int ystride, float *zp=NULL) {
				pbase.b = bp; dbase.f = zp;
				rowStride = ystride;
				if (NCSAMP > 3) {
					dtyp = RDTscolr; primp = NULL;
				} else {
					dtyp = RDTrgbe; primp = stdprims;
				}
				if (zp) dtyp |= RDTdfloat;
			}
	void		Init(COLRV *bp, int ystride, short *dp) {
				pbase.b = bp; dbase.s = dp;
				rowStride = ystride;
				if (NCSAMP > 3) {
					dtyp = RDTscolr; primp = NULL;
				} else {
					dtyp = RDTrgbe; primp = stdprims;
				}
				if (dp) dtyp |= RDTdshort;
			}
			/// Set color space after non-empty initialization
	bool		SetColorSpace(RenderDataType cs, RGBPRIMP pr=NULL);
			/// Get color space
	RenderDataType	ColorSpace() const {
				return RDTcolorT(dtyp);
			}
			/// Get color primaries (NULL if spectral)
	RGBPRIMP	Primaries() const {
				return primp;
			}
			/// Number of represented color/spectral components
	int		NC() const {
				switch (ColorSpace()) {
				case RDTrgb: case RDTrgbe:
				case RDTxyz: case RDTxyze:
					return 3;
				case RDTscolor: case RDTscolr:
					return NCSAMP;
				default:;
				}
				return 0;
			}
			/// Get depth type
	RenderDataType	DepthType() const {
				return RDTdepthT(dtyp);
			}
			/// Get row stride
	long		GetRowStride() const {
				return rowStride;
			}
			/// Assign a pixel value (& depth) from rendered ray value
	bool		SetPixel(int x, int y, const RAY *rp);
			/// Assign pixel color (& depth) -- may re-represent either/both
	bool		SetPixel(int x, int y, const COLORV *pv, float z=0) {
				if (NC() == 3) {
					if (RDTcommonE(dtyp))
						setcolr(CB3(x,y), pv[RED], pv[GRN], pv[BLU]);
					else
						copycolor(CF3(x,y), pv);
				} else if (RDTcommonE(dtyp))
					scolor_scolr(SCB(x,y), pv);
				else
					copyscolor(SCF(x,y), pv);
				if (RDTdepthT(dtyp) == RDTdfloat)
					dbase.f[rowStride*y + x] = z;
				else if (RDTdepthT(dtyp) == RDTdshort)
					dbase.s[rowStride*y + x] = depth2code(z, refDepth);
				return true;
			}
			/// Retrieve pixel color (& depth) -- may convert either/both
	bool		GetPixel(int x, int y, COLORV *pv, float *zp=NULL) const {
				if (NC() == 3) {
					if (RDTcommonE(dtyp))
						colr_color(pv, GetCB3(x,y));
					else
						copycolor(pv, GetCF3(x,y));
				} else if (RDTcommonE(dtyp))
					scolr_scolor(pv, GetSCB(x,y));
				else
					copyscolor(pv, GetSCF(x,y));
				if (!zp) return true;
				if (RDTdepthT(dtyp) == RDTdfloat)
					*zp = dbase.f[rowStride*y + x];
				else if (RDTdepthT(dtyp) == RDTdshort)
					*zp = code2depth(dbase.s[rowStride*y + x], refDepth);
				else
					*zp = .0f;
				return true;
			}
			/// Copy pixel from one location to another (no conversion)
	bool		CopyPixel(int dx, int dy, int sx, int sy) {
				if ((dx==sx) & (dy==sy)) return true;
				const int	nc = NC();
				if (nc == 3) {
					if (RDTcommonE(dtyp))
						copycolr(CB3(dx,dy), GetCB3(sx,sy));
					else
						copycolor(CF3(dx,dy), GetCF3(sx,sy));
				} else if (RDTcommonE(dtyp))
					copyscolr(SCB(dx,dy), GetSCB(sx,sy));
				else
					copyscolor(SCF(dx,dy), GetSCF(sx,sy));
				switch (RDTdepthT(dtyp)) {
				case RDTdfloat:
					dbase.f[rowStride*dy + dx] =
							dbase.f[rowStride*sy + sx];
					break;
				case RDTdshort:
					dbase.s[rowStride*dy + dx] =
							dbase.s[rowStride*sy + sx];
					break;
				default:;
				}
				return true;
			}
};

/// Call-back function for progress reporting
typedef void	ProgReportCB(double pct);

/// rpict-like simulation manager (at most one such object)
class RpictSimulManager : protected RtraceSimulManager {
	static RayReportCall	RtCall;			// our callback for cooked rays
	VIEW			vw;			// frame view
	VIEW			pvw;			// previous view
	int			hvres[2];		// overall picture dimensions
	int			tgsize[2];		// tile grid size
	int			thvres[2];		// tile dimensions
	VIEW			tvw, ptvw;		// this tile's view (& previous)
	PixelAccess		pacc;			// pixel accessor
	char			dunit[32];		// depth with units (if any)
	ABitMap2		doneMap;		// which tile pixels are done
	COLORV *		barPix;			// current render bar pixels
	float *			barDepth;		// current render bar depths
	bool			SetTile(const int ti[2]);
	bool			RenderRect();
	bool			ComputePixel(int x, int y);
	bool			BelowSampThresh(int x, int y, const int noff[4][2]) const;
	void			FillSquare(int x, int y, const int noff[4][2]);
	void			NewBar(int ht = 0);
	bool			LowerBar(int v);
	bool			RenderBelow(int ytop, const int vstep, FILE *pfp,
							const int dt, FILE *dfp=NULL);
public:
	ProgReportCB *		prCB;			// progress report call-back
	RGBPRIMP		prims;			// output primaries (NULL if spectral)
	int			frameNo;		// frame number (0 if not sequence)
				RpictSimulManager(const char *octn = NULL) :
							RtraceSimulManager(RtCall, this, octn) {
					tvw.type = vw.type = 0;
					hvres[0] = hvres[1] = 0;
					thvres[0] = thvres[1] = 0;
					pacc.refDepth = 1.; dunit[0] = '1'; dunit[1] = '\0';
					barPix = NULL; barDepth = NULL;
					prCB = NULL;
					prims = NULL;
					frameNo = 0;
				}
				~RpictSimulManager() {
					NewBar();
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
				/// Get header lines if any
	const char *		GetHeader() const {
					return RtraceSimulManager::GetHeader();
				}
				/// Set number of computation threads (0 => #cores)
	int			SetThreadCount(int nt = 0) {
					return RtraceSimulManager::SetThreadCount(nt);
				}
				/// Check thread count (1 means no multi-threading)
	int			NThreads() const {
					return RtraceSimulManager::NThreads();
				}
				/// How many threads are currently unoccupied?
	int			ThreadsAvailable() const {
					return RtraceSimulManager::ThreadsAvailable();
				}
				/// Are we ready?
	bool			Ready() const {
					return RtraceSimulManager::Ready();
				}
				/// Assign reference depth string (e.g., "2.5/meter")
	bool			SetReferenceDepth(const char *dstr) {
					double	dref = atof(dstr);
					if (dref <= .0) return false;
					strlcpy(dunit, dstr, sizeof(dunit));
					pacc.refDepth = dref;
					return true;
				}
	bool			SetReferenceDepth(double dref, const char *unit=NULL) {
					if (dref <= .0) return false;
					if (unit) sprintf(dunit, "%g/%s", dref, unit);
					else sprintf(dunit, "%g", dref);
					pacc.refDepth = dref;
					return true;
				}
				/// Return reference depth
	double			GetReferenceDepth(char *du=NULL) const {
					if (du) strcpy(du, dunit);
					return pacc.refDepth;
				}
				/// Set up rendering frame (call after octree loaded)
				/// Overall dimensions may be adjusted for view,
				/// optional pixel aspect ratio and tile grid
				/// Increments frameNo if >0
	bool			NewFrame(const VIEW &v, int xydim[2], double *ap=NULL,
						const int *tgrid=NULL);
				/// Get current picture width
	int			GetWidth() const {
					return hvres[0];
				}
				/// Get current picture height
	int			GetHeight() const {
					return hvres[1];
				}
				/// Tile width
	int			TWidth() const {
					return thvres[0];
				}
				/// Tile height
	int			THeight() const {
					return thvres[1];
				}
				/// Render the specified tile in frame
				/// Tile pixels are contiguous unless ystride != 0
				/// Tiles numbered from lower-left at (0,0)
				/// Pixel type influenced by this->prims assignment
	bool			RenderTile(COLORV *rp, int ystride=0, float *zp=NULL,
						const int *tile=NULL);
				/// Same but store as common-exponent COLR or SCOLR
	bool			RenderTile(COLRV *bp, int ystride=0, float *zp=NULL,
						const int *tile=NULL);
				/// Same but also use 16-bit encoded depth buffer
	bool			RenderTile(COLRV *bp, int ystride, short *dp,
						const int *tile=NULL);
				/// Render and write a frame to the named file
				/// Include any header lines set prior to call
				/// Picture file must not already exist
				/// Picture to stdout if pfname==NULL
				/// Depth written to a command if dfname[0]=='!'
	RenderDataType		RenderFrame(const char *pfname,
						RenderDataType dt=RDTrgbe,
						const char *dfname=NULL);
				/// Resume partially finished rendering
				/// Picture file must exist with valid header
	RenderDataType		ResumeFrame(const char *pfname,
						const char *dfname=NULL);
				/// Close octree, free data, return status
	int			Cleanup(bool everything = false) {
					NewBar();
					tvw.type = vw.type = 0;
					hvres[0] = hvres[1] = 0;
					thvres[0] = thvres[1] = 0;
					prims = NULL;
					frameNo = 0;
					return RtraceSimulManager::Cleanup(everything);
				}
};

#endif /* RpictSimulManager_h */
