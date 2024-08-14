#ifndef lint
static const char RCSid[] = "$Id: RpictSimulManager.cpp,v 2.1 2024/08/14 20:05:23 greg Exp $";
#endif
/*
 *  RpictSimulManager.cpp
 *
 *	Rpict simulation manager implementation
 *
 *  Created by Greg Ward on 07/11/2024.
 */

#include <ctype.h>
#include "platform.h"
#include "RpictSimulManager.h"
#include "depthcodec.h"
#include "random.h"

/************* Imported globals from rxpmain.c *************/

extern VIEW  ourview;			/* viewing parameters */
extern int  hres, vres;			/* current image resolution */

extern int  psample;			/* pixel sample size */
extern double  maxdiff;			/* max. sample difference */
extern double  dstrpix;			/* square pixel distribution */

extern double  mblur;			/* motion blur parameter */

extern double  dblur;			/* depth-of-field blur parameter */

// Assign a pixel value (& depth) from rendered ray value
bool
PixelAccess::SetPixel(int x, int y, const RAY *rp)
{
	if (!rp) return false;

	COLOR	col;
	float	zv = 0;			// get depth if needed
	if (DepthType())
		zv = raydistance(rp);

	switch (ColorSpace()) {
	case RDTscolor:			// keeping rendered spectrum?
	case RDTscolr:
		return SetPixel(x, y, rp->rcol, zv);
	case RDTrgb:
	case RDTrgbe:
	case RDTxyz:
	case RDTxyze:
		scolor_out(col, primp, rp->rcol);
		return SetPixel(x, y, col, zv);
	default:
		error(INTERNAL, "missing color space type in SetPixel()");
	}
	return false;
}

// Set color space after non-empty initialization
bool
PixelAccess::SetColorSpace(RenderDataType cs, RGBPRIMP pr)
{
	if (!dtyp) return false;

	if (!(cs = RDTcolorT(cs)))
		cs = RDTcolorT(dtyp);
	else if (RDTcommonE(cs) ^ RDTcommonE(dtyp))
		return false;

	if (NCSAMP == 3) {
		if (cs == RDTscolr) cs = RDTrgbe;
		else if (cs == RDTscolor) cs = RDTrgb;
	}
	switch (cs) {
	case RDTxyze:
	case RDTxyz:
		primp = xyzprims;
		break;
	case RDTrgbe:
	case RDTrgb:
		primp = pr ? pr : stdprims;
		break;
	default:			// RDTscolr | RDTscolor
		primp = NULL;
		break;
	}
	dtyp = RDTnewCT(dtyp, cs);
	if (primp)
		xyz2myrgbmat[0][0] = 0;
	return true;
}

/*
 * Set up rendering frame (call after octree loaded)
 * Overall dimensions may be adjusted for view,
 * optional pixel aspect ratio and tile grid
 * Increments frameNo if >0
 */
bool
RpictSimulManager::NewFrame(const VIEW &v, int xydim[2], double *ap, const int *tgrid)
{
	double	pasp = 1.;

	if (!xydim) return false;
	if (!ap) ap = &pasp;
	pvw = vw;			// save previous view for motion blur
	vw = v;
	const char *	verr = setview(&vw);
	if (verr) {
		error(WARNING, verr);
		vw = pvw;
		return false;
	}
	const double	va = viewaspect(&vw);
	normaspect(va, ap, &xydim[0], &xydim[1]);
					// set up tiling?
	if (tgrid && (tgrid[0] > 0) & (tgrid[1] > 0) & (tgrid[0]*tgrid[1] > 1)) {
		if ((8*tgrid[0] >= xydim[0]) | (8*tgrid[1] >= xydim[1])) {
			error(WARNING, "Excessive tiling for image size");
			return false;
		}
		xydim[0] -= xydim[0] % (tgsize[0] = tgrid[0]);
		xydim[1] -= xydim[1] % (tgsize[1] = tgrid[1]);
		*ap = va * xydim[0] / xydim[1];
	} else
		tgsize[0] = tgsize[1] = 1;

	if (vw.vaft > FTINY) rtFlags |= RTlimDist;
	else rtFlags &= ~RTlimDist;
	hvres[0] = xydim[0]; hvres[1] = xydim[1];
	thvres[0] = hvres[0]/tgsize[0];		// presumed tile width
	thvres[1] = hvres[1]/tgsize[1];		// ...and height
	frameNo += (frameNo > 0);		// caller may override after
	return true;
}

// Call-back for rendered pixel
int
RpictSimulManager::RtCall(RAY *r, void *cd)
{
	RpictSimulManager *	rsp = (RpictSimulManager *)cd;
	const int		ty = r->rno / rsp->TWidth();
	const int		tx = r->rno - (RNUMBER)ty*rsp->TWidth();

	if (ty >= rsp->THeight()) {
		error(INTERNAL, "bad pixel calculation position in RtCall()");
		return -1;
	}
	if (!rsp->doneMap.TestAndSet(tx, ty)) {
		error(WARNING, "duplicate pixel calculation");
		return 0;
	}
	return rsp->pacc.SetPixel(tx, ty, r);
}

// Set up the specified tile (or entire image if NULL)
bool
RpictSimulManager::SetTile(const int ti[2])
{
	tvw = vw; ptvw = pvw;

 	if (ti) {
		if ((ti[0] < 0) | (ti[0] >= tgsize[0]) |
				(ti[1] < 0) | (ti[1] >= tgsize[1])) {
			error(INTERNAL, "illegal tile specification in SetTile()");
			return false;
		}
		const char *	verr = cropview(&tvw,
						(double)ti[0]/tgsize[0],
						(double)ti[1]/tgsize[1],
						(ti[0]+1.)/tgsize[0],
						(ti[1]+1.)/tgsize[1]);
		if (verr) {
			sprintf(errmsg, "crop failure @ tile (%d,%d)/(%d,%d): %s",
					ti[0], ti[1], tgsize[0], tgsize[1], verr);
			error(USER, errmsg);
			return false;
		}			// previous tile view for blur
		if (!ptvw.type | (mblur <= FTINY) ||
				cropview(&ptvw, (double)ti[0]/tgsize[0],
						(double)ti[1]/tgsize[1],
						(ti[0]+1.)/tgsize[0],
						(ti[1]+1.)/tgsize[1]))
			ptvw.type = 0;

	} else if ((tgsize[0] > 1) | (tgsize[1] > 1)) {
		error(INTERNAL, "missing tile specification in SetTile()");
		return false;
	}
	return doneMap.NewBitMap(TWidth(), THeight());
}

#define	 pixjitter()	(.5+dstrpix*(.5-frandom()))

// Send the indicated pixel to ray tracer
bool
RpictSimulManager::ComputePixel(int x, int y)
{
	static const SCOLOR	scBlack = {0};
	int	i;
	FVECT	rodir[2];
	double	hpos = (x+pixjitter())/TWidth();
	double	vpos = (y+pixjitter())/THeight();
	double	dlim = viewray(rodir[0], rodir[1], &tvw, hpos, vpos);
	if (dlim < -FTINY) {	// off view?
		pacc.SetPixel(x, y, scBlack);
		doneMap.Set(x, y);
		return true;
	}
	if (ptvw.type) {	// add motion blur if requested
		FVECT	rorg2, rdir2;
		double	dlim2 = viewray(rorg2, rdir2, &ptvw, hpos, vpos);
		if (dlim2 >= -FTINY) {
			const double	d = mblur*(.5-frandom());
			dlim = (1.-d)*dlim + d*dlim2;
			for (i = 3; i--; ) {
				rodir[0][i] = (1.-d)*rodir[0][i] + d*rorg2[i];
				rodir[1][i] = (1.-d)*rodir[1][i] + d*rdir2[i];
			}
			if (normalize(rodir[1]) == 0)
				return false;
		}
	}
				// depth-of-field blur if any
	if (!jitteraperture(rodir[0], rodir[1], &tvw, dblur))
		return false;
				// include aft clipping distance?
	for (i = (dlim > FTINY)*3; i--; )
		rodir[1][i] *= dlim;

	return EnqueueRay(rodir[0], rodir[1], (RNUMBER)y*TWidth()+x);
}

// Check if neighbor differences are below pixel sampling threshold
bool
RpictSimulManager::BelowSampThresh(int x, int y, const int noff[4][2]) const
{
	SCOLOR	pval[4];
	float	dist[4];
	int	i, j;

	for (i = 4; i--; ) {		// get pixels from tile store
		int	px = x + noff[i][0];
		int	py = y + noff[i][1];
		if (!doneMap.Check(px, py) ||
				!pacc.GetPixel(px, py, pval[i], &dist[i]))
			return false;
	}
	const bool	spectr = (pacc.NC() > 3);
	for (i = 4; --i; )		// do pairwise comparisons
	    for (j = i; j--; ) {
	        if (pacc.DepthType() &&
	        		(dist[i] - dist[j] > maxdiff*dist[j]) |
	        		(dist[j] - dist[i] > maxdiff*dist[i]))
			return false;
		if (spectr ? sbigsdiff(pval[i], pval[j], maxdiff) :
				bigdiff(pval[i], pval[j], maxdiff))
			return false;
	    }
	return true;			// linear interpolation OK
}

// Fill an interior square patch with interpolated values
void
RpictSimulManager::FillSquare(int x, int y, const int noff[4][2])
{
	SCOLOR	pval[4];
	float	dist[4];
	int	i, j;
					// assumes 4 corners are valid!
	for (i = 4; i--; )
		pacc.GetPixel(x+noff[i][0], y+noff[i][1], pval[i], &dist[i]);

	i = abs(noff[1][0]-noff[0][0]);
	j = abs(noff[1][1]-noff[0][1]);
	const int	slen =  i > j ? i : j;
	const double	sf = 1./slen;
	const bool	spectr = (pacc.NC() > 3);
	for (i = 0; i <= slen; i++) {	// bilinear interpolant
	    const double	c1 = i*sf;
	    for (j = 0; j <= slen; j++) {
	    	const double	c2 = j*sf;
		const int	px = int(x + (1.-c1)*(1.-c2)*noff[0][0] +
						c1*(1.-c2)*noff[1][0] +
						(1.-c1)*c2*noff[2][0] +
						c1*c2*noff[3][0] + .5);
		const int	py = int(y + (1.-c1)*(1.-c2)*noff[0][1] +
						c1*(1.-c2)*noff[1][1] +
						(1.-c1)*c2*noff[2][1] +
						c1*c2*noff[3][1] + .5);
		if (!doneMap.TestAndSet(px, py))
			continue;
		float		zval = 0;
		if (pacc.DepthType())
			zval = (1.-c1)*(1.-c2)*dist[0] + c1*(1.-c2)*dist[1] +
					(1.-c1)*c2*dist[2] + c1*c2*dist[3];
		if (spectr) {		// XXX assumes pacc.NC() == NCSAMP
			SCOLOR	ipval, tpval;
			copyscolor(ipval, pval[0]);
			scalescolor(ipval, (1.-c1)*(1.-c2));
			copyscolor(tpval, pval[1]);
			scalescolor(tpval, c1*(1.-c2));
			saddscolor(ipval, tpval);
			copyscolor(tpval, pval[2]);
			scalescolor(tpval, (1.-c1)*c2);
			saddscolor(ipval, tpval);
			copyscolor(tpval, pval[3]);
			scalescolor(tpval, c1*c2);
			saddscolor(ipval, tpval);
			pacc.SetPixel(px, py, ipval, zval);
		} else {		// tristimulus interpolation
			COLOR	ipval, tpval;
			copycolor(ipval, pval[0]);
			scalecolor(ipval, (1.-c1)*(1.-c2));
			copycolor(tpval, pval[1]);
			scalecolor(tpval, c1*(1.-c2));
			addcolor(ipval, tpval);
			copycolor(tpval, pval[2]);
			scalecolor(tpval, (1.-c1)*c2);
			addcolor(ipval, tpval);
			copycolor(tpval, pval[3]);
			scalecolor(tpval, c1*c2);
			addcolor(ipval, tpval);
			pacc.SetPixel(px, py, ipval, zval);
		}
	    }
	}
}

// helper function to set up quincunx sampling
static void
SetQuincunx(ABitMap2 *bmp2, int noff[4][2], const int spc, const bool odd)
{
	for (int y = 0; y < bmp2->Height(); y += spc>>1)
	    for (int x = odd*(spc>>1); x < bmp2->Width(); x += spc)
	    	bmp2->Set(x, y);
					// order neighbors CCW
	if (odd) {
		noff[0][0] = spc>>1; noff[0][1] = 0;
		noff[1][0] = 0; noff[1][1] = spc>>1;
		noff[2][0] = -(spc>>1); noff[2][1] = 0;
		noff[3][0] = 0; noff[3][1] = -(spc>>1);
	} else {
		noff[0][0] = spc>>1; noff[0][1] = spc>>1;
		noff[1][0] = -(spc>>1); noff[1][1] = spc>>1;
		noff[2][0] = -(spc>>1); noff[2][1] = -(spc>>1);
		noff[3][0] = spc>>1; noff[3][1] = -(spc>>1);
	}
}

// Render (or finish rendering) current tile
bool
RpictSimulManager::RenderRect()
{
	if (!tvw.type || !Ready()) {
		error(INTERNAL, "need octree and view for RenderRect()");
		return false;
	}
	ABitMap2	doneSamples = doneMap;
	int		sp2 = ceil(log2((TWidth()>THeight() ? TWidth() : THeight()) - 1.));
	int		layer = 0;
	int		x, y;
fprintf(stderr, "Rendering %dx%d tile with psample=%d, maxdiff=%.3f ...\n",
TWidth(), THeight(), psample, maxdiff);
	while (sp2 > 0) {
		ABitMap2	sampMap(TWidth(), THeight());
		int		noff[4][2];
		if ((prCB != NULL) & (barPix == NULL))
			(*prCB)(100.*doneMap.SumTotal()/doneMap.Width()/doneMap.Height());
		SetQuincunx(&sampMap, noff, 1<<sp2, layer&1);
		sampMap -= doneSamples;	// avoid resampling pixels
		// Are we into adaptive sampling realm?
		if (noff[0][0]*noff[0][0] + noff[0][1]*noff[0][1] < psample*psample) {
			if (FlushQueue() < 0)	// need results to check threshold
				return false;
			ABitMap2	fillMap = sampMap;
			for (x = y = 0; sampMap.Find(&x, &y); x++)
				if (BelowSampThresh(x, y, noff))
					sampMap.Reset(x, y);
					// spread sampling to neighbors...
			const ABitMap2	origSampMap = sampMap;
			for (int yoff = -(1<<(sp2-1));
					yoff <= 1<<(sp2-1); yoff += 1<<sp2)
			    for (int xoff = -(1<<(sp2-1));
			    		xoff <= 1<<(sp2-1); xoff += 1<<sp2) {
				ABitMap2	stamp = origSampMap;
				stamp.Shift(xoff, yoff);
				sampMap |= stamp;
			}		// ...but don't resample what's done
			sampMap -= doneSamples;
					// interpolate smooth regions
			fillMap -= sampMap;
			for (x = y = 0; fillMap.Find(&x, &y); x++)
				FillSquare(x, y, noff);
			doneSamples |= doneMap;
		}			// compute required ray samples
		for (x = y = 0; sampMap.Find(&x, &y); x++)
			if (!ComputePixel(x, y))
				return false;
		doneSamples |= sampMap;	// samples now done or at least queued
fprintf(stderr, "Sampled %ld pixels at (sp2,layer)=(%d,%d)\n",
(long)sampMap.SumTotal(), sp2, layer);
fprintf(stderr, "\t%ld pixels (%.3f%%) completed (+%ld in process)\n",
(long)doneMap.SumTotal(), 100.*doneMap.SumTotal()/doneMap.Width()/doneMap.Height(),
(long)(doneSamples.SumTotal()-doneMap.SumTotal()));
		sp2 -= layer++ & 1;	// next denser sampling
	}
	if (FlushQueue() < 0)		// make sure we got everyone
		return false;
	x = y = 0;
	if (doneMap.Find(&x, &y, false)) {
		sprintf(errmsg, "missed %ld tile pixels, e.g. (%d,%d)",
				(long)doneMap.Width()*doneMap.Height() -
					doneMap.SumTotal(), x, y);
		error(WARNING, errmsg);
	}
	if ((prCB != NULL) & (barPix == NULL))
		(*prCB)(100.);
	return true;
}

/*
 * Render the specified tile in frame
 * Tile pixels are contiguous unless ystride != 0
 * Tiles numbered from upper-left at (0,0)
 * Pixel type influenced by this->prims assignment
 */
bool
RpictSimulManager::RenderTile(COLORV *rp, int ystride, float *zp, const int *tile)
{
	if (!rp | (GetWidth() <= 0) | (GetHeight() <= 0) | !vw.type)
		return false;
	if (!ystride)			// contiguous rows?
		ystride = TWidth();
	pacc.Init(rp, ystride, zp);
	if (prims == xyzprims)
		pacc.SetColorSpace(RDTxyz);
	else if (prims)
		pacc.SetColorSpace(RDTrgb, prims);
		
	return SetTile(tile) && RenderRect();
}

// Same but store as common-exponent COLR or SCOLR
bool
RpictSimulManager::RenderTile(COLRV *bp, int ystride, float *zp, const int *tile)
{
	if (!bp | (GetWidth() <= 0) | (GetHeight() <= 0) | !vw.type)
		return false;
	if (!ystride)			// contiguous rows?
		ystride = TWidth();
	pacc.Init(bp, ystride, zp);
	if (prims == xyzprims)
		pacc.SetColorSpace(RDTxyze);
	else if (prims)
		pacc.SetColorSpace(RDTrgbe, prims);

	return SetTile(tile) && RenderRect();
}

// Same but also use 16-bit encoded depth buffer
bool
RpictSimulManager::RenderTile(COLRV *bp, int ystride, short *dp, const int *tile)
{
	if (!bp | (GetWidth() <= 0) | (GetHeight() <= 0) | !vw.type)
		return false;
	if (!ystride)			// contiguous rows?
		ystride = TWidth();
	pacc.Init(bp, ystride, dp);
	if (prims == xyzprims)
		pacc.SetColorSpace(RDTxyze);
	else if (prims)
		pacc.SetColorSpace(RDTrgbe, prims);

	return SetTile(tile) && RenderRect();
}

// Allocate a new render bar
void
RpictSimulManager::NewBar(int ht)
{
	delete [] barPix;
	delete [] barDepth;
	if (ht > GetHeight()) ht = GetHeight();
	if ((ht <= 0) | (GetWidth() <= 0)) {
		doneMap.NewBitMap(0,0);
		pacc.Init();
		barPix = NULL; barDepth = NULL;
		return;
	}
	thvres[0] = GetWidth();
	thvres[1] = ht;
	const int	NC = prims ? 3 : NCSAMP;
	barPix = new COLORV [ht*thvres[0]*NC];
	barDepth = new float [ht*thvres[0]];
	pacc.Init(barPix + (ht-1)*thvres[0]*NC,
			-thvres[0], barDepth + (ht-1)*thvres[0]);
	if (prims == xyzprims)
		pacc.SetColorSpace(RDTxyz);
	else if (prims)
		pacc.SetColorSpace(RDTrgb, prims);

	doneMap.NewBitMap(TWidth(), THeight());
}

// Shift render bar area the specified amount down the frame
bool
RpictSimulManager::LowerBar(int v)
{
	if (v <= 0) return !v;
	if (!barPix | !barDepth | (v > THeight()) | !tvw.type)
		return false;
	tvw.voff -= double(v)/GetHeight();
	ptvw.voff -= double(v)/GetHeight();
	if (v == THeight()) {
		doneMap.ClearBitMap();
		return true;
	}
	const int	NC = pacc.NC();
	doneMap.Shift(0, v, false);		// lift finished pixel samples
	memmove(barPix, barPix + NC*TWidth()*v,
			sizeof(COLORV)*NC*TWidth()*(THeight()-v));
	memmove(barDepth, barDepth + TWidth()*v,
			sizeof(float)*TWidth()*(THeight()-v));
	return true;
}

// Continue rendering from the specified position
bool
RpictSimulManager::RenderBelow(int ytop, const int vstep, FILE *pfp, const int dt, FILE *dfp)
{
	if (ytop <= 0)
		return true;
	ptvw = pvw;				// set starting bar's view
	tvw = vw;
	const char *	verr = cropview(&tvw, 0., double(ytop-THeight())/GetHeight(),
						1., double(ytop)/GetHeight());
	if (verr) {
		sprintf(errmsg, "illegal render bar below y=%d: %s", ytop, verr);
		error(INTERNAL, errmsg);
		return false;
	}
	if (!ptvw.type | (mblur <= FTINY) ||
			cropview(&ptvw, 0., double(ytop-THeight())/GetHeight(),
						1., double(ytop)/GetHeight()))
		ptvw.type = 0;
						// set up spectral sampling
	if (setspectrsamp(CNDX, WLPART) <= 0) {
		error(USER, "unsupported spectral sampling");
		return false;
	}
	COLORV **	parr = NULL;		// set up tiny source drawing
	float **	zarr = NULL;
	if (!ptvw.type && directvis && dblur <= FTINY) {
		parr = new COLORV * [THeight()];
		zarr = new float * [THeight()];
		for (int n = THeight(); n-- > 0; ) {
			parr[THeight()-1-n] = barPix + pacc.NC()*TWidth()*n;
			zarr[THeight()-1-n] = barDepth + TWidth()*n;
		}
		ourview = vw; hres = GetWidth(); vres = GetHeight();
		init_drawsources(psample);
	}
	int		lastOut = ytop;		// render down frame
	while (ytop > 0) {
		if (ytop < THeight())		// mark what we won't do as finished
			doneMap.ClearRect(0, 0, TWidth(), THeight()-ytop, true);
		if (prCB)
			(*prCB)(100.*(GetHeight()-ytop)/GetHeight());
		if (!RenderRect())		// render this bar
			return false;
		int	nlines = lastOut - ytop + THeight();
		if (nlines > ytop)
			nlines = ytop;
		else if (parr)			// drawing sources?
			drawsources(parr, prims, zarr,
					0, hres, lastOut-nlines, nlines);

		if (dfp) {			// write out depth scanlines?
			const float *	dp = barDepth + TWidth()*(ytop-lastOut);
			if (RDTdepthT(dt) == RDTdshort) {
				for (int n = TWidth()*nlines; n-- > 0; dp++)
					if (putint(depth2code(*dp, pacc.refDepth), 2, dfp) == EOF)
						error(SYSTEM, "cannot write 16-bit depth buffer");
			} else if (putbinary(dp, sizeof(float), TWidth()*nlines, dfp)
						!= TWidth()*nlines)
				error(SYSTEM, "cannot write raw depth buffer");
		}
		COLORV *	bpos = barPix + pacc.NC()*TWidth()*(ytop-lastOut);
		while (nlines-- > 0) {		// write pixel scanlines
			switch (RDTcolorT(dt)) {
			case RDTrgbe:
			case RDTxyze:
				if (fwritescan((COLOR *)bpos, TWidth(), pfp) < 0)
					error(SYSTEM, "cannot write RGBE/XYZE output");
				break;
			case RDTscolr:
				if (fwritesscan(bpos, pacc.NC(), TWidth(), pfp) < 0)
					error(SYSTEM, "cannot write SCOLOR output");
				break;
			case RDTrgb:
			case RDTxyz:
			case RDTscolor:
				if (putbinary(bpos, sizeof(COLORV)*pacc.NC(), TWidth(), pfp)
						!= TWidth())
					error(SYSTEM, "cannot write SCOLOR output");
				break;
			default:
				error(INTERNAL, "missing output color type in RenderBelow()");
				break;
			}
			bpos += pacc.NC()*TWidth();
			--lastOut;
		}				// flush each scan bar
		if (fflush(pfp) == EOF || (dfp && fflush(dfp) == EOF))
			error(SYSTEM, "output write error");
						// advance down the frame
		if (lastOut > 0 && !LowerBar(vstep))
			return false;
		ytop -= vstep;
	}
	delete [] parr;
	delete [] zarr;
	if (prCB)
		(*prCB)(100.);
	return true;
}

/*
 * Render and write a frame to the named file
 * Include any header lines set prior to call
 * Picture file must not already exist
 * Picture to stdout if pfname==NULL
 * Depth written to a command if dfname[0]=='!'
 */
RenderDataType
RpictSimulManager::RenderFrame(const char *pfname, RenderDataType dt, const char *dfname)
{
	int	fd = 1;
	FILE *	pfp = NULL;
	FILE *	dfp = NULL;

	if (!RDTcolorT(dt))
		error(INTERNAL, "missing pixel output type in RenderFrame()");
	if (NCSAMP == 3) {
		if (RDTcolorT(dt) == RDTscolr)
			dt = RDTnewCT(dt, prims==xyzprims ? RDTxyze : RDTrgbe);
		else if (RDTcolorT(dt) == RDTscolor)
			dt = RDTnewCT(dt, prims==xyzprims ? RDTxyz : RDTrgb);
	}
	if (!RDTdepthT(dt) ^ !dfname)
		error(INTERNAL, "depth output requires file name and type in RenderFrame()");
	if (pfname) {				// open picture output file
		if (pfname[0] == '!') {
			error(INTERNAL, "writing picture to a command not supported");
			return RDTnone;
		}
		fd = open(pfname, O_WRONLY|O_CREAT|O_EXCL, 0666);
	}
	if (fd < 0) {
		if ((frameNo <= 0) | (errno != EEXIST)) {
			sprintf(errmsg, "cannot open picture file '%s'", pfname);
			error(SYSTEM, errmsg);
		}
		return RDTnone;			// expected in parallel sequence
	}
	if (fd == 1)
		pfp = stdout;
	else if (!(pfp = fdopen(fd, "w")))
		error(SYSTEM, "failure calling fdopen()");
	SET_FILE_BINARY(pfp);			// write picture header
	if ((pfp != stdout) | (frameNo <= 1)) {
		newheader("RADIANCE", pfp);
		fputs(GetHeader(), pfp);
	}
	fputs(VIEWSTR, pfp); fprintview(&vw, pfp); fputc('\n', pfp);
	if (frameNo > 0)
		fprintf(pfp, "FRAME=%d\n", frameNo);
	double	pasp = viewaspect(&vw) * GetWidth() / GetHeight();
	if (!FABSEQ(pasp, 1.0))
		fputaspect(pasp, pfp);
	fputnow(pfp);
	switch (RDTcolorT(dt)) {		// set primaries and picture format
	case RDTrgbe:
		if (!prims | (prims == xyzprims)) prims = stdprims;
		fputprims(prims, pfp);
		fputformat(COLRFMT, pfp);
		break;
	case RDTxyze:
		prims = xyzprims;
		fputformat(CIEFMT, pfp);
		break;
	case RDTscolr:
		prims = NULL;
		fputwlsplit(WLPART, pfp);
		fputncomp(NCSAMP, pfp);
		fputformat(SPECFMT, pfp);
		break;
	case RDTrgb:
		if (!prims | (prims == xyzprims)) prims = stdprims;
		fputprims(prims, pfp);
		fputncomp(3, pfp);
		fputendian(pfp);
		fputformat("float", pfp);
		break;
	case RDTxyz:
		prims = xyzprims;
		fputprims(prims, pfp);
		fputncomp(3, pfp);
		fputendian(pfp);
		fputformat("float", pfp);
		break;
	case RDTscolor:
		prims = NULL;
		fputwlsplit(WLPART, pfp);
		fputncomp(NCSAMP, pfp);
		fputendian(pfp);
		fputformat("float", pfp);
		break;
	default:;
	}
	fputc('\n', pfp);			// end picture header
	fprtresolu(GetWidth(), GetHeight(), pfp);
	if (dfname) {
		if (dfname[0] == '!')
			dfp = popen(dfname+1, "w");
		else
			dfp = fopen(dfname, "w");
		if (!dfp) {
			sprintf(errmsg, "cannot open depth output '%s'", dfname);
			error(SYSTEM, errmsg);
			return RDTnone;
		}
		SET_FILE_BINARY(dfp);
	}
	if (RDTdepthT(dt) == RDTdshort) {	// write header for 16-bit depth?
		newheader("RADIANCE", dfp);
		fputs(GetHeader(), dfp);
		fputs(VIEWSTR, dfp); fprintview(&vw, dfp); fputc('\n', dfp);
		fputs(DEPTHSTR, dfp); fputs(dunit, dfp); fputc('\n', dfp);
		fputformat(DEPTH16FMT, dfp);
		fputc('\n', dfp);		// end-of-info
		fprtresolu(GetWidth(), GetHeight(), dfp);
	}
	const int	bheight = (psample > 1) ? int(2*psample+.99) : 4;
	const int	vstep =  bheight >> (psample > 1);

	NewBar(bheight);			// render frame if we can
	if (!RenderBelow(GetHeight(), vstep, pfp, dt, dfp)) {
		fclose(pfp);
		if (dfp) (dfname[0] == '!') ? pclose(dfp) : fclose(dfp);
		Cleanup();
		return RDTnone;
	}
	NewBar();				// clean up and return
	if (pfp != stdout)
		fclose(pfp);
	if (dfp) {
		if (dfname[0] == '!') {
			int	status = pclose(dfp);
			if (status) {
				sprintf(errmsg, "depth output (%s) error status: %d",
						dfname, status);
				error(USER, errmsg);
				return RDTnone;
			}
		} else
			fclose(dfp);
	}
	return dt;
}

// passed struct for header line callback
struct HeaderInfo {
	char		fmt[MAXFMTLEN];
	char		depth_unit[32];
	int		ncomp;
	RGBPRIMS	prims;
	VIEW		vw;
	bool		gotWL;
	bool		gotprims;
	bool		gotview;
	bool		endianMatch;

			HeaderInfo() {
				strcpy(fmt, "MISSING");
				depth_unit[0] = '\0';
				ncomp = 3;
				vw = stdview;
				gotWL = false;
				gotprims = false;
				gotview = false;
				endianMatch = true;
			}
};

// helper function checks header line and records req. info.
static int
head_check(char *s, void *p)
{
	HeaderInfo *	hp = (HeaderInfo *)p;
	int		rval;

	if (isncomp(s)) {
		hp->ncomp = ncompval(s);
		return 1;
	}
	if (iswlsplit(s)) {
		hp->gotWL = wlsplitval(WLPART, s);
		return 1;
	}
	if (isprims(s)) {
		hp->gotprims = primsval(hp->prims, s);
		return 1;
	}
	if (isview(s)) {
		hp->gotview |= (sscanview(&hp->vw, s) > 0);
		return 1;
	}
	if (!strncmp(s, DEPTHSTR, LDEPTHSTR)) {
		strlcpy(hp->depth_unit, s+LDEPTHSTR, sizeof(hp->depth_unit));
		char *	cp = hp->depth_unit;
		while (*cp) cp++;
		while (cp > hp->depth_unit && isspace(cp[-1])) cp--;
		*cp = '\0';
		return 1;
	}
	if ((rval = isbigendian(s)) >= 0) {
		hp->endianMatch = (rval == nativebigendian());
		return 1;
	}
	if (formatval(hp->fmt, s))
		return 1;
	return 0;
}

// Resume partially finished rendering
// Picture file must exist
RenderDataType
RpictSimulManager::ResumeFrame(const char *pfname, const char *dfname)
{
	if (!pfname || pfname[0] == '!')
		return RDTnone;

	RenderDataType	dt = RDTnone;
	FILE *		dfp = NULL;
	FILE *		pfp = fopen(pfname, "r+");
	if (!pfp) {
		sprintf(errmsg, "cannot reopen output picture '%s'", pfname);
		error(SYSTEM, errmsg);
		return RDTnone;
	}
	SET_FILE_BINARY(pfp);
	HeaderInfo	hinfo;		// read header information & dimensions
	RESOLU		res;
	if (getheader(pfp, head_check, &hinfo) < 0) {
		fclose(pfp);
		return RDTnone;
	}
	if (!fgetsresolu(&res, pfp) || res.rt != PIXSTANDARD) {
		sprintf(errmsg, "missing/bad resolution for '%s'", pfname);
		error(USER, errmsg);
		fclose(pfp);
		return RDTnone;
	}
	if (!hinfo.gotview) {
		sprintf(errmsg, "missing view for '%s'", pfname);
		error(USER, errmsg);
		fclose(pfp);
		return RDTnone;
	}
	if (hinfo.ncomp < 3) {
		sprintf(errmsg, "bad # components (%d) in '%s'", hinfo.ncomp, pfname);
		error(USER, errmsg);
		fclose(pfp);
		return RDTnone;
	}
	int	bytesPer = 0;		// complicated part to set rendering/output space
	if (!strcmp(hinfo.fmt, COLRFMT)) {
		prims = hinfo.prims;
		int	n = 8*hinfo.gotprims;
		while (n--)
			if (!FABSEQ(hinfo.prims[0][n], stdprims[0][n]))
				break;
		if (n < 0)
			prims = stdprims;
		dt = RDTnewCT(dt, RDTrgbe);
	} else if (!strcmp(hinfo.fmt, CIEFMT)) {
		prims = xyzprims;
		dt = RDTnewCT(dt, RDTxyze);
	} else if (!strcmp(hinfo.fmt, SPECFMT)) {
		if ((hinfo.ncomp <= 3) | (hinfo.ncomp > MAXCSAMP)) {
			sprintf(errmsg, "incompatible sample count (%d) in '%s'",
					hinfo.ncomp, pfname);
			error(USER, errmsg);
			fclose(pfp);
			return RDTnone;
		}
		NCSAMP = hinfo.ncomp;	// overrides global setting
		prims = NULL;
		dt = RDTnewCT(dt, RDTscolr);
		bytesPer = hinfo.ncomp + 1;	// XXX assumes no compression
	} else if (!strcmp(hinfo.fmt, "float")) {
		if (!hinfo.endianMatch) {
			sprintf(errmsg, "incompatible byte ordering in '%s'", pfname);
			error(USER, errmsg);
			fclose(pfp);
			return RDTnone;
		}
		if (hinfo.ncomp == 3) {
			prims = hinfo.prims;		// custom primaries?
			int	n = 8*hinfo.gotprims;
			while (n--)
				if (!FABSEQ(hinfo.prims[0][n], stdprims[0][n]))
					break;
			if (n < 0)			// standard primaries?
				prims = stdprims;
			else if (hinfo.gotprims) {	// or check if XYZ
				for (n = 8; n--; )
					if (!FABSEQ(prims[0][n], xyzprims[0][n]))
						break;
				if (n < 0)
					prims = xyzprims;
			}
			if (prims == xyzprims)
				dt = RDTnewCT(dt, RDTxyz);
			else
				dt = RDTnewCT(dt, RDTrgb);
		} else {
			NCSAMP = hinfo.ncomp;	// overrides global setting
			prims = NULL;
			dt = RDTnewCT(dt, RDTscolor);
		}
		bytesPer = sizeof(float)*hinfo.ncomp;
	} else {
		sprintf(errmsg, "unknown format (%s) for '%s'", hinfo.fmt, pfname);
		error(USER, errmsg);
		fclose(pfp);
		return RDTnone;
	}
	vw.type = 0;				// set up new (unreferenced) frame
	frameNo = 0;
	int	hvdim[2] = {res.xr, res.yr};
	double	noAdj = 0;
	if (!NewFrame(hinfo.vw, hvdim, &noAdj) ||
			(hvdim[0] != res.xr) | (hvdim[1] != res.yr)) {
		fclose(pfp);
		return RDTnone;
	}
	long	dataStart = ftell(pfp);		// picture starting point
	if (dataStart < 0) {
		sprintf(errmsg, "cannot seek on '%s'", pfname);
		error(SYSTEM, errmsg);
		fclose(pfp);
		return RDTnone;
	}
	long	doneScans = 0;
	if (bytesPer) {				// fixed-width records?
		fseek(pfp, 0, SEEK_END);
		long	dataEnd = ftell(pfp);
		doneScans = (dataEnd - dataStart)/(bytesPer*GetWidth());
		if (dataEnd-dataStart > bytesPer*GetWidth()*doneScans)
			fseek(pfp, dataStart + bytesPer*GetWidth()*doneScans, SEEK_SET);
	} else {				// else get compressed scanlines
		COLR *	scan = (COLR *)tempbuffer(sizeof(COLR)*GetWidth());
		while (freadcolrs(scan, GetWidth(), pfp) >= 0)
			++doneScans;
		if (!feof(pfp)) {
			sprintf(errmsg, "error reading compressed scanline from '%s'", pfname);
			error(USER, errmsg);
			fclose(pfp);
			return RDTnone;
		}
	}
	if (doneScans >= GetHeight()) {		// nothing left to do?
		sprintf(errmsg, "output file '%s' is already complete", pfname);
		error(WARNING, errmsg);
		fclose(pfp);
		return dt;
	}
	if (!doneScans) {
		sprintf(errmsg, "restarting empty frame '%s'", pfname);
		error(WARNING, errmsg);
	}
	if (dfname) {				// append depth file, too?
		if (dfname[0] == '!') {
			error(USER, "depth data cannot be reloaded from command");
			fclose(pfp);
			return RDTnone;
		}
		dfp = fopen(dfname, "a");
		if (!dfp) {
			sprintf(errmsg, "cannot reopen depth file '%s'", dfname);
			error(SYSTEM, errmsg);
			fclose(pfp);
			return RDTnone;
		}
		SET_FILE_BINARY(dfp);
		const long	dflen = ftell(dfp);
		if (dflen != sizeof(float)*GetWidth()*doneScans) {
			fclose(dfp);
			dfp = fopen(dfname, "r+");
			if (!dfp) return RDTnone;	// WTH?
			SET_FILE_BINARY(dfp);
		}
		if (dflen < sizeof(float)*GetWidth()*doneScans) {
			HeaderInfo	dinfo;
			if (getheader(dfp, head_check, &dinfo) < 0)
				sprintf(errmsg, "bad header in encoded depth file '%s'",
						dfname);
			else if (strcmp(dinfo.fmt, DEPTH16FMT))
				sprintf(errmsg, "wrong format (%s) for depth file '%s'",
						dinfo.fmt, dfname);
			else if (!SetReferenceDepth(dinfo.depth_unit))
				sprintf(errmsg, "bad/missing reference depth (%s) in '%s'",
						dinfo.depth_unit, dfname);
			else if (!fscnresolu(hvdim, hvdim+1, dfp) ||
					(hvdim[0] != GetWidth()) | (hvdim[1] != GetHeight()))
				sprintf(errmsg, "bad/mismatched resolution in '%s'",
						dfname);
			else
				errmsg[0] = '\0';

			if (errmsg[0]) {
				error(USER, errmsg);
				fclose(dfp);
				fclose(pfp);
				return RDTnone;
			}
			const long	dStart = ftell(dfp);
			if (dflen-dStart < 2*GetWidth()*doneScans) {
				sprintf(errmsg, "missing %ld depths in '%s'",
					(long)GetWidth()*doneScans - (dflen-dStart)/2,
					dfname);
				error(WARNING, errmsg);
			}
			fseek(dfp, dStart + 2*GetWidth()*doneScans, SEEK_SET);
			dt = RDTnewDT(dt, RDTdshort);
		} else {
			if (dflen > sizeof(float)*GetWidth()*doneScans)
				fseek(dfp, sizeof(float)*GetWidth()*doneScans, SEEK_SET);
			dt = RDTnewDT(dt, RDTdfloat);
		}
	}
	int	bheight = (psample > 1) ? int(2*psample+.99) : 4;
	if (bheight > GetHeight()-doneScans)
		bheight = GetHeight()-doneScans;
	int	vstep =  bheight >> (psample > 1);
	vstep += !vstep;

	NewBar(bheight);			// render remainder if we can
	if (!RenderBelow(GetHeight()-doneScans, vstep, pfp, dt, dfp)) {
		fclose(pfp);
		if (dfp) fclose(dfp);
		Cleanup();
		return RDTnone;
	}
	NewBar();				// close up and return success
	fclose(pfp);
	if (dfp) fclose(dfp);
	return dt;
}
