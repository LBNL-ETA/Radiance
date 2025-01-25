#ifndef lint
static const char RCSid[] = "$Id: RpictSimulManager.cpp,v 2.15 2025/01/25 04:57:27 greg Exp $";
#endif
/*
 *  RpictSimulManager.cpp
 *
 *	Rpict simulation manager implementation
 *
 *  Created by Greg Ward on 07/11/2024.
 */

#define DEBUG	1	// XXX temporary!

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
		error(INTERNAL, "botched color space type in SetPixel()");
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
	case RDTscolr:
	case RDTscolor:
		primp = NULL;
		break;
	default:
		error(INTERNAL, "botched color space type in SetColorSpace()");
	}
	dtyp = RDTnewCT(dtyp, cs);
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
	if (&v == &vw) {
		pvw.type = 0;
	} else {
		pvw = vw;		// save previous view for motion blur
		vw = v;
	}
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
	const int		ty = (r->rno-1) / rsp->TWidth();
	const int		tx = r->rno-1 - (RNUMBER)ty*rsp->TWidth();

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
	DCHECK(doneMap.OffBitMap(x,y),
			CONSISTENCY, "illegal pixel index in ComputPixel()");
	int	i;
	FVECT	rodir[2];
	double	hpos = (x+pixjitter())/TWidth();
	double	vpos = (y+pixjitter())/THeight();
	double	dlim = viewray(rodir[0], rodir[1], &tvw, hpos, vpos);
	if (dlim < -FTINY) {	// off view?
		pacc.SetPixel(x, y, scblack);
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

	return EnqueueRay(rodir[0], rodir[1], (RNUMBER)y*TWidth()+x+1);
}

// Check if neighbor differences are below pixel sampling threshold
bool
RpictSimulManager::BelowSampThresh(const int x, const int y, const int noff[4][2]) const
{
	SCOLOR	pval[4];
	float	dist[4];
	int	i, j;

	for (i = 4; i--; ) {		// get pixels from tile store
		const int	px = x + noff[i][0];
		const int	py = y + noff[i][1];
		if (!doneMap.Check(px, py) ||
				!pacc.GetPixel(px, py, pval[i], &dist[i]))
			return false;
	}
					// do pairwise comparisons
	for (i = (pacc.DepthType() != RDTnone)*4; --i > 0; )
	    for (j = i; j--; )
	        if ((dist[i] - dist[j] > maxdiff*dist[j]) |
	        		(dist[j] - dist[i] > maxdiff*dist[i]))
			return false;
	if (pacc.NC() > 3) {
	    for (i = 4; --i; )
	    	for (j = i; j--; )
		    if (sbigsdiff(pval[i], pval[j], maxdiff))
			return false;
	} else {
	    for (i = 4; --i; )
	    	for (j = i; j--; )
		    if (bigdiff(pval[i], pval[j], maxdiff))
			return false;
	}
	return true;
}

// Fill an interior square patch with interpolated values
void
RpictSimulManager::FillSquare(const int x, const int y, const int noff[4][2])
{
	SCOLOR	pval[4];
	float	dist[4];
	int	i, j;
					// assumes 4 corners are valid!
	for (i = 4; i--; ) {
		DCHECK(!doneMap.Check(x+noff[i][0], y+noff[i][1]),
			CONSISTENCY, "inclusion of bad pixel in FillSquare()");
		pacc.GetPixel(x+noff[i][0], y+noff[i][1], pval[i], &dist[i]);
	}
	i = abs(noff[1][0]-noff[0][0]);
	j = abs(noff[1][1]-noff[0][1]);	// i==j for diamond fill
	const int	slen =  (i > j) ? i : j;
	const bool	spectr = (pacc.NC() > 3);
	for (i = slen+1 + (i==j)*slen; i--; ) {
	    const double	c1 = (i>slen ? i-slen-.5 : (double)i)/slen;
	    for (j = slen + (i<=slen); j--; ) {
	    	const double	c2 = (j + (i>slen)*.5)/slen;
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
SetQuincunx(ABitMap2 *bmp2, int noff[4][2], const int spc, bool odd, int x0, int y)
{
	if (odd) {			// order neighbors CCW
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
	int	nsteps;			// non-negative range
	if (x0 < -(spc>>1)) {
		nsteps = (spc-1 - x0 - (spc>>1))/spc;
		x0 += nsteps*spc;
	}
	if (y < 0) {			// get past y==0
		nsteps = ((spc>>1)-1 - y)/(spc>>1);
		y += nsteps*(spc>>1);
		odd ^= nsteps&1;
	}
	while (y < bmp2->Height()) {
	    for (int x = x0 + odd*(spc>>1); x < bmp2->Width(); x += spc)
	    	bmp2->Set(x, y);
	    y += spc>>1;
	    odd = !odd;
	}
}

// Render (or finish rendering) current tile
bool
RpictSimulManager::RenderRect(const int x0, const int y0)
{
	if (!tvw.type || !Ready()) {
		error(INTERNAL, "need octree and view for RenderRect()");
		return false;
	}
	ABitMap2	doneSamples = doneMap;
	int		sp2 = ceil(log2((TWidth()>THeight() ? TWidth() : THeight()) - 1.));
	int		layer = 0;
	int		x, y;
	while (sp2 > 0) {
		ABitMap2	sampMap(TWidth(), THeight());
		int		noff[4][2];
		if ((prCB != NULL) & (barPix == NULL))
			(*prCB)(100.*doneMap.SumTotal()/doneMap.Width()/doneMap.Height());
		SetQuincunx(&sampMap, noff, 1<<sp2, layer&1, x0, y0);
		sampMap -= doneSamples;	// avoid resampling pixels
		// Are we into adaptive sampling realm?
		if (noff[0][0]*noff[0][0] + noff[0][1]*noff[0][1] < psample*psample) {
			if (FlushQueue() < 0)	// need results to check thresholds
				return false;
			ABitMap2	fillMap = sampMap;
			for (x = y = 0; sampMap.Find(&x, &y); x++)
				if (BelowSampThresh(x, y, noff))
					sampMap.Reset(x, y);
					// spread sampling to neighbors...
			const ABitMap2	origSampMap = sampMap;
			for (x = 4; x--; ) {
				ABitMap2	stamp = origSampMap;
				stamp.Shift(noff[x][0] + noff[(x+1)&3][0],
						noff[x][1] + noff[(x+1)&3][1]);
				sampMap |= stamp;
			}		// ...but don't resample what's done
			sampMap -= doneSamples;
					// interpolate smooth regions
			fillMap -= sampMap;
			fillMap -= doneSamples;
			for (x = y = 0; fillMap.Find(&x, &y); x++)
				FillSquare(x, y, noff);
			doneSamples |= doneMap;
		}			// compute required ray samples
		for (x = y = 0; sampMap.Find(&x, &y); x++)
			if (!ComputePixel(x, y)) {
				sprintf(errmsg, "ComputePixel(%d,%d) failed", x, y);
				error(WARNING, errmsg);
				return false;
			}
		doneSamples |= sampMap;	// samples now done or at least queued
		sp2 -= layer++ & 1;	// next denser sampling
	}
	if (FlushQueue() < 0)		// compute stragglers
		return false;
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

// Back to float color with 16-bit depth
bool
RpictSimulManager::RenderTile(COLORV *rp, int ystride, short *dp, const int *tile)
{
	if (!rp | (GetWidth() <= 0) | (GetHeight() <= 0) | !vw.type)
		return false;
	if (!ystride)			// contiguous rows?
		ystride = TWidth();
	pacc.Init(rp, ystride, dp);
	if (prims == xyzprims)
		pacc.SetColorSpace(RDTxyz);
	else if (prims)
		pacc.SetColorSpace(RDTrgb, prims);

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
RpictSimulManager::LowerBar(int v, int ytop)
{
	if (!barPix | !barDepth | (v > THeight()) | !tvw.type)
		return false;
	if (v <= 0) return !v;
	if ((ytop -= v) <= 0)
		return true;
	tvw.voff -= double(v)/THeight();
	ptvw.voff -= double(v)/THeight();
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
	if (ytop < THeight())			// mark what we won't do as finished
		doneMap.ClearRect(0, 0, TWidth(), THeight()-ytop, true);
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
						// update spectral sampling
	int	rv = setspectrsamp(CNDX, WLPART);
	if (rv < 0) {
		error(USER, "unsupported spectral sampling");
		return false;
	}
	if (!rv & (RDTcolorT(dt) != RDTscolor) & (RDTcolorT(dt) != RDTscolr))
		error(WARNING, "spectral range incompatible with color output");
	COLORV **	parr = NULL;		// set up tiny source drawing
	float **	zarr = NULL;
	if (!ptvw.type && directvis && (dblur <= FTINY) & (mblur <= FTINY)) {
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
		if (prCB)
			(*prCB)(100.*(GetHeight()-ytop)/GetHeight());
		if (!RenderRect(0, THeight()-ytop))	// render this bar
			return false;
		int	nlines = lastOut - ytop + vstep;
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
				error(INTERNAL, "botched output color type in RenderBelow()");
				break;
			}
			bpos += pacc.NC()*TWidth();
			--lastOut;
		}				// flush each scan bar
		if (fflush(pfp) == EOF || (dfp && fflush(dfp) == EOF))
			error(SYSTEM, "output write error");
						// advance down the frame
		if (lastOut > 0 && !LowerBar(vstep, ytop))
			return false;
		ytop -= vstep;
	}
	delete [] parr;
	delete [] zarr;
	if (prCB)
		(*prCB)(100.);
	return true;
}

// Open new output picture file (and optional depth file)
RenderDataType
RpictSimulManager::NewOutput(FILE *pdfp[2], const char *pfname,
				RenderDataType dt, const char *dfname)
{
	pdfp[0] = pdfp[1] = NULL;
	if (!RDTcolorT(dt))
		error(INTERNAL, "missing color output type in NewOutput()");
	if (NCSAMP == 3) {
		if (RDTcolorT(dt) == RDTscolr)
			dt = RDTnewCT(dt, prims==xyzprims ? RDTxyze : RDTrgbe);
		else if (RDTcolorT(dt) == RDTscolor)
			dt = RDTnewCT(dt, prims==xyzprims ? RDTxyz : RDTrgb);
	}
	if (!RDTdepthT(dt) ^ !dfname)
		error(INTERNAL, "depth output requires file name and type in NewOutput()");
	int	fd = 1;
	if (pfname) {				// open picture output file
		if (pfname[0] == '!') {
			error(INTERNAL, "writing picture to a command not supported");
			return RDTnone;
		}
		fd = open(pfname, O_RDWR|O_CREAT|O_EXCL, 0666);
	}
	if (fd < 0) {
		if ((frameNo <= 0) | (errno != EEXIST)) {
			sprintf(errmsg, "cannot open picture file '%s'", pfname);
			error(SYSTEM, errmsg);
		}
		return RDTnone;			// may be expected in sequence run
	}
	if (fd == 1)
		pdfp[0] = stdout;
	else if (!(pdfp[0] = fdopen(fd, "w+")))
		error(SYSTEM, "failure calling fdopen()");
	SET_FILE_BINARY(pdfp[0]);		// write picture header
	if ((pdfp[0] != stdout) | (frameNo <= 1)) {
		newheader("RADIANCE", pdfp[0]);
		fputs(GetHeadStr(), pdfp[0]);
	}
	fputs(VIEWSTR, pdfp[0]); fprintview(&vw, pdfp[0]); fputc('\n', pdfp[0]);
	if (frameNo > 0)
		fprintf(pdfp[0], "FRAME=%d\n", frameNo);
	double	pasp = viewaspect(&vw) * GetWidth() / GetHeight();
	if ((0.99 > pasp) | (pasp > 1.01))
		fputaspect(pasp, pdfp[0]);
	fputnow(pdfp[0]);
	switch (RDTcolorT(dt)) {		// set primaries and picture format
	case RDTrgbe:
		if (!prims | (prims == xyzprims)) prims = stdprims;
		fputprims(prims, pdfp[0]);
		fputformat(COLRFMT, pdfp[0]);
		break;
	case RDTxyze:
		prims = xyzprims;
		fputformat(CIEFMT, pdfp[0]);
		break;
	case RDTscolr:
		prims = NULL;
		fputwlsplit(WLPART, pdfp[0]);
		fputncomp(NCSAMP, pdfp[0]);
		fputformat(SPECFMT, pdfp[0]);
		break;
	case RDTrgb:
		if (!prims | (prims == xyzprims)) prims = stdprims;
		fputprims(prims, pdfp[0]);
		fputncomp(3, pdfp[0]);
		fputendian(pdfp[0]);
		fputformat("float", pdfp[0]);
		break;
	case RDTxyz:
		prims = xyzprims;
		fputprims(prims, pdfp[0]);
		fputncomp(3, pdfp[0]);
		fputendian(pdfp[0]);
		fputformat("float", pdfp[0]);
		break;
	case RDTscolor:
		prims = NULL;
		fputwlsplit(WLPART, pdfp[0]);
		fputncomp(NCSAMP, pdfp[0]);
		fputendian(pdfp[0]);
		fputformat("float", pdfp[0]);
		break;
	default:;	// pro forma - caught this above
	}
	fputc('\n', pdfp[0]);			// flush picture header
	if (fflush(pdfp[0]) == EOF) {
		sprintf(errmsg, "cannot write header to picture '%s'", pfname);
		error(SYSTEM, errmsg);
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	if (dfname) {				// open depth output
		if (dfname[0] == '!')
			pdfp[1] = popen(dfname+1, "w");
		else
			pdfp[1] = fopen(dfname, "w+");
		if (!pdfp[1]) {
			sprintf(errmsg, "cannot open depth output '%s'", dfname);
			error(SYSTEM, errmsg);
			fclose(pdfp[0]);
			pdfp[0] = NULL;
			return RDTnone;
		}
		SET_FILE_BINARY(pdfp[1]);
	}
	if (RDTdepthT(dt) == RDTdshort) {	// write header for 16-bit depth?
		newheader("RADIANCE", pdfp[1]);
		fputs(GetHeadStr(), pdfp[1]);
		fputs(VIEWSTR, pdfp[1]); fprintview(&vw, pdfp[1]); fputc('\n', pdfp[1]);
		fputs(DEPTHSTR, pdfp[1]); fputs(dunit, pdfp[1]); fputc('\n', pdfp[1]);
		fputformat(DEPTH16FMT, pdfp[1]);
		fputc('\n', pdfp[1]);		// end-of-info
		if (fflush(pdfp[1]) == EOF) {
			sprintf(errmsg, "cannot write header to '%s'", dfname);
			error(SYSTEM, errmsg);
			fclose(pdfp[0]); fclose(pdfp[1]);
			pdfp[0] = pdfp[1] = NULL;
			return RDTnone;
		}
	}
	return dt;				// ready to roll
}

/*
 * Render and write a frame to the named file
 * Include any header lines set prior to call
 * Picture file must not exist
 * Write pixels to stdout if !pfname
 * Write depth to a command if dfname[0]=='!'
 */
RenderDataType
RpictSimulManager::RenderFrame(const char *pfname, RenderDataType dt, const char *dfname)
{
	FILE	*pdfp[2];
						// prepare output file(s)
	dt = NewOutput(pdfp, pfname, dt, dfname);
	if (dt == RDTnone)
		return RDTnone;
						// add resolution string(s)
	fprtresolu(GetWidth(), GetHeight(), pdfp[0]);
	if (RDTdepthT(dt) == RDTdshort)
		fprtresolu(GetWidth(), GetHeight(), pdfp[1]);

	const int	bheight = (psample > 1) ? int(8*psample+.99) : 16;
	const int	vstep = bheight >> (psample > 1);

	NewBar(bheight);			// render frame if we can
	if (!RenderBelow(GetHeight(), vstep, pdfp[0], dt, pdfp[1])) {
		fclose(pdfp[0]);
		if (pdfp[1]) (dfname[0] == '!') ? pclose(pdfp[1]) : fclose(pdfp[1]);
		return RDTnone;
	}
	NewBar();				// clean up and return
	if (pdfp[0] != stdout)
		fclose(pdfp[0]);
	if (pdfp[1]) {
		if (dfname[0] == '!') {
			int	status = pclose(pdfp[1]);
			if (status) {
				sprintf(errmsg, "depth output (%s) error status: %d",
						dfname, status);
				error(USER, errmsg);
				return RDTnone;
			}
		} else
			fclose(pdfp[1]);
	}
	return dt;
}

// passed struct for header line callback
static struct HeaderInfo {
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
}	hinfo;		// XXX single copy to hold custom primitives

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

// Reopen output file(s), leaving pointers at end of (each) header
RenderDataType
RpictSimulManager::ReopenOutput(FILE *pdfp[2], const char *pfname, const char *dfname)
{
	extern const char	HDRSTR[];

	if (!pfname || pfname[0] == '!') {
		pdfp[0] = pdfp[1] = NULL;
		return RDTnone;
	}
	RenderDataType	dt = RDTnone;
	pdfp[1] = NULL;
	pdfp[0] = fopen(pfname, "r+");
	if (!pdfp[0]) {
		sprintf(errmsg, "cannot reopen output picture '%s'", pfname);
		error(SYSTEM, errmsg);
		return RDTnone;
	}
	SET_FILE_BINARY(pdfp[0]);	// read header information
	if (getheader(pdfp[0], head_check, &hinfo) < 0) {
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	if (!hinfo.gotview) {
		sprintf(errmsg, "missing view for '%s'", pfname);
		error(USER, errmsg);
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	if (hinfo.ncomp < 3) {
		sprintf(errmsg, "bad # components (%d) in '%s'", hinfo.ncomp, pfname);
		error(USER, errmsg);
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
					// set rendering/output space
	if (!strcmp(hinfo.fmt, COLRFMT)) {
		prims = hinfo.prims;	// XXX static array
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
			fclose(pdfp[0]);
			pdfp[0] = NULL;
			return RDTnone;
		}
		NCSAMP = hinfo.ncomp;		// overrides global setting
		prims = NULL;
		dt = RDTnewCT(dt, RDTscolr);
	} else if (!strcmp(hinfo.fmt, "float")) {
		if (!hinfo.endianMatch) {
			sprintf(errmsg, "incompatible byte ordering in '%s'", pfname);
			error(USER, errmsg);
			fclose(pdfp[0]);
			pdfp[0] = NULL;
			return RDTnone;
		}
		if (hinfo.ncomp == 3) {
			prims = hinfo.prims;	// custom primaries?
			int	n = 8*hinfo.gotprims;
			while (n--)
				if (!FABSEQ(hinfo.prims[0][n], stdprims[0][n]))
					break;
			if (n < 0)		// standard primaries?
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
	} else {
		sprintf(errmsg, "unknown format (%s) for '%s'", hinfo.fmt, pfname);
		error(USER, errmsg);
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	if (hinfo.gotview) {			// header view overrides
		pvw = vw;
		vw = hinfo.vw;
	}
	if (!dfname)				// no depth file?
		return dt;

	if (dfname[0] == '!') {
		error(USER, "depth data cannot be reloaded from command");
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	pdfp[1] = fopen(dfname, "r+");
	if (!pdfp[1]) {
		sprintf(errmsg, "cannot reopen depth file '%s'", dfname);
		error(SYSTEM, errmsg);
		fclose(pdfp[0]);
		pdfp[0] = NULL;
		return RDTnone;
	}
	SET_FILE_BINARY(pdfp[1]);
	int	n, len = strlen(HDRSTR);
	char	buf[32];			// sniff for 16-bit header
	if (getbinary(buf, 1, len+1, pdfp[1]) < len+1) {
		sprintf(errmsg, "empty depth file '%s'", dfname);
		error(SYSTEM, errmsg);
		fclose(pdfp[0]); fclose(pdfp[1]);
		pdfp[0] = pdfp[1] = NULL;
		return RDTnone;
	}
	for (n = 0; n < len; n++)
		if (buf[n] != HDRSTR[n])
			break;			// not a Radiance header
	rewind(pdfp[1]);
	if ((n < len) | !isprint(buf[len]))
		return RDTnewDT(dt, RDTdfloat);

	HeaderInfo	dinfo;			// thinking it's 16-bit encoded
	if (getheader(pdfp[1], head_check, &dinfo) < 0)
		sprintf(errmsg, "bad header in encoded depth file '%s'",
				dfname);
	else if (strcmp(dinfo.fmt, DEPTH16FMT))
		sprintf(errmsg, "wrong format (%s) for depth file '%s'",
				dinfo.fmt, dfname);
	else if (!SetReferenceDepth(dinfo.depth_unit))
		sprintf(errmsg, "bad/missing reference depth (%s) in '%s'",
				dinfo.depth_unit, dfname);
	else
		errmsg[0] = '\0';

	if (errmsg[0]) {
		error(USER, errmsg);
		fclose(pdfp[1]); fclose(pdfp[0]);
		pdfp[0] = pdfp[1] = NULL;
		return RDTnone;
	}
	return RDTnewDT(dt, RDTdshort);
}

// Resume partially finished rendering
// Picture file must exist
RenderDataType
RpictSimulManager::ResumeFrame(const char *pfname, const char *dfname)
{
	FILE		*pdfp[2];

	RenderDataType	dt = ReopenOutput(pdfp, pfname, dfname);
	if (dt == RDTnone)
		return RDTnone;

	int	bytesPer = 0;		// figure out how far we got...
	switch (RDTcolorT(dt)) {
	case RDTrgbe:
	case RDTxyze:
		break;
	case RDTscolr:
		bytesPer = NCSAMP + 1;	// XXX assumes no compression
		break;
	case RDTrgb:
	case RDTxyz:
		bytesPer = sizeof(float)*3;
		break;
	case RDTscolor:
		bytesPer = sizeof(float)*NCSAMP;
		break;
	default:
		sprintf(errmsg, "unknown format for '%s'", pfname);
		error(USER, errmsg);
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return RDTnone;
	}
	RESOLU	res;
	if (!fgetsresolu(&res, pdfp[0]) || res.rt != PIXSTANDARD) {
		sprintf(errmsg, "missing/bad resolution for '%s'", pfname);
		error(USER, errmsg);
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return RDTnone;
	}
	frameNo = 0;				// set up unreferenced frame
	int	hvdim[2] = {res.xr, res.yr};
	double	noAdj = 0;
	if (!NewFrame(vw, hvdim, &noAdj) ||
			(hvdim[0] != res.xr) | (hvdim[1] != res.yr)) {
		error(CONSISTENCY, "unexpected resolution change in ResumeFrame()");
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return RDTnone;
	}
	long	dataStart = ftell(pdfp[0]);	// picture starting point
	if (dataStart < 0) {
		sprintf(errmsg, "cannot seek on '%s'", pfname);
		error(SYSTEM, errmsg);
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return RDTnone;
	}
	long	doneScans = 0;
	if (bytesPer) {				// fixed-width records?
		fseek(pdfp[0], 0, SEEK_END);
		long	dataEnd = ftell(pdfp[0]);
		doneScans = (dataEnd - dataStart)/(bytesPer*GetWidth());
		if (dataEnd-dataStart > bytesPer*GetWidth()*doneScans)
			fseek(pdfp[0], dataStart + bytesPer*GetWidth()*doneScans, SEEK_SET);
	} else {				// else get compressed scanlines
		COLR *	scan = (COLR *)tempbuffer(sizeof(COLR)*GetWidth());
		while (freadcolrs(scan, GetWidth(), pdfp[0]) >= 0)
			++doneScans;
		if (!feof(pdfp[0])) {
			sprintf(errmsg, "error reading compressed scanline from '%s'", pfname);
			error(USER, errmsg);
			fclose(pdfp[0]);
			if (pdfp[1]) fclose(pdfp[1]);
			return RDTnone;
		}
	}
	if (doneScans >= GetHeight()) {		// nothing left to do?
		sprintf(errmsg, "output file '%s' is already complete", pfname);
		error(WARNING, errmsg);
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return dt;
	}
	if (!doneScans) {
		sprintf(errmsg, "restarting empty frame '%s'", pfname);
		error(WARNING, errmsg);
	}
	long	toSkip = 0;
	switch (RDTdepthT(dt)) {		// append depth file, too?
	case RDTdfloat:
		toSkip = sizeof(float)*GetWidth()*doneScans;
		break;
	case RDTdshort:
		if (!fgetsresolu(&res, pdfp[1]) || (res.rt != PIXSTANDARD) |
				(res.xr != GetWidth()) | (res.yr != GetHeight())) {
			sprintf(errmsg, "missing/bad resolution for '%s'", dfname);
			error(USER, errmsg);
			fclose(pdfp[0]); fclose(pdfp[0]);
			return RDTnone;
		}
		toSkip = 2L*GetWidth()*doneScans;
		break;
	default:;
	}					// fseek() needed for output
	if (pdfp[1] && fseek(pdfp[1], toSkip, SEEK_CUR) < 0) {
		sprintf(errmsg, "cannot seek on depth file '%s'", dfname);
		error(SYSTEM, errmsg);
		fclose(pdfp[0]); fclose(pdfp[1]);
		return RDTnone;
	}
	int	bheight = (psample > 1) ? int(8*psample+.99) : 16;
	if (bheight > GetHeight()-doneScans)
		bheight = GetHeight()-doneScans;
	int	vstep =  bheight >> (psample > 1);
	vstep += !vstep;

	NewBar(bheight);			// render remainder if we can
	if (!RenderBelow(GetHeight()-doneScans, vstep, pdfp[0], dt, pdfp[1])) {
		fclose(pdfp[0]);
		if (pdfp[1]) fclose(pdfp[1]);
		return RDTnone;
	}
	NewBar();				// close up and return success
	fclose(pdfp[0]);
	if (pdfp[1]) fclose(pdfp[1]);
	return dt;
}
