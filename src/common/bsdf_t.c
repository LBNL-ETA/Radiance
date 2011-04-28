#ifndef lint
static const char RCSid[] = "$Id: bsdf_t.c,v 3.11 2011/04/28 17:46:25 greg Exp $";
#endif
/*
 *  bsdf_t.c
 *  
 *  Definitions for variable-resolution BSDF trees
 *
 *  Created by Greg Ward on 2/2/11.
 *
 */

#include "rtio.h"
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "ezxml.h"
#include "bsdf.h"
#include "bsdf_t.h"
#include "hilbert.h"

/* Callback function type for SDtraverseTre() */
typedef int	SDtreCallback(float val, const double *cmin,
					double csiz, void *cptr);

					/* reference width maximum (1.0) */
static const unsigned	iwbits = sizeof(unsigned)*4;
static const unsigned	iwmax = (1<<(sizeof(unsigned)*4))-1;
					/* maximum cumulative value */
static const unsigned	cumlmax = ~0;

/* Struct used for our distribution-building callback */
typedef struct {
	int		nic;		/* number of input coordinates */
	unsigned	alen;		/* current array length */
	unsigned	nall;		/* number of allocated entries */
	unsigned	wmin;		/* minimum square size so far */
	unsigned	wmax;		/* maximum square size */
	struct outdir_s {
		unsigned	hent;		/* entering Hilbert index */
		int		wid;		/* this square size */
		float		bsdf;		/* BSDF for this square */
	}		*darr;		/* output direction array */
} SDdistScaffold;

/* Allocate a new scattering distribution node */
static SDNode *
SDnewNode(int nd, int lg)
{
	SDNode	*st;
	
	if (nd <= 0) {
		strcpy(SDerrorDetail, "Zero dimension BSDF node request");
		return NULL;
	}
	if (nd > SD_MAXDIM) {
		sprintf(SDerrorDetail, "Illegal BSDF dimension (%d > %d)",
				nd, SD_MAXDIM);
		return NULL;
	}
	if (lg < 0) {
		st = (SDNode *)malloc(sizeof(SDNode) +
				sizeof(st->u.t[0])*((1<<nd) - 1));
		if (st != NULL)
			memset(st->u.t, 0, sizeof(st->u.t[0])<<nd);
	} else
		st = (SDNode *)malloc(sizeof(SDNode) +
				sizeof(st->u.v[0])*((1 << nd*lg) - 1));
		
	if (st == NULL) {
		if (lg < 0)
			sprintf(SDerrorDetail,
				"Cannot allocate %d branch BSDF tree", 1<<nd);
		else
			sprintf(SDerrorDetail,
				"Cannot allocate %d BSDF leaves", 1 << nd*lg);
		return NULL;
	}
	st->ndim = nd;
	st->log2GR = lg;
	return st;
}

/* Free an SD tree */
static void
SDfreeTre(SDNode *st)
{
	int	i;

	if (st == NULL)
		return;
	for (i = (st->log2GR < 0) << st->ndim; i--; )
		SDfreeTre(st->u.t[i]);
	free((void *)st);
}

/* Free a variable-resolution BSDF */
static void
SDFreeBTre(void *p)
{
	SDTre	*sdt = (SDTre *)p;

	if (sdt == NULL)
		return;
	SDfreeTre(sdt->st);
	free(sdt);
}

/* Fill branch's worth of grid values from subtree */
static void
fill_grid_branch(float *dptr, const float *sptr, int nd, int shft)
{
	unsigned	n = 1 << (shft-1);

	if (!--nd) {			/* end on the line */
		memcpy(dptr, sptr, sizeof(*dptr)*n);
		return;
	}
	while (n--)			/* recurse on each slice */
		fill_grid_branch(dptr + (n << shft*nd),
				sptr + (n << (shft-1)*nd), nd, shft);
}

/* Get pointer at appropriate offset for the given branch */
static float *
grid_branch_start(SDNode *st, int n)
{
	unsigned	skipsiz = 1 << st->log2GR;
	float		*vptr = st->u.v;
	int		i;

	for (i = 0; i < st->ndim; skipsiz <<= st->log2GR)
		if (1<<i++ & n)
			vptr += skipsiz >> 1;
	return vptr;
}

/* Simplify (consolidate) a tree by flattening uniform depth regions */
static SDNode *
SDsimplifyTre(SDNode *st)
{
	int		match, n;

	if (st == NULL)			/* check for invalid tree */
		return NULL;
	if (st->log2GR >= 0)		/* grid just returns unaltered */
		return st;
	match = 1;			/* check if grids below match */
	for (n = 0; n < 1<<st->ndim; n++) {
		if ((st->u.t[n] = SDsimplifyTre(st->u.t[n])) == NULL)
			return NULL;	/* propogate error up call stack */
		match &= (st->u.t[n]->log2GR == st->u.t[0]->log2GR);
	}
	if (match && (match = st->u.t[0]->log2GR) >= 0) {
		SDNode	*stn = SDnewNode(st->ndim, match + 1);
		if (stn == NULL)	/* out of memory? */
			return st;
					/* transfer values to new grid */
		for (n = 1 << st->ndim; n--; )
			fill_grid_branch(grid_branch_start(stn, n),
					st->u.t[n]->u.v, stn->ndim, stn->log2GR);
		SDfreeTre(st);		/* free old tree */
		st = stn;		/* return new one */
	}
	return st;
}

/* Find smallest leaf in tree */
static double
SDsmallestLeaf(const SDNode *st)
{
	if (st->log2GR < 0) {		/* tree branches */
		double	lmin = 1.;
		int	n;
		for (n = 1<<st->ndim; n--; ) {
			double	lsiz = SDsmallestLeaf(st->u.t[n]);
			if (lsiz < lmin)
				lmin = lsiz;
		}
		return .5*lmin;
	}
					/* leaf grid width */
	return 1. / (double)(1 << st->log2GR);
}

/* Add up N-dimensional hypercube array values over the given box */
static double
SDiterSum(const float *va, int nd, int shft, const int *imin, const int *imax)
{
	const unsigned	skipsiz = 1 << --nd*shft;
	double		sum = .0;
	int		i;
	
	if (skipsiz == 1)
		for (i = *imin; i < *imax; i++)
			sum += va[i];
	else
		for (i = *imin; i < *imax; i++)
			sum += SDiterSum(va + i*skipsiz, nd, shft, imin+1, imax+1);
	return sum;
}

/* Average BSDF leaves over an orthotope defined by the unit hypercube */
static double
SDavgTreBox(const SDNode *st, const double *bmin, const double *bmax)
{
	int		imin[SD_MAXDIM], imax[SD_MAXDIM];
	unsigned	n;
	int		i;

	if (!st)
		return .0;
					/* check box limits */
	for (i = st->ndim; i--; ) {
		if (bmin[i] >= 1.)
			return .0;
		if (bmax[i] <= .0)
			return .0;
		if (bmin[i] >= bmax[i])
			return .0;
	}
	if (st->log2GR < 0) {		/* iterate on subtree */
		double		sum = .0, wsum = 1e-20;
		double		sbmin[SD_MAXDIM], sbmax[SD_MAXDIM], w;

		for (n = 1 << st->ndim; n--; ) {
			w = 1.;
			for (i = st->ndim; i--; ) {
				sbmin[i] = 2.*bmin[i];
				sbmax[i] = 2.*bmax[i];
				if (n & 1<<i) {
					sbmin[i] -= 1.;
					sbmax[i] -= 1.;
				}
				if (sbmin[i] < .0) sbmin[i] = .0;
				if (sbmax[i] > 1.) sbmax[i] = 1.;
				w *= sbmax[i] - sbmin[i];
			}
			if (w > 1e-10) {
				sum += w * SDavgTreBox(st->u.t[n], sbmin, sbmax);
				wsum += w;
			}
		}
		return sum / wsum;
	}
	n = 1;				/* iterate over leaves */
	for (i = st->ndim; i--; ) {
		imin[i] = (bmin[i] <= 0) ? 0 
				: (int)((1 << st->log2GR)*bmin[i]);
		imax[i] = (bmax[i] >= 1.) ? (1 << st->log2GR)
				: (int)((1 << st->log2GR)*bmax[i] + .999999);
		n *= imax[i] - imin[i];
	}
	if (!n)
		return .0;
	
	return SDiterSum(st->u.v, st->ndim, st->log2GR, imin, imax) / (double)n;
}

/* Recursive call for SDtraverseTre() */
static int
SDdotravTre(const SDNode *st, const double *pos, int cmask,
				SDtreCallback *cf, void *cptr,
				const double *cmin, double csiz)
{
	int	rv, rval = 0;
	double	bmin[SD_MAXDIM];
	int	i, n;
					/* in branches? */
	if (st->log2GR < 0) {
		unsigned	skipmask = 0;

		csiz *= .5;
		for (i = st->ndim; i--; )
			if (1<<i & cmask)
				if (pos[i] < cmin[i] + csiz)
					for (n = 1 << st->ndim; n--; )
						if (n & 1<<i)
							skipmask |= 1<<n;
				else
					for (n = 1 << st->ndim; n--; )
						if (!(n & 1<<i))
							skipmask |= 1<<n;
		for (n = 1 << st->ndim; n--; ) {
			if (1<<n & skipmask)
				continue;
			for (i = st->ndim; i--; )
				if (1<<i & n)
					bmin[i] = cmin[i] + csiz;
				else
					bmin[i] = cmin[i];

			rval += rv = SDdotravTre(st->u.t[n], pos, cmask,
							cf, cptr, bmin, csiz);
			if (rv < 0)
				return rv;
		}
	} else {			/* else traverse leaves */
		int	clim[SD_MAXDIM][2];
		int	cpos[SD_MAXDIM];

		if (st->log2GR == 0)	/* short cut */
			return (*cf)(st->u.v[0], cmin, csiz, cptr);

		csiz /= (double)(1 << st->log2GR);
					/* assign coord. ranges */
		for (i = st->ndim; i--; )
			if (1<<i & cmask) {
				clim[i][0] = (pos[i] - cmin[i])/csiz;
					/* check overflow from f.p. error */
				clim[i][0] -= clim[i][0] >> st->log2GR;
				clim[i][1] = clim[i][0] + 1;
			} else {
				clim[i][0] = 0;
				clim[i][1] = 1 << st->log2GR;
			}
					/* fill in unused dimensions */
		for (i = SD_MAXDIM; i-- > st->ndim; ) {
			clim[i][0] = 0; clim[i][1] = 1;
		}
#if (SD_MAXDIM == 4)
		bmin[0] = cmin[0] + csiz*clim[0][0];
		for (cpos[0] = clim[0][0]; cpos[0] < clim[0][1]; cpos[0]++) {
		    bmin[1] = cmin[1] + csiz*clim[1][0];
		    for (cpos[1] = clim[1][0]; cpos[1] < clim[1][1]; cpos[1]++) {
			bmin[2] = cmin[2] + csiz*clim[2][0];
			for (cpos[2] = clim[2][0]; cpos[2] < clim[2][1]; cpos[2]++) {
			    bmin[3] = cmin[3] + csiz*(cpos[3] = clim[3][0]);
			    n = cpos[0];
			    for (i = 1; i < st->ndim; i++)
				n = (n << st->log2GR) + cpos[i];
			    for ( ; cpos[3] < clim[3][1]; cpos[3]++) {
				rval += rv = (*cf)(st->u.v[n++], bmin, csiz, cptr);
				if (rv < 0)
				    return rv;
				bmin[3] += csiz;
			    }
			    bmin[2] += csiz;
			}
			bmin[1] += csiz;
		    }
		    bmin[0] += csiz;
		}
#else
	_!_ "broken code segment!"
#endif
	}
	return rval;
}

/* Traverse a tree, visiting nodes in a slice that fits partial position */
static int
SDtraverseTre(const SDNode *st, const double *pos, int cmask,
				SDtreCallback *cf, void *cptr)
{
	static double	czero[SD_MAXDIM];
	int		i;
					/* check arguments */
	if ((st == NULL) | (cf == NULL))
		return -1;
	for (i = st->ndim; i--; )
		if (1<<i & cmask && (pos[i] < 0) | (pos[i] >= 1.))
			return -1;

	return SDdotravTre(st, pos, cmask, cf, cptr, czero, 1.);
}

/* Look up tree value at the given grid position */
static float
SDlookupTre(const SDNode *st, const double *pos, double *hcube)
{
	double	spos[SD_MAXDIM];
	int	i, n, t;
					/* initialize voxel return */
	if (hcube) {
		hcube[i = st->ndim] = 1.;
		while (i--)
			hcube[i] = .0;
	}
					/* climb the tree */
	while (st->log2GR < 0) {
		n = 0;			/* move to appropriate branch */
		if (hcube) hcube[st->ndim] *= .5;
		for (i = st->ndim; i--; ) {
			spos[i] = 2.*pos[i];
			t = (spos[i] >= 1.);
			n |= t<<i;
			spos[i] -= (double)t;
			if (hcube) hcube[i] += (double)t * hcube[st->ndim];
		}
		st = st->u.t[n];	/* avoids tail recursion */
		pos = spos;
	}
	if (st->log2GR == 0)		/* short cut */
		return st->u.v[0];
	n = t = 0;			/* find grid array index */
	for (i = st->ndim; i--; ) {
		n += (int)((1<<st->log2GR)*pos[i]) << t;
		t += st->log2GR;
	}
	if (hcube) {			/* compute final hypercube */
		hcube[st->ndim] /= (double)(1<<st->log2GR);
		for (i = st->ndim; i--; )
			hcube[i] += floor((1<<st->log2GR)*pos[i])*hcube[st->ndim];
	}
	return st->u.v[n];		/* no interpolation */
}

/* Query BSDF value and sample hypercube for the given vectors */
static float
SDqueryTre(const SDTre *sdt, const FVECT outVec, const FVECT inVec, double *hc)
{
	static const FVECT	zvec = {.0, .0, 1.};
	FVECT			rOutVec;
	double			gridPos[4];

	switch (sdt->sidef) {		/* whose side are you on? */
	case SD_UFRONT:
		if ((outVec[2] < 0) | (inVec[2] < 0))
			return -1.;
		break;
	case SD_UBACK:
		if ((outVec[2] > 0) | (inVec[2] > 0))
			return -1.;
		break;
	case SD_XMIT:
		if ((outVec[2] > 0) == (inVec[2] > 0))
			return -1.;
		break;
	default:
		return -1.;
	}
					/* convert vector coordinates */
	if (sdt->st->ndim == 3) {
		spinvector(rOutVec, outVec, zvec, -atan2(inVec[1],inVec[0]));
		gridPos[0] = .5 - .5*sqrt(inVec[0]*inVec[0] + inVec[1]*inVec[1]);
		SDdisk2square(gridPos+1, rOutVec[0], rOutVec[1]);
	} else if (sdt->st->ndim == 4) {
		SDdisk2square(gridPos, -inVec[0], -inVec[1]);
		SDdisk2square(gridPos+2, outVec[0], outVec[1]);
	} else
		return -1.;		/* should be internal error */

	return SDlookupTre(sdt->st, gridPos, hc);
}

/* Compute non-diffuse component for variable-resolution BSDF */
static int
SDgetTreBSDF(float coef[SDmaxCh], const FVECT outVec,
				const FVECT inVec, SDComponent *sdc)
{
					/* check arguments */
	if ((coef == NULL) | (outVec == NULL) | (inVec == NULL) | (sdc == NULL)
				|| sdc->dist == NULL)
		return 0;
					/* get nearest BSDF value */
	coef[0] = SDqueryTre((SDTre *)sdc->dist, outVec, inVec, NULL);
	return (coef[0] >= 0);		/* monochromatic for now */
}

/* Callback to build cumulative distribution using SDtraverseTre() */
static int
build_scaffold(float val, const double *cmin, double csiz, void *cptr)
{
	SDdistScaffold	*sp = (SDdistScaffold *)cptr;
	int		wid = csiz*(double)iwmax + .5;
	bitmask_t	bmin[2], bmax[2];

	cmin += sp->nic;		/* skip to output coords */
	if (wid < sp->wmin)		/* new minimum width? */
		sp->wmin = wid;
	if (wid > sp->wmax)		/* new maximum? */
		sp->wmax = wid;
	if (sp->alen >= sp->nall) {	/* need more space? */
		struct outdir_s	*ndarr;
		sp->nall += 8192;
		ndarr = (struct outdir_s *)realloc(sp->darr,
					sizeof(struct outdir_s)*sp->nall);
		if (ndarr == NULL)
			return -1;	/* abort build */
		sp->darr = ndarr;
	}
					/* find Hilbert entry index */
	bmin[0] = cmin[0]*(double)iwmax + .5;
	bmin[1] = cmin[1]*(double)iwmax + .5;
	bmax[0] = bmin[0] + wid-1;
	bmax[1] = bmin[1] + wid-1;
	hilbert_box_vtx(2, sizeof(bitmask_t), iwbits, 1, bmin, bmax);
	sp->darr[sp->alen].hent = hilbert_c2i(2, iwbits, bmin);
	sp->darr[sp->alen].wid = wid;
	sp->darr[sp->alen].bsdf = val;
	sp->alen++;			/* on to the next entry */
	return 0;
}

/* Scaffold comparison function for qsort -- ascending Hilbert index */
static int
sscmp(const void *p1, const void *p2)
{
	unsigned	h1 = (*(const struct outdir_s *)p1).hent;
	unsigned	h2 = (*(const struct outdir_s *)p2).hent;

	if (h1 > h2)
		return 1;
	if (h1 < h2)
		return -1;
	return 0;
}

/* Create a new cumulative distribution for the given input direction */
static SDTreCDst *
make_cdist(const SDTre *sdt, const double *pos)
{
	SDdistScaffold	myScaffold;
	SDTreCDst	*cd;
	struct outdir_s	*sp;
	double		scale, cursum;
	int		i;
					/* initialize scaffold */
	myScaffold.wmin = iwmax;
	myScaffold.wmax = 0;
	myScaffold.nic = sdt->st->ndim - 2;
	myScaffold.alen = 0;
	myScaffold.nall = 8192;
	myScaffold.darr = (struct outdir_s *)malloc(sizeof(struct outdir_s) *
							myScaffold.nall);
	if (myScaffold.darr == NULL)
		return NULL;
					/* grow the distribution */
	if (SDtraverseTre(sdt->st, pos, (1<<myScaffold.nic)-1,
				&build_scaffold, &myScaffold) < 0) {
		free(myScaffold.darr);
		return NULL;
	}
					/* allocate result holder */
	cd = (SDTreCDst *)malloc(sizeof(SDTreCDst) +
				sizeof(cd->carr[0])*myScaffold.alen);
	if (cd == NULL) {
		free(myScaffold.darr);
		return NULL;
	}
					/* sort the distribution */
	qsort(myScaffold.darr, cd->calen = myScaffold.alen,
				sizeof(struct outdir_s), &sscmp);

					/* record input range */
	scale = myScaffold.wmin / (double)iwmax;
	for (i = myScaffold.nic; i--; ) {
		cd->clim[i][0] = floor(pos[i]/scale) * scale;
		cd->clim[i][1] = cd->clim[i][0] + scale;
	}
	cd->max_psa = myScaffold.wmax / (double)iwmax;
	cd->max_psa *= cd->max_psa * M_PI;
	cd->sidef = sdt->sidef;
	cd->cTotal = 1e-20;		/* compute directional total */
	sp = myScaffold.darr;
	for (i = myScaffold.alen; i--; sp++)
		cd->cTotal += sp->bsdf * (double)sp->wid * sp->wid;
	cursum = .0;			/* go back and get cumulative values */
	scale = (double)cumlmax / cd->cTotal;
	sp = myScaffold.darr;
	for (i = 0; i < cd->calen; i++, sp++) {
		cd->carr[i].hndx = sp->hent;
		cd->carr[i].cuml = scale*cursum + .5;
		cursum += sp->bsdf * (double)sp->wid * sp->wid;
	}
	cd->carr[i].hndx = ~0;		/* make final entry */
	cd->carr[i].cuml = cumlmax;
	cd->cTotal *= M_PI/(double)iwmax/iwmax;
					/* all done, clean up and return */
	free(myScaffold.darr);
	return cd;
}

/* Find or allocate a cumulative distribution for the given incoming vector */
const SDCDst *
SDgetTreCDist(const FVECT inVec, SDComponent *sdc)
{
	const SDTre	*sdt;
	double		inCoord[2];
	int		vflags;
	int		i;
	SDTreCDst	*cd, *cdlast;
					/* check arguments */
	if ((inVec == NULL) | (sdc == NULL) ||
			(sdt = (SDTre *)sdc->dist) == NULL)
		return NULL;
	if (sdt->st->ndim == 3)		/* isotropic BSDF? */
		inCoord[0] = .5 - .5*sqrt(inVec[0]*inVec[0] + inVec[1]*inVec[1]);
	else if (sdt->st->ndim == 4)
		SDdisk2square(inCoord, -inVec[0], -inVec[1]);
	else
		return NULL;		/* should be internal error */
	cdlast = NULL;			/* check for direction in cache list */
	for (cd = (SDTreCDst *)sdc->cdList; cd != NULL;
				cdlast = cd, cd = (SDTreCDst *)cd->next) {
		for (i = sdt->st->ndim - 2; i--; )
			if ((cd->clim[i][0] > inCoord[i]) |
					(inCoord[i] >= cd->clim[i][1]))
				break;
		if (i < 0)
			break;		/* means we have a match */
	}
	if (cd == NULL)			/* need to create new entry? */
		cdlast = cd = make_cdist(sdt, inCoord);
	if (cdlast != NULL) {		/* move entry to head of cache list */
		cdlast->next = cd->next;
		cd->next = sdc->cdList;
		sdc->cdList = (SDCDst *)cd;
	}
	return (SDCDst *)cd;		/* ready to go */
}

/* Query solid angle for vector(s) */
static SDError
SDqueryTreProjSA(double *psa, const FVECT v1, const RREAL *v2,
					int qflags, SDComponent *sdc)
{
	double		myPSA[2];
					/* check arguments */
	if ((psa == NULL) | (v1 == NULL) | (sdc == NULL) ||
				sdc->dist == NULL)
		return SDEargument;
					/* get projected solid angle(s) */
	if (v2 != NULL) {
		const SDTre	*sdt = (SDTre *)sdc->dist;
		double		hcube[SD_MAXDIM];
		if (SDqueryTre(sdt, v1, v2, hcube) < 0) {
			strcpy(SDerrorDetail, "Bad call to SDqueryTreProjSA");
			return SDEinternal;
		}
		myPSA[0] = hcube[sdt->st->ndim];
		myPSA[1] = myPSA[0] *= myPSA[0] * M_PI;
	} else {
		const SDTreCDst	*cd = (const SDTreCDst *)SDgetTreCDist(v1, sdc);
		if (cd == NULL)
			return SDEmemory;
		myPSA[0] = M_PI * (cd->clim[0][1] - cd->clim[0][0]) *
				(cd->clim[1][1] - cd->clim[1][0]);
		myPSA[1] = cd->max_psa;
	}
	switch (qflags) {		/* record based on flag settings */
	case SDqueryVal:
		*psa = myPSA[0];
		break;
	case SDqueryMax:
		if (myPSA[1] > *psa)
			*psa = myPSA[1];
		break;
	case SDqueryMin+SDqueryMax:
		if (myPSA[1] > psa[1])
			psa[1] = myPSA[1];
		/* fall through */
	case SDqueryMin:
		if (myPSA[0] < psa[0])
			psa[0] = myPSA[0];
		break;
	}
	return SDEnone;
}

/* Sample cumulative distribution */
static SDError
SDsampTreCDist(FVECT ioVec, double randX, const SDCDst *cdp)
{
	const unsigned	nBitsC = 4*sizeof(bitmask_t);
	const unsigned	nExtraBits = 8*(sizeof(bitmask_t)-sizeof(unsigned));
	const SDTreCDst	*cd = (const SDTreCDst *)cdp;
	const unsigned	target = randX*cumlmax;
	bitmask_t	hndx, hcoord[2];
	double		gpos[3];
	int		i, iupper, ilower;
					/* check arguments */
	if ((ioVec == NULL) | (cd == NULL))
		return SDEargument;
	if (ioVec[2] > 0) {
		if (!(cd->sidef & SD_UFRONT))
			return SDEargument;
	} else if (!(cd->sidef & SD_UBACK))
		return SDEargument;
					/* binary search to find position */
	ilower = 0; iupper = cd->calen;
	while ((i = (iupper + ilower) >> 1) != ilower)
		if ((long)target >= (long)cd->carr[i].cuml)
			ilower = i;
		else
			iupper = i;
					/* localize random position */
	randX = (randX*cumlmax - cd->carr[ilower].cuml) /
		    (double)(cd->carr[iupper].cuml - cd->carr[ilower].cuml);
					/* index in longer Hilbert curve */
	hndx = (randX*cd->carr[iupper].hndx + (1.-randX)*cd->carr[ilower].hndx)
				* (double)((bitmask_t)1 << nExtraBits);
					/* convert Hilbert index to vector */
	hilbert_i2c(2, nBitsC, hndx, hcoord);
	for (i = 2; i--; )
		gpos[i] = ((double)hcoord[i] + rand()*(1./(RAND_MAX+.5))) /
				(double)((bitmask_t)1 << nBitsC);
	SDsquare2disk(gpos, gpos[0], gpos[1]);
					/* compute Z-coordinate */
	gpos[2] = 1. - gpos[0]*gpos[0] - gpos[1]*gpos[1];
	if (gpos[2] > 0)		/* paranoia, I hope */
		gpos[2] = sqrt(gpos[2]);
					/* emit from back? */
	if (ioVec[2] > 0 ^ cd->sidef != SD_XMIT)
		gpos[2] = -gpos[2];
	VCOPY(ioVec, gpos);
	return SDEnone;
}

/* Advance pointer to the next non-white character in the string (or nul) */
static int
next_token(char **spp)
{
	while (isspace(**spp))
		++*spp;
	return **spp;
}

#define eat_token(spp,c)	(next_token(spp)==(c) ? *(*(spp))++ : 0)

/* Count words from this point in string to '}' */
static int
count_values(char *cp)
{
	int	n = 0;

	while (next_token(&cp) != '}' && *cp) {
		while (!isspace(*cp) & (*cp != ',') & (*cp != '}'))
			if (!*++cp)
				break;
		++n;
		eat_token(&cp, ',');
	}
	return n;
}

/* Load an array of real numbers, returning total */
static int
load_values(char **spp, float *va, int n)
{
	float	*v = va;
	char	*svnext;

	while (n-- > 0 && (svnext = fskip(*spp)) != NULL) {
		*v++ = atof(*spp);
		*spp = svnext;
		eat_token(spp, ',');
	}
	return v - va;
}

/* Load BSDF tree data */
static SDNode *
load_tree_data(char **spp, int nd)
{
	SDNode	*st;
	int	n;

	if (!eat_token(spp, '{')) {
		strcpy(SDerrorDetail, "Missing '{' in tensor tree");
		return NULL;
	}
	if (next_token(spp) == '{') {	/* tree branches */
		st = SDnewNode(nd, -1);
		if (st == NULL)
			return NULL;
		for (n = 0; n < 1<<nd; n++)
			if ((st->u.t[n] = load_tree_data(spp, nd)) == NULL) {
				SDfreeTre(st);
				return NULL;
			}
	} else {			/* else load value grid */
		int	bsiz;
		n = count_values(*spp);	/* see how big the grid is */
		for (bsiz = 0; bsiz < 8*sizeof(size_t)-1; bsiz += nd)
			if (1<<bsiz == n)
				break;
		if (bsiz >= 8*sizeof(size_t)) {
			strcpy(SDerrorDetail, "Illegal value count in tensor tree");
			return NULL;
		}
		st = SDnewNode(nd, bsiz/nd);
		if (st == NULL)
			return NULL;
		if (load_values(spp, st->u.v, n) != n) {
			strcpy(SDerrorDetail, "Real format error in tensor tree");
			SDfreeTre(st);
			return NULL;
		}
	}
	if (!eat_token(spp, '}')) {
		strcpy(SDerrorDetail, "Missing '}' in tensor tree");
		SDfreeTre(st);
		return NULL;
	}
	eat_token(spp, ',');
	return st;
}

/* Compute min. proj. solid angle and max. direct hemispherical scattering */
static SDError
get_extrema(SDSpectralDF *df)
{
	SDNode	*st = (*(SDTre *)df->comp[0].dist).st;
	double	stepWidth, dhemi, bmin[4], bmax[4];

	stepWidth = SDsmallestLeaf(st);
	df->minProjSA = M_PI*stepWidth*stepWidth;
	if (stepWidth < .03125)
		stepWidth = .03125;	/* 1/32 resolution good enough */
	df->maxHemi = .0;
	if (st->ndim == 3) {		/* isotropic BSDF */
		bmin[1] = bmin[2] = .0;
		bmax[1] = bmax[2] = 1.;
		for (bmin[0] = .0; bmin[0] < .5-FTINY; bmin[0] += stepWidth) {
			bmax[0] = bmin[0] + stepWidth;
			dhemi = SDavgTreBox(st, bmin, bmax);
			if (dhemi > df->maxHemi)
				df->maxHemi = dhemi;
		}
	} else if (st->ndim == 4) {	/* anisotropic BSDF */
		bmin[2] = bmin[3] = .0;
		bmax[2] = bmax[3] = 1.;
		for (bmin[0] = .0; bmin[0] < 1.-FTINY; bmin[0] += stepWidth) {
			bmax[0] = bmin[0] + stepWidth;
			for (bmin[1] = .0; bmin[1] < 1.-FTINY; bmin[1] += stepWidth) {
				bmax[1] = bmin[1] + stepWidth;
				dhemi = SDavgTreBox(st, bmin, bmax);
				if (dhemi > df->maxHemi)
					df->maxHemi = dhemi;
			}
		}
	} else
		return SDEinternal;
					/* correct hemispherical value */
	df->maxHemi *= M_PI;
	return SDEnone;
}

/* Load BSDF distribution for this wavelength */
static SDError
load_bsdf_data(SDData *sd, ezxml_t wdb, int ndim)
{
	SDSpectralDF	*df;
	SDTre		*sdt;
	char		*sdata;
	int		i;
					/* allocate BSDF component */
	sdata = ezxml_txt(ezxml_child(wdb, "WavelengthDataDirection"));
	if (!sdata)
		return SDEnone;
	/*
	 * Remember that front and back are reversed from WINDOW 6 orientations
	 */
	if (!strcasecmp(sdata, "Transmission")) {
		if (sd->tf != NULL)
			SDfreeSpectralDF(sd->tf);
		if ((sd->tf = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->tf;
	} else if (!strcasecmp(sdata, "Reflection Front")) {
		if (sd->rb != NULL)	/* note back-front reversal */
			SDfreeSpectralDF(sd->rb);
		if ((sd->rb = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->rb;
	} else if (!strcasecmp(sdata, "Reflection Back")) {
		if (sd->rf != NULL)	/* note front-back reversal */
			SDfreeSpectralDF(sd->rf);
		if ((sd->rf = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->rf;
	} else
		return SDEnone;
	/* XXX should also check "ScatteringDataType" for consistency? */
					/* get angle bases */
	sdata = ezxml_txt(ezxml_child(wdb,"AngleBasis"));
	if (!sdata || strcasecmp(sdata, "LBNL/Shirley-Chiu")) {
		sprintf(SDerrorDetail, "%s angle basis for BSDF '%s'",
				!sdata ? "Missing" : "Unsupported", sd->name);
		return !sdata ? SDEformat : SDEsupport;
	}
					/* allocate BSDF tree */
	sdt = (SDTre *)malloc(sizeof(SDTre));
	if (sdt == NULL)
		return SDEmemory;
	if (df == sd->rf)
		sdt->sidef = SD_UFRONT;
	else if (df == sd->rb)
		sdt->sidef = SD_UBACK;
	else
		sdt->sidef = SD_XMIT;
	sdt->st = NULL;
	df->comp[0].cspec[0] = c_dfcolor; /* XXX monochrome for now */
	df->comp[0].dist = sdt;
	df->comp[0].func = &SDhandleTre;
					/* read BSDF data */
	sdata = ezxml_txt(ezxml_child(wdb, "ScatteringData"));
	if (!sdata || !next_token(&sdata)) {
		sprintf(SDerrorDetail, "Missing BSDF ScatteringData in '%s'",
				sd->name);
		return SDEformat;
	}
	sdt->st = load_tree_data(&sdata, ndim);
	if (sdt->st == NULL)
		return SDEformat;
	if (next_token(&sdata)) {	/* check for unconsumed characters */
		sprintf(SDerrorDetail,
			"Extra characters at end of ScatteringData in '%s'",
				sd->name);
		return SDEformat;
	}
					/* flatten branches where possible */
	sdt->st = SDsimplifyTre(sdt->st);
	if (sdt->st == NULL)
		return SDEinternal;
	return get_extrema(df);		/* compute global quantities */
}

/* Find minimum value in tree */
static float
SDgetTreMin(const SDNode *st)
{
	float	vmin = FHUGE;
	int	n;

	if (st->log2GR < 0) {
		for (n = 1<<st->ndim; n--; ) {
			float	v = SDgetTreMin(st->u.t[n]);
			if (v < vmin)
				vmin = v;
		}
	} else {
		for (n = 1<<(st->ndim*st->log2GR); n--; )
			if (st->u.v[n] < vmin)
				vmin = st->u.v[n];
	}
	return vmin;
}

/* Subtract the given value from all tree nodes */
static void
SDsubtractTreVal(SDNode *st, float val)
{
	int	n;

	if (st->log2GR < 0) {
		for (n = 1<<st->ndim; n--; )
			SDsubtractTreVal(st->u.t[n], val);
	} else {
		for (n = 1<<(st->ndim*st->log2GR); n--; )
			st->u.v[n] -= val;
	}
}

/* Subtract minimum value from BSDF */
static double
subtract_min(SDNode *st)
{
	float	vmin;
					/* be sure to skip unused portion */
	if (st->ndim == 3) {
		int	n;
		vmin = 1./M_PI;
		if (st->log2GR < 0) {
			for (n = 0; n < 4; n++) {
				float	v = SDgetTreMin(st->u.t[n]);
				if (v < vmin)
					vmin = v;
			}
		} else if (st->log2GR) {
			for (n = 1 << (3*st->log2GR - 1); n--; )
				if (st->u.v[n] < vmin)
					vmin = st->u.v[n];
		} else
			vmin = st->u.v[0];
	} else				/* anisotropic covers entire tree */
		vmin = SDgetTreMin(st);

	if (vmin <= FTINY)
		return .0;

	SDsubtractTreVal(st, vmin);

	return M_PI * vmin;		/* return hemispherical value */
}

/* Extract and separate diffuse portion of BSDF */
static void
extract_diffuse(SDValue *dv, SDSpectralDF *df)
{
	int	n;

	if (df == NULL || df->ncomp <= 0) {
		dv->spec = c_dfcolor;
		dv->cieY = .0;
		return;
	}
	dv->spec = df->comp[0].cspec[0];
	dv->cieY = subtract_min((*(SDTre *)df->comp[0].dist).st);
					/* in case of multiple components */
	for (n = df->ncomp; --n; ) {
		double	ymin = subtract_min((*(SDTre *)df->comp[n].dist).st);
		c_cmix(&dv->spec, dv->cieY, &dv->spec, ymin, &df->comp[n].cspec[0]);
		dv->cieY += ymin;
	}
	df->maxHemi -= dv->cieY;	/* adjust maximum hemispherical */
					/* make sure everything is set */
	c_ccvt(&dv->spec, C_CSXY+C_CSSPEC);
}

/* Load a variable-resolution BSDF tree from an open XML file */
SDError
SDloadTre(SDData *sd, ezxml_t wtl)
{
	SDError		ec;
	ezxml_t		wld, wdb;
	int		rank;
	char		*txt;
					/* basic checks and tensor rank */
	txt = ezxml_txt(ezxml_child(ezxml_child(wtl,
			"DataDefinition"), "IncidentDataStructure"));
	if (txt == NULL || !*txt) {
		sprintf(SDerrorDetail,
			"BSDF \"%s\": missing IncidentDataStructure",
				sd->name);
		return SDEformat;
	}
	if (!strcasecmp(txt, "TensorTree3"))
		rank = 3;
	else if (!strcasecmp(txt, "TensorTree4"))
		rank = 4;
	else {
		sprintf(SDerrorDetail,
			"BSDF \"%s\": unsupported IncidentDataStructure",
				sd->name);
		return SDEsupport;
	}
					/* load BSDF components */
	for (wld = ezxml_child(wtl, "WavelengthData");
				wld != NULL; wld = wld->next) {
		if (strcasecmp(ezxml_txt(ezxml_child(wld,"Wavelength")),
				"Visible"))
			continue;	/* just visible for now */
		for (wdb = ezxml_child(wld, "WavelengthDataBlock");
					wdb != NULL; wdb = wdb->next)
			if ((ec = load_bsdf_data(sd, wdb, rank)) != SDEnone)
				return ec;
	}
					/* separate diffuse components */
	extract_diffuse(&sd->rLambFront, sd->rf);
	extract_diffuse(&sd->rLambBack, sd->rb);
	extract_diffuse(&sd->tLamb, sd->tf);
					/* return success */
	return SDEnone;
}

/* Variable resolution BSDF methods */
SDFunc SDhandleTre = {
	&SDgetTreBSDF,
	&SDqueryTreProjSA,
	&SDgetTreCDist,
	&SDsampTreCDist,
	&SDFreeBTre,
};
