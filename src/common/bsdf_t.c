#ifndef lint
static const char RCSid[] = "$Id: bsdf_t.c,v 3.51 2021/03/27 20:08:30 greg Exp $";
#endif
/*
 *  bsdf_t.c
 *  
 *  Definitions for variable-resolution BSDF trees
 *
 *  Created by Greg Ward on 2/2/11.
 *
 */

#define	_USE_MATH_DEFINES
#include "rtio.h"
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "ezxml.h"
#include "bsdf.h"
#include "bsdf_t.h"
#include "hilbert.h"

/* Callback function type for SDtraverseTre() */
typedef int	SDtreCallback(float val, const double *cmin, double csiz,
						void *cptr);
					/* reference width maximum (1.0) */
static const unsigned	iwbits = sizeof(unsigned)*4;
static const unsigned	iwmax = 1<<(sizeof(unsigned)*4);
					/* maximum cumulative value */
static const unsigned	cumlmax = ~0;
					/* constant z-vector */
static const FVECT	zvec = {.0, .0, 1.};
					/* quantization value */
static double		quantum = 1./256.;
					/* our RGB primaries */
static C_COLOR		tt_RGB_prim[3];
static float		tt_RGB_coef[3];

static const double	czero[SD_MAXDIM];

enum {tt_Y, tt_u, tt_v};		/* tree components (tt_Y==0) */

/* Struct used for our distribution-building callback */
typedef struct {
	short		nic;		/* number of input coordinates */
	short		rev;		/* reversing query */
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
		if (st == NULL) {
			sprintf(SDerrorDetail,
				"Cannot allocate %d branch BSDF tree", 1<<nd);
			return NULL;
		}
		memset(st->u.t, 0, sizeof(st->u.t[0])<<nd);
	} else {
		st = (SDNode *)malloc(sizeof(SDNode) +
				sizeof(st->u.v[0])*((1 << nd*lg) - 1));		
		if (st == NULL) {
			sprintf(SDerrorDetail,
				"Cannot allocate %d BSDF leaves", 1 << nd*lg);
			return NULL;
		}
	}
	st->ndim = nd;
	st->log2GR = lg;
	return st;
}

/* Free an SD tree */
static void
SDfreeTre(SDNode *st)
{
	int	n;

	if (st == NULL)
		return;
	for (n = (st->log2GR < 0) << st->ndim; n--; )
		SDfreeTre(st->u.t[n]);
	free(st);
}

/* Free a variable-resolution BSDF */
static void
SDFreeBTre(void *p)
{
	SDTre	*sdt = (SDTre *)p;

	if (sdt == NULL)
		return;
	SDfreeTre(sdt->stc[tt_Y]);
	SDfreeTre(sdt->stc[tt_u]);
	SDfreeTre(sdt->stc[tt_v]);
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
	unsigned	skipsiz = 1 << (st->log2GR - 1);
	float		*vptr = st->u.v;
	int		i;

	for (i = st->ndim; i--; skipsiz <<= st->log2GR)
		if (1<<i & n)
			vptr += skipsiz;
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

/* Assign the given voxel in tree (produces no grid nodes) */
static SDNode *
SDsetVoxel(SDNode *sroot, int nd, const double *tmin, const double tsiz, float val)
{
	double	ctrk[SD_MAXDIM];
	double	csiz = 1.;
	SDNode	*st;
	int	i, n;
					/* check arguments */
	for (i = nd; i-- > 0; )
		if ((tmin[i] < .0) | (tmin[i] >= 1.-FTINY))
			break;
	if ((i >= 0) | (nd <= 0) | (tsiz <= FTINY) | (tsiz > 1.+FTINY) |
			(sroot != NULL && sroot->ndim != nd)) {
		SDfreeTre(sroot);
		return NULL;
	}
	if (tsiz >= 1.-FTINY) {		/* special case when tree is a leaf */
		SDfreeTre(sroot);
		if ((sroot = SDnewNode(nd, 0)) != NULL)
			sroot->u.v[0] = val;
		return sroot;
	}
					/* make sure we have branching root */
	if (sroot != NULL && sroot->log2GR >= 0) {
		SDfreeTre(sroot); sroot = NULL;
	}
	if (sroot == NULL && (sroot = SDnewNode(nd, -1)) == NULL)
		return NULL;
	st = sroot;			/* climb/grow tree */
	memset(ctrk, 0, sizeof(ctrk));
	for ( ; ; ) {
		csiz *= .5;		/* find appropriate branch */
		n = 0;
		for (i = nd; i--; )
			if (ctrk[i]+csiz <= tmin[i]+FTINY) {
				ctrk[i] += csiz;
				n |= 1 << i;
			}
					/* reached desired voxel? */
		if (csiz <= tsiz+FTINY) {
			SDfreeTre(st->u.t[n]);
			st = st->u.t[n] = SDnewNode(nd, 0);
			break;
		}
					/* else grow tree as needed */
		if (st->u.t[n] != NULL && st->u.t[n]->log2GR >= 0) {
			SDfreeTre(st->u.t[n]); st->u.t[n] = NULL;
		}
		if (st->u.t[n] == NULL)
			st->u.t[n] = SDnewNode(nd, -1);
		if ((st = st->u.t[n]) == NULL)
			break;
	}
	if (st == NULL) {
		SDfreeTre(sroot);
		return NULL;
	}
	st->u.v[0] = val;		/* assign leaf and return root */
	return sroot;
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

	va += *imin * skipsiz;

	if (skipsiz == 1)
		for (i = *imin; i < *imax; i++)
			sum += *va++;
	else
		for (i = *imin; i < *imax; i++, va += skipsiz)
			sum += SDiterSum(va, nd, shft, imin+1, imax+1);
	return sum;
}

/* Average BSDF leaves over an orthotope defined by the unit hypercube */
static double
SDavgTreBox(const SDNode *st, const double *bmin, const double *bmax)
{
	unsigned	n;
	int		i;

	if (!st)
		return .0;
					/* check box limits */
	for (i = st->ndim; i--; ) {
		if (bmin[i] >= 1.)
			return .0;
		if (bmax[i] <= 0)
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
				if (sbmin[i] >= sbmax[i]) {
					w = .0;
					break;
				}
				w *= sbmax[i] - sbmin[i];
			}
			if (w > 1e-10) {
				sum += w * SDavgTreBox(st->u.t[n], sbmin, sbmax);
				wsum += w;
			}
		}
		return sum / wsum;
	} else {			/* iterate over leaves */
		int		imin[SD_MAXDIM], imax[SD_MAXDIM];

		n = 1;
		for (i = st->ndim; i--; ) {
			imin[i] = (bmin[i] <= 0) ? 0 :
					(int)((1 << st->log2GR)*bmin[i]);
			imax[i] = (bmax[i] >= 1.) ? (1 << st->log2GR) :
				(int)((1 << st->log2GR)*bmax[i] + .999999);
			n *= imax[i] - imin[i];
		}
		if (n)
			return SDiterSum(st->u.v, st->ndim,
					st->log2GR, imin, imax) / (double)n;
	}
	return .0;
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
					/* paranoia */
	if (st == NULL)
		return 0;
					/* in branches? */
	if (st->log2GR < 0) {
		unsigned	skipmask = 0;
		csiz *= .5;
		for (i = st->ndim; i--; )
			if (1<<i & cmask) {
				if (pos[i] < cmin[i] + csiz)
					for (n = 1 << st->ndim; n--; ) {
						if (n & 1<<i)
							skipmask |= 1<<n;
					}
				else
					for (n = 1 << st->ndim; n--; ) {
						if (!(n & 1<<i))
							skipmask |= 1<<n;
					}
			}
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
#if (SD_MAXDIM == 4)
		bmin[0] = cmin[0] + csiz*clim[0][0];
		for (cpos[0] = clim[0][0]; cpos[0] < clim[0][1]; cpos[0]++) {
		    bmin[1] = cmin[1] + csiz*clim[1][0];
		    for (cpos[1] = clim[1][0]; cpos[1] < clim[1][1]; cpos[1]++) {
			bmin[2] = cmin[2] + csiz*clim[2][0];
			if (st->ndim == 3) {
			    cpos[2] = clim[2][0];
			    n = cpos[0];
			    for (i = 1; i < 3; i++)
				n = (n << st->log2GR) + cpos[i];
			    for ( ; cpos[2] < clim[2][1]; cpos[2]++) {
				rval += rv = (*cf)(st->u.v[n++], bmin, csiz, cptr);
				if (rv < 0)
				    return rv;
				bmin[2] += csiz;
			    }
			} else {
			    for (cpos[2] = clim[2][0]; cpos[2] < clim[2][1]; cpos[2]++) {
				bmin[3] = cmin[3] + csiz*(cpos[3] = clim[3][0]);
				n = cpos[0];
				for (i = 1; i < 4; i++)
				    n = (n << st->log2GR) + cpos[i];
				for ( ; cpos[3] < clim[3][1]; cpos[3]++) {
				    rval += rv = (*cf)(st->u.v[n++], bmin, csiz, cptr);
				    if (rv < 0)
					return rv;
				    bmin[3] += csiz;
				}
				bmin[2] += csiz;
			    }
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
	while (st != NULL && st->log2GR < 0) {
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
	if (st == NULL)			/* should never happen? */
		return .0;
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

/* Convert CIE (Y,u',v') color to our RGB */
static void
SDyuv2rgb(double yval, double uprime, double vprime, float rgb[3])
{
	const double	dfact = 1./(6.*uprime - 16.*vprime + 12.);
	C_COLOR		cxy;

	c_cset(&cxy, 9.*uprime*dfact, 4.*vprime*dfact);
	c_toSharpRGB(&cxy, yval, rgb);
}

/* Query BSDF value and sample hypercube for the given vectors */
static int
SDqueryTre(const SDTre *sdt, float *coef,
		const FVECT outVec, const FVECT inVec, double *hc)
{
	const RREAL	*vtmp;
	float		yval;
	FVECT		rOutVec;
	double		gridPos[4];

	if (sdt->stc[tt_Y] == NULL)	/* paranoia, I hope */
		return 0;

	switch (sdt->sidef) {		/* whose side are you on? */
	case SD_FREFL:
		if ((outVec[2] < 0) | (inVec[2] < 0))
			return 0;
		break;
	case SD_BREFL:
		if ((outVec[2] > 0) | (inVec[2] > 0))
			return 0;
		break;
	case SD_FXMIT:
		if (outVec[2] > 0) {
			if (inVec[2] > 0)
				return 0;
			vtmp = outVec; outVec = inVec; inVec = vtmp;
		} else if (inVec[2] < 0)
			return 0;
		break;
	case SD_BXMIT:
		if (inVec[2] > 0) {
			if (outVec[2] > 0)
				return 0;
			vtmp = outVec; outVec = inVec; inVec = vtmp;
		} else if (outVec[2] < 0)
			return 0;
		break;
	default:
		return 0;
	}
					/* convert vector coordinates */
	if (sdt->stc[tt_Y]->ndim == 3) {
		spinvector(rOutVec, outVec, zvec, -atan2(-inVec[1],-inVec[0]));
		gridPos[0] = (.5-FTINY) -
				.5*sqrt(inVec[0]*inVec[0] + inVec[1]*inVec[1]);
		SDdisk2square(gridPos+1, rOutVec[0], rOutVec[1]);
	} else if (sdt->stc[tt_Y]->ndim == 4) {
		SDdisk2square(gridPos, -inVec[0], -inVec[1]);
		SDdisk2square(gridPos+2, outVec[0], outVec[1]);
	} else
		return 0;		/* should be internal error */
					/* get BSDF value */
	yval = SDlookupTre(sdt->stc[tt_Y], gridPos, hc);
	if (coef == NULL)		/* just getting hypercube? */
		return 1;
	if (sdt->stc[tt_u] == NULL || sdt->stc[tt_v] == NULL) {
		*coef = yval;
		return 1;		/* no color */
	}
					/* else decode color */
	SDyuv2rgb(yval, SDlookupTre(sdt->stc[tt_u], gridPos, NULL),
			SDlookupTre(sdt->stc[tt_v], gridPos, NULL), coef);
	coef[0] *= tt_RGB_coef[0];
	coef[1] *= tt_RGB_coef[1];
	coef[2] *= tt_RGB_coef[2];
	return 3;
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
	return SDqueryTre((SDTre *)sdc->dist, coef, outVec, inVec, NULL);
}

/* Callback to build cumulative distribution using SDtraverseTre() */
static int
build_scaffold(float val, const double *cmin, double csiz, void *cptr)
{
	SDdistScaffold	*sp = (SDdistScaffold *)cptr;
	int		wid = csiz*(double)iwmax + .5;
	double		revcmin[2];
	bitmask_t	bmin[2], bmax[2];

	if (sp->rev) {			/* need to reverse sense? */
		revcmin[0] = 1. - cmin[0] - csiz;
		revcmin[1] = 1. - cmin[1] - csiz;
		cmin = revcmin;
	} else {
		cmin += sp->nic;	/* else skip to output coords */
	}
	if (wid < sp->wmin)		/* new minimum width? */
		sp->wmin = wid;
	if (wid > sp->wmax)		/* new maximum? */
		sp->wmax = wid;
	if (sp->alen >= sp->nall) {	/* need more space? */
		struct outdir_s	*ndarr;
		sp->nall = (int)(1.5*sp->nall) + 256;
		ndarr = (struct outdir_s *)realloc(sp->darr,
					sizeof(struct outdir_s)*sp->nall);
		if (ndarr == NULL) {
			sprintf(SDerrorDetail,
				"Cannot grow scaffold to %u entries", sp->nall);
			return -1;	/* abort build */
		}
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
make_cdist(const SDTre *sdt, const double *invec, int rev)
{
	SDdistScaffold	myScaffold;
	double		pos[4];
	int		cmask;
	SDTreCDst	*cd;
	struct outdir_s	*sp;
	double		scale, cursum;
	int		i;
					/* initialize scaffold */
	myScaffold.wmin = iwmax;
	myScaffold.wmax = 0;
	myScaffold.nic = sdt->stc[tt_Y]->ndim - 2;
	myScaffold.rev = rev;
	myScaffold.alen = 0;
	myScaffold.nall = 512;
	myScaffold.darr = (struct outdir_s *)malloc(sizeof(struct outdir_s) *
							myScaffold.nall);
	if (myScaffold.darr == NULL)
		return NULL;
					/* set up traversal */
	cmask = (1<<myScaffold.nic) - 1;
	for (i = myScaffold.nic; i--; )
			pos[i+2*rev] = invec[i];
	cmask <<= 2*rev;
					/* grow the distribution */
	if (SDtraverseTre(sdt->stc[tt_Y], pos, cmask,
				build_scaffold, &myScaffold) < 0) {
		free(myScaffold.darr);
		return NULL;
	}
					/* allocate result holder */
	cd = (SDTreCDst *)malloc(sizeof(SDTreCDst) +
				sizeof(cd->carr[0])*myScaffold.alen);
	if (cd == NULL) {
		sprintf(SDerrorDetail,
			"Cannot allocate %u entry cumulative distribution",
				myScaffold.alen);
		free(myScaffold.darr);
		return NULL;
	}
	cd->isodist = (myScaffold.nic == 1);
					/* sort the distribution */
	qsort(myScaffold.darr, cd->calen = myScaffold.alen,
				sizeof(struct outdir_s), sscmp);

					/* record input range */
	scale = myScaffold.wmin / (double)iwmax;
	for (i = myScaffold.nic; i--; ) {
		cd->clim[i][0] = floor(pos[i+2*rev]/scale) * scale;
		cd->clim[i][1] = cd->clim[i][0] + scale;
	}
	if (cd->isodist) {		/* avoid issue in SDqueryTreProjSA() */
		cd->clim[1][0] = cd->clim[0][0];
		cd->clim[1][1] = cd->clim[0][1];
	}
	cd->max_psa = myScaffold.wmax / (double)iwmax;
	cd->max_psa *= cd->max_psa * M_PI;
	if (rev)
		cd->sidef = (sdt->sidef==SD_BXMIT) ? SD_FXMIT : SD_BXMIT;
	else
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
	unsigned long	cacheLeft = SDmaxCache;
	const SDTre	*sdt;
	double		inCoord[2];
	int		i;
	int		mode;
	SDTreCDst	*cd, *cdlast, *cdlimit;
					/* check arguments */
	if ((inVec == NULL) | (sdc == NULL) ||
			(sdt = (SDTre *)sdc->dist) == NULL)
		return NULL;
	switch (mode = sdt->sidef) {	/* check direction */
	case SD_FREFL:
		if (inVec[2] < 0)
			return NULL;
		break;
	case SD_BREFL:
		if (inVec[2] > 0)
			return NULL;
		break;
	case SD_FXMIT:
		if (inVec[2] < 0)
			mode = SD_BXMIT;
		break;
	case SD_BXMIT:
		if (inVec[2] > 0)
			mode = SD_FXMIT;
		break;
	default:
		return NULL;
	}
	if (sdt->stc[tt_Y]->ndim == 3) {	/* isotropic BSDF? */
		if (mode != sdt->sidef)	/* XXX unhandled reciprocity */
			return &SDemptyCD;
		inCoord[0] = (.5-FTINY) -
				.5*sqrt(inVec[0]*inVec[0] + inVec[1]*inVec[1]);
	} else if (sdt->stc[tt_Y]->ndim == 4) {
		if (mode != sdt->sidef)	/* use reciprocity? */
			SDdisk2square(inCoord, inVec[0], inVec[1]);
		else
			SDdisk2square(inCoord, -inVec[0], -inVec[1]);
	} else
		return NULL;		/* should be internal error */
					/* quantize to avoid f.p. errors */
	for (i = sdt->stc[tt_Y]->ndim - 2; i--; )
		inCoord[i] = floor(inCoord[i]/quantum)*quantum + .5*quantum;
	cdlast = cdlimit = NULL;	/* check for direction in cache list */
	/* PLACE MUTEX LOCK HERE FOR THREAD-SAFE */
	for (cd = (SDTreCDst *)sdc->cdList; cd != NULL;
					cdlast = cd, cd = cd->next) {
		if (cacheLeft) {	/* check cache size limit */
			long	csiz = sizeof(SDTreCDst) +
					sizeof(cd->carr[0])*cd->calen;
			if (cacheLeft > csiz)
				cacheLeft -= csiz;
			else {
				cdlimit = cdlast;
				cacheLeft = 0;
			}
		}
		if (cd->sidef != mode)
			continue;
		for (i = sdt->stc[tt_Y]->ndim - 2; i--; )
			if ((cd->clim[i][0] > inCoord[i]) |
					(inCoord[i] >= cd->clim[i][1]))
				break;
		if (i < 0)
			break;		/* means we have a match */
	}
	if (cd == NULL) {		/* need to create new entry? */
		if (cdlimit != NULL)	/* exceeded cache size limit? */
			while ((cd = cdlimit->next) != NULL) {
				cdlimit->next = cd->next;
				free(cd);
			}
		cdlast = cd = make_cdist(sdt, inCoord, mode != sdt->sidef);
	}
	if (cdlast != NULL) {		/* move entry to head of cache list */
		cdlast->next = cd->next;
		cd->next = (SDTreCDst *)sdc->cdList;
		sdc->cdList = (SDCDst *)cd;
	}
	/* END MUTEX LOCK */
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
		double		hcube[SD_MAXDIM+1];
		if (!SDqueryTre(sdt, NULL, v1, v2, hcube)) {
			strcpy(SDerrorDetail, "Bad call to SDqueryTreProjSA");
			return SDEinternal;
		}
		myPSA[0] = hcube[sdt->stc[tt_Y]->ndim];
		myPSA[1] = myPSA[0] *= myPSA[0] * M_PI;
	} else {
		const SDTreCDst	*cd = (const SDTreCDst *)SDgetTreCDist(v1, sdc);
		if (cd == NULL)
			myPSA[0] = myPSA[1] = 0;
		else {
			myPSA[0] = M_PI * (cd->clim[0][1] - cd->clim[0][0]) *
					(cd->clim[1][1] - cd->clim[1][0]);
			myPSA[1] = cd->max_psa;
		}
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
		if ((myPSA[0] > 0) & (myPSA[0] < psa[0]))
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
	double		gpos[3], rotangle;
	int		i, iupper, ilower;
					/* check arguments */
	if ((ioVec == NULL) | (cd == NULL))
		return SDEargument;
	if (!cd->sidef)
		return SDEnone;		/* XXX should never happen */
	if (ioVec[2] > 0) {
		if ((cd->sidef != SD_FREFL) & (cd->sidef != SD_FXMIT))
			return SDEargument;
	} else if ((cd->sidef != SD_BREFL) & (cd->sidef != SD_BXMIT))
		return SDEargument;
					/* binary search to find position */
	ilower = 0; iupper = cd->calen;
	while ((i = (iupper + ilower) >> 1) != ilower)
		if (target >= cd->carr[i].cuml)
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
	gpos[2] = sqrt(gpos[2]*(gpos[2]>0));
					/* emit from back? */
	if ((cd->sidef == SD_BREFL) | (cd->sidef == SD_FXMIT))
		gpos[2] = -gpos[2];
	if (cd->isodist) {		/* rotate isotropic sample */
		rotangle = atan2(-ioVec[1],-ioVec[0]);
		spinvector(ioVec, gpos, zvec, rotangle);
	} else
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

/* Advance pointer past matching token (or any token if c==0) */
#define eat_token(spp,c)	((next_token(spp)==(c)) ^ !(c) ? *(*(spp))++ : 0)

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
		if ((*v++ = atof(*spp)) < 0)
			v[-1] = 0;
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
		for (bsiz = 0; bsiz < 8*sizeof(size_t); bsiz += nd)
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
	SDNode	*st = (*(SDTre *)df->comp[0].dist).stc[tt_Y];
	double	stepWidth, dhemi, bmin[4], bmax[4];

	stepWidth = SDsmallestLeaf(st);
	if (quantum > stepWidth)	/* adjust quantization factor */
		quantum = stepWidth;
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
load_bsdf_data(SDData *sd, ezxml_t wdb, int ct, int ndim)
{
	SDSpectralDF	*df;
	SDTre		*sdt;
	char		*sdata;
					/* allocate BSDF component */
	sdata = ezxml_txt(ezxml_child(wdb, "WavelengthDataDirection"));
	if (!sdata)
		return SDEnone;
	/*
	 * Remember that front and back are reversed from WINDOW 6 orientations
	 */
	if (!strcasecmp(sdata, "Transmission Front")) {
		if (sd->tb == NULL && (sd->tb = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->tb;
	} else if (!strcasecmp(sdata, "Transmission Back")) {
		if (sd->tf == NULL && (sd->tf = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->tf;
	} else if (!strcasecmp(sdata, "Reflection Front")) {
		if (sd->rb == NULL && (sd->rb = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->rb;
	} else if (!strcasecmp(sdata, "Reflection Back")) {
		if (sd->rf == NULL && (sd->rf = SDnewSpectralDF(1)) == NULL)
			return SDEmemory;
		df = sd->rf;
	} else
		return SDEnone;
					/* get angle bases */
	sdata = ezxml_txt(ezxml_child(wdb,"AngleBasis"));
	if (!sdata || strcasecmp(sdata, "LBNL/Shirley-Chiu")) {
		sprintf(SDerrorDetail, "%s angle basis for BSDF '%s'",
				!sdata ? "Missing" : "Unsupported", sd->name);
		return !sdata ? SDEformat : SDEsupport;
	}
	if (df->comp[0].dist == NULL) {	/* need to allocate BSDF tree? */
		sdt = (SDTre *)malloc(sizeof(SDTre));
		if (sdt == NULL)
			return SDEmemory;
		if (df == sd->rf)
			sdt->sidef = SD_FREFL;
		else if (df == sd->rb)
			sdt->sidef = SD_BREFL;
		else if (df == sd->tf)
			sdt->sidef = SD_FXMIT;
		else /* df == sd->tb */
			sdt->sidef = SD_BXMIT;
		sdt->stc[tt_Y] = sdt->stc[tt_u] = sdt->stc[tt_v] = NULL;
		df->comp[0].dist = sdt;
		df->comp[0].func = &SDhandleTre;
	} else {
		sdt = (SDTre *)df->comp[0].dist;
		if (sdt->stc[ct] != NULL) {
			SDfreeTre(sdt->stc[ct]);
			sdt->stc[ct] = NULL;
		}
	}
					/* read BSDF data */
	sdata = ezxml_txt(ezxml_child(wdb, "ScatteringData"));
	if (!sdata || !next_token(&sdata)) {
		sprintf(SDerrorDetail, "Missing BSDF ScatteringData in '%s'",
				sd->name);
		return SDEformat;
	}
	sdt->stc[ct] = load_tree_data(&sdata, ndim);
	if (sdt->stc[ct] == NULL)
		return SDEformat;
	if (next_token(&sdata)) {	/* check for unconsumed characters */
		sprintf(SDerrorDetail,
			"Extra characters at end of ScatteringData in '%s'",
				sd->name);
		return SDEformat;
	}
					/* flatten branches where possible */
	sdt->stc[ct] = SDsimplifyTre(sdt->stc[ct]);
	if (sdt->stc[ct] == NULL)
		return SDEinternal;
					/* compute global quantities for Y */
	return (ct == tt_Y) ? get_extrema(df) : SDEnone;
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
			if ((st->u.v[n] -= val) < 0)
				st->u.v[n] = .0f;
	}
}

/* Subtract minimum Y value from BSDF */
static double
subtract_min_Y(SDNode *st)
{
	const float	vmaxmin = 1.5/M_PI;
	float		vmin;
					/* be sure to skip unused portion */
	if (st->ndim == 3) {
		int	n;
		vmin = vmaxmin;
		if (st->log2GR < 0) {
			for (n = 0; n < 8; n += 2) {
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

	if ((vmin >= vmaxmin) | (vmin <= .01/M_PI))
		return .0;		/* not worth bothering about */

	SDsubtractTreVal(st, vmin);

	return M_PI * vmin;		/* return hemispherical value */
}

/* Struct used in callback to find RGB extrema */
typedef struct {
	SDNode	**stc;			/* original Y, u' & v' trees */
	float	rgb[3];			/* RGB value */
	SDNode	*new_stu, *new_stv;	/* replacement u' & v' trees */
} SDextRGBs;

/* Callback to find minimum RGB from Y value plus CIE (u',v') trees */
static int
get_min_RGB(float yval, const double *cmin, double csiz, void *cptr)
{
	SDextRGBs	*mp = (SDextRGBs *)cptr;
	double		cmax[SD_MAXDIM];
	float		rgb[3];

	if (mp->stc[tt_Y]->ndim == 3) {
		if (cmin[0] + .5*csiz >= .5)
			return 0;	/* ignore dead half of isotropic */
	} else
		cmax[3] = cmin[3] + csiz;
	cmax[0] = cmin[0] + csiz;
	cmax[1] = cmin[1] + csiz;
	cmax[2] = cmin[2] + csiz;
					/* average RGB color over voxel */
	SDyuv2rgb(yval, SDavgTreBox(mp->stc[tt_u], cmin, cmax),
			SDavgTreBox(mp->stc[tt_v], cmin, cmax), rgb);
					/* track smallest components */
	if (rgb[0] < mp->rgb[0]) mp->rgb[0] = rgb[0];
	if (rgb[1] < mp->rgb[1]) mp->rgb[1] = rgb[1];
	if (rgb[2] < mp->rgb[2]) mp->rgb[2] = rgb[2];
	return 0;
}

/* Callback to build adjusted u' tree */
static int
adjust_utree(float uprime, const double *cmin, double csiz, void *cptr)
{
	SDextRGBs	*mp = (SDextRGBs *)cptr;
	double		cmax[SD_MAXDIM];
	double		yval;
	float		rgb[3];
	C_COLOR		clr;

	if (mp->stc[tt_Y]->ndim == 3) {
		if (cmin[0] + .5*csiz >= .5)
			return 0;	/* ignore dead half of isotropic */
	} else
		cmax[3] = cmin[3] + csiz;
	cmax[0] = cmin[0] + csiz;
	cmax[1] = cmin[1] + csiz;
	cmax[2] = cmin[2] + csiz;
					/* average RGB color over voxel */
	SDyuv2rgb(yval=SDavgTreBox(mp->stc[tt_Y], cmin, cmax), uprime,
			SDavgTreBox(mp->stc[tt_v], cmin, cmax), rgb);
					/* subtract minimum (& clamp) */
	if ((rgb[0] -= mp->rgb[0]) < 1e-5*yval) rgb[0] = 1e-5*yval;
	if ((rgb[1] -= mp->rgb[1]) < 1e-5*yval) rgb[1] = 1e-5*yval;
	if ((rgb[2] -= mp->rgb[2]) < 1e-5*yval) rgb[2] = 1e-5*yval;
	c_fromSharpRGB(rgb, &clr);	/* compute new u' for adj. RGB */
	uprime = 4.*clr.cx/(-2.*clr.cx + 12.*clr.cy + 3.);
					/* assign in new u' tree */
	mp->new_stu = SDsetVoxel(mp->new_stu, mp->stc[tt_Y]->ndim,
					cmin, csiz, uprime);
	return -(mp->new_stu == NULL);
}

/* Callback to build adjusted v' tree */
static int
adjust_vtree(float vprime, const double *cmin, double csiz, void *cptr)
{
	SDextRGBs	*mp = (SDextRGBs *)cptr;
	double		cmax[SD_MAXDIM];
	double		yval;
	float		rgb[3];
	C_COLOR		clr;

	if (mp->stc[tt_Y]->ndim == 3) {
		if (cmin[0] + .5*csiz >= .5)
			return 0;	/* ignore dead half of isotropic */
	} else
		cmax[3] = cmin[3] + csiz;
	cmax[0] = cmin[0] + csiz;
	cmax[1] = cmin[1] + csiz;
	cmax[2] = cmin[2] + csiz;
					/* average RGB color over voxel */
	SDyuv2rgb(yval=SDavgTreBox(mp->stc[tt_Y], cmin, cmax),
			SDavgTreBox(mp->stc[tt_u], cmin, cmax),
			vprime, rgb);
					/* subtract minimum (& clamp) */
	if ((rgb[0] -= mp->rgb[0]) < 1e-5*yval) rgb[0] = 1e-5*yval;
	if ((rgb[1] -= mp->rgb[1]) < 1e-5*yval) rgb[1] = 1e-5*yval;
	if ((rgb[2] -= mp->rgb[2]) < 1e-5*yval) rgb[2] = 1e-5*yval;
	c_fromSharpRGB(rgb, &clr);	/* compute new v' for adj. RGB */
	vprime = 9.*clr.cy/(-2.*clr.cx + 12.*clr.cy + 3.);
					/* assign in new v' tree */
	mp->new_stv = SDsetVoxel(mp->new_stv, mp->stc[tt_Y]->ndim,
					cmin, csiz, vprime);
	return -(mp->new_stv == NULL);
}

/* Subtract minimum (diffuse) color and return luminance & CIE (x,y) */
static double
subtract_min_RGB(C_COLOR *cs, SDNode *stc[])
{
	SDextRGBs	my_min;
	double		ymin;

	my_min.stc = stc;
	my_min.rgb[0] = my_min.rgb[1] = my_min.rgb[2] = FHUGE;
	my_min.new_stu = my_min.new_stv = NULL;
					/* get minimum RGB value */
	SDtraverseTre(stc[tt_Y], NULL, 0, get_min_RGB, &my_min);
					/* convert to C_COLOR */
	ymin = 	c_fromSharpRGB(my_min.rgb, cs);
	if ((ymin >= .5*FHUGE) | (ymin <= .01/M_PI))
		return .0;		/* close to zero or no tree */
					/* adjust u' & v' trees */
	SDtraverseTre(stc[tt_u], NULL, 0, adjust_utree, &my_min);
	SDtraverseTre(stc[tt_v], NULL, 0, adjust_vtree, &my_min);
	SDfreeTre(stc[tt_u]); SDfreeTre(stc[tt_v]);
	stc[tt_u] = SDsimplifyTre(my_min.new_stu);
	stc[tt_v] = SDsimplifyTre(my_min.new_stv);
					/* subtract Y & return hemispherical */
	SDsubtractTreVal(stc[tt_Y], ymin);

	return M_PI * ymin;
}

/* Extract and separate diffuse portion of BSDF */
static void
extract_diffuse(SDValue *dv, SDSpectralDF *df)
{
	int	n;
	SDTre	*sdt;

	if (df == NULL || df->ncomp <= 0) {
		dv->spec = c_dfcolor;
		dv->cieY = .0;
		return;
	}
	sdt = (SDTre *)df->comp[0].dist;
					/* subtract minimum color/grayscale */
	if (sdt->stc[tt_u] != NULL && sdt->stc[tt_v] != NULL) {
		int	i = 3*(tt_RGB_coef[1] < .001);
		while (i--) {		/* initialize on first call */
			float	rgb[3];
			rgb[0] = rgb[1] = rgb[2] = .0f; rgb[i] = 1.f;
			tt_RGB_coef[i] = c_fromSharpRGB(rgb, &tt_RGB_prim[i]);
		}
		memcpy(df->comp[0].cspec, tt_RGB_prim, sizeof(tt_RGB_prim));
		dv->cieY = subtract_min_RGB(&dv->spec, sdt->stc);
	} else {
		df->comp[0].cspec[0] = dv->spec = c_dfcolor;
		dv->cieY = subtract_min_Y(sdt->stc[tt_Y]);
	}
	df->maxHemi -= dv->cieY;	/* adjust maximum hemispherical */
				
	c_ccvt(&dv->spec, C_CSXY);	/* make sure (x,y) is set */
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
		const char	*cnm = ezxml_txt(ezxml_child(wld,"Wavelength"));
		int		ct = -1;
		if (!strcasecmp(cnm, "Visible"))
			ct = tt_Y;
		else if (!strcasecmp(cnm, "CIE-u"))
			ct = tt_u;
		else if (!strcasecmp(cnm, "CIE-v"))
			ct = tt_v;
		else
			continue;
		for (wdb = ezxml_child(wld, "WavelengthDataBlock");
					wdb != NULL; wdb = wdb->next)
			if ((ec = load_bsdf_data(sd, wdb, ct, rank)) != SDEnone)
				return ec;
	}
					/* separate diffuse components */
	extract_diffuse(&sd->rLambFront, sd->rf);
	extract_diffuse(&sd->rLambBack, sd->rb);
	extract_diffuse(&sd->tLambFront, sd->tf);
	if (sd->tb != NULL) {
		extract_diffuse(&sd->tLambBack, sd->tb);
		if (sd->tf == NULL)
			sd->tLambFront = sd->tLambBack;
	} else if (sd->tf != NULL)
		sd->tLambBack = sd->tLambFront;
					/* return success */
	return SDEnone;
}

/* Variable resolution BSDF methods */
const SDFunc SDhandleTre = {
	&SDgetTreBSDF,
	&SDqueryTreProjSA,
	&SDgetTreCDist,
	&SDsampTreCDist,
	&SDFreeBTre,
};
