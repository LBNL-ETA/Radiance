#ifndef lint
static const char RCSid[] = "$Id: bsdf_t.c,v 3.6 2011/04/24 19:39:21 greg Exp $";
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
static const unsigned	iwmax = (1<<(sizeof(unsigned)*4))-1;

/* Struct used for our distribution-building callback */
typedef struct {
	int		wmin;		/* minimum square size so far */
	int		wmax;		/* maximum square size */
	int		nic;		/* number of input coordinates */
	int		alen;		/* current array length */
	int		nall;		/* number of allocated entries */
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
				((1<<nd) - 1)*sizeof(st->u.t[0]));
		if (st != NULL)
			memset(st->u.t, 0, sizeof(st->u.t[0])<<nd);
	} else
		st = (SDNode *)malloc(sizeof(SDNode) +
				((1 << nd*lg) - 1)*sizeof(st->u.v[0]));
		
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

/* Add up N-dimensional hypercube array values over the given box */
static double
SDiterSum(const float *va, int nd, int siz, const int *imin, const int *imax)
{
	double		sum = .0;
	unsigned	skipsiz = 1;
	int		i;
	
	for (i = nd; --i > 0; )
		skipsiz *= siz;
	if (skipsiz == 1)
		for (i = *imin; i < *imax; i++)
			sum += va[i];
	else
		for (i = *imin; i < *imax; i++)
			sum += SDiterSum(va + i*skipsiz,
					nd-1, siz, imin+1, imax+1);
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
		imin[i] = (bmin[i] <= .0) ? 0 
				: (int)((1 << st->log2GR)*bmin[i]);
		imax[i] = (bmax[i] >= 1.) ? (1 << st->log2GR)
				: (int)((1 << st->log2GR)*bmax[i] + .999999);
		n *= imax[i] - imin[i];
	}
	if (!n)
		return .0;
	
	return SDiterSum(st->u.v, st->ndim, 1 << st->log2GR, imin, imax) /
			(double)n;
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
					/* check transmission */
	if (!sdt->isxmit ^ outVec[2] > 0 ^ inVec[2] > 0)
		return -1.;
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
	bmax[0] = bmin[0] + wid;
	bmax[1] = bmin[1] + wid;
	hilbert_box_vtx(2, sizeof(bitmask_t), sizeof(unsigned)*4,
							1, bmin, bmax);
	sp->darr[sp->alen].hent = hilbert_c2i(2, sizeof(unsigned)*4, bmin);
	sp->darr[sp->alen].wid = wid;
	sp->darr[sp->alen].bsdf = val;
	sp->alen++;			/* on to the next entry */
	return 0;
}

/* Scaffold comparison function for qsort -- ascending Hilbert index */
static int
sscmp(const void *p1, const void *p2)
{
	return (int)((*(const struct outdir_s *)p1).hent -
			(*(const struct outdir_s *)p2).hent);
}

/* Create a new cumulative distribution for the given input direction */
static SDTreCDst *
make_cdist(const SDTre *sdt, const double *pos)
{
	const unsigned	cumlmax = ~0;
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
	scale = (double)myScaffold.wmin / iwmax;
	for (i = myScaffold.nic; i--; ) {
		cd->clim[i][0] = floor(pos[i]/scale + .5) * scale;
		cd->clim[i][1] = cd->clim[i][0] + scale;
	}
	cd->max_psa = myScaffold.wmax / (double)iwmax;
	cd->max_psa *= cd->max_psa * M_PI;
	cd->isxmit = sdt->isxmit;
	cd->cTotal = 1e-20;		/* compute directional total */
	sp = myScaffold.darr;
	for (i = myScaffold.alen; i--; sp++)
		cd->cTotal += sp->bsdf * (double)sp->wid * sp->wid;
	cursum = .0;			/* go back and get cumulative values */
	scale = (double)cumlmax / cd->cTotal;
	sp = myScaffold.darr;
	for (i = 0; i < cd->calen; i++, sp++) {
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
			if (qflags == SDqueryVal)
				*psa = M_PI;
			return SDEnone;
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
	const unsigned	maxval = ~0;
	const SDTreCDst	*cd = (const SDTreCDst *)cdp;
	const unsigned	target = randX*maxval;
	bitmask_t	hndx, hcoord[2];
	double		gpos[3];
	int		i, iupper, ilower;
					/* check arguments */
	if ((ioVec == NULL) | (cd == NULL))
		return SDEargument;
					/* binary search to find position */
	ilower = 0; iupper = cd->calen;
	while ((i = (iupper + ilower) >> 1) != ilower)
		if ((long)target >= (long)cd->carr[i].cuml)
			ilower = i;
		else
			iupper = i;
					/* localize random position */
	randX = (randX*maxval - cd->carr[ilower].cuml) /
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
	gpos[2] = 1. - gpos[0]*gpos[0] - gpos[1]*gpos[1];
	if (gpos[2] > 0)		/* paranoia, I hope */
		gpos[2] = sqrt(gpos[2]);
	if (ioVec[2] > 0 ^ !cd->isxmit)
		gpos[2] = -gpos[2];
	VCOPY(ioVec, gpos);
	return SDEnone;
}

/* Load a variable-resolution BSDF tree from an open XML file */
SDError
SDloadTre(SDData *sd, ezxml_t wtl)
{
	return SDEsupport;
}

/* Variable resolution BSDF methods */
SDFunc SDhandleTre = {
	&SDgetTreBSDF,
	&SDqueryTreProjSA,
	&SDgetTreCDist,
	&SDsampTreCDist,
	&SDFreeBTre,
};
