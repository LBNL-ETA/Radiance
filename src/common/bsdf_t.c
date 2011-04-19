#ifndef lint
static const char RCSid[] = "$Id: bsdf_t.c,v 3.5 2011/04/19 21:31:22 greg Exp $";
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
				"Cannot allocate %d branch BSDF tree", nd);
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
SDfreeTre(void *p)
{
	SDNode	*st = (SDNode *)p;
	int	i;

	if (st == NULL)
		return;
	for (i = (st->log2GR < 0) << st->ndim; i--; )
		SDfreeTre(st->u.t[i]);
	free((void *)st);
}

#if 0	/* Unused and untested routines */

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
SDavgBox(const SDNode *st, const double *bmin, const double *bmax)
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
				sum += w * SDavgBox(st->u.t[n], sbmin, sbmax);
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

#endif	/* 0 */

/* Look up tree value at the given grid position */
static float
SDlookupTre(const SDNode *st, const double *pos)
{
	double	spos[SD_MAXDIM];
	int	i, n, t;
					/* climb the tree */
	while (st->log2GR < 0) {
		n = 0;			/* move to appropriate branch */
		for (i = st->ndim; i--; ) {
			spos[i] = 2.*pos[i];
			t = (spos[i] >= 1.);
			n |= t<<i;
			spos[i] -= (double)t;
		}
		st = st->u.t[n];	/* avoids tail recursion */
		pos = spos;
	}
	n = t = 0;			/* find grid array index */
	for (i = st->ndim; i--; ) {
		n += (int)((1<<st->log2GR)*pos[i]) << t;
		t += st->log2GR;
	}
	return st->u.v[n];		/* XXX no interpolation */
}

/* Compute non-diffuse component for variable-resolution BSDF */
static int
SDgetTreBSDF(float coef[SDmaxCh], const FVECT outVec,
				const FVECT inVec, const void *dist)
{
	const SDNode	*st = (const SDNode *)dist;
	double		gridPos[4];
					/* convert vector coordinates */
	if (st->ndim == 3) {		/* reduce for isotropic BSDF? */
		static const FVECT	zvec = {.0, .0, 1.};
		FVECT			rOutVec;
		spinvector(rOutVec, outVec, zvec, -atan2(inVec[1],inVec[0]));
		SDdisk2square(gridPos, rOutVec[0], rOutVec[1]);
		gridPos[2] = .5 - .5*sqrt(inVec[0]*inVec[0] + inVec[1]*inVec[1]);
	} else if (st->ndim == 4) {	/* XXX no component identity checks */
		SDdisk2square(gridPos, outVec[0], outVec[1]);
		SDdisk2square(gridPos+2, -inVec[0], -inVec[1]);
	} else
		return 0;		/* should be internal error */
					/* get nearest BSDF value */
	coef[0] = SDlookupTre(st, gridPos);
	return 1;			/* monochromatic for now */
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
	NULL,
	NULL,
	NULL,
	&SDfreeTre,
};
