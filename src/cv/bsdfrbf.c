#ifndef lint
static const char RCSid[] = "$Id: bsdfrbf.c,v 2.23 2014/03/21 00:42:46 greg Exp $";
#endif
/*
 * Radial basis function representation for BSDF data.
 *
 *	G. Ward
 */

/****************************************************************
1) Collect samples into a grid using the Shirley-Chiu
	angular mapping from a hemisphere to a square.

2) Compute an adaptive quadtree by subdividing the grid so that
	each leaf node has at least one sample up to as many
	samples as fit nicely on a plane to within a certain
	MSE tolerance.

3) Place one Gaussian lobe at each leaf node in the quadtree,
	sizing it to have a radius equal to the leaf size and
	a volume equal to the energy in that node.
*****************************************************************/

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bsdfrep.h"

#ifndef RSCA
#define RSCA		2.2		/* radius scaling factor (empirical) */
#endif
#ifndef SMOOTH_MSE
#define SMOOTH_MSE	5e-5		/* acceptable mean squared error */
#endif
#ifndef SMOOTH_MSER
#define SMOOTH_MSER	0.03		/* acceptable relative MSE */
#endif
#define MAX_RAD		(GRIDRES/8)	/* maximum lobe radius */

#define	RBFALLOCB	10		/* RBF allocation block size */

				/* our loaded grid or comparison DSFs */
GRIDVAL			dsf_grid[GRIDRES][GRIDRES];

/* Start new DSF input grid */
void
new_bsdf_data(double new_theta, double new_phi)
{
	if (!new_input_direction(new_theta, new_phi))
		exit(1);
	memset(dsf_grid, 0, sizeof(dsf_grid));
}

/* Add BSDF data point */
void
add_bsdf_data(double theta_out, double phi_out, double val, int isDSF)
{
	FVECT	ovec;
	int	pos[2];

	if (!output_orient)		/* check output orientation */
		output_orient = 1 - 2*(theta_out > 90.);
	else if (output_orient > 0 ^ theta_out < 90.) {
		fputs("Cannot handle output angles on both sides of surface\n",
				stderr);
		exit(1);
	}
	ovec[2] = sin((M_PI/180.)*theta_out);
	ovec[0] = cos((M_PI/180.)*phi_out) * ovec[2];
	ovec[1] = sin((M_PI/180.)*phi_out) * ovec[2];
	ovec[2] = sqrt(1. - ovec[2]*ovec[2]);

	if (!isDSF)
		val *= ovec[2];		/* convert from BSDF to DSF */

					/* update BSDF histogram */
	if (val < BSDF2BIG*ovec[2] && val > BSDF2SML*ovec[2])
		++bsdf_hist[histndx(val/ovec[2])];

	pos_from_vec(pos, ovec);

	dsf_grid[pos[0]][pos[1]].sum.v += val;
	dsf_grid[pos[0]][pos[1]].sum.n++;
}

/* Compute minimum BSDF from histogram (does not clear) */
static void
comp_bsdf_min()
{
	unsigned long	cnt, target;
	int		i;

	cnt = 0;
	for (i = HISTLEN; i--; )
		cnt += bsdf_hist[i];
	if (!cnt) {				/* shouldn't happen */
		bsdf_min = 0;
		return;
	}
	target = cnt/100;			/* ignore bottom 1% */
	cnt = 0;
	for (i = 0; cnt <= target; i++)
		cnt += bsdf_hist[i];
	bsdf_min = histval(i-1);
}

/* Determine if the given region is empty of grid samples */
static int
empty_region(int x0, int x1, int y0, int y1)
{
	int	x, y;

	for (x = x0; x < x1; x++)
	    for (y = y0; y < y1; y++)
		if (dsf_grid[x][y].sum.n)
			return(0);
	return(1);
}

/* Determine if the given region is smooth enough to be a single lobe */
static int
smooth_region(int x0, int x1, int y0, int y1)
{
	RREAL	rMtx[3][3];
	FVECT	xvec;
	double	A, B, C, nvs, sqerr;
	int	x, y, n;
					/* compute planar regression */
	memset(rMtx, 0, sizeof(rMtx));
	memset(xvec, 0, sizeof(xvec));
	for (x = x0; x < x1; x++)
	    for (y = y0; y < y1; y++)
		if ((n = dsf_grid[x][y].sum.n) > 0) {
			double	z = dsf_grid[x][y].sum.v;
			rMtx[0][0] += x*x*(double)n;
			rMtx[0][1] += x*y*(double)n;
			rMtx[0][2] += x*(double)n;
			rMtx[1][1] += y*y*(double)n;
			rMtx[1][2] += y*(double)n;
			rMtx[2][2] += (double)n;
			xvec[0] += x*z;
			xvec[1] += y*z;
			xvec[2] += z;
		}
	rMtx[1][0] = rMtx[0][1];
	rMtx[2][0] = rMtx[0][2];
	rMtx[2][1] = rMtx[1][2];
	nvs = rMtx[2][2];
	if (SDinvXform(rMtx, rMtx) != SDEnone)
		return(1);		/* colinear values */
	A = DOT(rMtx[0], xvec);
	B = DOT(rMtx[1], xvec);
	C = DOT(rMtx[2], xvec);
	sqerr = 0.0;			/* compute mean squared error */
	for (x = x0; x < x1; x++)
	    for (y = y0; y < y1; y++)
		if ((n = dsf_grid[x][y].sum.n) > 0) {
			double	d = A*x + B*y + C - dsf_grid[x][y].sum.v/n;
			sqerr += n*d*d;
		}
	if (sqerr <= nvs*SMOOTH_MSE)	/* below absolute MSE threshold? */
		return(1);
					/* OR below relative MSE threshold? */
	return(sqerr*nvs <= xvec[2]*xvec[2]*SMOOTH_MSER);
}

/* Create new lobe based on integrated samples in region */
static int
create_lobe(RBFVAL *rvp, int x0, int x1, int y0, int y1)
{
	double	vtot = 0.0;
	int	nv = 0;
	double	rad;
	int	x, y;
					/* compute average for region */
	for (x = x0; x < x1; x++)
	    for (y = y0; y < y1; y++) {
		vtot += dsf_grid[x][y].sum.v;
		nv += dsf_grid[x][y].sum.n;
	    }
	if (!nv) {
		fprintf(stderr, "%s: internal - missing samples in create_lobe\n",
				progname);
		exit(1);
	}
	if (vtot <= 0)			/* only create positive lobes */
		return(0);
					/* peak value based on integral */
	vtot *= (x1-x0)*(y1-y0)*(2.*M_PI/GRIDRES/GRIDRES)/(double)nv;
	rad = (RSCA/(double)GRIDRES)*(x1-x0);
	rvp->peak =  vtot / ((2.*M_PI) * rad*rad);
	rvp->crad = ANG2R(rad);
	rvp->gx = (x0+x1)>>1;
	rvp->gy = (y0+y1)>>1;
	return(1);
}

/* Recursive function to build radial basis function representation */
static int
build_rbfrep(RBFVAL **arp, int *np, int x0, int x1, int y0, int y1)
{
	int	xmid = (x0+x1)>>1;
	int	ymid = (y0+y1)>>1;
	int	branched[4];
	int	nadded, nleaves;
					/* need to make this a leaf? */
	if (empty_region(x0, xmid, y0, ymid) ||
			empty_region(xmid, x1, y0, ymid) ||
			empty_region(x0, xmid, ymid, y1) ||
			empty_region(xmid, x1, ymid, y1))
		return(0);
					/* add children (branches+leaves) */
	if ((branched[0] = build_rbfrep(arp, np, x0, xmid, y0, ymid)) < 0)
		return(-1);
	if ((branched[1] = build_rbfrep(arp, np, xmid, x1, y0, ymid)) < 0)
		return(-1);
	if ((branched[2] = build_rbfrep(arp, np, x0, xmid, ymid, y1)) < 0)
		return(-1);
	if ((branched[3] = build_rbfrep(arp, np, xmid, x1, ymid, y1)) < 0)
		return(-1);
	nadded = branched[0] + branched[1] + branched[2] + branched[3];
	nleaves = !branched[0] + !branched[1] + !branched[2] + !branched[3];
	if (!nleaves)			/* nothing but branches? */
		return(nadded);
					/* combine 4 leaves into 1? */
	if ((nleaves == 4) & (x1-x0 <= MAX_RAD) &&
			smooth_region(x0, x1, y0, y1))
		return(0);
					/* need more array space? */
	if ((*np+nleaves-1)>>RBFALLOCB != (*np-1)>>RBFALLOCB) {
		*arp = (RBFVAL *)realloc(*arp,
				sizeof(RBFVAL)*(*np+nleaves-1+(1<<RBFALLOCB)));
		if (*arp == NULL)
			return(-1);
	}
					/* create lobes for leaves */
	if (!branched[0] && create_lobe(*arp+*np, x0, xmid, y0, ymid)) {
		++(*np); ++nadded;
	}
	if (!branched[1] && create_lobe(*arp+*np, xmid, x1, y0, ymid)) {
		++(*np); ++nadded;
	}
	if (!branched[2] && create_lobe(*arp+*np, x0, xmid, ymid, y1)) {
		++(*np); ++nadded;
	}
	if (!branched[3] && create_lobe(*arp+*np, xmid, x1, ymid, y1)) {
		++(*np); ++nadded;
	}
	return(nadded);
}

/* Count up filled nodes and build RBF representation from current grid */ 
RBFNODE *
make_rbfrep()
{
	RBFNODE	*newnode;
	RBFVAL	*rbfarr;
	int	nn;
				/* compute minimum BSDF */
	comp_bsdf_min();
				/* create RBF node list */
	rbfarr = NULL; nn = 0;
	if (build_rbfrep(&rbfarr, &nn, 0, GRIDRES, 0, GRIDRES) <= 0) {
		if (nn)
			goto memerr;
		fprintf(stderr, "%s: no usable data in make_rbfrep()\n",
				progname);
		exit(1);
	}
				/* (re)allocate RBF array */
	newnode = (RBFNODE *)realloc(rbfarr,
			sizeof(RBFNODE) + sizeof(RBFVAL)*(nn-1));
	if (newnode == NULL)
		goto memerr;
				/* copy computed lobes into RBF node */
	memmove(newnode->rbfa, newnode, sizeof(RBFVAL)*nn);
	newnode->ord = -1;
	newnode->next = NULL;
	newnode->ejl = NULL;
	newnode->invec[2] = sin((M_PI/180.)*theta_in_deg);
	newnode->invec[0] = cos((M_PI/180.)*phi_in_deg)*newnode->invec[2];
	newnode->invec[1] = sin((M_PI/180.)*phi_in_deg)*newnode->invec[2];
	newnode->invec[2] = input_orient*sqrt(1. - newnode->invec[2]*newnode->invec[2]);
	newnode->vtotal = .0;
	newnode->nrbf = nn;
				/* compute sum for normalization */
	while (nn-- > 0)
		newnode->vtotal += rbf_volume(&newnode->rbfa[nn]);
#ifdef DEBUG
	fprintf(stderr, "Built RBF with %d lobes\n", newnode->nrbf);
	fprintf(stderr, "Integrated DSF at (%.1f,%.1f) deg. is %.2f\n",
			get_theta180(newnode->invec), get_phi360(newnode->invec),
			newnode->vtotal);
#endif
	insert_dsf(newnode);

	return(newnode);
memerr:
	fprintf(stderr, "%s: Out of memory in make_rbfrep()\n", progname);
	exit(1);
}
