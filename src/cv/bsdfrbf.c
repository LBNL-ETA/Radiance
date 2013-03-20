#ifndef lint
static const char RCSid[] = "$Id: bsdfrbf.c,v 2.4 2013/03/20 01:00:22 greg Exp $";
#endif
/*
 * Radial basis function representation for BSDF data.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bsdfrep.h"

#ifndef RSCA
#define RSCA		2.7		/* radius scaling factor (empirical) */
#endif
				/* our loaded grid for this incident angle */
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

	dsf_grid[pos[0]][pos[1]].vsum += val;
	dsf_grid[pos[0]][pos[1]].nval++;
}

/* Compute radii for non-empty bins */
/* (distance to furthest empty bin for which non-empty bin is the closest) */
static void
compute_radii(void)
{
	unsigned int	fill_grid[GRIDRES][GRIDRES];
	unsigned short	fill_cnt[GRIDRES][GRIDRES];
	FVECT		ovec0, ovec1;
	double		ang2, lastang2;
	int		r, i, j, jn, ii, jj, inear, jnear;

	r = GRIDRES/2;				/* proceed in zig-zag */
	for (i = 0; i < GRIDRES; i++)
	    for (jn = 0; jn < GRIDRES; jn++) {
		j = (i&1) ? jn : GRIDRES-1-jn;
		if (dsf_grid[i][j].nval)	/* find empty grid pos. */
			continue;
		ovec_from_pos(ovec0, i, j);
		inear = jnear = -1;		/* find nearest non-empty */
		lastang2 = M_PI*M_PI;
		for (ii = i-r; ii <= i+r; ii++) {
		    if (ii < 0) continue;
		    if (ii >= GRIDRES) break;
		    for (jj = j-r; jj <= j+r; jj++) {
			if (jj < 0) continue;
			if (jj >= GRIDRES) break;
			if (!dsf_grid[ii][jj].nval)
				continue;
			ovec_from_pos(ovec1, ii, jj);
			ang2 = 2. - 2.*DOT(ovec0,ovec1);
			if (ang2 >= lastang2)
				continue;
			lastang2 = ang2;
			inear = ii; jnear = jj;
		    }
		}
		if (inear < 0) {
			fprintf(stderr,
				"%s: Could not find non-empty neighbor!\n",
					progname);
			exit(1);
		}
		ang2 = sqrt(lastang2);
		r = ANG2R(ang2);		/* record if > previous */
		if (r > dsf_grid[inear][jnear].crad)
			dsf_grid[inear][jnear].crad = r;
						/* next search radius */
		r = ang2*(2.*GRIDRES/M_PI) + 3;
	    }
						/* blur radii over hemisphere */
	memset(fill_grid, 0, sizeof(fill_grid));
	memset(fill_cnt, 0, sizeof(fill_cnt));
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++) {
		if (!dsf_grid[i][j].crad)
			continue;		/* missing distance */
		r = R2ANG(dsf_grid[i][j].crad)*(2.*RSCA*GRIDRES/M_PI);
		for (ii = i-r; ii <= i+r; ii++) {
		    if (ii < 0) continue;
		    if (ii >= GRIDRES) break;
		    for (jj = j-r; jj <= j+r; jj++) {
			if (jj < 0) continue;
			if (jj >= GRIDRES) break;
			if ((ii-i)*(ii-i) + (jj-j)*(jj-j) > r*r)
				continue;
			fill_grid[ii][jj] += dsf_grid[i][j].crad;
			fill_cnt[ii][jj]++;
		    }
		}
	    }
						/* copy back blurred radii */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (fill_cnt[i][j])
			dsf_grid[i][j].crad = fill_grid[i][j]/fill_cnt[i][j];
}

/* Cull points for more uniform distribution, leave all nval 0 or 1 */
static void
cull_values(void)
{
	FVECT	ovec0, ovec1;
	double	maxang, maxang2;
	int	i, j, ii, jj, r;
						/* simple greedy algorithm */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++) {
		if (!dsf_grid[i][j].nval)
			continue;
		if (!dsf_grid[i][j].crad)
			continue;		/* shouldn't happen */
		ovec_from_pos(ovec0, i, j);
		maxang = 2.*R2ANG(dsf_grid[i][j].crad);
		if (maxang > ovec0[2])		/* clamp near horizon */
			maxang = ovec0[2];
		r = maxang*(2.*GRIDRES/M_PI) + 1;
		maxang2 = maxang*maxang;
		for (ii = i-r; ii <= i+r; ii++) {
		    if (ii < 0) continue;
		    if (ii >= GRIDRES) break;
		    for (jj = j-r; jj <= j+r; jj++) {
			if (jj < 0) continue;
			if (jj >= GRIDRES) break;
			if (!dsf_grid[ii][jj].nval)
				continue;
			if ((ii == i) & (jj == j))
				continue;	/* don't get self-absorbed */
			ovec_from_pos(ovec1, ii, jj);
			if (2. - 2.*DOT(ovec0,ovec1) >= maxang2)
				continue;
						/* absorb sum */
			dsf_grid[i][j].vsum += dsf_grid[ii][jj].vsum;
			dsf_grid[i][j].nval += dsf_grid[ii][jj].nval;
						/* keep value, though */
			dsf_grid[ii][jj].vsum /= (float)dsf_grid[ii][jj].nval;
			dsf_grid[ii][jj].nval = 0;
		    }
		}
	    }
						/* final averaging pass */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (dsf_grid[i][j].nval > 1) {
			dsf_grid[i][j].vsum /= (float)dsf_grid[i][j].nval;
			dsf_grid[i][j].nval = 1;
		}
}

/* Count up filled nodes and build RBF representation from current grid */ 
RBFNODE *
make_rbfrep(void)
{
	int	niter = 16;
	double	lastVar, thisVar = 100.;
	int	nn;
	RBFNODE	*newnode;
	RBFVAL	*itera;
	int	i, j;
				/* compute RBF radii */
	compute_radii();
				/* coagulate lobes */
	cull_values();
	nn = 0;			/* count selected bins */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		nn += dsf_grid[i][j].nval;
				/* allocate RBF array */
	newnode = (RBFNODE *)malloc(sizeof(RBFNODE) + sizeof(RBFVAL)*(nn-1));
	if (newnode == NULL)
		goto memerr;
	newnode->ord = -1;
	newnode->next = NULL;
	newnode->ejl = NULL;
	newnode->invec[2] = sin((M_PI/180.)*theta_in_deg);
	newnode->invec[0] = cos((M_PI/180.)*phi_in_deg)*newnode->invec[2];
	newnode->invec[1] = sin((M_PI/180.)*phi_in_deg)*newnode->invec[2];
	newnode->invec[2] = input_orient*sqrt(1. - newnode->invec[2]*newnode->invec[2]);
	newnode->vtotal = 0;
	newnode->nrbf = nn;
	nn = 0;			/* fill RBF array */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (dsf_grid[i][j].nval) {
			newnode->rbfa[nn].peak = dsf_grid[i][j].vsum;
			newnode->rbfa[nn].crad = RSCA*dsf_grid[i][j].crad + .5;
			newnode->rbfa[nn].gx = i;
			newnode->rbfa[nn].gy = j;
			++nn;
		}
				/* iterate to improve interpolation accuracy */
	itera = (RBFVAL *)malloc(sizeof(RBFVAL)*newnode->nrbf);
	if (itera == NULL)
		goto memerr;
	memcpy(itera, newnode->rbfa, sizeof(RBFVAL)*newnode->nrbf);
	do {
		double	dsum = 0, dsum2 = 0;
		nn = 0;
		for (i = 0; i < GRIDRES; i++)
		    for (j = 0; j < GRIDRES; j++)
			if (dsf_grid[i][j].nval) {
				FVECT	odir;
				double	corr;
				ovec_from_pos(odir, i, j);
				itera[nn++].peak *= corr =
					dsf_grid[i][j].vsum /
						eval_rbfrep(newnode, odir);
				dsum += 1. - corr;
				dsum2 += (1.-corr)*(1.-corr);
			}
		memcpy(newnode->rbfa, itera, sizeof(RBFVAL)*newnode->nrbf);
		lastVar = thisVar;
		thisVar = dsum2/(double)nn;
#ifdef DEBUG
		fprintf(stderr, "Avg., RMS error: %.1f%%  %.1f%%\n",
					100.*dsum/(double)nn,
					100.*sqrt(thisVar));
#endif
	} while (--niter > 0 && lastVar-thisVar > 0.02*lastVar);

	free(itera);
	nn = 0;			/* compute sum for normalization */
	while (nn < newnode->nrbf)
		newnode->vtotal += rbf_volume(&newnode->rbfa[nn++]);
#ifdef DEBUG
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
