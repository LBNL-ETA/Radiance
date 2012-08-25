#ifndef lint
static const char RCSid[] = "$Id: pabopto2xml.c,v 2.5 2012/08/25 22:39:03 greg Exp $";
#endif
/*
 * Convert PAB-Opto measurements to XML format using tensor tree representation
 * Employs Bonneel et al. Earth Mover's Distance interpolant.
 *
 *	G.Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "bsdf.h"

#ifndef GRIDRES
#define GRIDRES		200		/* max. grid resolution per side */
#endif

#define RSCA		2.7		/* radius scaling factor (empirical) */

#define R2ANG(c)	(((c)+.5)*(M_PI/(1<<16)))
#define ANG2R(r)	(int)((r)*((1<<16)/M_PI))

typedef struct {
	float		vsum;		/* DSF sum */
	unsigned short	nval;		/* number of values in sum */
	unsigned short	crad;		/* radius (coded angle) */
} GRIDVAL;			/* grid value */

typedef struct {
	float		peak;		/* lobe value at peak */
	unsigned short	crad;		/* radius (coded angle) */
	unsigned char	gx, gy;		/* grid position */
} RBFVAL;			/* radial basis function value */

typedef struct s_rbflist {
	struct s_rbflist	*next;		/* next in our RBF list */
	FVECT			invec;		/* incident vector direction */
	int			nrbf;		/* number of RBFs */
	RBFVAL			rbfa[1];	/* RBF array (extends struct) */
} RBFLIST;			/* RBF representation of DSF @ 1 incidence */

				/* our loaded grid for this incident angle */
static double	theta_in_deg, phi_in_deg;
static GRIDVAL	dsf_grid[GRIDRES][GRIDRES];

				/* processed incident DSF measurements */
static RBFLIST	*dsf_list = NULL;

/* Compute outgoing vector from grid position */
static void
vec_from_pos(FVECT vec, int xpos, int ypos)
{
	double	uv[2];
	double	r2;
	
	SDsquare2disk(uv, (1./GRIDRES)*(xpos+.5), (1./GRIDRES)*(ypos+.5));
				/* uniform hemispherical projection */
	r2 = uv[0]*uv[0] + uv[1]*uv[1];
	vec[0] = vec[1] = sqrt(2. - r2);
	vec[0] *= uv[0];
	vec[1] *= uv[1];
	vec[2] = 1. - r2;
}

/* Compute grid position from normalized outgoing vector */
static void
pos_from_vec(int pos[2], const FVECT vec)
{
	double	sq[2];		/* uniform hemispherical projection */
	double	norm = 1./sqrt(1. + vec[2]);

	SDdisk2square(sq, vec[0]*norm, vec[1]*norm);

	pos[0] = (int)(sq[0]*GRIDRES);
	pos[1] = (int)(sq[1]*GRIDRES);
}

/* Evaluate RBF for DSF at the given normalized outgoing direction */
static double
eval_rbfrep(const RBFLIST *rp, const FVECT outvec)
{
	double		res = .0;
	const RBFVAL	*rbfp;
	FVECT		odir;
	double		sig2;
	int		n;

	rbfp = rp->rbfa;
	for (n = rp->nrbf; n--; rbfp++) {
		vec_from_pos(odir, rbfp->gx, rbfp->gy);
		sig2 = R2ANG(rbfp->crad);
		sig2 = (DOT(odir,outvec) - 1.) / (sig2*sig2);
		if (sig2 > -19.)
			res += rbfp->peak * exp(sig2);
	}
	return(res);
}

/* Count up filled nodes and build RBF representation from current grid */ 
static RBFLIST *
make_rbfrep(void)
{
	int	niter = 6;
	int	nn;
	RBFLIST	*newnode;
	int	i, j;
	
	nn = 0;			/* count selected bins */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		nn += (dsf_grid[i][j].nval > 0);
				/* allocate RBF array */
	newnode = (RBFLIST *)malloc(sizeof(RBFLIST) + sizeof(RBFVAL)*(nn-1));
	if (newnode == NULL) {
		fputs("Out of memory in make_rbfrep\n", stderr);
		exit(1);
	}
	newnode->next = NULL;
	newnode->invec[2] = sin(M_PI/180.*theta_in_deg);
	newnode->invec[0] = cos(M_PI/180.*phi_in_deg)*newnode->invec[2];
	newnode->invec[1] = sin(M_PI/180.*phi_in_deg)*newnode->invec[2];
	newnode->invec[2] = sqrt(1. - newnode->invec[2]*newnode->invec[2]);
	newnode->nrbf = nn;
	nn = 0;			/* fill RBF array */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (dsf_grid[i][j].nval) {
			newnode->rbfa[nn].peak =
					dsf_grid[i][j].vsum /=
						(double)dsf_grid[i][j].nval;
			dsf_grid[i][j].nval = 1;
			newnode->rbfa[nn].crad = RSCA*dsf_grid[i][j].crad + .5;
			newnode->rbfa[nn].gx = i;
			newnode->rbfa[nn].gy = j;
			++nn;
		}
				/* iterate for better convergence */
	while (niter--) {
		double	dsum = .0, dsum2 = .0;
		nn = 0;
		for (i = 0; i < GRIDRES; i++)
		    for (j = 0; j < GRIDRES; j++)
			if (dsf_grid[i][j].nval) {
				FVECT	odir;
				/* double	corr; */
				vec_from_pos(odir, i, j);
				newnode->rbfa[nn++].peak *= /* corr = */
					dsf_grid[i][j].vsum /
						eval_rbfrep(newnode, odir);
				/*
				dsum += corr - 1.;
				dsum2 += (corr-1.)*(corr-1.);
				*/
			}
		/*
		fprintf(stderr, "Avg., RMS error: %.1f%%  %.1f%%\n",
					100.*dsum/(double)nn,
					100.*sqrt(dsum2/(double)nn));
		*/
	}
	newnode->next = dsf_list;
	return(dsf_list = newnode);
}

/* Load a set of measurements corresponding to a particular incident angle */
static int
load_bsdf_meas(const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	int	inp_is_DSF = -1;
	double	theta_out, phi_out, val;
	char	buf[2048];
	int	n, c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
	memset(dsf_grid, 0, sizeof(dsf_grid));
				/* read header information */
	while ((c = getc(fp)) == '#' || c == EOF) {
		if (fgets(buf, sizeof(buf), fp) == NULL) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			fclose(fp);
			return(0);
		}
		if (!strcmp(buf, "format: theta phi DSF\n")) {
			inp_is_DSF = 1;
			continue;
		}
		if (!strcmp(buf, "format: theta phi BSDF\n")) {
			inp_is_DSF = 0;
			continue;
		}
		if (sscanf(buf, "intheta %lf", &theta_in_deg) == 1)
			continue;
		if (sscanf(buf, "inphi %lf", &phi_in_deg) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&theta_in_deg, &phi_in_deg) == 2)
			continue;
	}
	if (inp_is_DSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		fclose(fp);
		return(0);
	}
	ungetc(c, fp);		/* read actual data */
	while (fscanf(fp, "%lf %lf %lf\n", &theta_out, &phi_out, &val) == 3) {
		FVECT	ovec;
		int	pos[2];

		ovec[2] = sin(M_PI/180.*theta_out);
		ovec[0] = cos(M_PI/180.*phi_out) * ovec[2];
		ovec[1] = sin(M_PI/180.*phi_out) * ovec[2];
		ovec[2] = sqrt(1. - ovec[2]*ovec[2]);

		if (!inp_is_DSF)
			val *= ovec[2];	/* convert from BSDF to DSF */

		pos_from_vec(pos, ovec);

		dsf_grid[pos[0]][pos[1]].vsum += val;
		dsf_grid[pos[0]][pos[1]].nval++;
	}
	n = 0;
	while ((c = getc(fp)) != EOF)
		n += !isspace(c);
	if (n) 
		fprintf(stderr,
			"%s: warning: %d unexpected characters past EOD\n",
				fname, n);
	fclose(fp);
	return(1);
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
		vec_from_pos(ovec0, i, j);
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
			vec_from_pos(ovec1, ii, jj);
			ang2 = 2. - 2.*DOT(ovec0,ovec1);
			if (ang2 >= lastang2)
				continue;
			lastang2 = ang2;
			inear = ii; jnear = jj;
		    }
		}
		if (inear < 0) {
			fputs("Could not find non-empty neighbor!\n", stderr);
			exit(1);
		}
		ang2 = sqrt(lastang2);
		r = ANG2R(ang2);		/* record if > previous */
		if (r > dsf_grid[inear][jnear].crad)
			dsf_grid[inear][jnear].crad = r;
						/* next search radius */
		r = ang2*(2.*GRIDRES/M_PI) + 1;
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
						/* copy back averaged radii */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (fill_cnt[i][j])
			dsf_grid[i][j].crad = fill_grid[i][j]/fill_cnt[i][j];
}

/* Cull points for more uniform distribution */
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
		vec_from_pos(ovec0, i, j);
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
			vec_from_pos(ovec1, ii, jj);
			if (2. - 2.*DOT(ovec0,ovec1) >= maxang2)
				continue;
						/* absorb sum */
			dsf_grid[i][j].vsum += dsf_grid[ii][jj].vsum;
			dsf_grid[i][j].nval += dsf_grid[ii][jj].nval;
						/* keep value, though */
			dsf_grid[ii][jj].vsum /= (double)dsf_grid[ii][jj].nval;
			dsf_grid[ii][jj].nval = 0;
		    }
		}
	    }
}


#if 1
/* Test main produces a Radiance model from the given input file */
int
main(int argc, char *argv[])
{
	char	buf[128];
	FILE	*pfp;
	double	bsdf;
	FVECT	dir;
	int	i, j, n;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s input.dat > output.rad\n", argv[0]);
		return(1);
	}
	if (!load_bsdf_meas(argv[1]))
		return(1);

	compute_radii();
	cull_values();
	make_rbfrep();
						/* produce spheres at meas. */
	puts("void plastic yellow\n0\n0\n5 .6 .4 .01 .04 .08\n");
	puts("void plastic pink\n0\n0\n5 .5 .05 .9 .04 .08\n");
	n = 0;
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (dsf_grid[i][j].vsum > .0f) {
			vec_from_pos(dir, i, j);
			bsdf = dsf_grid[i][j].vsum / dir[2];
			if (dsf_grid[i][j].nval) {
				printf("pink cone c%04d\n0\n0\n8\n", ++n);
				printf("\t%.6g %.6g %.6g\n",
					dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
				printf("\t%.6g %.6g %.6g\n",
					dir[0]*(bsdf+.005), dir[1]*(bsdf+.005),
					dir[2]*(bsdf+.005));
				puts("\t.003\t0\n");
			} else {
				vec_from_pos(dir, i, j);
				printf("yellow sphere s%04d\n0\n0\n", ++n);
				printf("4 %.6g %.6g %.6g .0015\n\n",
					dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
			}
		}
						/* output continuous surface */
	puts("void trans tgreen\n0\n0\n7 .7 1 .7 .04 .04 .9 .9\n");
	fflush(stdout);
	sprintf(buf, "gensurf tgreen bsdf - - - %d %d", GRIDRES-1, GRIDRES-1);
	pfp = popen(buf, "w");
	if (pfp == NULL) {
		fputs(buf, stderr);
		fputs(": cannot start command\n", stderr);
		return(1);
	}
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++) {
		vec_from_pos(dir, i, j);
		bsdf = eval_rbfrep(dsf_list, dir) / dir[2];
		fprintf(pfp, "%.8e %.8e %.8e\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
	    }
	return(pclose(pfp)==0 ? 0 : 1);
}
#endif
