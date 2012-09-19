#ifndef lint
static const char RCSid[] = "$Id: pabopto2xml.c,v 2.11 2012/09/19 22:03:37 greg Exp $";
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

#define DEBUG		1

#ifndef GRIDRES
#define GRIDRES		200		/* grid resolution per side */
#endif

#define RSCA		2.7		/* radius scaling factor (empirical) */

					/* convert to/from coded radians */
#define ANG2R(r)	(int)((r)*((1<<16)/M_PI))
#define R2ANG(c)	(((c)+.5)*(M_PI/(1<<16)))

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

struct s_rbfnode;		/* forward declaration of RBF struct */

typedef struct s_migration {
	struct s_migration	*next;		/* next in global edge list */
	struct s_rbfnode	*rbfv[2];	/* from,to vertex */
	struct s_migration	*enxt[2];	/* next from,to sibling */
	float			mtx[1];		/* matrix (extends struct) */
} MIGRATION;			/* migration link (winged edge structure) */

typedef struct s_rbfnode {
	struct s_rbfnode	*next;		/* next in global RBF list */
	MIGRATION		*ejl;		/* edge list for this vertex */
	FVECT			invec;		/* incident vector direction */
	double			vtotal;		/* volume for normalization */
	int			nrbf;		/* number of RBFs */
	RBFVAL			rbfa[1];	/* RBF array (extends struct) */
} RBFNODE;			/* RBF representation of DSF @ 1 incidence */

				/* our loaded grid for this incident angle */
static double		theta_in_deg, phi_in_deg;
static GRIDVAL		dsf_grid[GRIDRES][GRIDRES];

				/* all incident angles in-plane so far? */
static int		single_plane_incident = -1;

				/* input/output orientations */
static int		input_orient = 0;
static int		output_orient = 0;

				/* processed incident DSF measurements */
static RBFNODE		*dsf_list = NULL;

				/* RBF-linking matrices (edges) */
static MIGRATION	*mig_list = NULL;

				/* migration edges drawn in raster fashion */
static MIGRATION	*mig_grid[GRIDRES][GRIDRES];

#define mtx_nrows(m)	((m)->rbfv[0]->nrbf)
#define mtx_ncols(m)	((m)->rbfv[1]->nrbf)
#define mtx_ndx(m,i,j)	((i)*mtx_ncols(m) + (j))
#define is_src(rbf,m)	((rbf) == (m)->rbfv[0])
#define is_dest(rbf,m)	((rbf) == (m)->rbfv[1])
#define nextedge(rbf,m)	(m)->enxt[is_dest(rbf,m)]
#define opp_rbf(rbf,m)	(m)->rbfv[is_src(rbf,m)]

#define	round(v)	(int)((v) + .5 - ((v) < -.5))

/* Compute volume associated with Gaussian lobe */
static double
rbf_volume(const RBFVAL *rbfp)
{
	double	rad = R2ANG(rbfp->crad);

	return((2.*M_PI) * rbfp->peak * rad*rad);
}

/* Compute outgoing vector from grid position */
static void
ovec_from_pos(FVECT vec, int xpos, int ypos)
{
	double	uv[2];
	double	r2;
	
	SDsquare2disk(uv, (1./GRIDRES)*(xpos+.5), (1./GRIDRES)*(ypos+.5));
				/* uniform hemispherical projection */
	r2 = uv[0]*uv[0] + uv[1]*uv[1];
	vec[0] = vec[1] = sqrt(2. - r2);
	vec[0] *= uv[0];
	vec[1] *= uv[1];
	vec[2] = output_orient*(1. - r2);
}

/* Compute grid position from normalized input/output vector */
static void
pos_from_vec(int pos[2], const FVECT vec)
{
	double	sq[2];		/* uniform hemispherical projection */
	double	norm = 1./sqrt(1. + fabs(vec[2]));

	SDdisk2square(sq, vec[0]*norm, vec[1]*norm);

	pos[0] = (int)(sq[0]*GRIDRES);
	pos[1] = (int)(sq[1]*GRIDRES);
}

/* Evaluate RBF for DSF at the given normalized outgoing direction */
static double
eval_rbfrep(const RBFNODE *rp, const FVECT outvec)
{
	double		res = .0;
	const RBFVAL	*rbfp;
	FVECT		odir;
	double		sig2;
	int		n;

	rbfp = rp->rbfa;
	for (n = rp->nrbf; n--; rbfp++) {
		ovec_from_pos(odir, rbfp->gx, rbfp->gy);
		sig2 = R2ANG(rbfp->crad);
		sig2 = (DOT(odir,outvec) - 1.) / (sig2*sig2);
		if (sig2 > -19.)
			res += rbfp->peak * exp(sig2);
	}
	return(res);
}

/* Insert a new directional scattering function in our global list */
static void
insert_dsf(RBFNODE *newrbf)
{
	RBFNODE		*rbf, *rbf_last;

					/* keep in ascending theta order */
	for (rbf_last = NULL, rbf = dsf_list;
			single_plane_incident & (rbf != NULL);
					rbf_last = rbf, rbf = rbf->next)
		if (input_orient*rbf->invec[2] < input_orient*newrbf->invec[2])
			break;
	if (rbf_last == NULL) {
		newrbf->next = dsf_list;
		dsf_list = newrbf;
		return;
	}
	newrbf->next = rbf;
	rbf_last->next = newrbf;
}

/* Count up filled nodes and build RBF representation from current grid */ 
static RBFNODE *
make_rbfrep(void)
{
	int	niter = 16;
	double	lastVar, thisVar = 100.;
	int	nn;
	RBFNODE	*newnode;
	int	i, j;
	
	nn = 0;			/* count selected bins */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		nn += dsf_grid[i][j].nval;
				/* allocate RBF array */
	newnode = (RBFNODE *)malloc(sizeof(RBFNODE) + sizeof(RBFVAL)*(nn-1));
	if (newnode == NULL) {
		fputs("Out of memory in make_rbfrep()\n", stderr);
		exit(1);
	}
	newnode->next = NULL;
	newnode->ejl = NULL;
	newnode->invec[2] = sin(M_PI/180.*theta_in_deg);
	newnode->invec[0] = cos(M_PI/180.*phi_in_deg)*newnode->invec[2];
	newnode->invec[1] = sin(M_PI/180.*phi_in_deg)*newnode->invec[2];
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
	do {
		double	dsum = .0, dsum2 = .0;
		nn = 0;
		for (i = 0; i < GRIDRES; i++)
		    for (j = 0; j < GRIDRES; j++)
			if (dsf_grid[i][j].nval) {
				FVECT	odir;
				double	corr;
				ovec_from_pos(odir, i, j);
				newnode->rbfa[nn++].peak *= corr =
					dsf_grid[i][j].vsum /
						eval_rbfrep(newnode, odir);
				dsum += corr - 1.;
				dsum2 += (corr-1.)*(corr-1.);
			}
		lastVar = thisVar;
		thisVar = dsum2/(double)nn;
#ifdef DEBUG
		fprintf(stderr, "Avg., RMS error: %.1f%%  %.1f%%\n",
					100.*dsum/(double)nn,
					100.*sqrt(thisVar));
#endif
	} while (--niter > 0 && lastVar-thisVar > 0.02*lastVar);

	nn = 0;			/* compute sum for normalization */
	while (nn < newnode->nrbf)
		newnode->vtotal += rbf_volume(&newnode->rbfa[nn++]);

	insert_dsf(newnode);
	return(newnode);
}

/* Load a set of measurements corresponding to a particular incident angle */
static int
load_pabopto_meas(const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	int	inp_is_DSF = -1;
	double	new_phi, theta_out, phi_out, val;
	char	buf[2048];
	int	n, c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
	memset(dsf_grid, 0, sizeof(dsf_grid));
#ifdef DEBUG
	fprintf(stderr, "Loading measurement file '%s'...\n", fname);
#endif
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
		if (sscanf(buf, "inphi %lf", &new_phi) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&theta_in_deg, &new_phi) == 2)
			continue;
	}
	if (inp_is_DSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		fclose(fp);
		return(0);
	}
	if (!input_orient)		/* check input orientation */
		input_orient = 1 - 2*(theta_in_deg > 90.);
	else if (input_orient > 0 ^ theta_in_deg < 90.) {
		fputs("Cannot handle input angles on both sides of surface\n",
				stderr);
		exit(1);
	}
	if (single_plane_incident > 0)	/* check if still in plane */
		single_plane_incident = (round(new_phi) == round(phi_in_deg));
	else if (single_plane_incident < 0)
		single_plane_incident = 1;
	phi_in_deg = new_phi;
	ungetc(c, fp);			/* read actual data */
	while (fscanf(fp, "%lf %lf %lf\n", &theta_out, &phi_out, &val) == 3) {
		FVECT	ovec;
		int	pos[2];

		if (!output_orient)	/* check output orientation */
			output_orient = 1 - 2*(theta_out > 90.);
		else if (output_orient > 0 ^ theta_out < 90.) {
			fputs("Cannot handle output angles on both sides of surface\n",
					stderr);
			exit(1);
		}
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
			fputs("Could not find non-empty neighbor!\n", stderr);
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

/* Compute (and allocate) migration price matrix for optimization */
static float *
price_routes(const RBFNODE *from_rbf, const RBFNODE *to_rbf)
{
	float	*pmtx = (float *)malloc(sizeof(float) *
					from_rbf->nrbf * to_rbf->nrbf);
	FVECT	*vto = (FVECT *)malloc(sizeof(FVECT) * to_rbf->nrbf);
	int	i, j;

	if ((pmtx == NULL) | (vto == NULL)) {
		fputs("Out of memory in migration_costs()\n", stderr);
		exit(1);
	}
	for (j = to_rbf->nrbf; j--; )		/* save repetitive ops. */
		ovec_from_pos(vto[j], to_rbf->rbfa[j].gx, to_rbf->rbfa[j].gy);

	for (i = from_rbf->nrbf; i--; ) {
	    const double	from_ang = R2ANG(from_rbf->rbfa[i].crad);
	    FVECT		vfrom;
	    ovec_from_pos(vfrom, from_rbf->rbfa[i].gx, from_rbf->rbfa[i].gy);
	    for (j = to_rbf->nrbf; j--; )
		pmtx[i*to_rbf->nrbf + j] = acos(DOT(vfrom, vto[j])) +
				fabs(R2ANG(to_rbf->rbfa[j].crad) - from_ang);
	}
	free(vto);
	return(pmtx);
}

/* Comparison routine needed for sorting price row */
static const float	*price_arr;
static int
msrt_cmp(const void *p1, const void *p2)
{
	float	c1 = price_arr[*(const int *)p1];
	float	c2 = price_arr[*(const int *)p2];

	if (c1 > c2) return(1);
	if (c1 < c2) return(-1);
	return(0);
}

/* Compute minimum (optimistic) cost for moving the given source material */
static double
min_cost(double amt2move, const double *avail, const float *price, int n)
{
	static int	*price_sort = NULL;
	static int	n_alloc = 0;
	double		total_cost = 0;
	int		i;

	if (amt2move <= FTINY)			/* pre-emptive check */
		return(0.);
	if (n > n_alloc) {			/* (re)allocate sort array */
		if (n_alloc) free(price_sort);
		price_sort = (int *)malloc(sizeof(int)*n);
		if (price_sort == NULL) {
			fputs("Out of memory in min_cost()\n", stderr);
			exit(1);
		}
		n_alloc = n;
	}
	for (i = n; i--; )
		price_sort[i] = i;
	price_arr = price;
	qsort(price_sort, n, sizeof(int), &msrt_cmp);
						/* move cheapest first */
	for (i = 0; i < n && amt2move > FTINY; i++) {
		int	d = price_sort[i];
		double	amt = (amt2move < avail[d]) ? amt2move : avail[d];

		total_cost += amt * price[d];
		amt2move -= amt;
	}
	return(total_cost);
}

/* Take a step in migration by choosing optimal bucket to transfer */
static double
migration_step(MIGRATION *mig, double *src_rem, double *dst_rem, const float *pmtx)
{
	static double	*src_cost = NULL;
	int		n_alloc = 0;
	const double	maxamt = 2./(mtx_nrows(mig)*mtx_ncols(mig));
	double		amt = 0;
	struct {
		int	s, d;	/* source and destination */
		double	price;	/* price estimate per amount moved */
		double	amt;	/* amount we can move */
	} cur, best;
	int		i;

	if (mtx_nrows(mig) > n_alloc) {		/* allocate cost array */
		if (n_alloc)
			free(src_cost);
		src_cost = (double *)malloc(sizeof(double)*mtx_nrows(mig));
		if (src_cost == NULL) {
			fputs("Out of memory in migration_step()\n", stderr);
			exit(1);
		}
		n_alloc = mtx_nrows(mig);
	}
	for (i = mtx_nrows(mig); i--; )		/* starting costs for diff. */
		src_cost[i] = min_cost(src_rem[i], dst_rem,
					pmtx+i*mtx_ncols(mig), mtx_ncols(mig));

						/* find best source & dest. */
	best.s = best.d = -1; best.price = FHUGE; best.amt = 0;
	for (cur.s = mtx_nrows(mig); cur.s--; ) {
	    const float	*price = pmtx + cur.s*mtx_ncols(mig);
	    double	cost_others = 0;
	    if (src_rem[cur.s] <= FTINY)
		    continue;
	    cur.d = -1;				/* examine cheapest dest. */
	    for (i = mtx_ncols(mig); i--; )
		if (dst_rem[i] > FTINY &&
				(cur.d < 0 || price[i] < price[cur.d]))
			cur.d = i;
	    if (cur.d < 0)
		    return(.0);
	    if ((cur.price = price[cur.d]) >= best.price)
		    continue;			/* no point checking further */
	    cur.amt = (src_rem[cur.s] < dst_rem[cur.d]) ?
				src_rem[cur.s] : dst_rem[cur.d];
	    if (cur.amt > maxamt) cur.amt = maxamt;
	    dst_rem[cur.d] -= cur.amt;		/* add up differential costs */
	    for (i = mtx_nrows(mig); i--; ) {
		if (i == cur.s) continue;
		cost_others += min_cost(src_rem[i], dst_rem, price, mtx_ncols(mig))
					- src_cost[i];
	    }
	    dst_rem[cur.d] += cur.amt;		/* undo trial move */
	    cur.price += cost_others/cur.amt;	/* adjust effective price */
	    if (cur.price < best.price)		/* are we better than best? */
		    best = cur;
	}
	if ((best.s < 0) | (best.d < 0))
		return(.0);
						/* make the actual move */
	mig->mtx[mtx_ndx(mig,best.s,best.d)] += best.amt;
	src_rem[best.s] -= best.amt;
	dst_rem[best.d] -= best.amt;
	return(best.amt);
}

/* Compute (and insert) migration along directed edge */
static MIGRATION *
make_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	const double	end_thresh = 0.02/(from_rbf->nrbf*to_rbf->nrbf);
	float		*pmtx = price_routes(from_rbf, to_rbf);
	MIGRATION	*newmig = (MIGRATION *)malloc(sizeof(MIGRATION) +
							sizeof(float) *
					(from_rbf->nrbf*to_rbf->nrbf - 1));
	double		*src_rem = (double *)malloc(sizeof(double)*from_rbf->nrbf);
	double		*dst_rem = (double *)malloc(sizeof(double)*to_rbf->nrbf);
	double		total_rem = 1.;
	int		i;

	if ((newmig == NULL) | (src_rem == NULL) | (dst_rem == NULL)) {
		fputs("Out of memory in make_migration()\n", stderr);
		exit(1);
	}
#ifdef DEBUG
	{
		double	theta, phi;
		theta = acos(from_rbf->invec[2])*(180./M_PI);
		phi = atan2(from_rbf->invec[1],from_rbf->invec[0])*(180./M_PI);
		fprintf(stderr, "Building path from (theta,phi) (%d,%d) to ",
				round(theta), round(phi));
		theta = acos(to_rbf->invec[2])*(180./M_PI);
		phi = atan2(to_rbf->invec[1],to_rbf->invec[0])*(180./M_PI);
		fprintf(stderr, "(%d,%d)\n", round(theta), round(phi));
	}
#endif
	newmig->next = NULL;
	newmig->rbfv[0] = from_rbf;
	newmig->rbfv[1] = to_rbf;
	newmig->enxt[0] = newmig->enxt[1] = NULL;
	memset(newmig->mtx, 0, sizeof(float)*from_rbf->nrbf*to_rbf->nrbf);
						/* starting quantities */
	for (i = from_rbf->nrbf; i--; )
		src_rem[i] = rbf_volume(&from_rbf->rbfa[i]) / from_rbf->vtotal;
	for (i = to_rbf->nrbf; i--; )
		dst_rem[i] = rbf_volume(&to_rbf->rbfa[i]) / to_rbf->vtotal;
						/* move a bit at a time */
	while (total_rem > end_thresh)
		total_rem -= migration_step(newmig, src_rem, dst_rem, pmtx);

	free(pmtx);				/* free working arrays */
	free(src_rem);
	free(dst_rem);
	for (i = from_rbf->nrbf; i--; ) {	/* normalize final matrix */
	    float	nf = rbf_volume(&from_rbf->rbfa[i]);
	    int		j;
	    if (nf <= FTINY) continue;
	    nf = from_rbf->vtotal / nf;
	    for (j = to_rbf->nrbf; j--; )
		newmig->mtx[mtx_ndx(newmig,i,j)] *= nf;
	}
						/* insert in edge lists */
	newmig->enxt[0] = from_rbf->ejl;
	from_rbf->ejl = newmig;
	newmig->enxt[1] = to_rbf->ejl;
	to_rbf->ejl = newmig;
	newmig->next = mig_list;
	return(mig_list = newmig);
}

/* Get triangle surface orientation (unnormalized) */
static void
tri_orient(FVECT vres, const FVECT v1, const FVECT v2, const FVECT v3)
{
	FVECT	v2minus1, v3minus2;

	VSUB(v2minus1, v2, v1);
	VSUB(v3minus2, v3, v2);
	VCROSS(vres, v2minus1, v3minus2);
}

/* Determine if vertex order is reversed (inward normal) */
static int
is_rev_tri(const FVECT v1, const FVECT v2, const FVECT v3)
{
	FVECT	tor;

	tri_orient(tor, v1, v2, v3);

	return(DOT(tor, v2) < 0.);
}

/* Find vertices completing triangles on either side of the given edge */
static int
get_triangles(RBFNODE *rbfv[2], const MIGRATION *mig)
{
	const MIGRATION	*ej, *ej2;
	RBFNODE		*tv;

	rbfv[0] = rbfv[1] = NULL;
	for (ej = mig->rbfv[0]->ejl; ej != NULL;
				ej = nextedge(mig->rbfv[0],ej)) {
		if (ej == mig)
			continue;
		tv = opp_rbf(mig->rbfv[0],ej);
		for (ej2 = tv->ejl; ej2 != NULL; ej2 = nextedge(tv,ej2))
			if (opp_rbf(tv,ej2) == mig->rbfv[1]) {
				rbfv[is_rev_tri(mig->rbfv[0]->invec,
						mig->rbfv[1]->invec,
						tv->invec)] = tv;
				break;
			}
	}
	return((rbfv[0] != NULL) + (rbfv[1] != NULL));
}

/* Find context hull vertex to complete triangle (oriented call) */
static RBFNODE *
find_chull_vert(const RBFNODE *rbf0, const RBFNODE *rbf1)
{
	FVECT	vmid, vor;
	RBFNODE	*rbf, *rbfbest = NULL;
	double	dprod2, bestdprod2 = 0.5;

	VADD(vmid, rbf0->invec, rbf1->invec);
	if (normalize(vmid) == 0)
		return(NULL);
						/* XXX exhaustive search */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		if ((rbf == rbf0) | (rbf == rbf1))
			continue;
		tri_orient(vor, rbf0->invec, rbf1->invec, rbf->invec);
		dprod2 = DOT(vor, vmid);
		if (dprod2 <= FTINY)
			continue;		/* wrong orientation */
		dprod2 *= dprod2 / DOT(vor,vor);
		if (dprod2 > bestdprod2) {	/* more convex than prev? */
			rbfbest = rbf;
			bestdprod2 = dprod2;
		}
	}
	return(rbf);
}

/* Create new migration edge and grow mesh recursively around it */
static void
mesh_from_edge(RBFNODE *rbf0, RBFNODE *rbf1)
{
	MIGRATION	*newej;
	RBFNODE		*tvert[2];

	if (rbf0 < rbf1)			/* avoid migration loops */
		newej = make_migration(rbf0, rbf1);
	else
		newej = make_migration(rbf1, rbf0);
						/* triangle on either side? */
	get_triangles(tvert, newej);
	if (tvert[0] == NULL) {			/* recurse on new right edge */
		tvert[0] = find_chull_vert(newej->rbfv[0], newej->rbfv[1]);
		if (tvert[0] != NULL) {
			mesh_from_edge(rbf0, tvert[0]);
			mesh_from_edge(rbf1, tvert[0]);
		}
	}
	if (tvert[1] == NULL) {			/* recurse on new left edge */
		tvert[1] = find_chull_vert(newej->rbfv[1], newej->rbfv[0]);
		if (tvert[1] != NULL) {
			mesh_from_edge(rbf0, tvert[1]);
			mesh_from_edge(rbf1, tvert[1]);
		}
	}
}

/* Draw edge list into mig_grid array */
static void
draw_edges()
{
	int		nnull = 0, ntot = 0;
	MIGRATION	*ej;
	int		p0[2], p1[2];

	/* memset(mig_grid, 0, sizeof(mig_grid)); */
	for (ej = mig_list; ej != NULL; ej = ej->next) {
		++ntot;
		pos_from_vec(p0, ej->rbfv[0]->invec);
		pos_from_vec(p1, ej->rbfv[1]->invec);
		if ((p0[0] == p1[0]) & (p0[1] == p1[1])) {
			++nnull;
			mig_grid[p0[0]][p0[1]] = ej;
			continue;
		}
		if (abs(p1[0]-p0[0]) > abs(p1[1]-p0[1])) {
			const int	xstep = 2*(p1[0] > p0[0]) - 1;
			const double	ystep = (double)((p1[1]-p0[1])*xstep) /
							(double)(p1[0]-p0[0]);
			int		x;
			double		y;
			for (x = p0[0], y = p0[1]+.5; x != p1[0];
						x += xstep, y += ystep)
				mig_grid[x][(int)y] = ej;
			mig_grid[x][(int)y] = ej;
		} else {
			const int	ystep = 2*(p1[1] > p0[1]) - 1;
			const double	xstep = (double)((p1[0]-p0[0])*ystep) /
							(double)(p1[1]-p0[1]);
			int		y;
			double		x;
			for (y = p0[1], x = p0[0]+.5; y != p1[1];
						y += ystep, x += xstep)
				mig_grid[(int)x][y] = ej;
			mig_grid[(int)x][y] = ej;
		}
	}
	if (nnull)
		fprintf(stderr, "Warning: %d of %d edges are null\n",
				nnull, ntot);
}
	
/* Build our triangle mesh from recorded RBFs */
static void
build_mesh()
{
	double		best2 = M_PI*M_PI;
	RBFNODE		*rbf, *rbf_near = NULL;
						/* check if isotropic */
	if (single_plane_incident) {
		for (rbf = dsf_list; rbf != NULL; rbf = rbf->next)
			if (rbf->next != NULL)
				make_migration(rbf, rbf->next);
		return;
	}
						/* find RBF nearest to head */
	if (dsf_list == NULL)
		return;
	for (rbf = dsf_list->next; rbf != NULL; rbf = rbf->next) {
		double	dist2 = 2. - 2.*DOT(dsf_list->invec,rbf->invec);
		if (dist2 < best2) {
			rbf_near = rbf;
			best2 = dist2;
		}
	}
	if (rbf_near == NULL) {
		fputs("Cannot find nearest point for first edge\n", stderr);
		exit(1);
	}
						/* build mesh from this edge */
	mesh_from_edge(dsf_list, rbf_near);
						/* draw edge list into grid */
	draw_edges();
}

/* Identify enclosing triangle for this position (flood fill raster check) */
static int
identify_tri(MIGRATION *miga[3], unsigned char vmap[GRIDRES][(GRIDRES+7)/8],
			int px, int py)
{
	const int	btest = 1<<(py&07);

	if (vmap[px][py>>3] & btest)		/* already visited here? */
		return(1);
						/* else mark it */
	vmap[px][py>>3] |= btest;

	if (mig_grid[px][py] != NULL) {		/* are we on an edge? */
		int	i;
		for (i = 0; i < 3; i++) {
			if (miga[i] == mig_grid[px][py])
				return(1);
			if (miga[i] != NULL)
				continue;
			miga[i] = mig_grid[px][py];
			return(1);
		}
		return(0);			/* outside triangle! */
	}
						/* check neighbors (flood) */
	if (px > 0 && !identify_tri(miga, vmap, px-1, py))
		return(0);
	if (px < GRIDRES-1 && !identify_tri(miga, vmap, px+1, py))
		return(0);
	if (py > 0 && !identify_tri(miga, vmap, px, py-1))
		return(0);
	if (py < GRIDRES-1 && !identify_tri(miga, vmap, px, py+1))
		return(0);
	return(1);				/* this neighborhood done */
}

/* Find edge(s) for interpolating the given incident vector */
static int
get_interp(MIGRATION *miga[3], const FVECT invec)
{
	miga[0] = miga[1] = miga[2] = NULL;
	if (single_plane_incident) {		/* isotropic BSDF? */
		RBFNODE	*rbf;			/* find edge we're on */
		for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
			if (input_orient*rbf->invec[2] < input_orient*invec[2])
				break;
			if (rbf->next != NULL &&
					input_orient*rbf->next->invec[2] <
							input_orient*invec[2]) {
				for (miga[0] = rbf->ejl; miga[0] != NULL;
						miga[0] = nextedge(rbf,miga[0]))
					if (opp_rbf(rbf,miga[0]) == rbf->next)
						return(1);
				break;
			}
		}
		return(0);			/* outside range! */
	}
	{					/* else use triagnle mesh */
		unsigned char	floodmap[GRIDRES][(GRIDRES+7)/8];
		int		pstart[2];

		pos_from_vec(pstart, invec);
		memset(floodmap, 0, sizeof(floodmap));
						/* call flooding function */
		if (!identify_tri(miga, floodmap, pstart[0], pstart[1]))
			return(0);		/* outside mesh */
		if ((miga[0] == NULL) | (miga[2] == NULL))
			return(0);		/* should never happen */
		if (miga[1] == NULL)
			return(1);		/* on edge */
		return(3);			/* else in triangle */
	}
}

/* Advect and allocate new RBF along edge */
static RBFNODE *
e_advect_rbf(const MIGRATION *mig, const FVECT invec)
{
	RBFNODE		*rbf;
	int		n, i, j;
	double		t, full_dist;
						/* get relative position */
	t = acos(DOT(invec, mig->rbfv[0]->invec));
	if (t < M_PI/GRIDRES) {			/* near first DSF */
		n = sizeof(RBFNODE) + sizeof(RBFVAL)*(mig->rbfv[0]->nrbf-1);
		rbf = (RBFNODE *)malloc(n);
		if (rbf == NULL)
			goto memerr;
		memcpy(rbf, mig->rbfv[0], n);	/* just duplicate */
		return(rbf);
	}
	full_dist = acos(DOT(mig->rbfv[0]->invec, mig->rbfv[1]->invec));
	if (t > full_dist-M_PI/GRIDRES) {	/* near second DSF */
		n = sizeof(RBFNODE) + sizeof(RBFVAL)*(mig->rbfv[1]->nrbf-1);
		rbf = (RBFNODE *)malloc(n);
		if (rbf == NULL)
			goto memerr;
		memcpy(rbf, mig->rbfv[1], n);	/* just duplicate */
		return(rbf);
	}
	t /= full_dist;	
	n = 0;					/* count migrating particles */
	for (i = 0; i < mtx_nrows(mig); i++)
	    for (j = 0; j < mtx_ncols(mig); j++)
		n += (mig->mtx[mtx_ndx(mig,i,j)] > FTINY);
	
	rbf = (RBFNODE *)malloc(sizeof(RBFNODE) + sizeof(RBFVAL)*(n-1));
	if (rbf == NULL)
		goto memerr;
	rbf->next = NULL; rbf->ejl = NULL;
	VCOPY(rbf->invec, invec);
	rbf->nrbf = n;
	rbf->vtotal = 1.-t + t*mig->rbfv[1]->vtotal/mig->rbfv[0]->vtotal;
	n = 0;					/* advect RBF lobes */
	for (i = 0; i < mtx_nrows(mig); i++) {
	    const RBFVAL	*rbf0i = &mig->rbfv[0]->rbfa[i];
	    const float		peak0 = rbf0i->peak;
	    const double	rad0 = R2ANG(rbf0i->crad);
	    FVECT		v0;
	    float		mv;
	    ovec_from_pos(v0, rbf0i->gx, rbf0i->gy);
	    for (j = 0; j < mtx_ncols(mig); j++)
		if ((mv = mig->mtx[mtx_ndx(mig,i,j)]) > FTINY) {
			const RBFVAL	*rbf1j = &mig->rbfv[1]->rbfa[j];
			double		rad1 = R2ANG(rbf1j->crad);
			FVECT		v;
			int		pos[2];
			rbf->rbfa[n].peak = peak0 * mv * rbf->vtotal;
			rbf->rbfa[n].crad = ANG2R(sqrt(rad0*rad0*(1.-t) +
							rad1*rad1*t));
			ovec_from_pos(v, rbf1j->gx, rbf1j->gy);
			geodesic(v, v0, v, t, GEOD_REL);
			pos_from_vec(pos, v);
			rbf->rbfa[n].gx = pos[0];
			rbf->rbfa[n].gy = pos[1];
			++n;
		}
	}
	rbf->vtotal *= mig->rbfv[0]->vtotal;	/* turn ratio into actual */
	return(rbf);
memerr:
	fputs("Out of memory in e_advect_rbf()\n", stderr);
	exit(1);
	return(NULL);	/* pro forma return */
}

/* Insert vertex in ordered list */
static void
insert_vert(RBFNODE **vlist, RBFNODE *v)
{
	int	i, j;

	for (i = 0; vlist[i] != NULL; i++) {
		if (v == vlist[i])
			return;
		if (v < vlist[i])
			break;
	}
	for (j = i; vlist[j] != NULL; j++)
		;
	while (j > i) {
		vlist[j] = vlist[j-1];
		--j;
	}
	vlist[i] = v;
}

/* Sort triangle edges in standard order */
static void
order_triangle(MIGRATION *miga[3])
{
	RBFNODE		*vert[4];
	MIGRATION	*ord[3];
	int		i;
						/* order vertices, first */
	memset(vert, 0, sizeof(vert));
	for (i = 0; i < 3; i++) {
		insert_vert(vert, miga[i]->rbfv[0]);
		insert_vert(vert, miga[i]->rbfv[1]);
	}
						/* identify edge 0 */
	for (i = 0; i < 3; i++)
		if (miga[i]->rbfv[0] == vert[0] &&
				miga[i]->rbfv[1] == vert[1]) {
			ord[0] = miga[i];
			break;
		}
						/* identify edge 1 */
	for (i = 0; i < 3; i++)
		if (miga[i]->rbfv[0] == vert[1] &&
				miga[i]->rbfv[1] == vert[2]) {
			ord[1] = miga[i];
			break;
		}
						/* identify edge 2 */
	for (i = 0; i < 3; i++)
		if (miga[i]->rbfv[0] == vert[0] &&
				miga[i]->rbfv[1] == vert[2]) {
			ord[2] = miga[i];
			break;
		}
	miga[0] = ord[0]; miga[1] = ord[1]; miga[2] = ord[2];
}

/* Partially advect between recorded incident angles and allocate new RBF */
static RBFNODE *
advect_rbf(const FVECT invec)
{
	MIGRATION	*miga[3];
	RBFNODE		*rbf;
	float		mbfact, mcfact;
	int		n, i, j, k;
	FVECT		v0, v1, v2;
	double		s, t;

	if (!get_interp(miga, invec))		/* can't interpolate? */
		return(NULL);
	if (miga[1] == NULL)			/* along edge? */
		return(e_advect_rbf(miga[0], invec));
						/* put in standard order */
	order_triangle(miga);
						/* figure out position */
	fcross(v0, miga[2]->rbfv[0]->invec, miga[2]->rbfv[1]->invec);
	normalize(v0);
	fcross(v2, miga[1]->rbfv[0]->invec, miga[1]->rbfv[1]->invec);
	normalize(v2);
	fcross(v1, invec, miga[1]->rbfv[1]->invec);
	normalize(v1);
	s = acos(DOT(v0,v1)) / acos(DOT(v0,v2));
	geodesic(v1, miga[0]->rbfv[0]->invec, miga[0]->rbfv[1]->invec,
			s, GEOD_REL);
	t = acos(DOT(v1,invec)) / acos(DOT(v1,miga[1]->rbfv[1]->invec));
	n = 0;					/* count migrating particles */
	for (i = 0; i < mtx_nrows(miga[0]); i++)
	    for (j = 0; j < mtx_ncols(miga[0]); j++)
		for (k = (miga[0]->mtx[mtx_ndx(miga[0],i,j)] > FTINY) *
					mtx_ncols(miga[2]); k--; )
			n += (miga[2]->mtx[mtx_ndx(miga[2],i,k)] > FTINY &&
				miga[1]->mtx[mtx_ndx(miga[1],j,k)] > FTINY);
				
	rbf = (RBFNODE *)malloc(sizeof(RBFNODE) + sizeof(RBFVAL)*(n-1));
	if (rbf == NULL) {
		fputs("Out of memory in advect_rbf()\n", stderr);
		exit(1);
	}
	rbf->next = NULL; rbf->ejl = NULL;
	VCOPY(rbf->invec, invec);
	rbf->nrbf = n;
	n = 0;					/* compute RBF lobes */
	mbfact = s * miga[0]->rbfv[1]->vtotal/miga[0]->rbfv[0]->vtotal *
		(1.-t + t*miga[1]->rbfv[1]->vtotal/miga[1]->rbfv[0]->vtotal);
	mcfact = (1.-s) *
		(1.-t + t*miga[2]->rbfv[1]->vtotal/miga[2]->rbfv[0]->vtotal);
	for (i = 0; i < mtx_nrows(miga[0]); i++) {
	    const RBFVAL	*rbf0i = &miga[0]->rbfv[0]->rbfa[i];
	    const float		w0i = rbf0i->peak;
	    const double	rad0i = R2ANG(rbf0i->crad);
	    ovec_from_pos(v0, rbf0i->gx, rbf0i->gy);
	    for (j = 0; j < mtx_ncols(miga[0]); j++) {
		const float	ma = miga[0]->mtx[mtx_ndx(miga[0],i,j)];
		const RBFVAL	*rbf1j;
		double		rad1j, srad2;
		if (ma <= FTINY)
			continue;
		rbf1j = &miga[0]->rbfv[1]->rbfa[j];
		rad1j = R2ANG(rbf1j->crad);
		srad2 = (1.-s)*(1.-t)*rad0i*rad0i + s*(1.-t)*rad1j*rad1j;
		ovec_from_pos(v1, rbf1j->gx, rbf1j->gy);
		geodesic(v1, v0, v1, s, GEOD_REL);
		for (k = 0; k < mtx_ncols(miga[2]); k++) {
		    float		mb = miga[1]->mtx[mtx_ndx(miga[1],j,k)];
		    float		mc = miga[2]->mtx[mtx_ndx(miga[2],i,k)];
		    const RBFVAL	*rbf2k;
		    double		rad2k;
		    FVECT		vout;
		    int			pos[2];
		    if ((mb <= FTINY) | (mc <= FTINY))
			continue;
		    rbf2k = &miga[2]->rbfv[1]->rbfa[k];
		    rbf->rbfa[n].peak = w0i * ma * (mb*mbfact + mc*mcfact);
		    rad2k = R2ANG(rbf2k->crad);
		    rbf->rbfa[n].crad = RAD2R(sqrt(srad2 + t*rad2k*rad2k));
		    ovec_from_pos(v2, rbf2k->gx, rbf2k->gy);
		    geodesic(vout, v1, v2, t, GEOD_REL);
		    pos_from_vec(pos, vout);
		    rbf->rbfa[n].gx = pos[0];
		    rbf->rbfa[n].gy = pos[1];
		    ++n;
		}
	    }
	}
	rbf->vtotal = miga[0]->rbfv[0]->vtotal * (mbfact + mcfact);
	return(rbf);
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
	if (!load_pabopto_meas(argv[1]))
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
			ovec_from_pos(dir, i, j);
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
				ovec_from_pos(dir, i, j);
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
		ovec_from_pos(dir, i, j);
		bsdf = eval_rbfrep(dsf_list, dir) / dir[2];
		fprintf(pfp, "%.8e %.8e %.8e\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
	    }
	return(pclose(pfp)==0 ? 0 : 1);
}
#endif
