#ifndef lint
static const char RCSid[] = "$Id: pabopto2xml.c,v 2.9 2012/09/06 00:07:43 greg Exp $";
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
} RBFLIST;			/* RBF representation of DSF @ 1 incidence */

				/* our loaded grid for this incident angle */
static double	theta_in_deg, phi_in_deg;
static GRIDVAL	dsf_grid[GRIDRES][GRIDRES];

				/* processed incident DSF measurements */
static RBFLIST		*dsf_list = NULL;

				/* RBF-linking matrices (edges) */
static MIGRATION	*mig_list = NULL;

#define mtx_nrows(m)	((m)->rbfv[0]->nrbf)
#define mtx_ncols(m)	((m)->rbfv[1]->nrbf)
#define mtx_ndx(m,i,j)	((i)*mtx_ncols(m) + (j))
#define is_src(rbf,m)	((rbf) == (m)->rbfv[0])
#define is_dest(rbf,m)	((rbf) == (m)->rbfv[1])
#define nextedge(rbf,m)	(m)->enxt[is_dest(rbf,m)]

/* Compute volume associated with Gaussian lobe */
static double
rbf_volume(const RBFVAL *rbfp)
{
	double	rad = R2ANG(rbfp->crad);

	return((2.*M_PI) * rbfp->peak * rad*rad);
}

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
	int	niter = 16;
	double	lastVar, thisVar = 100.;
	int	nn;
	RBFLIST	*newnode;
	int	i, j;
	
	nn = 0;			/* count selected bins */
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		nn += dsf_grid[i][j].nval;
				/* allocate RBF array */
	newnode = (RBFLIST *)malloc(sizeof(RBFLIST) + sizeof(RBFVAL)*(nn-1));
	if (newnode == NULL) {
		fputs("Out of memory in make_rbfrep()\n", stderr);
		exit(1);
	}
	newnode->next = NULL;
	newnode->ejl = NULL;
	newnode->invec[2] = sin(M_PI/180.*theta_in_deg);
	newnode->invec[0] = cos(M_PI/180.*phi_in_deg)*newnode->invec[2];
	newnode->invec[1] = sin(M_PI/180.*phi_in_deg)*newnode->invec[2];
	newnode->invec[2] = sqrt(1. - newnode->invec[2]*newnode->invec[2]);
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
				vec_from_pos(odir, i, j);
				newnode->rbfa[nn++].peak *= corr =
					dsf_grid[i][j].vsum /
						eval_rbfrep(newnode, odir);
				dsum += corr - 1.;
				dsum2 += (corr-1.)*(corr-1.);
			}
		lastVar = thisVar;
		thisVar = dsum2/(double)nn;
		/*
		fprintf(stderr, "Avg., RMS error: %.1f%%  %.1f%%\n",
					100.*dsum/(double)nn,
					100.*sqrt(thisVar));
		*/
	} while (--niter > 0 && lastVar-thisVar > 0.02*lastVar);

	nn = 0;			/* compute sum for normalization */
	while (nn < newnode->nrbf)
		newnode->vtotal += rbf_volume(&newnode->rbfa[nn++]);

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
price_routes(const RBFLIST *from_rbf, const RBFLIST *to_rbf)
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
		vec_from_pos(vto[j], to_rbf->rbfa[j].gx, to_rbf->rbfa[j].gy);

	for (i = from_rbf->nrbf; i--; ) {
	    const double	from_ang = R2ANG(from_rbf->rbfa[i].crad);
	    FVECT		vfrom;
	    vec_from_pos(vfrom, from_rbf->rbfa[i].gx, from_rbf->rbfa[i].gy);
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
if (amt2move > 1e-5) fprintf(stderr, "%g leftover!\n", amt2move);
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
make_migration(RBFLIST *from_rbf, RBFLIST *to_rbf)
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

#if 0
/* Partially advect between the given RBFs to a newly allocated one */
static RBFLIST *
advect_rbf(const RBFLIST *from_rbf, const RBFLIST *to_rbf,
			const float *mtx, const FVECT invec)
{
	RBFLIST		*rbf;

	if (from_rbf->nrbf > to_rbf->nrbf) {
		fputs("Internal error: source RBF won't fit destination\n",
				stderr);
		exit(1);
	}
	rbf = (RBFLIST *)malloc(sizeof(RBFLIST) + sizeof(RBFVAL)*(to_rbf->nrbf-1));
	if (rbf == NULL) {
		fputs("Out of memory in advect_rbf()\n", stderr);
		exit(1);
	}
	rbf->next = NULL; rbf->ejl = NULL;
	VCOPY(rbf->invec, invec);
	rbf->vtotal = 0;
	rbf->nrbf = to_rbf->nrbf;
	
	return rbf;
}
#endif

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
