#ifndef lint
static const char RCSid[] = "$Id: bsdfmesh.c,v 2.21 2014/03/09 01:51:48 greg Exp $";
#endif
/*
 * Create BSDF advection mesh from radial basis functions.
 *
 *	G. Ward
 */

#ifndef _WIN32
#include <unistd.h>
#include <sys/wait.h>
#include <sys/mman.h>
#endif
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bsdfrep.h"

#ifndef NEIGH_FACT2
#define NEIGH_FACT2	0.1	/* empirical neighborhood distance weight */
#endif
				/* number of processes to run */
int			nprocs = 1;
				/* number of children (-1 in child) */
static int		nchild = 0;

typedef struct {
	int		nrows, ncols;	/* array size (matches migration) */
	float		*price;		/* migration prices */
	short		*sord;		/* sort for each row, low to high */
	float		*prow;		/* current price row */
} PRICEMAT;			/* sorted pricing matrix */

#define	pricerow(p,i)	((p)->price + (i)*(p)->ncols)
#define psortrow(p,i)	((p)->sord + (i)*(p)->ncols)

/* Create a new migration holder (sharing memory for multiprocessing) */
static MIGRATION *
new_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	size_t		memlen = sizeof(MIGRATION) +
				sizeof(float)*(from_rbf->nrbf*to_rbf->nrbf - 1);
	MIGRATION	*newmig;
#ifdef _WIN32
	if (nprocs > 1)
		fprintf(stderr, "%s: warning - multiprocessing not supported\n",
				progname);
	nprocs = 1;
	newmig = (MIGRATION *)malloc(memlen);
#else
	if (nprocs <= 1) {			/* single process? */
		newmig = (MIGRATION *)malloc(memlen);
	} else {				/* else need to share memory */
		newmig = (MIGRATION *)mmap(NULL, memlen, PROT_READ|PROT_WRITE,
						MAP_ANON|MAP_SHARED, -1, 0);
		if ((void *)newmig == MAP_FAILED)
			newmig = NULL;
	}
#endif
	if (newmig == NULL) {
		fprintf(stderr, "%s: cannot allocate new migration\n", progname);
		exit(1);
	}
	newmig->rbfv[0] = from_rbf;
	newmig->rbfv[1] = to_rbf;
						/* insert in edge lists */
	newmig->enxt[0] = from_rbf->ejl;
	from_rbf->ejl = newmig;
	newmig->enxt[1] = to_rbf->ejl;
	to_rbf->ejl = newmig;
	newmig->next = mig_list;		/* push onto global list */
	return(mig_list = newmig);
}

#ifdef _WIN32
#define await_children(n)	(void)(n)
#define run_subprocess()	0
#define end_subprocess()	(void)0
#else

/* Wait for the specified number of child processes to complete */
static void
await_children(int n)
{
	int	exit_status = 0;

	if (n > nchild)
		n = nchild;
	while (n-- > 0) {
		int	status;
		if (wait(&status) < 0) {
			fprintf(stderr, "%s: missing child(ren)!\n", progname);
			nchild = 0;
			break;
		}
		--nchild;
		if (status) {			/* something wrong */
			if ((status = WEXITSTATUS(status)))
				exit_status = status;
			else
				exit_status += !exit_status;
			fprintf(stderr, "%s: subprocess died\n", progname);
			n = nchild;		/* wait for the rest */
		}
	}
	if (exit_status)
		exit(exit_status);
}

/* Start child process if multiprocessing selected */
static pid_t
run_subprocess(void)
{
	int	status;
	pid_t	pid;

	if (nprocs <= 1)			/* any children requested? */
		return(0);
	await_children(nchild + 1 - nprocs);	/* free up child process */
	if ((pid = fork())) {
		if (pid < 0) {
			fprintf(stderr, "%s: cannot fork subprocess\n",
					progname);
			await_children(nchild);
			exit(1);
		}
		++nchild;			/* subprocess started */
		return(pid);
	}
	nchild = -1;
	return(0);				/* put child to work */
}

/* If we are in subprocess, call exit */
#define	end_subprocess()	if (nchild < 0) _exit(0); else

#endif	/* ! _WIN32 */

/* Compute normalized distribution scattering functions for comparison */
static void
compute_nDSFs(const RBFNODE *rbf0, const RBFNODE *rbf1)
{
	const double	nf0 = (GRIDRES*GRIDRES) / rbf0->vtotal;
	const double	nf1 = (GRIDRES*GRIDRES) / rbf1->vtotal;
	int		x, y;
	FVECT		dv;

	for (x = GRIDRES; x--; )
	    for (y = GRIDRES; y--; ) {
		ovec_from_pos(dv, x, y);	/* cube root (brightness) */
		dsf_grid[x][y].val[0] = pow(nf0*eval_rbfrep(rbf0, dv), .3333);
		dsf_grid[x][y].val[1] = pow(nf1*eval_rbfrep(rbf1, dv), .3333);
	    }
}	

/* Compute neighborhood distance-squared (dissimilarity) */
static double
neighborhood_dist2(int x0, int y0, int x1, int y1)
{
	int	rad = GRIDRES>>5;
	double	sum2 = 0.;
	double	d;
	int	p[4];
	int	i, j;

	if ((x0 == x1) & (y0 == y1))
		return(0.);
						/* check radius */
	p[0] = x0; p[1] = y0; p[2] = x1; p[3] = y1;
	for (i = 4; i--; ) {
		if (p[i] < rad) rad = p[i];
		if (GRIDRES-1-p[i] < rad) rad = GRIDRES-1-p[i];
	}
	for (i = -rad; i <= rad; i++)
	    for (j = -rad; j <= rad; j++) {
	    	d = dsf_grid[x0+i][y0+j].val[0] -
			dsf_grid[x1+i][y1+j].val[1];
		sum2 += d*d;
	    }
	return(sum2 / (4*rad*(rad+1) + 1));
}

/* Comparison routine needed for sorting price row */
static int
msrt_cmp(void *b, const void *p1, const void *p2)
{
	PRICEMAT	*pm = (PRICEMAT *)b;
	float		c1 = pm->prow[*(const short *)p1];
	float		c2 = pm->prow[*(const short *)p2];

	if (c1 > c2) return(1);
	if (c1 < c2) return(-1);
	return(0);
}

/* Compute (and allocate) migration price matrix for optimization */
static void
price_routes(PRICEMAT *pm, const RBFNODE *from_rbf, const RBFNODE *to_rbf)
{
	FVECT	*vto = (FVECT *)malloc(sizeof(FVECT) * to_rbf->nrbf);
	int	i, j;

	compute_nDSFs(from_rbf, to_rbf);
	pm->nrows = from_rbf->nrbf;
	pm->ncols = to_rbf->nrbf;
	pm->price = (float *)malloc(sizeof(float) * pm->nrows*pm->ncols);
	pm->sord = (short *)malloc(sizeof(short) * pm->nrows*pm->ncols);
	
	if ((pm->price == NULL) | (pm->sord == NULL) | (vto == NULL)) {
		fprintf(stderr, "%s: Out of memory in migration_costs()\n",
				progname);
		exit(1);
	}
	for (j = to_rbf->nrbf; j--; )		/* save repetitive ops. */
		ovec_from_pos(vto[j], to_rbf->rbfa[j].gx, to_rbf->rbfa[j].gy);

	for (i = from_rbf->nrbf; i--; ) {
	    const double	from_ang = R2ANG(from_rbf->rbfa[i].crad);
	    FVECT		vfrom;
	    short		*srow;
	    ovec_from_pos(vfrom, from_rbf->rbfa[i].gx, from_rbf->rbfa[i].gy);
	    pm->prow = pricerow(pm,i);
	    srow = psortrow(pm,i);
	    for (j = to_rbf->nrbf; j--; ) {
		double	d;			/* quadratic cost function */
		d = Acos(DOT(vfrom, vto[j]));
		pm->prow[j] = d*d;
		d = R2ANG(to_rbf->rbfa[j].crad) - from_ang;
		pm->prow[j] += d*d;
						/* neighborhood difference */
		pm->prow[j] += NEIGH_FACT2 * neighborhood_dist2(
				from_rbf->rbfa[i].gx, from_rbf->rbfa[i].gy,
				to_rbf->rbfa[j].gx, to_rbf->rbfa[j].gy );
		srow[j] = j;
	    }
	    qsort_r(srow, pm->ncols, sizeof(short), pm, &msrt_cmp);
	}
	free(vto);
}

/* Free price matrix */
static void
free_routes(PRICEMAT *pm)
{
	free(pm->price); pm->price = NULL;
	free(pm->sord); pm->sord = NULL;
}

/* Compute minimum (optimistic) cost for moving the given source material */
static double
min_cost(double amt2move, const double *avail, const PRICEMAT *pm, int s)
{
	const short	*srow = psortrow(pm,s);
	const float	*prow = pricerow(pm,s);
	double		total_cost = 0;
	int		j;
						/* move cheapest first */
	for (j = 0; (j < pm->ncols) & (amt2move > FTINY); j++) {
		int	d = srow[j];
		double	amt = (amt2move < avail[d]) ? amt2move : avail[d];

		total_cost += amt * prow[d];
		amt2move -= amt;
	}
	return(total_cost);
}

typedef struct {
	short	s, d;		/* source and destination */
	float	dc;		/* discount to push inventory */
} ROWSENT;		/* row sort entry */

/* Compare entries by discounted moving price */
static int
rmovcmp(void *b, const void *p1, const void *p2)
{
	PRICEMAT	*pm = (PRICEMAT *)b;
	const ROWSENT	*re1 = (const ROWSENT *)p1;
	const ROWSENT	*re2 = (const ROWSENT *)p2;
	double		price_diff;

	if (re1->d < 0) return(re2->d >= 0);
	if (re2->d < 0) return(-1);
	price_diff = re1->dc*pricerow(pm,re1->s)[re1->d] -
			re2->dc*pricerow(pm,re2->s)[re2->d];
	if (price_diff > 0) return(1);
	if (price_diff < 0) return(-1);
	return(0);
}

/* Take a step in migration by choosing reasonable bucket to transfer */
static double
migration_step(MIGRATION *mig, double *src_rem, double *dst_rem, PRICEMAT *pm)
{
	const int	max2check = 100;
	const double	maxamt = 1./(double)pm->ncols;
	const double	minamt = maxamt*1e-4;
	double		*src_cost;
	ROWSENT		*rord;
	struct {
		int	s, d;	/* source and destination */
		double	price;	/* cost per amount moved */
		double	amt;	/* amount we can move */
	} cur, best;
	int		r2check, i, ri;
	/*
	 * Check cheapest available routes only -- a higher adjusted
	 * destination price implies that another source is closer, so
	 * we can hold off considering more expensive options until
	 * some other (hopefully better) moves have been made.
	 * A discount based on source remaining is supposed to prioritize
	 * movement from large lobes, but it doesn't seem to do much,
	 * so we have it set to 1.0 at the moment.
	 */
#define discount(qr)	1.0
						/* most promising row order */
	rord = (ROWSENT *)malloc(sizeof(ROWSENT)*pm->nrows);
	if (rord == NULL)
		goto memerr;
	for (ri = pm->nrows; ri--; ) {
	    rord[ri].s = ri;
	    rord[ri].d = -1;
	    rord[ri].dc = 1.f;
	    if (src_rem[ri] <= minamt)		/* enough source material? */
		    continue;
	    for (i = 0; i < pm->ncols; i++)
		if (dst_rem[ rord[ri].d = psortrow(pm,ri)[i] ] > minamt)
			break;
	    if (i >= pm->ncols) {		/* moved all we can? */
		free(rord);
		return(.0);
	    }
	    rord[ri].dc = discount(src_rem[ri]);
	}
	if (pm->nrows > max2check)		/* sort if too many sources */
		qsort_r(rord, pm->nrows, sizeof(ROWSENT), pm, &rmovcmp);
						/* allocate cost array */
	src_cost = (double *)malloc(sizeof(double)*pm->nrows);
	if (src_cost == NULL)
		goto memerr;
	for (i = pm->nrows; i--; )		/* starting costs for diff. */
		src_cost[i] = min_cost(src_rem[i], dst_rem, pm, i);
						/* find best source & dest. */
	best.s = best.d = -1; best.price = FHUGE; best.amt = 0;
	if ((r2check = pm->nrows) > max2check)
		r2check = max2check;		/* put a limit on search */
	for (ri = 0; ri < r2check; ri++) {	/* check each source row */
	    double	cost_others = 0;
	    cur.s = rord[ri].s;
	    if ((cur.d = rord[ri].d) < 0 ||
			rord[ri].dc*pricerow(pm,cur.s)[cur.d] >= best.price) {
		if (pm->nrows > max2check) break;	/* sorted end */
		continue;			/* else skip this one */
	    }
	    cur.amt = (src_rem[cur.s] < dst_rem[cur.d]) ?
				src_rem[cur.s] : dst_rem[cur.d];
						/* don't just leave smidgen */
	    if (cur.amt > maxamt*1.02) cur.amt = maxamt;
	    dst_rem[cur.d] -= cur.amt;		/* add up opportunity costs */
	    for (i = pm->nrows; i--; )
		if (i != cur.s)
		    cost_others += min_cost(src_rem[i], dst_rem, pm, i)
					- src_cost[i];
	    dst_rem[cur.d] += cur.amt;		/* undo trial move */
						/* discount effective price */
	    cur.price = ( pricerow(pm,cur.s)[cur.d] + cost_others/cur.amt ) *
					rord[ri].dc;
	    if (cur.price < best.price)		/* are we better than best? */
		best = cur;
	}
	free(src_cost);				/* clean up */
	free(rord);
	if ((best.s < 0) | (best.d < 0))	/* nothing left to move? */
		return(.0);
						/* else make the actual move */
	mtx_coef(mig,best.s,best.d) += best.amt;
	src_rem[best.s] -= best.amt;
	dst_rem[best.d] -= best.amt;
	return(best.amt);
memerr:
	fprintf(stderr, "%s: Out of memory in migration_step()\n", progname);
	exit(1);
#undef discount
}

/* Compute and insert migration along directed edge (may fork child) */
static MIGRATION *
create_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	const double	end_thresh = 5e-6;
	PRICEMAT	pmtx;
	MIGRATION	*newmig;
	double		*src_rem, *dst_rem;
	double		total_rem = 1., move_amt;
	int		i, j;
						/* check if exists already */
	for (newmig = from_rbf->ejl; newmig != NULL;
			newmig = nextedge(from_rbf,newmig))
		if (newmig->rbfv[1] == to_rbf)
			return(NULL);
						/* else allocate */
#ifdef DEBUG
	fprintf(stderr, "Building path from (theta,phi) (%.1f,%.1f) ",
			get_theta180(from_rbf->invec),
			get_phi360(from_rbf->invec));
	fprintf(stderr, "to (%.1f,%.1f) with %d x %d matrix\n",
			get_theta180(to_rbf->invec),
			get_phi360(to_rbf->invec), 
			from_rbf->nrbf, to_rbf->nrbf);
#endif
	newmig = new_migration(from_rbf, to_rbf);
	if (run_subprocess())
		return(newmig);			/* child continues */
	price_routes(&pmtx, from_rbf, to_rbf);
	src_rem = (double *)malloc(sizeof(double)*from_rbf->nrbf);
	dst_rem = (double *)malloc(sizeof(double)*to_rbf->nrbf);
	if ((src_rem == NULL) | (dst_rem == NULL)) {
		fprintf(stderr, "%s: Out of memory in create_migration()\n",
				progname);
		exit(1);
	}
						/* starting quantities */
	memset(newmig->mtx, 0, sizeof(float)*from_rbf->nrbf*to_rbf->nrbf);
	for (i = from_rbf->nrbf; i--; )
		src_rem[i] = rbf_volume(&from_rbf->rbfa[i]) / from_rbf->vtotal;
	for (j = to_rbf->nrbf; j--; )
		dst_rem[j] = rbf_volume(&to_rbf->rbfa[j]) / to_rbf->vtotal;

	do {					/* move a bit at a time */
		move_amt = migration_step(newmig, src_rem, dst_rem, &pmtx);
		total_rem -= move_amt;
	} while ((total_rem > end_thresh) & (move_amt > 0));

	for (i = from_rbf->nrbf; i--; ) {	/* normalize final matrix */
	    double	nf = rbf_volume(&from_rbf->rbfa[i]);
	    if (nf <= FTINY) continue;
	    nf = from_rbf->vtotal / nf;
	    for (j = to_rbf->nrbf; j--; )
		mtx_coef(newmig,i,j) *= nf;	/* row now sums to 1.0 */
	}
	end_subprocess();			/* exit here if subprocess */
	free_routes(&pmtx);			/* free working arrays */
	free(src_rem);
	free(dst_rem);
	return(newmig);
}

/* Check if prospective vertex would create overlapping triangle */
static int
overlaps_tri(const RBFNODE *bv0, const RBFNODE *bv1, const RBFNODE *pv)
{
	const MIGRATION	*ej;
	RBFNODE		*vother[2];
	int		im_rev;
						/* find shared edge in mesh */
	for (ej = pv->ejl; ej != NULL; ej = nextedge(pv,ej)) {
		const RBFNODE	*tv = opp_rbf(pv,ej);
		if (tv == bv0) {
			im_rev = is_rev_tri(ej->rbfv[0]->invec,
					ej->rbfv[1]->invec, bv1->invec);
			break;
		}
		if (tv == bv1) {
			im_rev = is_rev_tri(ej->rbfv[0]->invec,
					ej->rbfv[1]->invec, bv0->invec);
			break;
		}
	}
	if (!get_triangles(vother, ej))		/* triangle on same side? */
		return(0);
	return(vother[im_rev] != NULL);
}

/* Find convex hull vertex to complete triangle (oriented call) */
static RBFNODE *
find_chull_vert(const RBFNODE *rbf0, const RBFNODE *rbf1)
{
	FVECT	vmid, vejn, vp;
	RBFNODE	*rbf, *rbfbest = NULL;
	double	dprod, area2, bestarea2 = FHUGE, bestdprod = -.5;

	VSUB(vejn, rbf1->invec, rbf0->invec);
	VADD(vmid, rbf0->invec, rbf1->invec);
	if (normalize(vejn) == 0 || normalize(vmid) == 0)
		return(NULL);
						/* XXX exhaustive search */
	/* Find triangle with minimum rotation from perpendicular */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		if ((rbf == rbf0) | (rbf == rbf1))
			continue;
		tri_orient(vp, rbf0->invec, rbf1->invec, rbf->invec);
		if (DOT(vp, vmid) <= FTINY)
			continue;		/* wrong orientation */
		area2 = .25*DOT(vp,vp);
		VSUB(vp, rbf->invec, vmid);
		dprod = -DOT(vp, vejn);
		VSUM(vp, vp, vejn, dprod);	/* above guarantees non-zero */
		dprod = DOT(vp, vmid) / VLEN(vp);
		if (dprod <= bestdprod + FTINY*(1 - 2*(area2 < bestarea2)))
			continue;		/* found better already */
		if (overlaps_tri(rbf0, rbf1, rbf))
			continue;		/* overlaps another triangle */
		rbfbest = rbf;
		bestdprod = dprod;		/* new one to beat */
		bestarea2 = area2;
	}
	return(rbfbest);
}

/* Create new migration edge and grow mesh recursively around it */
static void
mesh_from_edge(MIGRATION *edge)
{
	MIGRATION	*ej0, *ej1;
	RBFNODE		*tvert[2];
	
	if (edge == NULL)
		return;
						/* triangle on either side? */
	get_triangles(tvert, edge);
	if (tvert[0] == NULL) {			/* grow mesh on right */
		tvert[0] = find_chull_vert(edge->rbfv[0], edge->rbfv[1]);
		if (tvert[0] != NULL) {
			if (tvert[0]->ord > edge->rbfv[0]->ord)
				ej0 = create_migration(edge->rbfv[0], tvert[0]);
			else
				ej0 = create_migration(tvert[0], edge->rbfv[0]);
			if (tvert[0]->ord > edge->rbfv[1]->ord)
				ej1 = create_migration(edge->rbfv[1], tvert[0]);
			else
				ej1 = create_migration(tvert[0], edge->rbfv[1]);
			mesh_from_edge(ej0);
			mesh_from_edge(ej1);
		}
	} else if (tvert[1] == NULL) {		/* grow mesh on left */
		tvert[1] = find_chull_vert(edge->rbfv[1], edge->rbfv[0]);
		if (tvert[1] != NULL) {
			if (tvert[1]->ord > edge->rbfv[0]->ord)
				ej0 = create_migration(edge->rbfv[0], tvert[1]);
			else
				ej0 = create_migration(tvert[1], edge->rbfv[0]);
			if (tvert[1]->ord > edge->rbfv[1]->ord)
				ej1 = create_migration(edge->rbfv[1], tvert[1]);
			else
				ej1 = create_migration(tvert[1], edge->rbfv[1]);
			mesh_from_edge(ej0);
			mesh_from_edge(ej1);
		}
	}
}

/* Add normal direction if missing */
static void
check_normal_incidence(void)
{
	static const FVECT	norm_vec = {.0, .0, 1.};
	const int		saved_nprocs = nprocs;
	RBFNODE			*near_rbf, *mir_rbf, *rbf;
	double			bestd;
	int			n;

	if (dsf_list == NULL)
		return;				/* XXX should be error? */
	near_rbf = dsf_list;
	bestd = input_orient*near_rbf->invec[2];
	if (single_plane_incident) {		/* ordered plane incidence? */
		if (bestd >= 1.-2.*FTINY)
			return;			/* already have normal */
	} else {
		switch (inp_coverage) {
		case INP_QUAD1:
		case INP_QUAD2:
		case INP_QUAD3:
		case INP_QUAD4:
			break;			/* quadrilateral symmetry? */
		default:
			return;			/* else we can interpolate */
		}
		for (rbf = near_rbf->next; rbf != NULL; rbf = rbf->next) {
			const double	d = input_orient*rbf->invec[2];
			if (d >= 1.-2.*FTINY)
				return;		/* seems we have normal */
			if (d > bestd) {
				near_rbf = rbf;
				bestd = d;
			}
		}
	}
	if (mig_list != NULL) {			/* need to be called first */
		fprintf(stderr, "%s: Late call to check_normal_incidence()\n",
				progname);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "Interpolating normal incidence by mirroring (%.1f,%.1f)\n",
			get_theta180(near_rbf->invec), get_phi360(near_rbf->invec));
#endif
						/* mirror nearest incidence */
	n = sizeof(RBFNODE) + sizeof(RBFVAL)*(near_rbf->nrbf-1);
	mir_rbf = (RBFNODE *)malloc(n);
	if (mir_rbf == NULL)
		goto memerr;
	memcpy(mir_rbf, near_rbf, n);
	mir_rbf->ord = near_rbf->ord - 1;	/* not used, I think */
	mir_rbf->next = NULL;
	rev_rbf_symmetry(mir_rbf, MIRROR_X|MIRROR_Y);
	nprocs = 1;				/* compute migration matrix */
	if (mig_list != create_migration(mir_rbf, near_rbf))
		exit(1);			/* XXX should never happen! */
						/* interpolate normal dist. */
	rbf = e_advect_rbf(mig_list, norm_vec, 2*near_rbf->nrbf);
	nprocs = saved_nprocs;			/* final clean-up */
	free(mir_rbf);
	free(mig_list);
	mig_list = near_rbf->ejl = NULL;
	insert_dsf(rbf);			/* insert interpolated normal */
	return;
memerr:
	fprintf(stderr, "%s: Out of memory in check_normal_incidence()\n",
				progname);
	exit(1);
}
	
/* Build our triangle mesh from recorded RBFs */
void
build_mesh(void)
{
	double		best2 = M_PI*M_PI;
	RBFNODE		*shrt_edj[2];
	RBFNODE		*rbf0, *rbf1;
						/* add normal if needed */
	check_normal_incidence();
						/* check if isotropic */
	if (single_plane_incident) {
		for (rbf0 = dsf_list; rbf0 != NULL; rbf0 = rbf0->next)
			if (rbf0->next != NULL)
				create_migration(rbf0, rbf0->next);
		await_children(nchild);
		return;
	}
	shrt_edj[0] = shrt_edj[1] = NULL;	/* start w/ shortest edge */
	for (rbf0 = dsf_list; rbf0 != NULL; rbf0 = rbf0->next)
	    for (rbf1 = rbf0->next; rbf1 != NULL; rbf1 = rbf1->next) {
		double	dist2 = 2. - 2.*DOT(rbf0->invec,rbf1->invec);
		if (dist2 < best2) {
			shrt_edj[0] = rbf0;
			shrt_edj[1] = rbf1;
			best2 = dist2;
		}
	}
	if (shrt_edj[0] == NULL) {
		fprintf(stderr, "%s: Cannot find shortest edge\n", progname);
		exit(1);
	}
						/* build mesh from this edge */
	if (shrt_edj[0]->ord < shrt_edj[1]->ord)
		mesh_from_edge(create_migration(shrt_edj[0], shrt_edj[1]));
	else
		mesh_from_edge(create_migration(shrt_edj[1], shrt_edj[0]));
						/* complete migrations */
	await_children(nchild);
}
