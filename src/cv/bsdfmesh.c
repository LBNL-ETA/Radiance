#ifndef lint
static const char RCSid[] = "$Id: bsdfmesh.c,v 2.2 2012/10/20 07:02:00 greg Exp $";
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
				/* number of processes to run */
int			nprocs = 1;
				/* number of children (-1 in child) */
static int		nchild = 0;

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

/* Compute (and allocate) migration price matrix for optimization */
static float *
price_routes(const RBFNODE *from_rbf, const RBFNODE *to_rbf)
{
	float	*pmtx = (float *)malloc(sizeof(float) *
					from_rbf->nrbf * to_rbf->nrbf);
	FVECT	*vto = (FVECT *)malloc(sizeof(FVECT) * to_rbf->nrbf);
	int	i, j;

	if ((pmtx == NULL) | (vto == NULL)) {
		fprintf(stderr, "%s: Out of memory in migration_costs()\n",
				progname);
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
			fprintf(stderr, "%s: Out of memory in min_cost()\n",
					progname);
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
	const double	maxamt = .1;
	const double	minamt = maxamt*5e-6;
	static double	*src_cost = NULL;
	static int	n_alloc = 0;
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
			fprintf(stderr, "%s: Out of memory in migration_step()\n",
					progname);
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
	    if (src_rem[cur.s] <= minamt)
		    continue;
	    cur.d = -1;				/* examine cheapest dest. */
	    for (i = mtx_ncols(mig); i--; )
		if (dst_rem[i] > minamt &&
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
	    for (i = mtx_nrows(mig); i--; )
		if (i != cur.s)
			cost_others += min_cost(src_rem[i], dst_rem,
						price, mtx_ncols(mig))
					- src_cost[i];
	    dst_rem[cur.d] += cur.amt;		/* undo trial move */
	    cur.price += cost_others/cur.amt;	/* adjust effective price */
	    if (cur.price < best.price)		/* are we better than best? */
		    best = cur;
	}
	if ((best.s < 0) | (best.d < 0))
		return(.0);
						/* make the actual move */
	mtx_coef(mig,best.s,best.d) += best.amt;
	src_rem[best.s] -= best.amt;
	dst_rem[best.d] -= best.amt;
	return(best.amt);
}

#ifdef DEBUG
static char *
thetaphi(const FVECT v)
{
	static char	buf[128];
	double		theta, phi;

	theta = 180./M_PI*acos(v[2]);
	phi = 180./M_PI*atan2(v[1],v[0]);
	sprintf(buf, "(%.0f,%.0f)", theta, phi);

	return(buf);
}
#endif

/* Compute and insert migration along directed edge (may fork child) */
static MIGRATION *
create_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	const double	end_thresh = 5e-6;
	float		*pmtx;
	MIGRATION	*newmig;
	double		*src_rem, *dst_rem;
	double		total_rem = 1., move_amt;
	int		i;
						/* check if exists already */
	for (newmig = from_rbf->ejl; newmig != NULL;
			newmig = nextedge(from_rbf,newmig))
		if (newmig->rbfv[1] == to_rbf)
			return(NULL);
						/* else allocate */
	newmig = new_migration(from_rbf, to_rbf);
	if (run_subprocess())
		return(newmig);			/* child continues */
	pmtx = price_routes(from_rbf, to_rbf);
	src_rem = (double *)malloc(sizeof(double)*from_rbf->nrbf);
	dst_rem = (double *)malloc(sizeof(double)*to_rbf->nrbf);
	if ((src_rem == NULL) | (dst_rem == NULL)) {
		fprintf(stderr, "%s: Out of memory in create_migration()\n",
				progname);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "Building path from (theta,phi) %s ",
			thetaphi(from_rbf->invec));
	fprintf(stderr, "to %s with %d x %d matrix\n",
			thetaphi(to_rbf->invec), 
			from_rbf->nrbf, to_rbf->nrbf);
#endif
						/* starting quantities */
	memset(newmig->mtx, 0, sizeof(float)*from_rbf->nrbf*to_rbf->nrbf);
	for (i = from_rbf->nrbf; i--; )
		src_rem[i] = rbf_volume(&from_rbf->rbfa[i]) / from_rbf->vtotal;
	for (i = to_rbf->nrbf; i--; )
		dst_rem[i] = rbf_volume(&to_rbf->rbfa[i]) / to_rbf->vtotal;
	do {					/* move a bit at a time */
		move_amt = migration_step(newmig, src_rem, dst_rem, pmtx);
		total_rem -= move_amt;
#ifdef DEBUG
		if (!nchild)
			fprintf(stderr, "\r%.9f remaining...", total_rem);
#endif
	} while ((total_rem > end_thresh) & (move_amt > 0));
#ifdef DEBUG
	if (!nchild) fputs("done.\n", stderr);
	else fprintf(stderr, "finished with %.9f remaining\n", total_rem);
#endif
	for (i = from_rbf->nrbf; i--; ) {	/* normalize final matrix */
	    float	nf = rbf_volume(&from_rbf->rbfa[i]);
	    int		j;
	    if (nf <= FTINY) continue;
	    nf = from_rbf->vtotal / nf;
	    for (j = to_rbf->nrbf; j--; )
		mtx_coef(newmig,i,j) *= nf;
	}
	end_subprocess();			/* exit here if subprocess */
	free(pmtx);				/* free working arrays */
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

/* Find context hull vertex to complete triangle (oriented call) */
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
		VSUB(vp, rbf->invec, rbf0->invec);
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
	
/* Build our triangle mesh from recorded RBFs */
void
build_mesh(void)
{
	double		best2 = M_PI*M_PI;
	RBFNODE		*shrt_edj[2];
	RBFNODE		*rbf0, *rbf1;
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
