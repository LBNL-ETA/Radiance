#ifndef lint
static const char RCSid[] = "$Id: bsdfmesh.c,v 2.39 2017/10/06 00:23:09 greg Exp $";
#endif
/*
 * Create BSDF advection mesh from radial basis functions.
 *
 *	G. Ward
 */

#if !defined(_WIN32) && !defined(_WIN64)
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

/* Compute average DSF value at the given radius from central vector */
static double
eval_DSFsurround(const RBFNODE *rbf, const FVECT outvec, const double rad)
{
	const int	ninc = 12;
	const double	phinc = 2.*M_PI/ninc;
	double		sum = 0;
	int		n = 0;
	FVECT		tvec;
	int		i;
						/* compute initial vector */
	if (output_orient*outvec[2] >= 1.-FTINY) {
		tvec[0] = tvec[2] = 0;
		tvec[1] = 1;
	} else {
		tvec[0] = tvec[1] = 0;
		tvec[2] = 1;
	}
	geodesic(tvec, outvec, tvec, rad, GEOD_RAD);
						/* average surrounding DSF */
	for (i = 0; i < ninc; i++) {
		if (i) spinvector(tvec, tvec, outvec, phinc);
		if (tvec[2] > 0 ^ output_orient > 0)
			continue;
		sum += eval_rbfrep(rbf, tvec) * COSF(tvec[2]);
		++n;
	}
	if (n < 2)				/* should never happen! */
		return(sum);
	return(sum/(double)n);
}

/* Estimate single-lobe radius for DSF at the given outgoing angle */
static double
est_DSFrad(const RBFNODE *rbf, const FVECT outvec)
{
	const double	rad_epsilon = 0.01;
	const double	DSFtarget = 0.60653066 * eval_rbfrep(rbf,outvec) *
							COSF(outvec[2]);
	double		inside_rad = rad_epsilon;
	double		outside_rad = 0.5;
	double		DSFinside = eval_DSFsurround(rbf, outvec, inside_rad);
	double		DSFoutside = eval_DSFsurround(rbf, outvec, outside_rad);
#define	interp_rad	inside_rad + (outside_rad-inside_rad) * \
				(DSFtarget-DSFinside) / (DSFoutside-DSFinside)
						/* Newton's method (sort of) */
	do {
		double	test_rad = interp_rad;
		double	DSFtest;
		if ((test_rad >= outside_rad) | (test_rad <= inside_rad))
			test_rad = .5*(inside_rad + outside_rad);
		DSFtest = eval_DSFsurround(rbf, outvec, test_rad);
		if (DSFtest > DSFtarget) {
			inside_rad = test_rad;
			DSFinside = DSFtest;
		} else {
			outside_rad = test_rad;
			DSFoutside = DSFtest;
		}
	} while (outside_rad-inside_rad > rad_epsilon);

	return(.5*(inside_rad + outside_rad));
#undef interp_rad
}

static int
dbl_cmp(const void *p1, const void *p2)
{
	double	d1 = *(const double *)p1;
	double	d2 = *(const double *)p2;

	if (d1 > d2) return(1);
	if (d1 < d2) return(-1);
	return(0);
}

/* Conservative estimate of average BSDF value from current DSF's */
static void
comp_bsdf_spec(void)
{
	double		vmod_sum = 0;
	double		rad_sum = 0;
	int		n = 0;
	double		*cost_list = NULL;
	double		max_cost = 1.;
	RBFNODE		*rbf;
	FVECT		sdv;
						/* sort by incident altitude */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next)
		n++;
	if (n >= 10)
		cost_list = (double *)malloc(sizeof(double)*n);
	if (cost_list == NULL) {
		bsdf_spec_val = 0;
		bsdf_spec_rad = 0;
		return;
	}
	n = 0;
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next)
		cost_list[n++] = rbf->invec[2]*input_orient;
	qsort(cost_list, n, sizeof(double), dbl_cmp);
	max_cost = cost_list[(n+3)/4];		/* accept 25% nearest grazing */
	free(cost_list);
	n = 0;
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		double	this_rad, cosfact, vest;
		if (rbf->invec[2]*input_orient > max_cost)
			continue;
		sdv[0] = -rbf->invec[0];
		sdv[1] = -rbf->invec[1];
		sdv[2] = rbf->invec[2]*(2*(input_orient==output_orient) - 1);
		cosfact = COSF(sdv[2]);
		this_rad = est_DSFrad(rbf, sdv);
		vest = eval_rbfrep(rbf, sdv) * cosfact *
				(2.*M_PI) * this_rad*this_rad;
		if (vest > rbf->vtotal)		/* don't over-estimate energy */
			vest = rbf->vtotal;
		vmod_sum += vest / cosfact;	/* remove cosine factor */
		rad_sum += this_rad;
		++n;
	}
	bsdf_spec_rad = rad_sum/(double)n;
	bsdf_spec_val = vmod_sum/(2.*M_PI*n*bsdf_spec_rad*bsdf_spec_rad);
}

/* Create a new migration holder (sharing memory for multiprocessing) */
static MIGRATION *
new_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	size_t		memlen = sizeof(MIGRATION) +
				sizeof(float)*(from_rbf->nrbf*to_rbf->nrbf - 1);
	MIGRATION	*newmig;
#if defined(_WIN32) || defined(_WIN64)
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

#if defined(_WIN32) || defined(_WIN64)
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

/* Compute distance between two RBF lobes */
double
lobe_distance(RBFVAL *rbf1, RBFVAL *rbf2)
{
	FVECT	vfrom, vto;
	double	d, res;
					/* quadratic cost function */
	ovec_from_pos(vfrom, rbf1->gx, rbf1->gy);
	ovec_from_pos(vto, rbf2->gx, rbf2->gy);
	d = Acos(DOT(vfrom, vto));
	res = d*d;
	d = R2ANG(rbf2->crad) - R2ANG(rbf1->crad);
	res += d*d;
					/* neighborhood difference */
	res += NEIGH_FACT2 * neighborhood_dist2( rbf1->gx, rbf1->gy,
						rbf2->gx, rbf2->gy );
	return(res);
}


/* Compute and insert migration along directed edge (may fork child) */
static MIGRATION *
create_migration(RBFNODE *from_rbf, RBFNODE *to_rbf)
{
	MIGRATION	*newmig;
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

						/* compute transport plan */
	compute_nDSFs(from_rbf, to_rbf);
	plan_transport(newmig);

	for (i = from_rbf->nrbf; i--; ) {	/* normalize final matrix */
	    double	nf = rbf_volume(&from_rbf->rbfa[i]);
	    if (nf <= FTINY) continue;
	    nf = from_rbf->vtotal / nf;
	    for (j = to_rbf->nrbf; j--; )
		mtx_coef(newmig,i,j) *= nf;	/* row now sums to 1.0 */
	}
	end_subprocess();			/* exit here if subprocess */
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
			return;
		}
	}
	if (tvert[1] == NULL) {			/* grow mesh on left */
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
	static FVECT		norm_vec = {.0, .0, 1.};
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
		for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
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
	mir_rbf->ejl = NULL;
	rev_rbf_symmetry(mir_rbf, MIRROR_X|MIRROR_Y);
	nprocs = 1;				/* compute migration matrix */
	if (create_migration(mir_rbf, near_rbf) == NULL)
		exit(1);			/* XXX should never happen! */
	norm_vec[2] = input_orient;		/* interpolate normal dist. */
	rbf = e_advect_rbf(mig_list, norm_vec, 0);
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
						/* average specular peak */
	comp_bsdf_spec();
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
