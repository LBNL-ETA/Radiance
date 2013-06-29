#ifndef lint
static const char RCSid[] = "$Id: bsdfinterp.c,v 2.11 2013/06/29 21:03:25 greg Exp $";
#endif
/*
 * Interpolate BSDF data from radial basis functions in advection mesh.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bsdfrep.h"

/* Insert vertex in ordered list */
static void
insert_vert(RBFNODE **vlist, RBFNODE *v)
{
	int	i, j;

	for (i = 0; vlist[i] != NULL; i++) {
		if (v == vlist[i])
			return;
		if (v->ord < vlist[i]->ord)
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
static int
order_triangle(MIGRATION *miga[3])
{
	RBFNODE		*vert[7];
	MIGRATION	*ord[3];
	int		i;
						/* order vertices, first */
	memset(vert, 0, sizeof(vert));
	for (i = 3; i--; ) {
		if (miga[i] == NULL)
			return(0);
		insert_vert(vert, miga[i]->rbfv[0]);
		insert_vert(vert, miga[i]->rbfv[1]);
	}
						/* should be just 3 vertices */
	if ((vert[2] == NULL) | (vert[3] != NULL))
		return(0);
						/* identify edge 0 */
	for (i = 3; i--; )
		if (miga[i]->rbfv[0] == vert[0] &&
				miga[i]->rbfv[1] == vert[1]) {
			ord[0] = miga[i];
			break;
		}
	if (i < 0)
		return(0);
						/* identify edge 1 */
	for (i = 3; i--; )
		if (miga[i]->rbfv[0] == vert[1] &&
				miga[i]->rbfv[1] == vert[2]) {
			ord[1] = miga[i];
			break;
		}
	if (i < 0)
		return(0);
						/* identify edge 2 */
	for (i = 3; i--; )
		if (miga[i]->rbfv[0] == vert[0] &&
				miga[i]->rbfv[1] == vert[2]) {
			ord[2] = miga[i];
			break;
		}
	if (i < 0)
		return(0);
						/* reassign order */
	miga[0] = ord[0]; miga[1] = ord[1]; miga[2] = ord[2];
	return(1);
}

/* Determine if we are close enough to an edge */
static int
on_edge(const MIGRATION *ej, const FVECT ivec)
{
	double	cos_a, cos_b, cos_c, cos_aplusb;
					/* use triangle inequality */
	cos_a = DOT(ej->rbfv[0]->invec, ivec);
	if (cos_a <= 0)
		return(0);

	cos_b = DOT(ej->rbfv[1]->invec, ivec);
	if (cos_b <= 0)
		return(0);

	cos_aplusb = cos_a*cos_b - sqrt((1.-cos_a*cos_a)*(1.-cos_b*cos_b));
	if (cos_aplusb <= 0)
		return(0);

	cos_c = DOT(ej->rbfv[0]->invec, ej->rbfv[1]->invec);

	return(cos_c - cos_aplusb < .001);
}

/* Determine if we are inside the given triangle */
static int
in_tri(const RBFNODE *v1, const RBFNODE *v2, const RBFNODE *v3, const FVECT p)
{
	FVECT	vc;
	int	sgn1, sgn2, sgn3;
					/* signed volume test */
	VCROSS(vc, v1->invec, v2->invec);
	sgn1 = (DOT(p, vc) > 0);
	VCROSS(vc, v2->invec, v3->invec);
	sgn2 = (DOT(p, vc) > 0);
	if (sgn1 != sgn2)
		return(0);
	VCROSS(vc, v3->invec, v1->invec);
	sgn3 = (DOT(p, vc) > 0);
	return(sgn2 == sgn3);
}

/* Test (and set) bitmap for edge */
static int
check_edge(unsigned char *emap, int nedges, const MIGRATION *mig, int mark)
{
	int	ejndx, bit2check;

	if (mig->rbfv[0]->ord > mig->rbfv[1]->ord)
		ejndx = mig->rbfv[1]->ord + (nedges-1)*mig->rbfv[0]->ord;
	else
		ejndx = mig->rbfv[0]->ord + (nedges-1)*mig->rbfv[1]->ord;

	bit2check = 1<<(ejndx&07);

	if (emap[ejndx>>3] & bit2check)
		return(0);
	if (mark)
		emap[ejndx>>3] |= bit2check;
	return(1);
}

/* Compute intersection with the given position over remaining mesh */
static int
in_mesh(MIGRATION *miga[3], unsigned char *emap, int nedges,
			const FVECT ivec, MIGRATION *mig)
{
	RBFNODE		*tv[2];
	MIGRATION	*sej[2], *dej[2];
	int		i;
						/* check visitation record */
	if (!check_edge(emap, nedges, mig, 1))
		return(0);
	if (on_edge(mig, ivec)) {
		miga[0] = mig;			/* close enough to edge */
		return(1);
	}
	if (!get_triangles(tv, mig))		/* do triangles either side? */
		return(0);
	for (i = 2; i--; ) {			/* identify edges to check */
		MIGRATION	*ej;
		sej[i] = dej[i] = NULL;
		if (tv[i] == NULL)
			continue;
		for (ej = tv[i]->ejl; ej != NULL; ej = nextedge(tv[i],ej)) {
			RBFNODE	*rbfop = opp_rbf(tv[i],ej);
			if (rbfop == mig->rbfv[0]) {
				if (check_edge(emap, nedges, ej, 0))
					sej[i] = ej;
			} else if (rbfop == mig->rbfv[1]) {
				if (check_edge(emap, nedges, ej, 0))
					dej[i] = ej;
			}
		}
	}
	for (i = 2; i--; ) {			/* check triangles just once */
		if (sej[i] != NULL && in_mesh(miga, emap, nedges, ivec, sej[i]))
			return(1);
		if (dej[i] != NULL && in_mesh(miga, emap, nedges, ivec, dej[i]))
			return(1);
		if ((sej[i] == NULL) | (dej[i] == NULL))
			continue;
		if (in_tri(mig->rbfv[0], mig->rbfv[1], tv[i], ivec)) {
			miga[0] = mig;
			miga[1] = sej[i];
			miga[2] = dej[i];
			return(1);
		}
	}
	return(0);				/* not near this edge */
}

/* Find edge(s) for interpolating the given vector, applying symmetry */
int
get_interp(MIGRATION *miga[3], FVECT invec)
{
	miga[0] = miga[1] = miga[2] = NULL;
	if (single_plane_incident) {		/* isotropic BSDF? */
	    RBFNODE	*rbf;			/* find edge we're on */
	    for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		if (input_orient*rbf->invec[2] < input_orient*invec[2])
			break;
		if (rbf->next != NULL && input_orient*rbf->next->invec[2] <
							input_orient*invec[2]) {
		    for (miga[0] = rbf->ejl; miga[0] != NULL;
					miga[0] = nextedge(rbf,miga[0]))
			if (opp_rbf(rbf,miga[0]) == rbf->next) {
				double	nf = 1. - rbf->invec[2]*rbf->invec[2];
				if (nf > FTINY) {	/* rotate to match */
					nf = sqrt((1.-invec[2]*invec[2])/nf);
					invec[0] = nf*rbf->invec[0];
					invec[1] = nf*rbf->invec[1];
				}
				return(0);
			}
		    break;
		}
	    }
	    return(-1);				/* outside range! */
	}
	{					/* else use triangle mesh */
		int		sym = use_symmetry(invec);
		int		nedges = 0;
		MIGRATION	*mep;
		unsigned char	*emap;
						/* clear visitation map */
		for (mep = mig_list; mep != NULL; mep = mep->next)
			++nedges;
		emap = (unsigned char *)calloc((nedges*(nedges-1) + 7)>>3, 1);
		if (emap == NULL) {
			fprintf(stderr, "%s: Out of memory in get_interp()\n",
					progname);
			exit(1);
		}
						/* identify intersection  */
		if (!in_mesh(miga, emap, nedges, invec, mig_list)) {
#ifdef DEBUG
			fprintf(stderr,
			"Incident angle (%.1f,%.1f) deg. outside mesh\n",
					get_theta180(invec), get_phi360(invec));
#endif
			sym = -1;		/* outside mesh */
		} else if (miga[1] != NULL &&
				(miga[2] == NULL || !order_triangle(miga))) {
#ifdef DEBUG
			fputs("Munged triangle in get_interp()\n", stderr);
#endif
			sym = -1;
		}
		free(emap);
		return(sym);			/* return in standard order */
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
	t = Acos(DOT(invec, mig->rbfv[0]->invec));
	if (t < M_PI/grid_res) {		/* near first DSF */
		n = sizeof(RBFNODE) + sizeof(RBFVAL)*(mig->rbfv[0]->nrbf-1);
		rbf = (RBFNODE *)malloc(n);
		if (rbf == NULL)
			goto memerr;
		memcpy(rbf, mig->rbfv[0], n);	/* just duplicate */
		rbf->next = NULL; rbf->ejl = NULL;
		return(rbf);
	}
	full_dist = acos(DOT(mig->rbfv[0]->invec, mig->rbfv[1]->invec));
	if (t > full_dist-M_PI/grid_res) {	/* near second DSF */
		n = sizeof(RBFNODE) + sizeof(RBFVAL)*(mig->rbfv[1]->nrbf-1);
		rbf = (RBFNODE *)malloc(n);
		if (rbf == NULL)
			goto memerr;
		memcpy(rbf, mig->rbfv[1], n);	/* just duplicate */
		rbf->next = NULL; rbf->ejl = NULL;
		return(rbf);
	}
	t /= full_dist;	
	n = 0;					/* count migrating particles */
	for (i = 0; i < mtx_nrows(mig); i++)
	    for (j = 0; j < mtx_ncols(mig); j++)
		n += (mtx_coef(mig,i,j) > FTINY);
#ifdef DEBUG
	fprintf(stderr, "Input RBFs have %d, %d nodes -> output has %d\n",
			mig->rbfv[0]->nrbf, mig->rbfv[1]->nrbf, n);
#endif
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
		if ((mv = mtx_coef(mig,i,j)) > FTINY) {
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
	fprintf(stderr, "%s: Out of memory in e_advect_rbf()\n", progname);
	exit(1);
	return(NULL);	/* pro forma return */
}

/* Partially advect between recorded incident angles and allocate new RBF */
RBFNODE *
advect_rbf(const FVECT invec)
{
	FVECT		sivec;
	MIGRATION	*miga[3];
	RBFNODE		*rbf;
	int		sym;
	float		mbfact, mcfact;
	int		n, i, j, k;
	FVECT		v0, v1, v2;
	double		s, t;

	VCOPY(sivec, invec);			/* find triangle/edge */
	sym = get_interp(miga, sivec);
	if (sym < 0)				/* can't interpolate? */
		return(NULL);
	if (miga[1] == NULL) {			/* advect along edge? */
		rbf = e_advect_rbf(miga[0], sivec);
		if (single_plane_incident)
			rotate_rbf(rbf, invec);
		else
			rev_rbf_symmetry(rbf, sym);
		return(rbf);
	}
#ifdef DEBUG
	if (miga[0]->rbfv[0] != miga[2]->rbfv[0] |
			miga[0]->rbfv[1] != miga[1]->rbfv[0] |
			miga[1]->rbfv[1] != miga[2]->rbfv[1]) {
		fprintf(stderr, "%s: Triangle vertex screw-up!\n", progname);
		exit(1);
	}
#endif
						/* figure out position */
	fcross(v0, miga[2]->rbfv[0]->invec, miga[2]->rbfv[1]->invec);
	normalize(v0);
	fcross(v2, miga[1]->rbfv[0]->invec, miga[1]->rbfv[1]->invec);
	normalize(v2);
	fcross(v1, sivec, miga[1]->rbfv[1]->invec);
	normalize(v1);
	s = acos(DOT(v0,v1)) / acos(DOT(v0,v2));
	geodesic(v1, miga[0]->rbfv[0]->invec, miga[0]->rbfv[1]->invec,
			s, GEOD_REL);
	t = acos(DOT(v1,sivec)) / acos(DOT(v1,miga[1]->rbfv[1]->invec));
	n = 0;					/* count migrating particles */
	for (i = 0; i < mtx_nrows(miga[0]); i++)
	    for (j = 0; j < mtx_ncols(miga[0]); j++)
		for (k = (mtx_coef(miga[0],i,j) > FTINY) *
					mtx_ncols(miga[2]); k--; )
			n += (mtx_coef(miga[2],i,k) > FTINY ||
				mtx_coef(miga[1],j,k) > FTINY);
#ifdef DEBUG
	fprintf(stderr, "Input RBFs have %d, %d, %d nodes -> output has %d\n",
			miga[0]->rbfv[0]->nrbf, miga[0]->rbfv[1]->nrbf,
			miga[2]->rbfv[1]->nrbf, n);
#endif
	rbf = (RBFNODE *)malloc(sizeof(RBFNODE) + sizeof(RBFVAL)*(n-1));
	if (rbf == NULL) {
		fprintf(stderr, "%s: Out of memory in advect_rbf()\n", progname);
		exit(1);
	}
	rbf->next = NULL; rbf->ejl = NULL;
	VCOPY(rbf->invec, sivec);
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
		const float	ma = mtx_coef(miga[0],i,j);
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
		    float		mb = mtx_coef(miga[1],j,k);
		    float		mc = mtx_coef(miga[2],i,k);
		    const RBFVAL	*rbf2k;
		    double		rad2k;
		    int			pos[2];
		    if ((mb <= FTINY) & (mc <= FTINY))
			continue;
		    rbf2k = &miga[2]->rbfv[1]->rbfa[k];
		    rbf->rbfa[n].peak = w0i * ma * (mb*mbfact + mc*mcfact);
		    rad2k = R2ANG(rbf2k->crad);
		    rbf->rbfa[n].crad = ANG2R(sqrt(srad2 + t*rad2k*rad2k));
		    ovec_from_pos(v2, rbf2k->gx, rbf2k->gy);
		    geodesic(v2, v1, v2, t, GEOD_REL);
		    pos_from_vec(pos, v2);
		    rbf->rbfa[n].gx = pos[0];
		    rbf->rbfa[n].gy = pos[1];
		    ++n;
		}
	    }
	}
	rbf->vtotal = miga[0]->rbfv[0]->vtotal * (mbfact + mcfact);
	rev_rbf_symmetry(rbf, sym);
	return(rbf);
}
