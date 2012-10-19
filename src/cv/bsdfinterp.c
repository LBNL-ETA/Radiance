#ifndef lint
static const char RCSid[] = "$Id: bsdfinterp.c,v 2.1 2012/10/19 04:14:29 greg Exp $";
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
				/* migration edges drawn in raster fashion */
MIGRATION		*mig_grid[GRIDRES][GRIDRES];

#ifdef DEBUG
#include "random.h"
#include "bmpfile.h"
/* Hash pointer to byte value (must return 0 for NULL) */
static int
byte_hash(const void *p)
{
	size_t	h = (size_t)p;
	h ^= (size_t)p >> 8;
	h ^= (size_t)p >> 16;
	h ^= (size_t)p >> 24;
	return(h & 0xff);
}
/* Write out BMP image showing edges */
static void
write_edge_image(const char *fname)
{
	BMPHeader	*hdr = BMPmappedHeader(GRIDRES, GRIDRES, 0, 256);
	BMPWriter	*wtr;
	int		i, j;

	fprintf(stderr, "Writing incident mesh drawing to '%s'\n", fname);
	hdr->compr = BI_RLE8;
	for (i = 256; --i; ) {			/* assign random color map */
		hdr->palette[i].r = random() & 0xff;
		hdr->palette[i].g = random() & 0xff;
		hdr->palette[i].b = random() & 0xff;
						/* reject dark colors */
		i += (hdr->palette[i].r + hdr->palette[i].g +
						hdr->palette[i].b < 128);
	}
	hdr->palette[0].r = hdr->palette[0].g = hdr->palette[0].b = 0;
						/* open output */
	wtr = BMPopenOutputFile(fname, hdr);
	if (wtr == NULL) {
		free(hdr);
		return;
	}
	for (i = 0; i < GRIDRES; i++) {		/* write scanlines */
		for (j = 0; j < GRIDRES; j++)
			wtr->scanline[j] = byte_hash(mig_grid[i][j]);
		if (BMPwriteScanline(wtr) != BIR_OK)
			break;
	}
	BMPcloseOutput(wtr);			/* close & clean up */
}
#endif

/* Draw edge list into mig_grid array */
void
draw_edges(void)
{
	int		nnull = 0, ntot = 0;
	MIGRATION	*ej;
	int		p0[2], p1[2];

	memset(mig_grid, 0, sizeof(mig_grid));
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
#ifdef DEBUG
	write_edge_image("bsdf_edges.bmp");
#endif
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
	if ((vert[3] == NULL) | (vert[4] != NULL))
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
			if (rbf->next != NULL &&
					input_orient*rbf->next->invec[2] <
							input_orient*invec[2]) {
				for (miga[0] = rbf->ejl; miga[0] != NULL;
						miga[0] = nextedge(rbf,miga[0]))
					if (opp_rbf(rbf,miga[0]) == rbf->next)
						return(0);
				break;
			}
		}
		return(-1);			/* outside range! */
	}
	{					/* else use triangle mesh */
		const int	sym = use_symmetry(invec);
		unsigned char	floodmap[GRIDRES][(GRIDRES+7)/8];
		int		pstart[2];
		RBFNODE		*vother;
		MIGRATION	*ej;
		int		i;

		pos_from_vec(pstart, invec);
		memset(floodmap, 0, sizeof(floodmap));
						/* call flooding function */
		if (!identify_tri(miga, floodmap, pstart[0], pstart[1]))
			return(-1);		/* outside mesh */
		if ((miga[0] == NULL) | (miga[2] == NULL))
			return(-1);		/* should never happen */
		if (miga[1] == NULL)
			return(sym);		/* on edge */
						/* verify triangle */
		if (!order_triangle(miga)) {
#ifdef DEBUG
			fputs("Munged triangle in get_interp()\n", stderr);
#endif
			vother = NULL;		/* find triangle from edge */
			for (i = 3; i--; ) {
			    RBFNODE	*tpair[2];
			    if (get_triangles(tpair, miga[i]) &&
					(vother = tpair[ is_rev_tri(
							miga[i]->rbfv[0]->invec,
							miga[i]->rbfv[1]->invec,
							invec) ]) != NULL)
					break;
			}
			if (vother == NULL) {	/* couldn't find 3rd vertex */
#ifdef DEBUG
				fputs("No triangle in get_interp()\n", stderr);
#endif
				return(-1);
			}
						/* reassign other two edges */
			for (ej = vother->ejl; ej != NULL;
						ej = nextedge(vother,ej)) {
				RBFNODE	*vorig = opp_rbf(vother,ej);
				if (vorig == miga[i]->rbfv[0])
					miga[(i+1)%3] = ej;
				else if (vorig == miga[i]->rbfv[1])
					miga[(i+2)%3] = ej;
			}
			if (!order_triangle(miga)) {
#ifdef DEBUG
				fputs("Bad triangle in get_interp()\n", stderr);
#endif
				return(-1);
			}
		}
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
		for (k = (miga[0]->mtx[mtx_ndx(miga[0],i,j)] > FTINY) *
					mtx_ncols(miga[2]); k--; )
			n += (miga[2]->mtx[mtx_ndx(miga[2],i,k)] > FTINY &&
				miga[1]->mtx[mtx_ndx(miga[1],j,k)] > FTINY);
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
		    rbf->rbfa[n].crad = ANG2R(sqrt(srad2 + t*rad2k*rad2k));
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
	rev_rbf_symmetry(rbf, sym);
	return(rbf);
}
