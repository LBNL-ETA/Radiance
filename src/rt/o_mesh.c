#ifndef lint
static const char RCSid[] = "$Id: o_mesh.c,v 2.5 2003/03/14 21:27:46 greg Exp $";
#endif
/*
 *  Routines for computing ray intersections with meshes.
 *
 *  Intersection with a triangle mesh is based on Segura and Feito's
 *  WSCG 2001 paper, "Algorithms to Test Ray-Triangle Intersection,
 *  Comparative Study."  This method avoids additional storage
 *  requirements, floating divides, and allows some savings by
 *  caching ray-edge comparisons that are otherwise repeated locally
 *  in typical mesh geometries.  (This is our own optimization.)
 *
 *  The code herein is quite similar to that in o_instance.c, the
 *  chief differences being the custom triangle intersection routines
 *  and the fact that an "OBJECT" in the mesh octree is not an index
 *  into the Radiance OBJREC list, but a mesh triangle index.  We still
 *  utilize the standard octree traversal code by setting the hitf
 *  function pointer in the RAY struct to our custom mesh_hit() call.
 */

#include  "copyright.h"

#include  "ray.h"

#include  "mesh.h"

#include  "tmesh.h"


#define  EDGE_CACHE_SIZ		251	/* length of mesh edge cache */

#define  curmi			(edge_cache.mi)
#define  curmsh			(curmi->msh)


/* Cache of signed volumes for this ray and this mesh */
struct EdgeCache {
	OBJREC		*o;	/* mesh object */
	MESHINST	*mi;	/* current mesh instance */
	struct EdgeSide {
		int4	v1i, v2i;	/* vertex indices (lowest first) */
		short	signum;		/* signed volume */
	}		cache[EDGE_CACHE_SIZ];
}	edge_cache;


static void
prep_edge_cache(o)		/* get instance and clear edge cache */
OBJREC	*o;
{
					/* get mesh instance */
	edge_cache.mi = getmeshinst(edge_cache.o = o, IO_ALL);
					/* clear edge cache */
	bzero((void *)edge_cache.cache, sizeof(edge_cache.cache));
}


static int
signed_volume(r, v1, v2)	/* get signed volume for ray and edge */
register RAY	*r;
int4		v1, v2;
{
	int				reversed = 0;
	register struct EdgeSide	*ecp;
	
	if (v1 > v2) {
		int4	t = v2; v2 = v1; v1 = t;
		reversed = 1;
	}
	ecp = &edge_cache.cache[((v2<<11 ^ v1) & 0x7fffffff) % EDGE_CACHE_SIZ];
	if (ecp->v1i != v1 || ecp->v2i != v2) {
		MESHVERT	tv1, tv2;	/* compute signed volume */
		FVECT		v2d;
		double		vol;
		if (!getmeshvert(&tv1, edge_cache.mi->msh, v1, MT_V) ||
			    !getmeshvert(&tv2, edge_cache.mi->msh, v2, MT_V))
			objerror(edge_cache.o, INTERNAL,
				"missing mesh vertex in signed_volume");
		VSUB(v2d, tv2.v, r->rorg);
		vol = (tv1.v[0] - r->rorg[0]) *
				(v2d[1]*r->rdir[2] - v2d[2]*r->rdir[1]);
		vol += (tv1.v[1] - r->rorg[1]) *
				(v2d[2]*r->rdir[0] - v2d[0]*r->rdir[2]);
		vol += (tv1.v[2] - r->rorg[2]) *
				(v2d[0]*r->rdir[1] - v2d[1]*r->rdir[0]);
		if (vol > .0)
			ecp->signum = 1;
		else if (vol < .0)
			ecp->signum = -1;
		else
			ecp->signum = 0;
		ecp->v1i = v1;
		ecp->v2i = v2;
	}
	return(reversed ? -ecp->signum : ecp->signum);
}


static void
mesh_hit(oset, r)		/* intersect ray with mesh triangle(s) */
OBJECT	*oset;
RAY	*r;
{
	int4		tvi[3];
	int		sv1, sv2, sv3;
	MESHVERT	tv[3];
	OBJECT		tmod;
	FVECT		va, vb, nrm;
	double		d;
	int		i;
					/* check each triangle */
	for (i = oset[0]; i > 0; i--) {
		if (!getmeshtrivid(tvi, &tmod, curmsh, oset[i]))
			objerror(edge_cache.o, INTERNAL,
				"missing triangle vertices in mesh_hit");
		sv1 = signed_volume(r, tvi[0], tvi[1]);
		sv2 = signed_volume(r, tvi[1], tvi[2]);
		sv3 = signed_volume(r, tvi[2], tvi[0]);
		if (sv1 != sv2 || sv2 != sv3)	/* compare volume signs */
			if (sv1 && sv2 && sv3)
				continue;
						/* compute intersection */
		getmeshvert(&tv[0], curmsh, tvi[0], MT_V);
		getmeshvert(&tv[1], curmsh, tvi[1], MT_V);
		getmeshvert(&tv[2], curmsh, tvi[2], MT_V);
		VSUB(va, tv[0].v, tv[2].v);
		VSUB(vb, tv[1].v, tv[0].v);
		VCROSS(nrm, va, vb);
		d = DOT(r->rdir, nrm);
		if (d == 0.0)
			continue;		/* ray is tangent */
		VSUB(va, tv[0].v, r->rorg);
		d = DOT(va, nrm) / d;
		if (d <= FTINY || d >= r->rot)
			continue;		/* not good enough */
		r->robj = oset[i];		/* else record hit */
		r->ro = edge_cache.o;
		r->rot = d;
		VSUM(r->rop, r->rorg, r->rdir, d);
		VCOPY(r->ron, nrm);
		/* normalize(r->ron) called & r->rod set in o_mesh() */
	}
}


int
o_mesh(o, r)			/* compute ray intersection with a mesh */
OBJREC		*o;
register RAY	*r;
{
	RAY		rcont;
	int		flags;
	MESHVERT	tv[3];
	OBJECT		tmod;
	FLOAT		wt[3];
	int		i;
					/* get the mesh instance */
	prep_edge_cache(o);
					/* copy and transform ray */
	copystruct(&rcont, r);
	multp3(rcont.rorg, r->rorg, curmi->x.b.xfm);
	multv3(rcont.rdir, r->rdir, curmi->x.b.xfm);
	for (i = 0; i < 3; i++)
		rcont.rdir[i] /= curmi->x.b.sca;
	rcont.rmax *= curmi->x.b.sca;
					/* clear and trace ray */
	rayclear(&rcont);
	rcont.hitf = mesh_hit;
	if (!localhit(&rcont, &curmi->msh->mcube))
		return(0);			/* missed */
	if (rcont.rot * curmi->x.f.sca >= r->rot)
		return(0);			/* not close enough */
					/* transform ray back */
	r->rot = rcont.rot * curmi->x.f.sca;
	multp3(r->rop, rcont.rop, curmi->x.f.xfm);
	multv3(r->ron, rcont.ron, curmi->x.f.xfm);
	normalize(r->ron);
	r->rod = -DOT(r->rdir, r->ron);
					/* get triangle */
	flags = getmeshtri(tv, &tmod, curmsh, rcont.robj, MT_ALL);
	if (!(flags & MT_V))
		objerror(o, INTERNAL, "missing mesh vertices in o_mesh");
	r->robj = objndx(o);		/* set object and material */
	if (o->omod == OVOID && tmod != OVOID) {
		r->ro = getmeshpseudo(curmsh, tmod);
		r->rox = &curmi->x;
	} else
		r->ro = o;
					/* compute barycentric weights */
	if (flags & (MT_N|MT_UV))
		if (get_baryc(wt, rcont.rop, tv[0].v, tv[1].v, tv[2].v) < 0) {
			objerror(o, WARNING, "bad triangle in o_mesh");
			flags &= ~(MT_N|MT_UV);
		}
	if (flags & MT_N) {		/* interpolate normal */
		for (i = 0; i < 3; i++)
			rcont.pert[i] = wt[0]*tv[0].n[i] +
					wt[1]*tv[1].n[i] +
					wt[2]*tv[2].n[i];
		multv3(r->pert, rcont.pert, curmi->x.f.xfm);
		if (normalize(r->pert) != 0.0)
			for (i = 0; i < 3; i++)
				r->pert[i] -= r->ron[i];
	} else
		r->pert[0] = r->pert[1] = r->pert[2] = .0;

	if (flags & MT_UV)		/* interpolate uv coordinates */
		for (i = 0; i < 2; i++)
			r->uv[i] = wt[0]*tv[0].uv[i] +
					wt[1]*tv[1].uv[i] +
					wt[2]*tv[2].uv[i];
	else
		r->uv[0] = r->uv[1] = .0;

					/* return hit */
	return(1);
}
