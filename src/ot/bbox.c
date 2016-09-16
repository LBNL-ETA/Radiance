#ifndef lint
static const char RCSid[] = "$Id: bbox.c,v 2.7 2016/09/16 15:09:21 greg Exp $";
#endif
/*
 *  bbox.c - routines for bounding box computation.
 */

#include  "copyright.h"

#include  "standard.h"
#include  "object.h"
#include  "octree.h"
#include  "otypes.h"
#include  "face.h"
#include  "cone.h"
#include  "instance.h"
#include  "mesh.h"
#include  "oconv.h"


static void point2bbox(FVECT  p, FVECT  bbmin, FVECT  bbmax);
static void circle2bbox(FVECT  cent, FVECT  norm, double  rad, FVECT  bbmin, FVECT  bbmax);


void
add2bbox(		/* expand bounding box to fit object */
	register OBJREC  *o,
	FVECT  bbmin,
	FVECT  bbmax
)
{
	CONE  *co;
	FACE  *fo;
	INSTANCE  *io;
	MESHINST  *mi;
	FVECT  v;
	register int  i, j;

	switch (o->otype) {
	case OBJ_SPHERE:
	case OBJ_BUBBLE:
		if (o->oargs.nfargs != 4)
			objerror(o, USER, "bad arguments");
		for (i = 0; i < 3; i++) {
			VCOPY(v, o->oargs.farg);
			v[i] -= o->oargs.farg[3];
			point2bbox(v, bbmin, bbmax);
			v[i] += 2.0 * o->oargs.farg[3];
			point2bbox(v, bbmin, bbmax);
		}
		break;
	case OBJ_FACE:
		fo = getface(o);
		j = fo->nv;
		while (j--)
			point2bbox(VERTEX(fo,j), bbmin, bbmax);
		break;
	case OBJ_CONE:
	case OBJ_CUP:
	case OBJ_CYLINDER:
	case OBJ_TUBE:
	case OBJ_RING:
		co = getcone(o, 0);
		if (co == NULL)
			break;
		if (o->otype != OBJ_RING)
			circle2bbox(CO_P0(co), co->ad, CO_R0(co), bbmin, bbmax);
		circle2bbox(CO_P1(co), co->ad, CO_R1(co), bbmin, bbmax);
		break;
	case OBJ_INSTANCE:
		io = getinstance(o, IO_BOUNDS);
		for (j = 0; j < 8; j++) {
			for (i = 0; i < 3; i++) {
				v[i] = io->obj->scube.cuorg[i];
				if (j & 1<<i)
					v[i] += io->obj->scube.cusize;
			}
			multp3(v, v, io->x.f.xfm);
			point2bbox(v, bbmin, bbmax);
		}
		break;
	case OBJ_MESH:
		mi = getmeshinst(o, IO_BOUNDS);
		for (j = 0; j < 8; j++) {
			for (i = 0; i < 3; i++) {
				v[i] = mi->msh->mcube.cuorg[i];
				if (j & 1<<i)
					v[i] += mi->msh->mcube.cusize;
			}
			multp3(v, v, mi->x.f.xfm);
			point2bbox(v, bbmin, bbmax);
		}
		break;
	}
}


static void
point2bbox(		/* expand bounding box to fit point */
	register FVECT  p,
	register FVECT  bbmin,
	register FVECT  bbmax
)
{
	register int  i;

	for (i = 0; i < 3; i++) {
		if (p[i] < bbmin[i])
			bbmin[i] = p[i];
		if (p[i] > bbmax[i])
			bbmax[i] = p[i];
	}
}


static void
circle2bbox(	/* expand bbox to fit circle */
	FVECT  cent,
	FVECT  norm,
	double  rad,
	FVECT  bbmin,
	FVECT  bbmax
)
{
	double	d, r;
	register int  i;

	for (i = 0; i < 3; i++) {
		r = sqrt(1. - norm[i]*norm[i]);
		d = cent[i] + r*rad;
		if (d > bbmax[i])
			bbmax[i] = d;
		d = cent[i] - r*rad;
		if (d < bbmin[i])
			bbmin[i] = d;
	}
}
