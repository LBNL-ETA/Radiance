#ifndef lint
static const char RCSid[] = "$Id: cvmesh.c,v 2.9 2004/03/27 12:41:45 schorsch Exp $";
#endif
/*
 *  Radiance triangle mesh conversion routines
 */

#include "copyright.h"
#include "standard.h"
#include "cvmesh.h"
#include "otypes.h"
#include "face.h"
#include "tmesh.h"

/*
 * We need to divide faces into triangles and record auxiliary information
 * if given (surface normal and uv coordinates).  We do this by extending
 * the face structure linked to the OBJREC os member and putting our
 * auxiliary after it -- a bit sly, but it works.
 */

/* Auxiliary data for triangle */
typedef struct {
	int		fl;		/* flags of what we're storing */
	OBJECT		obj;		/* mesh triangle ID */
	FVECT		vn[3];		/* normals */
	RREAL		vc[3][2];	/* uv coords. */
} TRIDATA;

#define tdsize(fl)	((fl)&MT_UV ? sizeof(TRIDATA) : \
				(fl)&MT_N ? sizeof(TRIDATA)-6*sizeof(RREAL) : \
				sizeof(int)+sizeof(OBJECT))

#define	 OMARGIN	(10*FTINY)	/* margin around global cube */

MESH	*ourmesh = NULL;		/* our global mesh data structure */

FVECT	meshbounds[2];			/* mesh bounding box */

static void add2bounds(FVECT vp, RREAL vc[2]);
static OBJECT cvmeshtri(OBJECT obj);
static OCTREE cvmeshoct(OCTREE ot);



MESH *
cvinit(			/* initialize empty mesh */
	char	*nm
)
{
				/* free old mesh, first */
	if (ourmesh != NULL) {
		freemesh(ourmesh);
		ourmesh = NULL;
		freeobjects(0, nobjects);
		donesets();
	}
	if (nm == NULL)
		return(NULL);
	ourmesh = (MESH *)calloc(1, sizeof(MESH));
	if (ourmesh == NULL)
		goto nomem;
	ourmesh->name = savestr(nm);
	ourmesh->nref = 1;
	ourmesh->ldflags = 0;
	ourmesh->mcube.cutree = EMPTY;
	ourmesh->uvlim[0][0] = ourmesh->uvlim[0][1] = FHUGE;
	ourmesh->uvlim[1][0] = ourmesh->uvlim[1][1] = -FHUGE;
	meshbounds[0][0] = meshbounds[0][1] = meshbounds[0][2] = FHUGE;
	meshbounds[1][0] = meshbounds[1][1] = meshbounds[1][2] = -FHUGE;
	return(ourmesh);
nomem:
	error(SYSTEM, "out of memory in cvinit");
	return(NULL);
}


int
cvpoly(	/* convert a polygon to extended triangles */
	OBJECT	mo,
	int	n,
	FVECT	*vp,
	FVECT	*vn,
	RREAL	(*vc)[2]
)
{
	int	tcnt = 0;
	int	flags;
	RREAL	*tn[3], *tc[3];
	int	*ord;
	int	i, j;

	if (n < 3)		/* degenerate face */
		return(0);
	flags = MT_V;
	if (vn != NULL) {
		tn[0] = vn[0]; tn[1] = vn[1]; tn[2] = vn[2];
		flags |= MT_N;
	} else {
		tn[0] = tn[1] = tn[2] = NULL;
	}
	if (vc != NULL) {
		tc[0] = vc[0]; tc[1] = vc[1]; tc[2] = vc[2];
		flags |= MT_UV;
	} else {
		tc[0] = tc[1] = tc[2] = NULL;
	}
	if (n == 3)		/* output single triangle */
		return(cvtri(mo, vp[0], vp[1], vp[2],
				tn[0], tn[1], tn[2],
				tc[0], tc[1], tc[2]));

				/* decimate polygon (assumes convex) */
	ord = (int *)malloc(n*sizeof(int));
	if (ord == NULL)
		error(SYSTEM, "out of memory in cvpoly");
	for (i = n; i--; )
		ord[i] = i;
	while (n >= 3) {
		if (flags & MT_N)
			for (i = 3; i--; )
				tn[i] = vn[ord[i]];
		if (flags & MT_UV)
			for (i = 3; i--; )
				tc[i] = vc[ord[i]];
		tcnt += cvtri(mo, vp[ord[0]], vp[ord[1]], vp[ord[2]],
				tn[0], tn[1], tn[2],
				tc[0], tc[1], tc[2]);
			/* remove vertex and rotate */
		n--;
		j = ord[0];
		for (i = 0; i < n-1; i++)
			ord[i] = ord[i+2];
		ord[i] = j;
	}
	free((void *)ord);
	return(tcnt);
}


static void
add2bounds(		/* add point and uv coordinate to bounds */
	FVECT	vp,
	RREAL	vc[2]
)
{
	register int	j;

	for (j = 3; j--; ) {
		if (vp[j] < meshbounds[0][j])
			meshbounds[0][j] = vp[j];
		if (vp[j] > meshbounds[1][j])
			meshbounds[1][j] = vp[j];
	}
	if (vc == NULL)
		return;
	for (j = 2; j--; ) {
		if (vc[j] < ourmesh->uvlim[0][j])
			ourmesh->uvlim[0][j] = vc[j];
		if (vc[j] > ourmesh->uvlim[1][j])
			ourmesh->uvlim[1][j] = vc[j];
	}
}


int				/* create an extended triangle */
cvtri(
	OBJECT	mo,
	FVECT	vp1,
	FVECT	vp2,
	FVECT	vp3,
	FVECT	vn1,
	FVECT	vn2,
	FVECT	vn3,
	RREAL	vc1[2],
	RREAL	vc2[2],
	RREAL	vc3[2]
)
{
	static OBJECT	fobj = OVOID;
	char		buf[32];
	int		flags;
	TRIDATA		*ts;
	FACE		*f;
	OBJREC		*fop;
	int		j;
	
	flags = MT_V;		/* check what we have */
	if (vn1 != NULL && vn2 != NULL && vn3 != NULL) {
		RREAL	*rp;
		switch (flat_tri(vp1, vp2, vp3, vn1, vn2, vn3)) {
		case ISBENT:
			flags |= MT_N;
			/* fall through */
		case ISFLAT:
			break;
		case RVBENT:
			flags |= MT_N;
			rp = vn1; vn1 = vn3; vn3 = rp;
			/* fall through */
		case RVFLAT:
			rp = vp1; vp1 = vp3; vp3 = rp;
			rp = vc1; vc1 = vc3; vc3 = rp;
			break;
		case DEGEN:
			error(WARNING, "degenerate triangle");
			return(0);
		default:
			error(INTERNAL, "bad return from flat_tri()");
		}
	}
	if (vc1 != NULL && vc2 != NULL && vc3 != NULL)
		flags |= MT_UV;
	if (fobj == OVOID) {	/* create new triangle object */
		fobj = newobject();
		if (fobj == OVOID)
			goto nomem;
		fop = objptr(fobj);
		fop->omod = mo;
		fop->otype = OBJ_FACE;
		sprintf(buf, "t%ld", fobj);
		fop->oname = savqstr(buf);
		fop->oargs.nfargs = 9;
		fop->oargs.farg = (RREAL *)malloc(9*sizeof(RREAL));
		if (fop->oargs.farg == NULL)
			goto nomem;
	} else {		/* else reuse failed one */
		fop = objptr(fobj);
		if (fop->otype != OBJ_FACE || fop->oargs.nfargs != 9)
			error(CONSISTENCY, "code error 1 in cvtri");
	}
	for (j = 3; j--; ) {
		fop->oargs.farg[j] = vp1[j];
		fop->oargs.farg[3+j] = vp2[j];
		fop->oargs.farg[6+j] = vp3[j];
	}
				/* create face record */
	f = getface(fop);
	if (f->area == 0.) {
		free_os(fop);
		return(0);
	}
	if (fop->os != (char *)f)
		error(CONSISTENCY, "code error 2 in cvtri");
				/* follow with auxliary data */
	f = (FACE *)realloc((void *)f, sizeof(FACE)+tdsize(flags));
	if (f == NULL)
		goto nomem;
	fop->os = (char *)f;
	ts = (TRIDATA *)(f+1);
	ts->fl = flags;
	ts->obj = OVOID;
	if (flags & MT_N)
		for (j = 3; j--; ) {
			ts->vn[0][j] = vn1[j];
			ts->vn[1][j] = vn2[j];
			ts->vn[2][j] = vn3[j];
		}
	if (flags & MT_UV)
		for (j = 2; j--; ) {
			ts->vc[0][j] = vc1[j];
			ts->vc[1][j] = vc2[j];
			ts->vc[2][j] = vc3[j];
		}
	else
		vc1 = vc2 = vc3 = NULL;
				/* update bounds */
	add2bounds(vp1, vc1);
	add2bounds(vp2, vc2);
	add2bounds(vp3, vc3);
	fobj = OVOID;		/* we used this one */
	return(1);
nomem:
	error(SYSTEM, "out of memory in cvtri");
	return(0);
}


static OBJECT
cvmeshtri(			/* add an extended triangle to our mesh */
	OBJECT	obj
)
{
	OBJREC		*o = objptr(obj);
	TRIDATA		*ts;
	MESHVERT	vert[3];
	int		i, j;
	
	if (o->otype != OBJ_FACE)
		error(CONSISTENCY, "non-face in mesh");
	if (o->oargs.nfargs != 9)
		error(CONSISTENCY, "non-triangle in mesh");
	if (o->os == NULL)
		error(CONSISTENCY, "missing face record in cvmeshtri");
	ts = (TRIDATA *)((FACE *)o->os + 1);
	if (ts->obj != OVOID)	/* already added? */
		return(ts->obj);
	vert[0].fl = vert[1].fl = vert[2].fl = ts->fl;
	for (i = 3; i--; )
		for (j = 3; j--; )
			vert[i].v[j] = o->oargs.farg[3*i+j];
	if (ts->fl & MT_N)
		for (i = 3; i--; )
			for (j = 3; j--; )
				vert[i].n[j] = ts->vn[i][j];
	if (ts->fl & MT_UV)
		for (i = 3; i--; )
			for (j = 2; j--; )
				vert[i].uv[j] = ts->vc[i][j];
	ts->obj = addmeshtri(ourmesh, vert, o->omod);
	if (ts->obj == OVOID)
		error(INTERNAL, "addmeshtri failed");
	return(ts->obj);
}


void
cvmeshbounds(void)			/* set mesh boundaries */
{
	int	i;

	if (ourmesh == NULL)
		return;
				/* fix coordinate bounds */
	for (i = 0; i < 3; i++) {
		if (meshbounds[0][i] > meshbounds[1][i])
			error(USER, "no polygons in mesh");
		meshbounds[0][i] -= OMARGIN;
		meshbounds[1][i] += OMARGIN;
		if (meshbounds[1][i]-meshbounds[0][i] > ourmesh->mcube.cusize)
			ourmesh->mcube.cusize = meshbounds[1][i] -
						meshbounds[0][i];
	}
	for (i = 0; i < 3; i++)
		ourmesh->mcube.cuorg[i] = (meshbounds[1][i]+meshbounds[0][i] -
						ourmesh->mcube.cusize)*.5;
	if (ourmesh->uvlim[0][0] > ourmesh->uvlim[1][0]) {
		ourmesh->uvlim[0][0] = ourmesh->uvlim[0][1] = 0.;
		ourmesh->uvlim[1][0] = ourmesh->uvlim[1][1] = 0.;
	} else {
		for (i = 0; i < 2; i++) {
			double	marg;		/* expand past endpoints */
			marg = (2./(1L<<(8*sizeof(uint16)))) *
					(ourmesh->uvlim[1][i] -
					 ourmesh->uvlim[0][i]);
			ourmesh->uvlim[0][i] -= marg;
			ourmesh->uvlim[1][i] += marg;
		}
	}
	ourmesh->ldflags |= IO_BOUNDS;
}


static OCTREE
cvmeshoct(			/* convert triangles in subtree */
	OCTREE	ot
)
{
	int	i;

	if (isempty(ot))
		return(EMPTY);

	if (isfull(ot)) {
		OBJECT	oset1[MAXSET+1];
		OBJECT	oset2[MAXSET+1];
		objset(oset1, ot);
		oset2[0] = 0;
		for (i = oset1[0]; i > 0; i--)
			insertelem(oset2, cvmeshtri(oset1[i]));
		return(fullnode(oset2));
	}

	for (i = 8; i--; )
		octkid(ot, i) = cvmeshoct(octkid(ot, i));
	return(ot);
}


MESH *
cvmesh(void)			/* convert mesh and octree leaf nodes */
{
	if (ourmesh == NULL)
		return(NULL);
				/* convert triangles in octree nodes */
	ourmesh->mcube.cutree = cvmeshoct(ourmesh->mcube.cutree);
	ourmesh->ldflags |= IO_SCENE|IO_TREE;

	return(ourmesh);
}
