#ifndef lint
static const char RCSid[] = "$Id: objtriangulate.c,v 2.2 2021/04/19 19:40:03 greg Exp $";
#endif
/*
 * Turn all faces with > 3 sides to triangles in Wavefront .OBJ scene
 */

#include "copyright.h"

#include <stdio.h>
#include <stdlib.h>
#include "rterror.h"
#include "objutil.h"
#include "triangulate.h"

typedef struct {
	Scene	*sc;
	Face	*f;
	int	rev;
} SceneFace;

/* Callback to create new triangle */
static int
addtriangle(const Vert2_list *tp, int a, int b, int c)
{
	SceneFace	*sf = (SceneFace *)tp->p;
	Face		*f = sf->f;
	VNDX		triv[3];

	if (sf->rev) {			/* need to reverse triangle? */
		int	t = a;
		a = c;
		c = t;
	}
	triv[0][0] = f->v[a].vid;
	triv[0][1] = f->v[a].tid;
	triv[0][2] = f->v[a].nid;
	triv[1][0] = f->v[b].vid;
	triv[1][1] = f->v[b].tid;
	triv[1][2] = f->v[b].nid;
	triv[2][0] = f->v[c].vid;
	triv[2][1] = f->v[c].tid;
	triv[2][2] = f->v[c].nid;

	f = addFace(sf->sc, triv, 3);
	if (f == NULL)
		return(0);

	f->flags |= sf->f->flags;
	f->grp = sf->f->grp;
	f->mat = sf->f->mat;
	return(1);
}

/* Callback to convert face to triangles */
static int
mktriangles(Scene *sc, Face *f, void *p)
{
	Normal		nrm;
	int		i, ax, ay;
	Vert2_list	*poly;
	SceneFace	mysf;

	if (f->nv <= 3)			/* already a triangle? */
		return(0);

	if (faceArea(sc, f, nrm) == 0)	/* hidden degenerate? */
		return(0);
					/* identify major axis from normal */
	i = (nrm[1]*nrm[1] > nrm[0]*nrm[0]);
	if (nrm[2]*nrm[2] > nrm[i]*nrm[i]) i = 2;
	ax = (i+1)%3;
	ay = (i+2)%3;
					/* assign 2-D vertices */
	poly = polyAlloc(f->nv);
	if (poly == NULL) {
		error(SYSTEM, "Out of memory in mktriangles");
		return(-1);
	}
	mysf.sc = sc;
	mysf.f = f;
	poly->p = &mysf;
	for (i = 0; i < f->nv; i++) {
		double	*pt = sc->vert[f->v[i].vid].p;
		poly->v[i].mX = pt[ax];
		poly->v[i].mY = pt[ay];
	}
	mysf.rev = (polyArea(poly) < .0);
	i = polyTriangulate(poly, addtriangle);
	polyFree(poly);
					/* flag face if replaced */
	f->flags |= (i > 0)*FACE_DUPLICATE;
	return(i);
}

/* Convert all faces with > 3 vertices to triangles */
int
triangulateScene(Scene *sc)
{
	int	orig_nfaces = sc->nfaces;
	int	nc = foreachFace(sc, mktriangles,
				0, FACE_DUPLICATE|FACE_DEGENERATE, NULL);

	if (nc < 0)			/* bad (memory) error? */
		return(-1);
	if (nc > 0)			/* else remove faces we split up */
		deleteFaces(sc, FACE_DUPLICATE, 0);

	return(sc->nfaces - orig_nfaces);
}
