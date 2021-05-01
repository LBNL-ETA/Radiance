#ifndef lint
static const char RCSid[] = "$Id: objutil.c,v 2.20 2021/05/01 00:58:18 greg Exp $";
#endif
/*
 *  Basic .OBJ scene handling routines.
 *
 *  Created by Greg Ward on Wed Feb 11 2004.
 */

#include <stdlib.h>
#include <ctype.h>
#include "rtio.h"
#include "rtmath.h"
#include "rterror.h"
#include "lookup.h"
#include "objutil.h"

/* Find an existing name in a list of names */
int
findName(const char *nm, const char **nmlist, int n)
{
	int	i;
	
	for (i = n; i-- > 0; )
		if (!strcmp(nmlist[i], nm))
			break;
	return(i);
}

/* Clear face selection */
void
clearSelection(Scene *sc, int set)
{
	Face	*f;
	
	for (f = sc->flist; f != NULL; f = f->next)
		if (set)
			f->flags |= FACE_SELECTED;
		else
			f->flags &= ~FACE_SELECTED;
}

/* Select faces by object name (modifies current) */
void
selectGroup(Scene *sc, const char *gname, int invert)
{
	int	gid = getGroupID(sc, gname);
	Face	*f;

	if (gid < 0)
		return;
	for (f = sc->flist; f != NULL; f = f->next)
		if (f->grp == gid) {
			if (invert)
				f->flags &= ~FACE_SELECTED;
			else
				f->flags |= FACE_SELECTED;
		}
}

/* Select faces by material name (modifies current) */
void
selectMaterial(Scene *sc, const char *mname, int invert)
{
	int	mid = getMaterialID(sc, mname);
	Face	*f;
	
	if (mid < 0)
		return;
	for (f = sc->flist; f != NULL; f = f->next)
		if (f->mat == mid) {
			if (invert)
				f->flags &= ~FACE_SELECTED;
			else
				f->flags |= FACE_SELECTED;
		}
}

/* Invert face selection */
void
invertSelection(Scene *sc)
{
	Face	*f;
	
	for (f = sc->flist; f != NULL; f = f->next)
		f->flags ^= FACE_SELECTED;
}

/* Count selected faces */
int
numberSelected(Scene *sc)
{
	int	nsel = 0;
	Face	*f;
	
	for (f = sc->flist; f != NULL; f = f->next)
		nsel += ((f->flags & FACE_SELECTED) != 0);
	return(nsel);
}

/* Execute callback on indicated faces */
int
foreachFace(Scene *sc, int (*cb)(Scene *, Face *, void *),
			int flreq, int flexc, void *c_data)
{
	const int	fltest = flreq | flexc;
	int		sum = 0;
	Face		*f;
	
	if ((sc == NULL) | (cb == NULL))
		return(0);
	for (f = sc->flist; f != NULL; f = f->next)
		if ((f->flags & fltest) == flreq) {
			int     res = (*cb)(sc, f, c_data);
			if (res < 0)
				return(res);
			sum += res;
		}
	return(sum);
}

/* Callback for removeTexture() */
static int
remFaceTexture(Scene *sc, Face *f, void *dummy)
{
	int	hadTexture = 0;
	int	i;
	
	for (i = f->nv; i-- > 0; ) {
		if (f->v[i].tid < 0)
			continue;
		f->v[i].tid = -1;
		hadTexture = 1;
	}
	return(hadTexture);
}

/* Remove texture coordinates from the indicated faces */
int
removeTexture(Scene *sc, int flreq, int flexc)
{
	return(foreachFace(sc, remFaceTexture, flreq, flexc, NULL));
}

/* Callback for removeNormals() */
static int
remFaceNormal(Scene *sc, Face *f, void *dummy)
{
	int	hadNormal = 0;
	int	i;
	
	for (i = f->nv; i-- > 0; ) {
		if (f->v[i].nid < 0)
			continue;
		f->v[i].nid = -1;
		hadNormal = 1;
	}
	return(hadNormal);
}

/* Remove surface normals from the indicated faces */
int
removeNormals(Scene *sc, int flreq, int flexc)
{
	return(foreachFace(sc, remFaceNormal, flreq, flexc, NULL));
}

/* Callback for changeGroup() */
static int
chngFaceGroup(Scene *sc, Face *f, void *ptr)
{
	int	grp = *(int *)ptr;

	if (f->grp == grp)
		return(0);
	f->grp = grp;
	return(1);
}

/* Change group for the indicated faces */
int
changeGroup(Scene *sc, const char *gname, int flreq, int flexc)
{
	int	grp = getGroupID(sc, gname);

	if (grp < 0) {
		sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
		sc->grpname[grp=sc->ngrps++] = savqstr((char *)gname);
	}
	return(foreachFace(sc, chngFaceGroup, flreq, flexc, (void *)&grp));
}

/* Callback for changeMaterial() */
static int
chngFaceMaterial(Scene *sc, Face *f, void *ptr)
{
	int	mat = *(int *)ptr;

	if (f->mat == mat)
		return(0);
	f->mat = mat;
	return(1);
}

/* Change material for the indicated faces */
int
changeMaterial(Scene *sc, const char *mname, int flreq, int flexc)
{
	int	mat = getMaterialID(sc, mname);

	if (mat < 0) {
		sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
		sc->matname[mat=sc->nmats++] = savqstr((char *)mname);
	}
	return(foreachFace(sc, chngFaceMaterial, flreq, flexc, (void *)&mat));
}

/* Compare floating point values for (near) equality */
static int
fdiffer(double f1, double f2, double eps)
{
	if (f2 != .0)
		f1 = f1/f2 - 1.;
	return((f1 > eps) | (f1 < -eps));
}

/* Compare two texture coordinates for (near) equality */
static int
tex_diff(const TexCoord *t1, const TexCoord *t2, double eps)
{
	if (fdiffer(t1->u, t2->u, eps))
		return(1);
	if (fdiffer(t1->v, t2->v, eps))
		return(1);
	return(0);
}

/* Compare two surface normals for (near) equality */
static int
norm_diff(const Normal n1, const Normal n2, double eps)
{
	if (fabs(n1[0]-n2[0]) > eps)
		return(1);
	if (fabs(n1[1]-n2[1]) > eps)
		return(1);
	if (fabs(n1[2]-n2[2]) > eps)
		return(1);
	return(0);
}

/* Replace the given vertex with an equivalent one */
static int
replace_vertex(Scene *sc, int prev, int repl, double eps)
{
	int	repl_tex[10];
	int	repl_norm[10];
	Face	*f, *flast=NULL;
	int	i, j=0;
					/* check to see if it's even used */
	if (sc->vert[prev].vflist == NULL)
		return(0);
					/* get replacement textures */
	repl_tex[0] = -1;
	for (f = sc->vert[repl].vflist; f != NULL; f = f->v[j].fnext) {
					/* make sure prev isn't in there */
		for (j = f->nv; j-- > 0; )
			if (f->v[j].vid == prev)
				return(0);
		for (j = f->nv; j-- > 0; )
			if (f->v[j].vid == repl)
				break;
		if (j < 0)
			goto linkerr;
		if (f->v[j].tid < 0)
			continue;
					/* see if it's new */
		for (i = 0; repl_tex[i] >= 0; i++) {
			if (repl_tex[i] == f->v[j].tid)
				break;
			if (!tex_diff(&sc->tex[repl_tex[i]],
					&sc->tex[f->v[j].tid], eps)) {
				f->v[j].tid = repl_tex[i];
				break;  /* consolidate */
			}
		}
		if (repl_tex[i] >= 0)
			continue;       /* have this one already */
					/* else add it */
		repl_tex[i++] = f->v[j].tid;
		repl_tex[i] = -1;
		if (i >= 9)
			break;		/* that's all we have room for */
	}
					/* get replacement normals */
	repl_norm[0] = -1;
	for (f = sc->vert[repl].vflist; f != NULL; f = f->v[j].fnext) {
		for (j = f->nv; j-- > 0; )
			if (f->v[j].vid == repl)
				break;
		if (f->v[j].nid < 0)
			continue;
					/* see if it's new */
		for (i = 0; repl_norm[i] >= 0; i++) {
			if (repl_norm[i] == f->v[j].nid)
				break;
			if (!norm_diff(sc->norm[repl_norm[i]],
					sc->norm[f->v[j].nid], eps)) {
				f->v[j].nid = repl_norm[i];
				break;  /* consolidate */
			}
		}
		if (repl_norm[i] >= 0)
			continue;       /* have this one already */
					/* else add it */
		repl_norm[i++] = f->v[j].nid;
		repl_norm[i] = -1;
		if (i >= 9)
			break;		/* that's all we have room for */
	}
					/* replace occurrences of vertex */
	for (f = sc->vert[prev].vflist; f != NULL; f = f->v[j].fnext) {
		for (j = f->nv; j-- > 0; )
			if (f->v[j].vid == prev)
				break;
		if (j < 0)
			goto linkerr;
		/* XXX doesn't allow for multiple references to prev in face */
		f->v[j].vid = repl;     /* replace vertex itself */
		if (faceArea(sc, f, NULL) <= FTINY*FTINY)
			f->flags |= FACE_DEGENERATE;
		if (f->v[j].tid >= 0)   /* replace texture if appropriate */
			for (i = 0; repl_tex[i] >= 0; i++) {
				if (repl_tex[i] == f->v[j].tid)
					break;
				if (!tex_diff(&sc->tex[repl_tex[i]],
						&sc->tex[f->v[j].tid], eps)) {
					f->v[j].tid = repl_tex[i];
					break;
				}
			}
		if (f->v[j].nid >= 0)   /* replace normal if appropriate */
			for (i = 0; repl_norm[i] >= 0; i++) {
				if (repl_norm[i] == f->v[j].nid)
					break;
				if (!norm_diff(sc->norm[repl_norm[i]],
						sc->norm[f->v[j].nid], eps)) {
					f->v[j].nid = repl_norm[i];
					break;
				}
			}
		flast = f;
	}
					/* transfer face list */
	flast->v[j].fnext = sc->vert[repl].vflist;
	sc->vert[repl].vflist = sc->vert[prev].vflist;
	sc->vert[prev].vflist = NULL;
	return(1);
linkerr:
	error(CONSISTENCY, "Link error in replace_vertex()");
	return(0);			/* shouldn't return */
}

/* Eliminate duplicate vertices, return # joined */
int
coalesceVertices(Scene *sc, double eps)
{
	int	nelim = 0;
	LUTAB	myLookup;
	LUENT	*le;
	char	vertfmt[32], vertbuf[64];
	double	d;
	int	i;
	
	if (eps >= 1.)
		return(0);
	if (sc->nverts <= 3)
		return(0);
					/* create hash table */
	myLookup.hashf = lu_shash;
	myLookup.keycmp = (lut_keycmpf_t *)strcmp;
	myLookup.freek = (lut_free_t *)freeqstr;
	myLookup.freed = NULL;
	if (!lu_init(&myLookup, sc->nverts))
		return(0);
	if (eps <= 5e-15)
		strcpy(vertfmt, "%.15e %.15e %.15e");
	else {
		int     nsigdig = 0;
		for (d = eps; d < 0.5; d *= 10.)
			++nsigdig;
		sprintf(vertfmt, "%%.%de %%.%de %%.%de",
				nsigdig, nsigdig, nsigdig);
	}
					/* find coicident vertices */
	for (i = 0; i < sc->nverts; i++) {
		if (verbose && !((i+1) & 0x3fff))
			fprintf(stderr, "%3d%% complete\r", 100*i/sc->nverts);
					/* check for match */
		sprintf(vertbuf, vertfmt, sc->vert[i].p[0],
				sc->vert[i].p[1], sc->vert[i].p[2]);
		le = lu_find(&myLookup, vertbuf);
		if (le->data != NULL) { /* coincident vertex */
			nelim += replace_vertex(sc, i,
					(Vertex *)le->data - sc->vert, eps);
			continue;
		}
		if (le->key == NULL)    /* else create new entry */
			le->key = savqstr(vertbuf);
		le->data = (char *)&sc->vert[i];
	}
	lu_done(&myLookup);		/* clean up */
	return(nelim);
}

/* Identify duplicate faces */
int
findDuplicateFaces(Scene *sc)
{
	int	nchecked = 0;
	int	nfound = 0;
	Face	*f, *f1;
	int	vid;
	int	i, j;
						/* start fresh */
	for (f = sc->flist; f != NULL; f = f->next)
		f->flags &= ~FACE_DUPLICATE;
						/* check each face */
	for (f = sc->flist; f != NULL; f = f->next) {
		nchecked++;
		if (verbose && !(nchecked & 0x3fff))
			fprintf(stderr, "%3d%% complete\r",
					100*nchecked/sc->nfaces);
		if (f->flags & FACE_DUPLICATE)
			continue;		/* already identified */
		vid = f->v[0].vid;
						/* look for duplicates */
		for (f1 = sc->vert[vid].vflist; f1 != NULL;
					f1 = f1->v[j].fnext) {
			for (j = f1->nv; j-- > 0; )
				if (f1->v[j].vid == vid)
					break;
			if (j < 0)
				break;		/* missing link! */
			if (f1 == f)
				continue;       /* shouldn't happen */
			if (f1->flags & FACE_DUPLICATE)
				continue;       /* already marked */
			if (f1->nv != f->nv)
				continue;       /* couldn't be dup. */
			for (i = f->nv; --i > 0; )
				if (f->v[i].vid != f1->v[(j+i)%f1->nv].vid)
					break;  /* vertex mismatch */
			if (i) {
#if DUP_CHECK_REVERSE				/* check reverse direction */
				for (i = f->nv; --i > 0; )
					if (f1->v[(j+f1->nv-i)%f1->nv].vid
							!= f->v[i].vid)
						break;
				if (i)		/* no match */
#endif
					continue;
			}
			f1->flags |= FACE_DUPLICATE;
			++nfound;
		}
	}
	if (verbose)
		fprintf(stderr, "Found %d duplicate faces\n", nfound);
	return(nfound);
}

/* Delete indicated faces */
int
deleteFaces(Scene *sc, int flreq, int flexc)
{
	const int	fltest = flreq | flexc;
	int		orig_nfaces = sc->nfaces;
	Face		fhead;
	Face		*f, *ftst;
	
	fhead.next = sc->flist;
	f = &fhead;
	while ((ftst = f->next) != NULL)
		if ((ftst->flags & fltest) == flreq) {
			Face    *vf;		/* remove from vertex lists */
			int     vid, i, j;
			for (i = 0; i < ftst->nv; i++) {
				vid = ftst->v[i].vid;
				vf = sc->vert[vid].vflist;
				if (vf == ftst) {
					sc->vert[vid].vflist = ftst->v[i].fnext;
					continue;
				}
				while (vf != NULL) {
					for (j = vf->nv; j-- > 0; )
						if (vf->v[j].vid == vid)
							break;
					if (j < 0)
						break;  /* error */
					if (vf->v[j].fnext == ftst) {
						vf->v[j].fnext =
							ftst->v[i].fnext;
						break;
					}
					vf = vf->v[j].fnext;
				}
			}
			f->next = ftst->next;   /* remove from scene list */
			efree((char *)ftst);
			sc->nfaces--;
		} else
			f = f->next;
	sc->flist = fhead.next;
	return(orig_nfaces - sc->nfaces);
}

/* Compute face area (and normal) */
double
faceArea(const Scene *sc, const Face *f, Normal nrm)
{
	FVECT   fnrm;
	double  area;
	FVECT   v1, v2, v3;
	double  *p0;
	int     i;
	
	if (f->flags & FACE_DEGENERATE)
		return(.0);			/* should we check this? */
	fnrm[0] = fnrm[1] = fnrm[2] = .0;
	p0 = sc->vert[f->v[0].vid].p;
	VSUB(v1, sc->vert[f->v[1].vid].p, p0);
	for (i = 2; i < f->nv; i++) {
		VSUB(v2, sc->vert[f->v[i].vid].p, p0);
		fcross(v3, v1, v2);
		VADD(fnrm, fnrm, v3);
		VCOPY(v1, v2);
	}
	area = 0.5*normalize(fnrm);
	if (nrm != NULL) {
		if (area != 0.)
			VCOPY(nrm, fnrm);
		else
			nrm[0] = nrm[1] = nrm[2] = 0.;
	}
	return(area);
}

/* Add a descriptive comment */
void
addComment(Scene *sc, const char *comment)
{
	sc->descr = chunk_alloc(char *, sc->descr, sc->ndescr);
	sc->descr[sc->ndescr++] = savqstr((char *)comment);
}

/* Find index for comment containing the given string (starting from n) */
int
findComment(Scene *sc, const char *match, int n)
{
	if (n >= sc->ndescr)
		return(-1);
	n *= (n > 0);
	while (n < sc->ndescr)
		if (strstr(sc->descr[n], match) != NULL)
			return(n);
	return(-1);
}

/* Clear comments */
void
clearComments(Scene *sc)
{
	while (sc->ndescr > 0)
		freeqstr(sc->descr[--sc->ndescr]);
	efree((char *)sc->descr);
	sc->ndescr = 0;
}

/* Add a face to a vertex face list */
static void
add2facelist(Scene *sc, Face *f, int i)
{
	int     vid = f->v[i].vid;
	Face    *fp = sc->vert[vid].vflist;
	int     j;
	
	f->v[i].fnext = NULL;			/* will put at end */
	if (fp == NULL) {			/* new list */
		sc->vert[vid].vflist = f;
		return;
	}
	for ( ; ; ) {				/* else find position */
		if (fp == f)
			return;			/* already in list */
		for (j = fp->nv; j-- > 0; )
			if (fp->v[j].vid == vid)
				break;
		if (j < 0)
			error(CONSISTENCY, "Link error in add2facelist()");
		if (fp->v[j].fnext == NULL)
			break;			/* reached the end */
		fp = fp->v[j].fnext;
	}
	fp->v[j].fnext = f;			/* append new face */
}

/* Add a vertex to a scene */
int
addVertex(Scene *sc, double x, double y, double z)
{
	sc->vert = chunk_alloc(Vertex, sc->vert, sc->nverts);
	sc->vert[sc->nverts].p[0] = x;
	sc->vert[sc->nverts].p[1] = y;
	sc->vert[sc->nverts].p[2] = z;
	sc->vert[sc->nverts].vflist = NULL;
	return(sc->nverts++);
}

/* Add a texture coordinate to a scene */
int
addTexture(Scene *sc, double u, double v)
{
	sc->tex = chunk_alloc(TexCoord, sc->tex, sc->ntex);
	sc->tex[sc->ntex].u = u;
	sc->tex[sc->ntex].v = v;
	return(sc->ntex++);
}

/* Add a surface normal to a scene */
int
addNormal(Scene *sc, double xn, double yn, double zn)
{
	FVECT	nrm;

	nrm[0] = xn; nrm[1] = yn; nrm[2] = zn;
	if (normalize(nrm) == .0)
		return(-1);
	sc->norm = chunk_alloc(Normal, sc->norm, sc->nnorms);
	VCOPY(sc->norm[sc->nnorms], nrm);
	return(sc->nnorms++);
}

/* Set current (last) group */
void
setGroup(Scene *sc, const char *nm)
{
	sc->lastgrp = getGroupID(sc, nm);
	if (sc->lastgrp >= 0)
		return;
	sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
	sc->grpname[sc->lastgrp=sc->ngrps++] = savqstr((char *)nm);
}

/* Set current (last) material */
void
setMaterial(Scene *sc, const char *nm)
{
	sc->lastmat = getMaterialID(sc, nm);
	if (sc->lastmat >= 0)
		return;
	sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
	sc->matname[sc->lastmat=sc->nmats++] = savqstr((char *)nm);
}

/* Add new face to a scene */
Face *
addFace(Scene *sc, VNDX vid[], int nv)
{
	Face	*f;
	int	i;
	
	if (nv < 3)
		return(NULL);
	f = (Face *)emalloc(sizeof(Face)+sizeof(VertEnt)*(nv-3));
	f->flags = 0;
	f->nv = nv;
	f->grp = sc->lastgrp;
	f->mat = sc->lastmat;
	for (i = 0; i < nv; i++) {		/* add each vertex */
		int     j;
		f->v[i].vid = vid[i][0];
		f->v[i].tid = vid[i][1];
		f->v[i].nid = vid[i][2];
		f->v[i].fnext = NULL;
		for (j = i; j-- > 0; )
			if (f->v[j].vid == vid[i][0])
				break;
		if (j < 0) {			/* first occurrence? */
			f->v[i].fnext = sc->vert[vid[i][0]].vflist;
			sc->vert[vid[i][0]].vflist = f;
		} else if (nv == 3)		/* degenerate triangle? */
			f->flags |= FACE_DEGENERATE;
	}
	f->next = sc->flist;			/* push onto face list */
	sc->flist = f;
	sc->nfaces++;
						/* check face area */
	if (!(f->flags & FACE_DEGENERATE) && faceArea(sc, f, NULL) <= FTINY*FTINY)
		f->flags |= FACE_DEGENERATE;
	return(f);
}

/* Get neighbor vertices: malloc array with valence in index[0] */
int *
getVertexNeighbors(Scene *sc, int vid)
{
	int	*varr;
	Face	*f;
	int	j=0;

	if (sc == NULL || (vid < 0) | (vid >= sc->nverts) ||
			(f = sc->vert[vid].vflist) == NULL)
		return(NULL);

	varr = (int *)emalloc(sizeof(int)*CHUNKSIZ);
	varr[0] = 0;				/* add unique neighbors */
	for ( ; f != NULL; f = f->v[j].fnext) {
		int	i, cvi;
		for (j = f->nv; j--; )		/* find ourself in poly */
			if (f->v[j].vid == vid)
				break;
		if (j < 0)			/* this is an error */
			break;
		if (f->nv < 3)			/* also never happens */
			continue;
						/* check previous neighbor */
		cvi = f->v[ (j > 0) ? j-1 : f->nv-1 ].vid;
		for (i = varr[0]+1; --i; )	/* making sure not in list */
			if (varr[i] == cvi)
				break;
		if (!i) {			/* new neighbor? */
			varr = chunk_alloc(int, varr, varr[0]+1);
			varr[++varr[0]] = cvi;
		}
						/* check next neighbor */
		cvi = f->v[ (j < f->nv-1) ? j+1 : 0 ].vid;
		for (i = varr[0]+1; --i; )
			if (varr[i] == cvi)
				break;
		if (!i) {
			varr = chunk_alloc(int, varr, varr[0]+1);
			varr[++varr[0]] = cvi;
		}
	}
	if (!varr[0]) {
		efree((char *)varr);		/* something went awry */
		return(NULL);
	}
						/* tighten array & return */
	return((int *)erealloc((char *)varr, sizeof(int)*(varr[0]+1)));
}

/* Callback for growBoundingBox() */
static int
addBBox(Scene *sc, Face *f, void *p)
{
	double	(*bbox)[3] = (double (*)[3])p;
	int	i, j;

	for (i = f->nv; i-- > 0; ) {
		double	*p3 = sc->vert[f->v[i].vid].p;
		for (j = 3; j--; ) {
			if (p3[j] < bbox[0][j])
				bbox[0][j] = p3[j];
			if (p3[j] > bbox[1][j])
				bbox[1][j] = p3[j];
		}
	}
	return(1);
}

/* Expand bounding box min & max (initialize bbox to all zeroes) */
int
growBoundingBox(Scene *sc, double bbox[2][3], int flreq, int flexc)
{
	if (sc == NULL || sc->nfaces <= 0 || bbox == NULL)
		return(0);

	if (VABSEQ(bbox[0], bbox[1])) {		/* first run */
		bbox[0][0] = bbox[0][1] = bbox[0][2] = FHUGE;
		bbox[1][0] = bbox[1][1] = bbox[1][2] = -FHUGE;
	}
	return(foreachFace(sc, addBBox, flreq, flexc, bbox));
}

/* Allocate an empty scene */
Scene *
newScene(void)
{
	Scene	*sc = (Scene *)ecalloc(1, sizeof(Scene));
						/* default group & material */
	sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
	sc->grpname[sc->ngrps++] = savqstr("DEFAULT_GROUP");
	sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
	sc->matname[sc->nmats++] = savqstr("DEFAULT_MATERIAL");

	return(sc);
}

/* Duplicate a scene, optionally selecting faces */
Scene *
dupScene(const Scene *osc, int flreq, int flexc)
{
	int     	fltest = flreq | flexc;
	Scene		*sc;
	const Face      *fo;
	Face		*f;
	int		i;

	if (osc == NULL)
		return(NULL);
	sc = newScene();
	for (i = 0; i < osc->ndescr; i++)
		addComment(sc, osc->descr[i]);
	if (osc->ngrps > 1) {
		sc->grpname = (char **)erealloc((char *)sc->grpname,
				sizeof(char *) * (osc->ngrps+(CHUNKSIZ-1)));
		for (i = 1; i < osc->ngrps; i++)
			sc->grpname[i] = savqstr(osc->grpname[i]);
		sc->ngrps = osc->ngrps;
	}
	if (osc->nmats > 1) {
		sc->matname = (char **)erealloc((char *)sc->matname,
				sizeof(char *) * (osc->nmats+(CHUNKSIZ-1)));
		for (i = 1; i < osc->nmats; i++)
			sc->matname[i] = savqstr(osc->matname[i]);
		sc->nmats = osc->nmats;
	}
	if (osc->nverts) {
		sc->vert = (Vertex *)emalloc(sizeof(Vertex) *
						(osc->nverts+(CHUNKSIZ-1)));
		memcpy(sc->vert, osc->vert, sizeof(Vertex)*osc->nverts);
		sc->nverts = osc->nverts;
		for (i = 0; i < sc->nverts; i++)
			sc->vert[i].vflist = NULL;
	}
	if (osc->ntex) {
		sc->tex = (TexCoord *)emalloc(sizeof(TexCoord) *
						(osc->ntex+(CHUNKSIZ-1)));
		memcpy(sc->tex, osc->tex, sizeof(TexCoord)*osc->ntex);
		sc->ntex = osc->ntex;
	}
	if (osc->nnorms) {
		sc->norm = (Normal *)emalloc(sizeof(Normal) *
						(osc->nnorms+CHUNKSIZ));
		memcpy(sc->norm, osc->norm, sizeof(Normal)*osc->nnorms);
		sc->nnorms = osc->nnorms;
	}
	for (fo = osc->flist; fo != NULL; fo = fo->next) {
		if ((fo->flags & fltest) != flreq)
			continue;
		f = (Face *)emalloc(sizeof(Face) + sizeof(VertEnt)*(fo->nv-3));
		memcpy(f, fo, sizeof(Face) + sizeof(VertEnt)*(fo->nv-3));
		for (i = 0; i < f->nv; i++)
			add2facelist(sc, f, i);
		f->next = sc->flist;
		sc->flist = f;
		sc->nfaces++;
	}
	deleteUnreferenced(sc);		/* jetsam */
	return(sc);
}

/* Add one scene to another, not checking for redundancies */
int
addScene(Scene *scdst, const Scene *scsrc)
{
	VNDX		my_vlist[4];
	int		*vert_map = NULL;
	int		tex_off = 0;
	int		norm_off = 0;
	VNDX		*vlist = my_vlist;
	int		vllen = sizeof(my_vlist)/sizeof(VNDX);
	int		cur_mat = 0;
	int		cur_grp = 0;
	int		fcnt = 0;
	const Face	*f;
	int		i;

	if ((scdst == NULL) | (scsrc == NULL))
		return(-1);
	if (scsrc->nfaces <= 0)
		return(0);
					/* map vertices */
	vert_map = (int *)emalloc(sizeof(int)*scsrc->nverts);
	for (i = 0; i < scsrc->nverts; i++) {
		const Vertex	*v = scsrc->vert + i;
		if (v->vflist == NULL) {
			vert_map[i] = -1;
			continue;
		}
		vert_map[i] = addVertex(scdst, v->p[0], v->p[1], v->p[2]);
	}
	tex_off = scdst->ntex;		/* append texture coords */
	if (scsrc->ntex > 0) {
		scdst->tex = (TexCoord *)erealloc((char *)scdst->tex,
			sizeof(TexCoord)*(tex_off+scsrc->ntex+(CHUNKSIZ-1)));
		memcpy(scdst->tex+tex_off, scsrc->tex,
				sizeof(TexCoord)*scsrc->ntex);
	}
	norm_off = scdst->nnorms;	/* append normals */
	if (scsrc->nnorms > 0) {
		scdst->norm = (Normal *)erealloc((char *)scdst->norm,
			sizeof(Normal)*(norm_off+scsrc->nnorms+(CHUNKSIZ-1)));
		memcpy(scdst->norm+norm_off, scsrc->norm,
				sizeof(Normal)*scsrc->nnorms);
	}
					/* add faces */
	scdst->lastgrp = scdst->lastmat = 0;
	for (f = scsrc->flist; f != NULL; f = f->next) {
		if (f->grp != cur_grp)
			setGroup(scdst, scsrc->grpname[cur_grp = f->grp]);
		if (f->mat != cur_mat)
			setMaterial(scdst, scsrc->matname[cur_mat = f->mat]);
		if (f->nv > vllen) {
			vlist = (VNDX *)( vlist == my_vlist ?
				emalloc(sizeof(VNDX)*f->nv) :
				erealloc((char *)vlist, sizeof(VNDX)*f->nv) );
			vllen = f->nv;
		}
		memset(vlist, 0xff, sizeof(VNDX)*f->nv);
		for (i = f->nv; i-- > 0; ) {
			if (f->v[i].vid >= 0)
				vlist[i][0] = vert_map[f->v[i].vid];
			if (f->v[i].tid >= 0)
				vlist[i][1] = f->v[i].tid + tex_off;
			if (f->v[i].nid >= 0)
				vlist[i][2] = f->v[i].nid + norm_off;
		}
		fcnt += (addFace(scdst, vlist, f->nv) != NULL);
	}
					/* clean up */
	if (vlist != my_vlist) efree((char *)vlist);
	efree((char *)vert_map);
	return(fcnt);
}

#define MAXAC	100

/* Transform entire scene */
int
xfScene(Scene *sc, int xac, char *xav[])
{
	char	comm[24+MAXAC*8];
	char	*cp;
	XF	myxf;
	FVECT	vec;
	int	i;

	if ((sc == NULL) | (xac <= 0) | (xav == NULL))
		return(0);
					/* compute matrix */
	if (xf(&myxf, xac, xav) < xac)
		return(0);
					/* transform vertices */
	for (i = 0; i < sc->nverts; i++) {
		VCOPY(vec, sc->vert[i].p);
		multp3(vec, vec, myxf.xfm);
		VCOPY(sc->vert[i].p, vec);
	}
					/* transform normals */
	for (i = 0; i < sc->nnorms; i++) {
		VCOPY(vec, sc->norm[i]);
		multv3(vec, vec, myxf.xfm);
		vec[0] /= myxf.sca; vec[1] /= myxf.sca; vec[2] /= myxf.sca;
		VCOPY(sc->norm[i], vec);
	}
					/* add comment */
	cp = strcpy(comm, "Transformed by:");
	for (i = 0; i < xac; i++) {
		while (*cp) cp++;
		*cp++ = ' ';
		strcpy(cp, xav[i]);
	}
	addComment(sc, comm);
	return(xac);			/* all done */
}

/* Ditto, using transform string rather than pre-parsed words */
int
xfmScene(Scene *sc, const char *xfm)
{
	char	*xav[MAXAC+1];
	int	xac, i;

	if ((sc == NULL) | (xfm == NULL))
		return(0);
					/* skip spaces at beginning */
	while (isspace(*xfm))
		xfm++;
	if (!*xfm)
		return(0);
					/* parse string into words */
	xav[0] = savqstr((char *)xfm);
	xac = 1; i = 0;
	for ( ; ; ) {
		while (!isspace(xfm[++i]))
			if (!xfm[i])
				break;
		while (isspace(xfm[i]))
			xav[0][i++] = '\0';
		if (!xfm[i])
			break;
		if (xac >= MAXAC-1) {
			freeqstr(xav[0]);
			return(0);
		}
		xav[xac++] = xav[0] + i;
	}
	xav[xac] = NULL;
	i = xfScene(sc, xac, xav);
	freeqstr(xav[0]);
	return(i);
}
#undef MAXAC

/* Delete unreferenced vertices, normals, texture coords */
void
deleteUnreferenced(Scene *sc)
{
	int     *vmap;
	Face    *f;
	int     nused, i;
						/* allocate index map */
	if (!sc->nverts)
		return;
	i = sc->nverts;
	if (sc->ntex > i)
		i = sc->ntex;
	if (sc->nnorms > i)
		i = sc->nnorms;
	vmap = (int *)emalloc(sizeof(int)*i);
						/* remap positions */
	for (i = nused = 0; i < sc->nverts; i++) {
		if (sc->vert[i].vflist == NULL) {
			vmap[i] = -1;
			continue;
		}
		if (nused != i)
			sc->vert[nused] = sc->vert[i];
		vmap[i] = nused++;
	}
	if (nused == sc->nverts)
		goto skip_pos;
	sc->vert = (Vertex *)erealloc((char *)sc->vert,
				sizeof(Vertex)*(nused+(CHUNKSIZ-1)));
	sc->nverts = nused;
	for (f = sc->flist; f != NULL; f = f->next)
		for (i = f->nv; i--; )
			if ((f->v[i].vid = vmap[f->v[i].vid]) < 0)
				error(CONSISTENCY,
					"Link error in del_unref_verts()");
skip_pos:
						/* remap texture coord's */
	if (!sc->ntex)
		goto skip_tex;
	memset((void *)vmap, 0, sizeof(int)*sc->ntex);
	for (f = sc->flist; f != NULL; f = f->next)
		for (i = f->nv; i--; )
			if (f->v[i].tid >= 0)
				vmap[f->v[i].tid] = 1;
	for (i = nused = 0; i < sc->ntex; i++) {
		if (!vmap[i])
			continue;
		if (nused != i)
			sc->tex[nused] = sc->tex[i];
		vmap[i] = nused++;
	}
	if (nused == sc->ntex)
		goto skip_tex;
	sc->tex = (TexCoord *)erealloc((char *)sc->tex,
				sizeof(TexCoord)*(nused+(CHUNKSIZ-1)));
	sc->ntex = nused;
	for (f = sc->flist; f != NULL; f = f->next)
		for (i = f->nv; i--; )
			if (f->v[i].tid >= 0)
				f->v[i].tid = vmap[f->v[i].tid];
skip_tex:
						/* remap normals */
	if (!sc->nnorms)
		goto skip_norms;
	memset((void *)vmap, 0, sizeof(int)*sc->nnorms);
	for (f = sc->flist; f != NULL; f = f->next)
		for (i = f->nv; i--; )
			if (f->v[i].nid >= 0)
				vmap[f->v[i].nid] = 1;
	for (i = nused = 0; i < sc->nnorms; i++) {
		if (!vmap[i])
			continue;
		if (nused != i)
			memcpy(sc->norm[nused], sc->norm[i], sizeof(Normal));
		vmap[i] = nused++;
	}
	if (nused == sc->nnorms)
		goto skip_norms;
	sc->norm = (Normal *)erealloc((char *)sc->norm,
				sizeof(Normal)*(nused+(CHUNKSIZ-1)));
	sc->nnorms = nused;
	for (f = sc->flist; f != NULL; f = f->next)
		for (i = f->nv; i--; )
			if (f->v[i].nid >= 0)
				f->v[i].nid = vmap[f->v[i].nid];
skip_norms:
						/* clean up */
	efree((char *)vmap);
}

/* Free a scene */
void
freeScene(Scene *sc)
{
	int     i;
	Face    *f;
	
	if (sc == NULL)
		return;
	clearComments(sc);
	for (i = sc->ngrps; i-- > 0; )
		freeqstr(sc->grpname[i]);
	efree((char *)sc->grpname);
	for (i = sc->nmats; i-- > 0; )
		freeqstr(sc->matname[i]);
	efree((char *)sc->matname);
	efree((char *)sc->vert);
	efree((char *)sc->tex);
	efree((char *)sc->norm);
	while ((f = sc->flist) != NULL) {
		sc->flist = f->next;
		efree((char *)f);
	}
	efree((char *)sc);
}
