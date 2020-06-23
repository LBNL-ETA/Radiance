#ifndef lint
static const char RCSid[] = "$Id: objutil.c,v 2.8 2020/06/23 19:29:40 greg Exp $";
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
	int    i;
	
	for (i = n; i-- > 0; )
		if (!strcmp(nmlist[i], nm))
			break;
	return(i);
}

/* Clear face selection */
void
clearSelection(Scene *sc, int set)
{
	Face    *f;
	
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
	int     gid = findName(gname, (const char **)sc->grpname, sc->ngrps);
	Face    *f;

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
	int     mid = findName(mname, (const char **)sc->matname, sc->nmats);
	Face    *f;
	
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
	Face    *f;
	
	for (f = sc->flist; f != NULL; f = f->next)
		f->flags ^= FACE_SELECTED;
}

/* Count selected faces */
int
numberSelected(Scene *sc)
{
	int     nsel = 0;
	Face    *f;
	
	for (f = sc->flist; f != NULL; f = f->next)
		nsel += ((f->flags & FACE_SELECTED) != 0);
	return(nsel);
}

/* Execute callback on indicated faces */
int
foreachFace(Scene *sc, int (*cb)(Scene *, Face *, void *),
			int flreq, int flexc, void *c_data)
{
	int     fltest = flreq | flexc;
	int     sum = 0;
	Face    *f;
	
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
	int     hadTexture = 0;
	int     i;
	
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
	int     hadNormal = 0;
	int     i;
	
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
	int     grp = *(int *)ptr;
	if (f->grp == grp)
		return(0);
	f->grp = grp;
	return(1);
}

/* Change group for the indicated faces */
int
changeGroup(Scene *sc, const char *gname, int flreq, int flexc)
{
	int     grp = findName(gname, (const char **)sc->grpname, sc->ngrps);
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
	int     mat = *(int *)ptr;
	if (f->mat == mat)
		return(0);
	f->mat = mat;
	return(1);
}

/* Change material for the indicated faces */
int
changeMaterial(Scene *sc, const char *mname, int flreq, int flexc)
{
	int     mat = findName(mname, (const char **)sc->matname, sc->nmats);
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
	int     repl_tex[10];
	int     repl_norm[10];
	Face    *f, *flast=NULL;
	int     i, j=0;
					/* check to see if it's even used */
	if (sc->vert[prev].vflist == NULL)
		return(0);
					/* get replacement textures */
	repl_tex[0] = -1;
	for (f = sc->vert[repl].vflist; f != NULL; f = f->v[j].fnext) {
					/* make sure prev isn't in there */
		for (j = 0; j < f->nv; j++)
			if (f->v[j].vid == prev)
				return(0);
		for (j = 0; j < f->nv; j++)
			if (f->v[j].vid == repl)
				break;
		if (j >= f->nv)
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
		for (j = 0; j < f->nv; j++)
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
		for (j = 0; j < f->nv; j++)
			if (f->v[j].vid == prev)
				break;
		if (j >= f->nv)
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
	int     nelim = 0;
	LUTAB   myLookup;
	LUENT   *le;
	char    vertfmt[32], vertbuf[64];
	double  d;
	int     i;
	
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
	int     nchecked = 0;
	int     nfound = 0;
	Face    *f, *f1;
	int     vid;
	int     i, j;
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
			for (j = 0; j < f1->nv; j++)
				if (f1->v[j].vid == vid)
					break;
			if (j >= f1->nv)
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
	return(nfound);
}

/* Delete indicated faces */
int
deleteFaces(Scene *sc, int flreq, int flexc)
{
	int     fltest = flreq | flexc;
	int     orig_nfaces = sc->nfaces;
	Face    fhead;
	Face    *f, *ftst;
	
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
					for (j = 0; j < vf->nv; j++)
						if (vf->v[j].vid == vid)
							break;
					if (j >= vf->nv)
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
		if (strcasestr(sc->descr[n], match) != NULL)
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
		for (j = 0; j < fp->nv; j++)
			if (fp->v[j].vid == vid)
				break;
		if (j >= fp->nv)
			error(CONSISTENCY, "Link error in add2facelist()");
		if (fp->v[j].fnext == NULL)
			break;			/* reached the end */
		fp = fp->v[j].fnext;
	}
	fp->v[j].fnext = f;			/* append new face */
}

/* Allocate an empty scene */
Scene *
newScene(void)
{
	Scene   *sc = (Scene *)ecalloc(1, sizeof(Scene));
						/* default group & material */
	sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
	sc->grpname[sc->ngrps++] = savqstr("DEFAULT_GROUP");
	sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
	sc->matname[sc->nmats++] = savqstr("DEFAULT_MATERIAL");

	return(sc);
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
	FVECT   nrm;

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
	sc->lastgrp = findName(nm, (const char **)sc->grpname, sc->ngrps);
	if (sc->lastgrp >= 0)
		return;
	sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
	sc->grpname[sc->lastgrp=sc->ngrps++] = savqstr((char *)nm);
}

/* Set current (last) material */
void
setMaterial(Scene *sc, const char *nm)
{
	sc->lastmat = findName(nm, (const char **)sc->matname, sc->nmats);
	if (sc->lastmat >= 0)
		return;
	sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
	sc->matname[sc->lastmat=sc->nmats++] = savqstr((char *)nm);
}

/* Add new face to a scene */
Face *
addFace(Scene *sc, VNDX vid[], int nv)
{
	Face    *f;
	int     i;
	
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

/* Duplicate a scene */
Scene *
dupScene(const Scene *osc)
{
	Scene		*sc;
	const Face      *fo;
	Face		*f;
	int		i;

	if (osc == NULL)
		return(NULL);
	sc = newScene();
	for (i = 0; i < osc->ndescr; i++)
		addComment(sc, osc->descr[i]);
	if (osc->ngrps) {
		sc->grpname = (char **)emalloc(sizeof(char *) *
						(osc->ngrps+(CHUNKSIZ-1)));
		for (i = 0; i < osc->ngrps; i++)
			sc->grpname[i] = savqstr(osc->grpname[i]);
		sc->ngrps = osc->ngrps;
	}
	if (osc->nmats) {
		sc->matname = (char **)emalloc(sizeof(char *) *
						(osc->nmats+(CHUNKSIZ-1)));
		for (i = 0; i < osc->nmats; i++)
			sc->matname[i] = savqstr(osc->matname[i]);
		sc->nmats = osc->nmats;
	}
	if (osc->nverts) {
		sc->vert = (Vertex *)emalloc(sizeof(Vertex) *
						(osc->nverts+(CHUNKSIZ-1)));
		memcpy((void *)sc->vert, (const void *)osc->vert,
				sizeof(Vertex)*osc->nverts);
		sc->nverts = osc->nverts;
		for (i = 0; i < sc->nverts; i++)
			sc->vert[i].vflist = NULL;
	}
	if (osc->ntex) {
		sc->tex = (TexCoord *)emalloc(sizeof(TexCoord) *
						(osc->ntex+(CHUNKSIZ-1)));
		memcpy((void *)sc->tex, (const void *)osc->tex,
				sizeof(TexCoord)*osc->ntex);
		sc->ntex = osc->ntex;
	}
	if (osc->nnorms) {
		sc->norm = (Normal *)emalloc(sizeof(Normal) *
						(osc->nnorms+CHUNKSIZ));
		memcpy((void *)sc->norm, (const void *)osc->norm,
				sizeof(Normal)*osc->nnorms);
		sc->nnorms = osc->nnorms;
	}
	for (fo = osc->flist; fo != NULL; fo = fo->next) {
		f = (Face *)emalloc(sizeof(Face) +
				sizeof(VertEnt)*(fo->nv-3));
		memcpy((void *)f, (const void *)fo, sizeof(Face) +
				sizeof(VertEnt)*(fo->nv-3));
		for (i = 0; i < f->nv; i++)
			add2facelist(sc, f, i);
		f->next = sc->flist;
		sc->flist = f;
		sc->nfaces++;
	}
	return(sc);
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
	xav[0] = strcpy((char *)malloc(strlen(xfm)+1), xfm);
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
			free(xav[0]);
			return(0);
		}
		xav[xac++] = xav[0] + i;
	}
	xav[xac] = NULL;
	i = xfScene(sc, xac, xav);
	free(xav[0]);
	return(i);
}
#undef MAXAC

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
