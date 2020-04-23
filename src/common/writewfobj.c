#ifndef lint
static const char RCSid[] = "$Id: writewfobj.c,v 2.2 2020/04/23 03:19:48 greg Exp $";
#endif
/*
 *  writewfobj.c
 *
 *  Routines for writing out Wavefront file.
 *
 *  Created by Greg Ward on Thu Feb 12 2004.
 */

#include "paths.h"
#include "rterror.h"
#include "objutil.h"

/* Delete unreferenced vertices */
static void
del_unref_verts(Scene *sc)
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
			memcpy((void *)sc->norm[nused],
					(void *)sc->norm[i],
					sizeof(Normal));
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

/* Write out header comments */
static void
write_header(Scene *sc, FILE *fp)
{
	int     i;
	
	for (i = 0; i < sc->ndescr; i++)
		fprintf(fp, "# %s\n", sc->descr[i]);
	fputs("#\n", fp);
	fprintf(fp, "# %d final faces\n", sc->nfaces);
	fprintf(fp,
		"#\t%d vertices, %d texture coordinates, %d surface normals\n",
			sc->nverts, sc->ntex, sc->nnorms);
	fputs("#\n\n", fp);
}

/* Write out vertex lists */
static void
write_verts(Scene *sc, FILE *fp)
{
	int     i;
	
	fputs("# Vertex positions\n", fp);
	for (i = 0; i < sc->nverts; i++)
		fprintf(fp, "v %.12g %.12g %.12g\n",
				sc->vert[i].p[0],
				sc->vert[i].p[1],
				sc->vert[i].p[2]);
	fputs("\n# Vertex texture coordinates\n", fp);
	for (i = 0; i < sc->ntex; i++)
		fprintf(fp, "vt %.6g %.6g\n",
				sc->tex[i].u, sc->tex[i].v);
	fputs("\n# Vertex normals\n", fp);
	for (i = 0; i < sc->nnorms; i++)
		fprintf(fp, "vn %.6f %.6f %.6f\n",
				sc->norm[i][0],
				sc->norm[i][1],
				sc->norm[i][2]);
	fputc('\n', fp);
}

/* Write out an object group */
static void
write_group(const Scene *sc, int gid, FILE *fp)
{
	int		mid = -1;
	const Face      *f;
	int		j;
	
	fprintf(fp, "# Face group\ng %s\n", sc->grpname[gid]);
	for (f = sc->flist; f != NULL; f = f->next) {
		if (f->grp != gid)
			continue;
		if (f->mat != mid)
			fprintf(fp, "usemtl %s\n", sc->matname[mid=f->mat]);
		fputc('f', fp);
		for (j = 0; j < f->nv; j++) {
			fprintf(fp, " %d/", f->v[j].vid - sc->nverts);
			if (f->v[j].tid >= 0)
				fprintf(fp, "%d/", f->v[j].tid - sc->ntex);
			else
				fputc('/', fp);
			if (f->v[j].nid >= 0)
				fprintf(fp, "%d", f->v[j].nid - sc->nnorms);
		}
		fputc('\n', fp);
	}
	fprintf(fp, "# End of face group %s\n\n", sc->grpname[gid]);
}

/* Write a .OBJ to stream, return # faces written */
int
toOBJ(Scene *sc, FILE *fp)
{
	int     i;

	if (sc == NULL || sc->nfaces <= 0 || fp == NULL)
		return(0);
	if (verbose)
		fputs(" Removing unreferenced vertices\r", stderr);
	del_unref_verts(sc);		/* clean up, first */
	if (verbose)
		fputs(" Writing vertices               \r", stderr);
	write_header(sc, fp);		/* write out header */
	write_verts(sc, fp);		/* write out vertex lists */
					/* write out each object group */
	for (i = 0; i < sc->ngrps; i++) {
		if (verbose)
			fprintf(stderr, " Writing faces %s               \r",
					sc->grpname[i]);
		write_group(sc, i, fp);
	}
	if (fflush(fp) < 0) {
		error(SYSTEM, "Error flushing .OBJ output");
		return(-1);
	}
	return(sc->nfaces);
}

/* Write a .OBJ file, return # faces written */
int
writeOBJ(Scene *sc, const char *fspec)
{
	extern char     *progname;
	FILE		*fp;
	int		n;

	if (sc == NULL || sc->nfaces <= 0 || fspec == NULL || !*fspec)
		return(0);
#if POPEN_SUPPORT
	if (fspec[0] == '!') {
		if ((fp = popen(fspec+1, "w")) == NULL) {
			sprintf(errmsg, "%s: cannot execute", fspec);
			error(SYSTEM, errmsg);
			return(0);
		}
	} else
#endif
	if ((fp = fopen(fspec, "w")) == NULL) {
		sprintf(errmsg, "%s: cannot open for writing", fspec);
		error(SYSTEM, errmsg);
		return(0);
	}
						/* start off header */
	fprintf(fp, "# Wavefront .OBJ file created by %s\n#\n", progname);
	n = toOBJ(sc, fp);			/* write scene */
#if POPEN_SUPPORT
	if (fspec[0] == '!')
		pclose(fp);
	else
#endif
	fclose(fp);				/* close it */
	return(n);
}
