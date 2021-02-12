/* RCSid $Id: objutil.h,v 2.9 2021/02/12 01:57:49 greg Exp $ */
/*
 *  Declarations for .OBJ file utility
 *
 *  Include after <stdio.h>
 *
 *  Created by Greg Ward on Wed Feb 11 2004.
 */

#ifndef _OBJUTIL_H_
#define _OBJUTIL_H_

#ifndef DUP_CHECK_REVERSE
#define DUP_CHECK_REVERSE       1       /* eliminate flipped duplicates */
#endif
#ifndef POPEN_SUPPORT
#define POPEN_SUPPORT		1       /* support "!command" i/o */
#endif
					/* face flags */
#define FACE_SELECTED		01
#define FACE_DEGENERATE		02
#define FACE_DUPLICATE		04

struct Face;				/* forward declaration */

typedef int		VNDX[3];	/* vertex indices (point,map,normal) */

/* Structure to hold vertex indices and link back to face list */
typedef struct {
	int		vid;		/* vertex id */
	int		tid;		/* texture id */
	int		nid;		/* normal id */
	struct Face     *fnext;		/* next in vertex face list */
} VertEnt;

/* Structure to hold face and its vertex references */
typedef struct Face {
	struct Face     *next;		/* next face in main list */
	short		flags;		/* face selected, etc */
	short		nv;		/* number of vertices */
	short		grp;		/* group/object index */
	short		mat;		/* material index */
	/* in cases where the same vertex appears twice, first has link */
	VertEnt		v[3];		/* vertex list (extends struct) */
} Face;

/* Structure to hold vertex */
typedef struct {
	double		p[3];		/* 3-D position */
	Face		*vflist;	/* linked face list (no repeats) */
} Vertex;

/* Structure to hold texture coordinate */
typedef struct {
	float		u, v;		/* 2-D local texture coordinate */
} TexCoord;

/* Array to hold surface normal */
typedef float		Normal[3];

/* Structure to hold a loaded .OBJ file */
typedef struct {
	char		**descr;	/* descriptive comments */
	int		ndescr;		/* number of comments */
	char		**grpname;      /* object/group name list */
	int		ngrps;		/* number of group names */
	int		lastgrp;	/* last group seen */
	char		**matname;      /* material name list */
	int		nmats;		/* number of materials */
	int		lastmat;	/* last material seen */
	Vertex		*vert;		/* vertex array */
	int		nverts;		/* number of vertices */
	TexCoord	*tex;		/* texture coordinate array */
	int		ntex;		/* number of texture coord's */
	Normal		*norm;		/* surface normal array */
	int		nnorms;		/* number of surface normals */
	Face		*flist;		/* linked face list */
	int		nfaces;		/* count of faces */
} Scene;

#ifdef __cplusplus
extern "C" {
#endif

/* Allocate a new scene holder */
Scene *		newScene(void);

/* Add a .OBJ file to a scene */
Scene *		loadOBJ(Scene *sc, const char *fspec);

/* Duplicate a scene, optionally selecting faces */
Scene *		dupScene(const Scene *sc, int flreq, int flexc);

/* Transform entire scene */
int		xfScene(Scene *sc, int xac, char *xav[]);
int		xfmScene(Scene *sc, const char *xfm);

/* Add a descriptive comment */
void		addComment(Scene *sc, const char *comment);

/* Find index for comment containing the given string (starting from n) */
int		findComment(Scene *sc, const char *match, int n);

/* Clear comments */
void		clearComments(Scene *sc);

/* Write a .OBJ file, return # faces written or -1 on error */
int		toOBJ(Scene *sc, FILE *fp);
int		writeOBJ(Scene *sc, const char *fspec);

/* Convert indicated faces to Radiance, return # written or -1 on error */
int		toRadiance(Scene *sc, FILE *fp, int flreq, int flexc);
int		writeRadiance(Scene *sc, const char *fspec,
				int flreq, int flexc);

/* Compute face area (and normal) */
double		faceArea(const Scene *sc, const Face *f, Normal nrm);

/* Eliminate duplicate vertices, return # joined */
int		coalesceVertices(Scene *sc, double eps);

/* Identify duplicate faces */
int		findDuplicateFaces(Scene *sc);

/* Delete indicated faces, return # deleted */
int		deleteFaces(Scene *sc, int flreq, int flexc);

/* Clear face selection */
void		clearSelection(Scene *sc, int set);

/* Invert face selection */
void		invertSelection(Scene *sc);

/* Count number of faces selected */
int		numberSelected(Scene *sc);

/* Select faces by object name (modifies current) */
void		selectGroup(Scene *sc, const char *gname, int invert);

/* Select faces by material name (modifies current) */
void		selectMaterial(Scene *sc, const char *mname, int invert);

/* Execute callback on indicated faces */
int		foreachFace(Scene *sc, int (*cb)(Scene *, Face *, void *),
					int flreq, int flexc, void *c_data);

/* Remove texture coordinates from the indicated faces */
int		removeTexture(Scene *sc, int flreq, int flexc);

/* Remove surface normals from the indicated faces */
int		removeNormals(Scene *sc, int flreq, int flexc);

/* Change group for the indicated faces */
int		changeGroup(Scene *sc, const char *gname,
					int flreq, int flexc);

/* Change material for the indicated faces */
int		changeMaterial(Scene *sc, const char *mname,
					int flreq, int flexc);

/* Add a vertex to our scene, returning index */
int		addVertex(Scene *sc, double x, double y, double z);

/* Add a texture coordinate to our scene, returning index */
int		addTexture(Scene *sc, double u, double v);

/* Add a surface normal to our scene, returning index */
int		addNormal(Scene *sc, double xn, double yn, double zn);

/* Set current group (sc->lastgrp) to given ID */
void		setGroup(Scene *sc, const char *nm);

/* Set current material (sc->lastmat) to given ID */
void		setMaterial(Scene *sc, const char *nm);

/* Add a new face to our scene, using current group and material */
Face *		addFace(Scene *sc, VNDX vid[], int nv);

/* Delete unreferenced vertices, normals, texture coords */
void		deleteUnreferenced(Scene *sc);

/* Free a scene */
void		freeScene(Scene *sc);

/* Find an existing name in a list of names */
int		findName(const char *nm, const char **nmlist, int n);

/* Verbose mode global */
extern int      verbose;

extern char     *emalloc(unsigned int n);
extern char     *ecalloc(unsigned int ne, unsigned int n);
extern char     *erealloc(char *cp, unsigned int n);
extern void     efree(char *cp); 

#define CHUNKSIZ	128	/* object allocation chunk size */

#define chunk_alloc(typ, arr, nold) \
		((nold)%CHUNKSIZ ? (arr) : \
		(typ *)erealloc((char *)(arr), sizeof(typ)*((nold)+CHUNKSIZ)))

#ifdef __cplusplus
}
#endif

#endif  /* ! _OBJUTIL_H_ */
