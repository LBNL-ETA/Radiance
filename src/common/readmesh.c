#ifndef lint
static const char RCSid[] = "$Id: readmesh.c,v 2.3 2003/06/05 19:29:34 schorsch Exp $";
#endif
/*
 *  Routines for reading a compiled mesh from a file
 */

#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "mesh.h"

static char	*meshfn;	/* input file name */
static FILE	*meshfp;	/* mesh file pointer */
static int	objsize;	/* sizeof(OBJECT) from writer */


static void
mesherror(etyp, msg)			/* mesh read error */
int  etyp;
char  *msg;
{
	char  msgbuf[128];

	sprintf(msgbuf, "(%s): %s", meshfn, msg);
	error(etyp, msgbuf);
}


static long
mgetint(siz)				/* get a siz-byte integer */
int  siz;
{
	register long  r;

	r = getint(siz, meshfp);
	if (feof(meshfp))
		mesherror(USER, "truncated mesh file");
	return(r);
}


static double
mgetflt()				/* get a floating point number */
{
	double	r;

	r = getflt(meshfp);
	if (feof(meshfp))
		mesherror(USER, "truncated mesh file");
	return(r);
}
	

static OCTREE
getfullnode()				/* get a set, return fullnode */
{
	OBJECT	set[MAXSET+1];
	register int  i;

	if ((set[0] = mgetint(objsize)) > MAXSET)
		mesherror(USER, "bad set in getfullnode");
	for (i = 1; i <= set[0]; i++)
		set[i] = mgetint(objsize);
	return(fullnode(set));
}	


static OCTREE
gettree()				/* get a pre-ordered octree */
{
	register OCTREE	 ot;
	register int  i;
	
	switch (getc(meshfp)) {
	case OT_EMPTY:
		return(EMPTY);
	case OT_FULL:
		return(getfullnode());
	case OT_TREE:
		if ((ot = octalloc()) == EMPTY)
			mesherror(SYSTEM, "out of tree space in readmesh");
		for (i = 0; i < 8; i++)
			octkid(ot, i) = gettree();
		return(ot);
	case EOF:
		mesherror(USER, "truncated mesh octree");
	default:
		mesherror(USER, "damaged mesh octree");
	}
}


static void
skiptree()				/* skip octree on input */
{
	register int  i;
	
	switch (getc(meshfp)) {
	case OT_EMPTY:
		return;
	case OT_FULL:
		for (i = mgetint(objsize)*objsize; i-- > 0; )
			if (getc(meshfp) == EOF)
				mesherror(USER, "truncated mesh octree");
		return;
	case OT_TREE:
		for (i = 0; i < 8; i++)
			skiptree();
		return;
	case EOF:
		mesherror(USER, "truncated mesh octree");
	default:
		mesherror(USER, "damaged mesh octree");
	}
}


static void
getpatch(pp)				/* load a mesh patch */
register MESHPATCH	*pp;
{
	int	flags;
	int	i, j;
					/* vertex flags */
	flags = mgetint(1);
	if (!(flags & MT_V) || flags & ~(MT_V|MT_N|MT_UV))
		mesherror(USER, "bad patch flags");
					/* allocate vertices */
	pp->nverts = mgetint(2);
	if (pp->nverts <= 0 || pp->nverts > 256)
		mesherror(USER, "bad number of patch vertices");
	pp->xyz = (uint4 (*)[3])malloc(pp->nverts*3*sizeof(uint4));
	if (pp->xyz == NULL)
		goto nomem;
	if (flags & MT_N) {
		pp->norm = (int4 *)calloc(pp->nverts, sizeof(int4));
		if (pp->norm == NULL)
			goto nomem;
	} else
		pp->norm = NULL;
	if (flags & MT_UV) {
		pp->uv = (uint4 (*)[2])calloc(pp->nverts, 2*sizeof(uint4));
		if (pp->uv == NULL)
			goto nomem;
	} else
		pp->uv = NULL;
					/* vertex xyz locations */
	for (i = 0; i < pp->nverts; i++)
		for (j = 0; j < 3; j++)
			pp->xyz[i][j] = mgetint(4);
					/* vertex normals */
	if (flags & MT_N)
		for (i = 0; i < pp->nverts; i++)
			pp->norm[i] = mgetint(4);
					/* uv coordinates */
	if (flags & MT_UV)
		for (i = 0; i < pp->nverts; i++)
			for (j = 0; j < 2; j++)
				pp->uv[i][j] = mgetint(4);
					/* local triangles */
	pp->ntris = mgetint(2);
	if (pp->ntris < 0 || pp->ntris > 512)
		mesherror(USER, "bad number of local triangles");
	if (pp->ntris) {
		pp->tri = (struct PTri *)malloc(pp->ntris *
						sizeof(struct PTri));
		if (pp->tri == NULL)
			goto nomem;
		for (i = 0; i < pp->ntris; i++) {
			pp->tri[i].v1 = mgetint(1);
			pp->tri[i].v2 = mgetint(1);
			pp->tri[i].v3 = mgetint(1);
		}
	} else
		pp->tri = NULL;
					/* local triangle material(s) */
	if (mgetint(2) > 1) {
		pp->trimat = (int2 *)malloc(pp->ntris*sizeof(int2));
		if (pp->trimat == NULL)
			goto nomem;
		for (i = 0; i < pp->ntris; i++)
			pp->trimat[i] = mgetint(2);
	} else {
		pp->solemat = mgetint(2);
		pp->trimat = NULL;
	}
					/* joiner triangles */
	pp->nj1tris = mgetint(2);
	if (pp->nj1tris < 0 || pp->nj1tris > 512)
		mesherror(USER, "bad number of joiner triangles");
	if (pp->nj1tris) {
		pp->j1tri = (struct PJoin1 *)malloc(pp->nj1tris *
						sizeof(struct PJoin1));
		if (pp->j1tri == NULL)
			goto nomem;
		for (i = 0; i < pp->nj1tris; i++) {
			pp->j1tri[i].v1j = mgetint(4);
			pp->j1tri[i].v2 = mgetint(1);
			pp->j1tri[i].v3 = mgetint(1);
			pp->j1tri[i].mat = mgetint(2);
		}
	} else
		pp->j1tri = NULL;
					/* double joiner triangles */
	pp->nj2tris = mgetint(2);
	if (pp->nj2tris < 0 || pp->nj2tris > 256)
		mesherror(USER, "bad number of double joiner triangles");
	if (pp->nj2tris) {
		pp->j2tri = (struct PJoin2 *)malloc(pp->nj2tris *
						sizeof(struct PJoin2));
		if (pp->j2tri == NULL)
			goto nomem;
		for (i = 0; i < pp->nj2tris; i++) {
			pp->j2tri[i].v1j = mgetint(4);
			pp->j2tri[i].v2j = mgetint(4);
			pp->j2tri[i].v3 = mgetint(1);
			pp->j2tri[i].mat = mgetint(2);
		}
	} else
		pp->j2tri = NULL;
	return;
nomem:
	error(SYSTEM, "out of mesh memory in getpatch");
}


void
readmesh(mp, path, flags)		/* read in mesh structures */
MESH	*mp;
char	*path;
int	flags;
{
	char	*err;
	char	sbuf[64];
	int	i;
					/* check what's loaded */
	flags &= (IO_INFO|IO_BOUNDS|IO_TREE|IO_SCENE) & ~mp->ldflags;
					/* open input file */
	if (path == NULL) {
		meshfn = "standard input";
		meshfp = stdin;
	} else if ((meshfp = fopen(meshfn=path, "r")) == NULL) {
		sprintf(errmsg, "cannot open mesh file \"%s\"", path);
		error(SYSTEM, errmsg);
	}
	SET_FILE_BINARY(meshfp);
					/* read header */
	checkheader(meshfp, MESHFMT, flags&IO_INFO ? stdout : (FILE *)NULL);
					/* read format number */
	objsize = getint(2, meshfp) - MESHMAGIC;
	if (objsize <= 0 || objsize > MAXOBJSIZ || objsize > sizeof(long))
		mesherror(USER, "incompatible mesh format");
					/* read boundaries */
	if (flags & IO_BOUNDS) {
		for (i = 0; i < 3; i++)
			mp->mcube.cuorg[i] = atof(getstr(sbuf, meshfp));
		mp->mcube.cusize = atof(getstr(sbuf, meshfp));
		for (i = 0; i < 2; i++) {
			mp->uvlim[0][i] = mgetflt();
			mp->uvlim[1][i] = mgetflt();
		}
	} else {
		for (i = 0; i < 4; i++)
			getstr(sbuf, meshfp);
		for (i = 0; i < 4; i++)
			mgetflt();
	}
					/* read the octree */
	if (flags & IO_TREE)
		mp->mcube.cutree = gettree();
	else if (flags & IO_SCENE)
		skiptree();
					/* read materials and patches */
	if (flags & IO_SCENE) {
		mp->mat0 = nobjects;
		readscene(meshfp, objsize);
		mp->nmats = nobjects - mp->mat0;
		mp->npatches = mgetint(4);
		mp->patch = (MESHPATCH *)calloc(mp->npatches,
					sizeof(MESHPATCH));
		if (mp->patch == NULL)
			mesherror(SYSTEM, "out of patch memory");
		for (i = 0; i < mp->npatches; i++)
			getpatch(&mp->patch[i]);
	}
					/* clean up */
	fclose(meshfp);
	mp->ldflags |= flags;
					/* verify data */
	if ((err = checkmesh(mp)) != NULL)
		mesherror(USER, err);
}
