#ifndef lint
static const char RCSid[] = "$Id: writemesh.c,v 2.5 2004/03/27 12:41:45 schorsch Exp $";
#endif
/*
 *  Routines for writing compiled mesh to a file stream
 */

#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "mesh.h"


static void putfullnode(OCTREE fn, FILE *fp);
static void puttree(OCTREE ot, FILE *fp);
static void putpatch(MESHPATCH *pp, FILE *fp);


static void
putfullnode(		/* write out a full node */
	OCTREE	fn,
	FILE	*fp
)
{
	OBJECT  oset[MAXSET+1];
	register int  i;

	objset(oset, fn);
	for (i = 0; i <= oset[0]; i++)
		putint((long)oset[i], sizeof(OBJECT), fp);
}


static void
puttree(			/* write octree to fp in pre-order form */
	register OCTREE	ot,
	FILE		*fp
)
{
	
	if (istree(ot)) {
		register int	i;
		putc(OT_TREE, fp);		/* indicate tree */
		for (i = 0; i < 8; i++)		/* write tree */
			puttree(octkid(ot, i), fp);
		return;
	}
	if (isfull(ot)) {
		putc(OT_FULL, fp);		/* indicate fullnode */
		putfullnode(ot, fp);		/* write fullnode */
		return;
	}
	putc(OT_EMPTY, fp);			/* indicate empty */
}


static void
putpatch(			/* write out a mesh patch */
	register MESHPATCH	*pp,
	FILE			*fp
)
{
	int	flags = MT_V;
	int	i, j;
					/* vertex flags */
	if (pp->norm != NULL)
		flags |= MT_N;
	if (pp->uv != NULL)
		flags |= MT_UV;
	putint((long)flags, 1, fp);
					/* number of vertices */
	putint((long)pp->nverts, 2, fp);
					/* vertex xyz locations */
	for (i = 0; i < pp->nverts; i++)
		for (j = 0; j < 3; j++)
			putint((long)pp->xyz[i][j], 4, fp);
					/* vertex normals */
	if (flags & MT_N)
		for (i = 0; i < pp->nverts; i++)
			putint((long)pp->norm[i], 4, fp);
					/* uv coordinates */
	if (flags & MT_UV)
		for (i = 0; i < pp->nverts; i++)
			for (j = 0; j < 2; j++)
				putint((long)pp->uv[i][j], 4, fp);
					/* local triangles */
	putint((long)pp->ntris, 2, fp);
	for (i = 0; i < pp->ntris; i++) {
		putint((long)pp->tri[i].v1, 1, fp);
		putint((long)pp->tri[i].v2, 1, fp);
		putint((long)pp->tri[i].v3, 1, fp);
	}
					/* local triangle material(s) */
	if (pp->trimat == NULL) {
		putint(1L, 2, fp);
		putint((long)pp->solemat, 2, fp);
	} else {
		putint((long)pp->ntris, 2, fp);
		for (i = 0; i < pp->ntris; i++)
			putint((long)pp->trimat[i], 2, fp);
	}
					/* joiner triangles */
	putint((long)pp->nj1tris, 2, fp);
	for (i = 0; i < pp->nj1tris; i++) {
		putint((long)pp->j1tri[i].v1j, 4, fp);
		putint((long)pp->j1tri[i].v2, 1, fp);
		putint((long)pp->j1tri[i].v3, 1, fp);
		putint((long)pp->j1tri[i].mat, 2, fp);
	}
					/* double joiner triangles */
	putint((long)pp->nj2tris, 2, fp);
	for (i = 0; i < pp->nj2tris; i++) {
		putint((long)pp->j2tri[i].v1j, 4, fp);
		putint((long)pp->j2tri[i].v2j, 4, fp);
		putint((long)pp->j2tri[i].v3, 1, fp);
		putint((long)pp->j2tri[i].mat, 2, fp);
	}
}


void
writemesh(mp, fp)			/* write mesh structures to fp */
MESH	*mp;
FILE	*fp;
{
	char	*err;
	char	sbuf[64];
	int	i;
					/* do we have everything? */
	if ((mp->ldflags & (IO_SCENE|IO_TREE|IO_BOUNDS)) !=
			(IO_SCENE|IO_TREE|IO_BOUNDS))
		error(INTERNAL, "missing data in writemesh");
					/* validate mesh data */
	if ((err = checkmesh(mp)) != NULL)
		error(USER, err);
					/* write format number */
	putint((long)(MESHMAGIC+sizeof(OBJECT)), 2, fp);
					/* write boundaries */
	for (i = 0; i < 3; i++) {
		sprintf(sbuf, "%.12g", mp->mcube.cuorg[i]);
		putstr(sbuf, fp);
	}
	sprintf(sbuf, "%.12g", mp->mcube.cusize);
	putstr(sbuf, fp);
	for (i = 0; i < 2; i++) {
		putflt(mp->uvlim[0][i], fp);
		putflt(mp->uvlim[1][i], fp);
	}
					/* write the octree */
	puttree(mp->mcube.cutree, fp);
					/* write the materials */
	writescene(mp->mat0, mp->nmats, fp);
					/* write the patches */
	putint((long)mp->npatches, 4, fp);
	for (i = 0; i < mp->npatches; i++)
		putpatch(&mp->patch[i], fp);
	if (ferror(fp))
		error(SYSTEM, "write error in writemesh");
}
