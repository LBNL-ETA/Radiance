#ifndef lint
static const char RCSid[] = "$Id: mesh.c,v 2.2 2003/03/12 04:59:04 greg Exp $";
#endif
/*
 * Mesh support routines
 */

#include "standard.h"
#include "octree.h"
#include "object.h"
#include "mesh.h"
#include "lookup.h"

/* An encoded mesh vertex */
typedef struct {
	int		fl;
	uint4		xyz[3];
	int4		norm;
	uint4		uv[2];
} MCVERT;

#define  MPATCHBLKSIZ	128		/* patch allocation block size */

#define  IO_LEGAL	(IO_BOUNDS|IO_TREE|IO_SCENE)

static MESH	*mlist = NULL;		/* list of loaded meshes */


static unsigned long
cvhash(cvp)				/* hash an encoded vertex */
MCVERT	*cvp;
{
	unsigned long	hval;
	
	if (!(cvp->fl & MT_V))
		return(0);
	hval = cvp->xyz[0] ^ cvp->xyz[1] << 11 ^ cvp->xyz[2] << 22;
	if (cvp->fl & MT_N)
		hval ^= cvp->norm;
	if (cvp->fl & MT_UV)
		hval ^= cvp->uv[0] ^ cvp->uv[1] << 16;
	return(hval);
}


static int
cvcmp(v1, v2)				/* compare encoded vertices */
register MCVERT	*v1, *v2;
{
	if (v1->fl != v2->fl)
		return(1);
	if (v1->xyz[0] != v2->xyz[0])
		return(1);
	if (v1->xyz[1] != v2->xyz[1])
		return(1);
	if (v1->xyz[2] != v2->xyz[2])
		return(1);
	if (v1->fl & MT_N && v1->norm != v2->norm)
		return(1);
	if (v1->fl & MT_UV) {
		if (v1->uv[0] != v2->uv[0])
			return(1);
		if (v1->uv[1] != v2->uv[1])
			return(1);
	}
	return(0);
}


MESH *
getmesh(mname, flags)			/* get mesh data */
char	*mname;
int	flags;
{
	char  *pathname;
	register MESH  *ms;

	flags &= IO_LEGAL;
	for (ms = mlist; ms != NULL; ms = ms->next)
		if (!strcmp(mname, ms->name)) {
			if ((ms->ldflags & flags) == flags) {
				ms->nref++;
				return(ms);		/* loaded */
			}
			break;			/* load the rest */
		}
	if (ms == NULL) {
		ms = (MESH *)calloc(1, sizeof(MESH));
		if (ms == NULL)
			error(SYSTEM, "out of memory in getmesh");
		ms->name = savestr(mname);
		ms->nref = 1;
		ms->mcube.cutree = EMPTY;
		ms->next = mlist;
		mlist = ms;
	}
	if ((pathname = getpath(mname, getlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find mesh file \"%s\"", mname);
		error(USER, errmsg);
	}
	flags &= ~ms->ldflags;
	if (flags)
		readmesh(ms, pathname, flags);
	return(ms);
}


MESHINST *
getmeshinst(o, flags)			/* create mesh instance */
OBJREC	*o;
int	flags;
{
	register MESHINST  *ins;

	flags &= IO_LEGAL;
	if ((ins = (MESHINST *)o->os) == NULL) {
		if ((ins = (MESHINST *)malloc(sizeof(MESHINST))) == NULL)
			error(SYSTEM, "out of memory in getmeshinst");
		if (o->oargs.nsargs < 1)
			objerror(o, USER, "bad # of arguments");
		if (fullxf(&ins->x, o->oargs.nsargs-1,
				o->oargs.sarg+1) != o->oargs.nsargs-1)
			objerror(o, USER, "bad transform");
		if (ins->x.f.sca < 0.0) {
			ins->x.f.sca = -ins->x.f.sca;
			ins->x.b.sca = -ins->x.b.sca;
		}
		ins->msh = NULL;
		o->os = (char *)ins;
	}
	if (ins->msh == NULL || (ins->msh->ldflags & flags) != flags)
		ins->msh = getmesh(o->oargs.sarg[0], flags);
	return(ins);
}


int
getmeshtrivid(tvid, mp, ti)		/* get triangle vertex ID's */
int4		tvid[3];
register MESH	*mp;
OBJECT		ti;
{
	int		pn = ti >> 10;

	if (pn >= mp->npatches)
		return(0);
	ti &= 0x3ff;
	if (!(ti & 0x200)) {		/* local triangle */
		struct PTri	*tp;
		if (ti >= mp->patch[pn].ntris)
			return(0);
		tp = &mp->patch[pn].tri[ti];
		tvid[0] = tvid[1] = tvid[2] = pn << 8;
		tvid[0] |= tp->v1;
		tvid[1] |= tp->v2;
		tvid[2] |= tp->v3;
		return(1);
	}
	ti &= ~0x200;
	if (!(ti & 0x100)) {		/* single link vertex */
		struct PJoin1	*tp1;
		if (ti >= mp->patch[pn].nj1tris)
			return(0);
		tp1 = &mp->patch[pn].j1tri[ti];
		tvid[0] = tp1->v1j;
		tvid[1] = tvid[2] = pn << 8;
		tvid[1] |= tp1->v2;
		tvid[2] |= tp1->v3;
		return(1);
	}
	ti &= ~0x100;
	{				/* double link vertex */
		struct PJoin2	*tp2;
		if (ti >= mp->patch[pn].nj2tris)
			return(0);
		tp2 = &mp->patch[pn].j2tri[ti];
		tvid[0] = tp2->v1j;
		tvid[1] = tp2->v2j;
		tvid[2] = pn << 8 | tp2->v3;
	}
	return(1);
}


int
getmeshvert(vp, mp, vid, what)	/* get triangle vertex from ID */
MESHVERT	*vp;
register MESH	*mp;
int4		vid;
int		what;
{
	int		pn = vid >> 8;
	MESHPATCH	*pp;
	double		vres;
	register int	i;
	
	vp->fl = 0;
	if (pn >= mp->npatches)
		return(0);
	pp = &mp->patch[pn];
	vid &= 0xff;
	if (vid >= pp->nverts)
		return(0);
					/* get location */
	if (what & MT_V) {
		vres = (1./4294967296.)*mp->mcube.cusize;
		for (i = 0; i < 3; i++)
			vp->v[i] = mp->mcube.cuorg[i] +
					(pp->xyz[vid][i] + .5)*vres;
		vp->fl |= MT_V;
	}
					/* get normal */
	if (what & MT_N && pp->norm != NULL && pp->norm[vid]) {
		decodedir(vp->n, pp->norm[vid]);
		vp->fl |= MT_N;
	}
					/* get (u,v) */
	if (what & MT_UV && pp->uv != NULL && pp->uv[vid][0]) {
		for (i = 0; i < 2; i++)
			vp->uv[i] = mp->uvlim[0][i] +
				(mp->uvlim[1][i] - mp->uvlim[0][i])*
				(pp->uv[vid][i] + .5)*(1./4294967296.);
		vp->fl |= MT_UV;
	}
	return(vp->fl);
}


int
getmeshtri(tv, mp, ti, what)	/* get triangle vertices */
MESHVERT	tv[3];
MESH		*mp;
OBJECT		ti;
int		what;
{
	int4	tvid[3];

	if (!getmeshtrivid(tvid, mp, ti))
		return(0);

	getmeshvert(&tv[0], mp, tvid[0], what);
	getmeshvert(&tv[1], mp, tvid[1], what);
	getmeshvert(&tv[2], mp, tvid[2], what);

	return(tv[0].fl & tv[1].fl & tv[2].fl);
}


int4
addmeshvert(mp, vp)		/* find/add a mesh vertex */
register MESH	*mp;
MESHVERT	*vp;
{
	LUTAB		*ltp;
	LUENT		*lvp;
	MCVERT		cv;
	register int	i;

	if (!(vp->fl & MT_V))
		return(-1);
					/* encode vertex */
	for (i = 0; i < 3; i++) {
		if (vp->v[i] < mp->mcube.cuorg[i])
			return(-1);
		if (vp->v[i] >= mp->mcube.cuorg[i] + mp->mcube.cusize)
			return(-1);
		cv.xyz[i] = (uint4)(4294967296. *
				(vp->v[i] - mp->mcube.cuorg[i]) /
				mp->mcube.cusize);
	}
	if (vp->fl & MT_N)
		cv.norm = encodedir(vp->n);
	if (vp->fl & MT_UV)
		for (i = 0; i < 2; i++) {
			if (vp->uv[i] <= mp->uvlim[0][i])
				return(-1);
			if (vp->uv[i] >= mp->uvlim[1][i])
				return(-1);
			cv.uv[i] = (uint4)(4294967296. *
					(vp->uv[i] - mp->uvlim[0][i]) /
					(mp->uvlim[1][i] - mp->uvlim[0][i]));
		}
	cv.fl = vp->fl;
	ltp = (LUTAB *)mp->cdata;	/* get lookup table */
	if (ltp == NULL) {
		ltp = (LUTAB *)calloc(1, sizeof(LUTAB));
		if (ltp == NULL)
			goto nomem;
		ltp->hashf = cvhash;
		ltp->keycmp = cvcmp;
		ltp->freek = free;
		if (!lu_init(ltp, 50000))
			goto nomem;
		mp->cdata = (char *)ltp;
	}
					/* find entry */
	lvp = lu_find(ltp, (char *)&cv);
	if (lvp == NULL)
		goto nomem;
	if (lvp->key == NULL) {
		lvp->key = (char *)malloc(sizeof(MCVERT)+sizeof(int4));
		bcopy((void *)&cv, (void *)lvp->key, sizeof(MCVERT));
	}
	if (lvp->data == NULL) {	/* new vertex */
		register MESHPATCH	*pp;
		if (mp->npatches <= 0) {
			mp->patch = (MESHPATCH *)calloc(MPATCHBLKSIZ,
					sizeof(MESHPATCH));
			if (mp->patch == NULL)
				goto nomem;
			mp->npatches = 1;
		} else if (mp->patch[mp->npatches-1].nverts >= 256) {
			if (mp->npatches % MPATCHBLKSIZ == 0) {
				mp->patch = (MESHPATCH *)realloc(
						(void *)mp->patch,
					(mp->npatches + MPATCHBLKSIZ)*
						sizeof(MESHPATCH));
				bzero((void *)(mp->patch + mp->npatches),
					MPATCHBLKSIZ*sizeof(MESHPATCH));
			}
			if (mp->npatches++ >= 1L<<22)
				error(INTERNAL, "too many mesh patches");
		}
		pp = &mp->patch[mp->npatches-1];
		if (pp->xyz == NULL) {
			pp->xyz = (uint4 (*)[3])calloc(256, 3*sizeof(int4));
			if (pp->xyz == NULL)
				goto nomem;
		}
		for (i = 0; i < 3; i++)
			pp->xyz[pp->nverts][i] = cv.xyz[i];
		if (cv.fl & MT_N) {
			if (pp->norm == NULL) {
				pp->norm = (int4 *)calloc(256, sizeof(int4));
				if (pp->norm == NULL)
					goto nomem;
			}
			pp->norm[pp->nverts] = cv.norm;
		}
		if (cv.fl & MT_UV) {
			if (pp->uv == NULL) {
				pp->uv = (uint4 (*)[2])calloc(256,
						2*sizeof(uint4));
				if (pp->uv == NULL)
					goto nomem;
			}
			for (i = 0; i < 2; i++)
				pp->uv[pp->nverts][i] = cv.uv[i];
		}
		pp->nverts++;
		lvp->data = lvp->key + sizeof(MCVERT);
		*(int4 *)lvp->data = (mp->npatches-1) << 8 | (pp->nverts-1);
	}
	return(*(int4 *)lvp->data);
nomem:
	error(SYSTEM, "out of memory in addmeshvert");
	return(-1);
}


OBJECT
addmeshtri(mp, tv)		/* add a new mesh triangle */
MESH		*mp;
MESHVERT	tv[3];
{
	int4			vid[3], t;
	int			pn[3], i;
	register MESHPATCH	*pp;

	if (!(tv[0].fl & tv[1].fl & tv[2].fl & MT_V))
		return(OVOID);
				/* find/allocate patch vertices */
	for (i = 0; i < 3; i++) {
		if ((vid[i] = addmeshvert(mp, &tv[i])) < 0)
			return(OVOID);
		pn[i] = vid[i] >> 8;
	}
				/* assign triangle */
	if (pn[0] == pn[1] && pn[1] == pn[2]) {	/* local case */
		pp = &mp->patch[pn[0]];
		if (pp->tri == NULL) {
			pp->tri = (struct PTri *)malloc(
					512*sizeof(struct PTri));
			if (pp->tri == NULL)
				goto nomem;
		}
		if (pp->ntris < 512) {
			pp->tri[pp->ntris].v1 = vid[0] & 0xff;
			pp->tri[pp->ntris].v2 = vid[1] & 0xff;
			pp->tri[pp->ntris].v3 = vid[2] & 0xff;
			return(pn[0] << 10 | pp->ntris++);
		}
	}
	if (pn[0] == pn[1]) {
		t = vid[2]; vid[2] = vid[1]; vid[1] = vid[0]; vid[0] = t;
		i = pn[2]; pn[2] = pn[1]; pn[1] = pn[0]; pn[0] = i;
	} else if (pn[0] == pn[2]) {
		t = vid[0]; vid[0] = vid[1]; vid[1] = vid[2]; vid[2] = t;
		i = pn[0]; pn[0] = pn[1]; pn[1] = pn[2]; pn[2] = i;
	}
	if (pn[1] == pn[2]) {			/* single link */
		pp = &mp->patch[pn[1]];
		if (pp->j1tri == NULL) {
			pp->j1tri = (struct PJoin1 *)malloc(
					256*sizeof(struct PJoin1));
			if (pp->j1tri == NULL)
				goto nomem;
		}
		if (pp->nj1tris < 256) {
			pp->j1tri[pp->nj1tris].v1j = vid[0];
			pp->j1tri[pp->nj1tris].v2 = vid[1] & 0xff;
			pp->j1tri[pp->nj1tris].v3 = vid[2] & 0xff;
			return(pn[1] << 10 | 0x200 | pp->nj1tris++);
		}
	}
						/* double link */
	pp = &mp->patch[pn[2]];
	if (pp->j2tri == NULL) {
		pp->j2tri = (struct PJoin2 *)malloc(
					256*sizeof(struct PJoin2));
		if (pp->j2tri == NULL)
			goto nomem;
	}
	if (pp->nj2tris >= 256)
		error(INTERNAL, "too many patch triangles in addmeshtri");
	pp->j2tri[pp->nj2tris].v1j = vid[0];
	pp->j2tri[pp->nj2tris].v2j = vid[1];
	pp->j2tri[pp->nj2tris].v3 = vid[2] & 0xff;
	return(pn[2] << 10 | 0x300 | pp->nj2tris++);
nomem:
	error(SYSTEM, "out of memory in addmeshtri");
	return(OVOID);
}


static void
tallyoctree(ot, ecp, lcp, ocp)	/* tally octree size */
OCTREE	ot;
int	*ecp, *lcp, *ocp;
{
	int	i;

	if (isempty(ot)) {
		(*ecp)++;
		return;
	}
	if (isfull(ot)) {
		OBJECT	oset[MAXSET+1];
		(*lcp)++;
		objset(oset, ot);
		*ocp += oset[0];
		return;
	}
	for (i = 0; i < 8; i++)
		tallyoctree(octkid(ot, i), ecp, lcp, ocp);
}


void
printmeshstats(ms, fp)		/* print out mesh statistics */
MESH	*ms;
FILE	*fp;
{
	int	lfcnt=0, lecnt=0, locnt=0;
	int	vcnt=0, ncnt=0, uvcnt=0;
	int	nscnt=0, uvscnt=0;
	int	tcnt=0, t1cnt=0, t2cnt=0;
	int	i, j;
	
	tallyoctree(ms->mcube.cutree, &lecnt, &lfcnt, &locnt);
	for (i = 0; i < ms->npatches; i++) {
		register MESHPATCH	*pp = &ms->patch[i];
		vcnt += pp->nverts;
		if (pp->norm != NULL) {
			for (j = pp->nverts; j--; )
				if (pp->norm[j])
					ncnt++;
			nscnt += pp->nverts;
		}
		if (pp->uv != NULL) {
			for (j = pp->nverts; j--; )
				if (pp->uv[j][0])
					uvcnt++;
			uvscnt += pp->nverts;
		}
		tcnt += pp->ntris;
		t1cnt += pp->nj1tris;
		t2cnt += pp->nj2tris;
	}
	fprintf(fp, "Mesh statistics:\n");
	fprintf(fp, "\t%d patches (%.2f MBytes)\n", ms->npatches,
			(ms->npatches*sizeof(MESHPATCH) +
			vcnt*3*sizeof(uint4) +
			nscnt*sizeof(int4) +
			uvscnt*2*sizeof(uint4) +
			tcnt*sizeof(struct PTri) +
			t1cnt*sizeof(struct PJoin1) +
			t2cnt*sizeof(struct PJoin2))/(1024.*1024.));
	fprintf(fp, "\t%d vertices (%.1f%% w/ normals, %.1f%% w/ uv)\n",
			vcnt, 100.*ncnt/vcnt, 100.*uvcnt/vcnt);
	fprintf(fp, "\t%d triangles (%.1f%% local, %.1f%% joiner)\n",
			tcnt+t1cnt+t2cnt,
			100.*tcnt/(tcnt+t1cnt+t2cnt),
			100.*t1cnt/(tcnt+t1cnt+t2cnt));
	fprintf(fp,
	"\t%d leaves in octree (%.1f%% empty, %.2f avg. set size)\n",
			lfcnt+lecnt, 100.*lecnt/(lfcnt+lecnt),
			(double)locnt/lfcnt);
}


void
freemesh(ms)			/* free mesh data */
MESH	*ms;
{
	MESH	mhead;
	MESH	*msp;
	
	if (ms == NULL)
		return;
	if (ms->nref <= 0)
		error(CONSISTENCY, "unreferenced mesh in freemesh");
	ms->nref--;
	if (ms->nref)		/* still in use */
		return;
				/* else remove from list */
	mhead.next = mlist;
	for (msp = &mhead; msp->next != NULL; msp = msp->next)
		if (msp->next == ms) {
			msp->next = ms->next;
			ms->next = NULL;
			break;
		}
	if (ms->next != NULL)	/* can't be in list anymore */
		error(CONSISTENCY, "unlisted mesh in freemesh");
	mlist = mhead.next;
				/* free mesh data */
	freestr(ms->name);
	octfree(ms->mcube.cutree);
	if (ms->cdata != NULL)
		lu_done((LUTAB *)ms->cdata);
	if (ms->npatches > 0) {
		register MESHPATCH	*pp = ms->patch + ms->npatches;
		while (pp-- > ms->patch) {
			if (pp->j2tri != NULL)
				free((void *)pp->j2tri);
			if (pp->j1tri != NULL)
				free((void *)pp->j1tri);
			if (pp->tri != NULL)
				free((void *)pp->tri);
			if (pp->uv != NULL)
				free((void *)pp->uv);
			if (pp->norm != NULL)
				free((void *)pp->norm);
			if (pp->xyz != NULL)
				free((void *)pp->xyz);
		}
		free((void *)ms->patch);
	}
	free((void *)ms);
}


void
freemeshinst(o)			/* free mesh instance */
OBJREC	*o;
{
	if (o->os == NULL)
		return;
	freemesh((*(MESHINST *)o->os).msh);
	free((void *)o->os);
	o->os = NULL;
}
