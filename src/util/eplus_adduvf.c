#ifndef lint
static const char RCSid[] = "$Id: eplus_adduvf.c,v 2.17 2014/04/04 17:56:45 greg Exp $";
#endif
/*
 * Add User View Factors to EnergyPlus Input Data File
 *
 *	G.Ward for LBNL
 */

#include <stdlib.h>
#include "rtio.h"
#include "rtmath.h"
#include "random.h"
#include "eplus_idf.h"
#include "triangulate.h"
#include "rtprocess.h"

#ifndef NSAMPLES
#define NSAMPLES	80000			/* default number of samples */
#endif

#define	SURF_EPS	0.0005			/* surface testing epsilon */

char		*progname;			/* global argv[0] */

int		nsamps = NSAMPLES;		/* number of samples to use */

char		temp_octree[128];		/* temporary octree */

const char	UVF_PNAME[] =
			"ZoneProperty:UserViewFactors:bySurfaceName";

const char	SUBSURF_PNAME[] =
			"FenestrationSurface:Detailed";

const char	ADD_HEADER[] =
			"\n!+++ User View Factors computed by Radiance +++!\n\n";

#define NAME_FLD	1			/* name field always first? */

#define SS_BASE_FLD	4			/* subsurface base surface */
#define SS_VERT_FLD	10			/* subsurface vertex count */

typedef struct {
	const char	*pname;			/* object type name */
	short		zone_fld;		/* zone field index */
	short		vert_fld;		/* vertex field index */
} SURF_PTYPE;		/* surface type we're interested in */

const SURF_PTYPE	surf_type[] = {
		{"BuildingSurface:Detailed", 4, 10},
		{"Floor:Detailed", 3, 9},
		{"RoofCeiling:Detailed", 3, 9},
		{"Wall:Detailed", 3, 9},
		{NULL}
	};				/* IDF surface types */

typedef struct s_zone {
	const char	*zname;			/* zone name */
	struct s_zone	*next;			/* next zone in list */
	int		nsurf;			/* surface count */
	int		ntotal;			/* surfaces+subsurfaces */
	IDF_OBJECT	*pfirst;		/* first matching object */
	IDF_OBJECT	*plast;			/* last before subsurfaces */
	float		*area_redu;		/* subsurface area per surf. */
} ZONE;			/* a list of collected zone surfaces */

ZONE		*zone_list = NULL;	/* our list of zones */

LUTAB		zonesurf_lut =		/* zone surface lookup table */
			LU_SINIT(NULL,NULL);

IDF_LOADED	*our_idf = NULL;	/* loaded/modified IDF */

typedef struct {
	FVECT		norm;		/* surface normal */
	double		area;		/* surface area */
	int		nv;		/* number of vertices */
	FVECT		vl[3];		/* vertex list (extends struct) */
} SURFACE;		/* a polygonal surface */

typedef struct {
	FVECT		sdir[3];	/* UVW unit sampling vectors */
	double		poff;		/* W-offset for plane of polygon */
	double		area_left;	/* area left to sample */
	int		samp_left;	/* remaining samples */
	int		wd;		/* output file descriptor */
} POLYSAMP;		/* structure for polygon sampling */

/* Create a new zone and push to top of our list */
static ZONE *
new_zone(const char *zname, IDF_OBJECT *param)
{
	ZONE	*znew = (ZONE *)malloc(sizeof(ZONE));

	if (znew == NULL)
		return(NULL);
	znew->zname = zname;			/* assumes static */
	znew->next = zone_list;
	znew->pfirst = znew->plast = param;
	znew->ntotal = znew->nsurf = 1;
	znew->area_redu = NULL;
	return(zone_list = znew);
}

/* Add the detailed surface (polygon) to the named zone */
static ZONE *
add2zone(IDF_OBJECT *param, const char *zname)
{
	IDF_FIELD	*nfp = idf_getfield(param, NAME_FLD);
	ZONE		*zptr;
	LUENT		*lep;

	if (nfp == NULL) {
		fputs(progname, stderr);
		fputs(": surface missing name field!\n", stderr);
		return(NULL);
	}
	for (zptr = zone_list; zptr != NULL; zptr = zptr->next)
		if (!strcmp(zptr->zname, zname))
			break;
	if (zptr == NULL) {
		zptr = new_zone(zname, param);
	} else {				/* keep surfaces together */
		if (!idf_movobject(our_idf, param, zptr->plast))
			return(NULL);
		zptr->plast = param;
		zptr->nsurf++;
		zptr->ntotal++;
	}
						/* add to lookup table */
	lep = lu_find(&zonesurf_lut, nfp->val);
	if (lep == NULL)
		return(NULL);
	if (lep->data != NULL) {
		fputs(progname, stderr);
		fputs(": duplicate surface name '", stderr);
		fputs(nfp->val, stderr);
		fputs("'\n", stderr);
		return(NULL);
	}
	lep->key = nfp->val;
	lep->data = (char *)zptr;
	return(zptr);
}

/* Add a subsurface by finding its base surface and the corresponding zone */
static ZONE *
add_subsurf(IDF_OBJECT *param)
{
	IDF_FIELD	*bfp = idf_getfield(param, SS_BASE_FLD);
	ZONE		*zptr;
	LUENT		*lep;

	if (bfp == NULL) {
		fputs(progname, stderr);
		fputs(": missing base field name in subsurface!\n", stderr);
		return(NULL);
	}
	lep = lu_find(&zonesurf_lut, bfp->val);	/* find base zone */
	if (lep == NULL || lep->data == NULL) {
		fputs(progname, stderr);
		fputs(": cannot find referenced base surface '", stderr);
		fputs(bfp->val, stderr);
		fputs("'\n", stderr);
		return(NULL);
	}
	zptr = (ZONE *)lep->data;		/* add this subsurface */
	if (!idf_movobject(our_idf, param, zptr->plast))
		return(NULL);
	zptr->ntotal++;
	return(zptr);
}

/* Return field for vertices in the given object */
static IDF_FIELD *
get_vlist(IDF_OBJECT *param, const char *zname)
{
#define	itm_len		(sizeof(IDF_FIELD)+6)
	static char		fld_buf[4*itm_len];
	static char		*next_fbp = fld_buf;
	int			i;
	IDF_FIELD		*res;
						/* check if subsurface */
	if (!strcmp(param->pname, SUBSURF_PNAME)) {
		if (zname != NULL) {
			LUENT	*lep = lu_find(&zonesurf_lut,
					idf_getfield(param,SS_BASE_FLD)->val);
			if (strcmp((*(ZONE *)lep->data).zname, zname))
				return(NULL);
		}
		res = idf_getfield(param, SS_VERT_FLD);
	} else {
		i = 0;				/* check for surface type */
		while (strcmp(surf_type[i].pname, param->pname))
			if (surf_type[++i].pname == NULL)
				return(NULL);

		if (zname != NULL) {		/* matches specified zone? */
			IDF_FIELD	*fptr = idf_getfield(param, surf_type[i].zone_fld);
			if (fptr == NULL || strcmp(fptr->val, zname))
				return(NULL);
		}
		res = idf_getfield(param, surf_type[i].vert_fld);
	}
	if (!res->val[0]) {			/* hack for missing #vert */
		IDF_FIELD	*fptr;
		if (next_fbp >= fld_buf+sizeof(fld_buf))
			next_fbp = fld_buf;
		i = 0;				/* count vertices */
		for (fptr = res->next; fptr != NULL; fptr = fptr->next)
			++i;
		if (i % 3)
			return(NULL);
		fptr = res->next;
		res = (IDF_FIELD *)next_fbp; next_fbp += itm_len;
		res->next = fptr;
		res->rem = "";
		sprintf(res->val, "%d", i/3);
	}
	return(res);
#undef itm_len
}

/* Get/allocate surface polygon */
static SURFACE *
get_surface(IDF_FIELD *fptr)
{
	SURFACE	*surf;
	int	nv;
	FVECT	e1, e2, vc;
	int	i, j;
					/* get number of vertices */
	if (fptr == NULL || (nv = atoi(fptr->val)) < 3) {
		fputs(progname, stderr);
		fputs(": bad vertex count for surface!\n", stderr);
		return(NULL);
	}
	surf = (SURFACE *)malloc(sizeof(SURFACE)+sizeof(FVECT)*(nv-3));
	if (surf == NULL)
		return(NULL);
	surf->nv = nv;
	for (i = nv; i--; )		/* reverse vertex order/orientation */
		for (j = 0; j < 3; j++) {
			fptr = fptr->next;
			if (fptr == NULL) {
				fputs(progname, stderr);
				fputs(": missing vertex in get_surface()\n", stderr);
				free(surf);
				return(NULL);
			}
			if ((surf->vl[i][j] = atof(fptr->val)) == 0 &&
							!isflt(fptr->val)) {
				fputs(progname, stderr);
				fputs(": bad vertex in get_surface()\n", stderr);
				free(surf);
				return(NULL);
			}
		}
					/* compute area and normal */
	surf->norm[0] = surf->norm[1] = surf->norm[2] = 0;
	VSUB(e1, surf->vl[1], surf->vl[0]);
	for (i = 2; i < nv; i++) {
		VSUB(e2, surf->vl[i], surf->vl[0]); 
		fcross(vc, e1, e2);
		surf->norm[0] += vc[0];
		surf->norm[1] += vc[1];
		surf->norm[2] += vc[2];
		VCOPY(e1, e2);
	}
	surf->area = .5 * normalize(surf->norm);
	if (surf->area == 0) {
		fputs(progname, stderr);
		fputs(": degenerate polygon in get_surface()\n", stderr);
		free(surf);
		return(NULL);
	}
	return(surf);
}

/* Convert surface to Radiance with modifier based on unique name */
static int
rad_surface(IDF_OBJECT *param, FILE *ofp)
{
	const char	*sname = idf_getfield(param,NAME_FLD)->val;
	IDF_FIELD	*fptr = get_vlist(param, NULL);
	const char	**fargs;
	int		nvert, i, j;

	if (fptr == NULL || (nvert = atoi(fptr->val)) < 3) {
		fprintf(stderr, "%s: bad surface '%s'\n", progname, sname);
		return(0);
	}
	fargs = (const char **)malloc(sizeof(const char *)*3*nvert);
	if (fargs == NULL)
		return(0);
	fprintf(ofp, "\nvoid glow '%s'\n0\n0\n4 1 1 1 0\n", sname);
	fprintf(ofp, "\n'%s' polygon 's_%s'\n0\n0\n%d\n", sname, sname, 3*nvert);
	for (j = nvert; j--; ) {		/* get arguments in reverse */
		for (i = 0; i < 3; i++) {
			fptr = fptr->next;
			if (fptr == NULL || !isflt(fptr->val)) {
				fprintf(stderr,
					"%s: missing/bad vertex for %s '%s'\n",
						progname, param->pname, sname);
				return(0);
			}
			fargs[3*j + i] = fptr->val;
		}
	}
	for (j = 0; j < nvert; j++) {		/* output reversed verts */
		for (i = 0; i < 3; i++) {
			fputc('\t', ofp);
			fputs(fargs[3*j + i], ofp);
		}
		fputc('\n', ofp);
	}
	free(fargs);
	return(!ferror(ofp));
}

/* Convert subsurface to Radiance with modifier based on unique name */
static double
rad_subsurface(IDF_OBJECT *param, FILE *ofp)
{
	const char	*sname = idf_getfield(param,NAME_FLD)->val;
	SURFACE		*surf = get_surface(get_vlist(param, NULL));
	double		area;
	int		i;

	if (surf == NULL) {
		fprintf(stderr, "%s: bad subsurface '%s'\n", progname, sname);
		return(-1.);
	}
	fprintf(ofp, "\nvoid glow '%s'\n0\n0\n4 1 1 1 0\n", sname);
	fprintf(ofp, "\n'%s' polygon 'ss_%s'\n0\n0\n%d\n",
						sname, sname, 3*surf->nv);
	for (i = 0; i < surf->nv; i++) {	/* offset vertices */
		FVECT	vert;
		VSUM(vert, surf->vl[i], surf->norm, 2.*SURF_EPS);
		fprintf(ofp, "\t%.12f %.12f %.12f\n", vert[0], vert[1], vert[2]);
	}
	area = surf->area;
	free(surf);
	if (ferror(ofp))
		return(-1.);
	return(area);
}

/* Start rcontrib process */
static int
start_rcontrib(SUBPROC *pd, ZONE *zp)
{
#define	BASE_AC		6
	static char	*base_av[BASE_AC] = {
				"rcontrib", "-V+", "-fff", "-h", "-x", "1"
			};
	char		cbuf[300];
	char		**av;
	FILE		*ofp;
	IDF_OBJECT	*pptr;
	IDF_FIELD	*fptr;
	int		i, j, n;
						/* start oconv command */
	sprintf(cbuf, "oconv - > \"%s\"", temp_octree);
	if ((ofp = popen(cbuf, "w")) == NULL) {
		fputs(progname, stderr);
		fputs(": cannot open oconv process\n", stderr);
		return(0);
	}
						/* allocate argument list */
	av = (char **)malloc(sizeof(char *)*(BASE_AC+4+2*(zp->ntotal)));
	if (av == NULL)
		return(0);
	for (i = 0; i < BASE_AC; i++)
		av[i] = base_av[i];
	sprintf(cbuf, "%d", nsamps);
	av[i++] = "-c";
	av[i++] = cbuf;				/* add modifiers & surfaces */
	for (n = 0, pptr = zp->pfirst; n < zp->nsurf; n++, pptr = pptr->dnext) {
		fptr = idf_getfield(pptr,NAME_FLD);
		if (fptr == NULL || !fptr->val[0]) {
			fputs(progname, stderr);
			fputs(": missing name for surface object\n", stderr);
			return(0);
		}
		if (!rad_surface(pptr, ofp))	/* add surface to octree */
			return(0);
		av[i++] = "-m";
		av[i++] = fptr->val;
	}
						/* now subsurfaces */
	if (zp->ntotal > zp->nsurf) {
		if (zp->area_redu != NULL)
			memset(zp->area_redu, 0, sizeof(float)*zp->ntotal);
		else if ((zp->area_redu = (float *)calloc(zp->ntotal,
						sizeof(float))) == NULL)
			return(0);
	}
	for ( ; n < zp->ntotal; n++, pptr = pptr->dnext) {
		double		ssarea;
		const char	*bname;
		IDF_OBJECT	*pptr1;
		fptr = idf_getfield(pptr,NAME_FLD);
		if (fptr == NULL || !fptr->val[0]) {
			fputs(progname, stderr);
			fputs(": missing name for subsurface object\n", stderr);
			return(0);
		}
						/* add subsurface to octree */
		if ((ssarea = rad_subsurface(pptr, ofp)) < 0)
			return(0);
						/* mark area for subtraction */
		bname = idf_getfield(pptr,SS_BASE_FLD)->val;
		for (j = 0, pptr1 = zp->pfirst;
				j < zp->nsurf; j++, pptr1 = pptr1->dnext)
			if (!strcmp(idf_getfield(pptr1,NAME_FLD)->val, bname)) {
				zp->area_redu[j] += ssarea;
				break;
			}
		av[i++] = "-m";
		av[i++] = fptr->val;
	}
	if (pclose(ofp) != 0) {			/* finish oconv */
		fputs(progname, stderr);
		fputs(": error running oconv process\n", stderr);
		return(0);
	}
	av[i++] = temp_octree;			/* add final octree argument */
	av[i] = NULL;
	if (open_process(pd, av) <= 0) {	/* start process */
		fputs(progname, stderr);
		fputs(": cannot start rcontrib process\n", stderr);
		return(0);
	}
	free(av);				/* all done -- clean up */
	return(1);
#undef BASE_AC
}

/* Initialize polygon sampling */
static Vert2_list *
init_poly(POLYSAMP *ps, IDF_FIELD *f0)
{
	SURFACE		*surf;
	Vert2_list	*vl2;
	int		i;
					/* get 3-D polygon vertices */
	if ((surf = get_surface(f0)) == NULL)
		return(NULL);
	vl2 = polyAlloc(surf->nv);
	if (vl2 == NULL)
		return(NULL);
					/* create X & Y axes */
	VCOPY(ps->sdir[2], surf->norm);
	VSUB(ps->sdir[0], surf->vl[1], surf->vl[0]);
	if (normalize(ps->sdir[0]) == 0)
		return(NULL);
	fcross(ps->sdir[1], ps->sdir[2], ps->sdir[0]);
					/* compute plane offset */
	ps->poff = DOT(surf->vl[0], ps->sdir[2]);
					/* assign 2-D vertices */
	for (i = 0; i < surf->nv; i++) {
		vl2->v[i].mX = DOT(surf->vl[i], ps->sdir[0]);
		vl2->v[i].mY = DOT(surf->vl[i], ps->sdir[1]);
	}
	ps->area_left = surf->area;
	free(surf);			/* it's ready! */
	vl2->p = (void *)ps;
	return(vl2);
}

/* Generate samples on 2-D triangle */
static int
sample_triangle(const Vert2_list *vl2, int a, int b, int c)
{
	POLYSAMP	*ps = (POLYSAMP *)vl2->p;
	float		*samp;
	FVECT		orig;
	FVECT		ab, ac;
	double		area;
	int		i, j, ns;
					/* compute sampling axes */
	for (i = 3; i--; ) {
		orig[i] = vl2->v[a].mX*ps->sdir[0][i] +
				vl2->v[a].mY*ps->sdir[1][i] +
				(ps->poff+SURF_EPS)*ps->sdir[2][i];
		ab[i] = (vl2->v[b].mX - vl2->v[a].mX)*ps->sdir[0][i] +
				(vl2->v[b].mY - vl2->v[a].mY)*ps->sdir[1][i];
		ac[i] = (vl2->v[c].mX - vl2->v[a].mX)*ps->sdir[0][i] +
				(vl2->v[c].mY - vl2->v[a].mY)*ps->sdir[1][i];
	}
					/* compute number of samples to take */
	area = .5*(vl2->v[a].mX*vl2->v[b].mY - vl2->v[b].mX*vl2->v[a].mY +
			vl2->v[b].mX*vl2->v[c].mY - vl2->v[c].mX*vl2->v[b].mY +
			vl2->v[c].mX*vl2->v[a].mY - vl2->v[a].mX*vl2->v[c].mY);
	if (area < .0) {
		fputs(progname, stderr);
		fputs(": negative triangle area in sample_triangle()\n", stderr);
		return(0);
	}
	if (area >= ps->area_left) {
		ns = ps->samp_left;
		ps->area_left = 0;
	} else {
		ns = (ps->samp_left*area/ps->area_left + .5);
		ps->samp_left -= ns;
		ps->area_left -= area;
	}
	if (ns <= 0)			/* XXX should be error? */
		return(1);
					/* buffer sample rays */
	samp = (float *)malloc(sizeof(float)*6*ns);
	if (samp == NULL)
		return(0);
	for (i = ns; i--; ) {		/* stratified Monte Carlo sampling */
		double	sv[4];
		FVECT	dv;
		multisamp(sv, 4, (i+frandom())/(double)ns);
		sv[0] *= sv[1] = sqrt(sv[1]);
		sv[1] = 1. - sv[1];
		for (j = 3; j--; )
			samp[i*6 + j] = orig[j] + sv[0]*ab[j] + sv[1]*ac[j];
		sv[2] = sqrt(sv[2]);
		sv[3] *= 2.*PI;
		dv[0] = tcos(sv[3]) * sv[2];
		dv[1] = tsin(sv[3]) * sv[2];
		dv[2] = sqrt(1. - sv[2]*sv[2]);
		for (j = 3; j--; )
			samp[i*6 + 3 + j] = dv[0]*ps->sdir[0][j] +
						dv[1]*ps->sdir[1][j] +
						dv[2]*ps->sdir[2][j] ;
	}
					/* send to our process */
	writebuf(ps->wd, (char *)samp, sizeof(float)*6*ns);
	free(samp);			/* that's it! */
	return(1);
}

/* Sample the given surface */
static double
sample_surface(IDF_OBJECT *param, int wd)
{
	POLYSAMP	psamp;
	double		area;
	int		nv;
	Vert2_list	*vlist2;
					/* set up our polygon sampler */
	if ((vlist2 = init_poly(&psamp, get_vlist(param, NULL))) == NULL) {
		fprintf(stderr, "%s: bad polygon %s '%s'\n",
				progname, param->pname,
				idf_getfield(param,NAME_FLD)->val);
		return(-1.);
	}
	psamp.samp_left = nsamps;	/* assign samples & destination */
	psamp.wd = wd;
					/* hack for subsurface sampling */
	psamp.poff += 2.*SURF_EPS * !strcmp(param->pname, SUBSURF_PNAME);

	area = psamp.area_left;		/* remember starting surface area */
					/* sample each subtriangle */
	if (!polyTriangulate(vlist2, &sample_triangle))
		return(-1.);
	polyFree(vlist2);		/* clean up and return */
	return(area);
}

/* Compute User View Factors using open rcontrib process */
static int
compute_uvfs(SUBPROC *pd, ZONE *zp)
{
	IDF_OBJECT	*pptr, *pout, *pptr1;
	float		*uvfa;
	char		uvfbuf[24];
	int		n, m;
						/* create output object */
	pout = idf_newobject(our_idf, UVF_PNAME,
			"    ! computed by Radiance\n        ", our_idf->plast);
	if (pout == NULL) {
		fputs(progname, stderr);
		fputs(": cannot create new IDF object\n", stderr);
		return(0);
	}
	if (!idf_addfield(pout, zp->zname,
			"    !- Zone Name\n        ")) {
		fputs(progname, stderr);
		fputs(": cannot add zone name field\n", stderr);
		return(0);
	}
						/* allocate read buffer */
	uvfa = (float *)malloc(sizeof(float)*3*zp->ntotal);
	if (uvfa == NULL)
		return(0);
						/* UVFs from each surface */
	for (n = 0, pptr = zp->pfirst; n < zp->ntotal; n++, pptr = pptr->dnext) {
		double	vfsum = 0;
		double	adj_factor;
						/* send samples to rcontrib */
		if ((adj_factor = sample_surface(pptr, pd->w)) < 0)
			return(0);
		if (zp->area_redu == NULL)
			adj_factor = 1.;
		else				/* comp. for subsurface area */
			adj_factor /= adj_factor - zp->area_redu[n];
						/* read results */
		if (readbuf(pd->r, (char *)uvfa, sizeof(float)*3*zp->ntotal) !=
				sizeof(float)*3*zp->ntotal) {
			fputs(progname, stderr);
			fputs(": read error from rcontrib process\n", stderr);
			return(0);
		}
						/* append UVF fields */
		for (m = 0, pptr1 = zp->pfirst;
				m < zp->ntotal; m++, pptr1 = pptr1->dnext) {
			const double	uvf = uvfa[3*m + 1] * adj_factor;
			vfsum += uvf;
			if (pptr1 == pptr && uvf > .001)
				fprintf(stderr,
		"%s: warning - non-zero self-VF (%.1f%%) for surface '%s'\n",
						progname, 100.*uvf,
						idf_getfield(pptr,NAME_FLD)->val);
			sprintf(uvfbuf, "%.4f", uvf);
			if (!idf_addfield(pout,
					idf_getfield(pptr,NAME_FLD)->val, NULL) ||
				!idf_addfield(pout,
					idf_getfield(pptr1,NAME_FLD)->val, NULL) ||
				!idf_addfield(pout, uvfbuf,
						(n+m < 2*zp->ntotal-2) ?
						"\n        " : "\n\n")) {
				fputs(progname, stderr);
				fputs(": error adding UVF fields\n", stderr);
				return(0);
			}
		}
		if (vfsum < 0.95)
			fprintf(stderr,
		"%s: warning - missing %.1f%% of energy from surface '%s'\n",
					progname, 100.*(1.-vfsum),
					idf_getfield(pptr,NAME_FLD)->val);
	}
	free(uvfa);				/* clean up and return */
	return(1);
}

/* Compute zone User View Factors */
static int
compute_zones(void)
{
	ZONE	*zptr;
						/* temporary octree name */
	mktemp(strcpy(temp_octree, TEMPLATE));
						/* compute each zone */
	for (zptr = zone_list; zptr != NULL; zptr = zptr->next) {
		SUBPROC	rcproc;
						/* start rcontrib process */
		if (!start_rcontrib(&rcproc, zptr))
			return(0);
						/* compute+add view factors */
		if (!compute_uvfs(&rcproc, zptr))
			return(0);
		if (close_process(&rcproc) != 0) {
			fputs(progname, stderr);
			fputs(": bad return status from rcontrib\n", stderr);
			return(0);
		}
	}
	unlink(temp_octree);			/* remove octree file */
	return(1);
}

/* Load IDF and compute User View Factors */
int
main(int argc, char *argv[])
{
	int		incl_comments = 1;
	char		*origIDF, *revIDF;
	IDF_OBJECT	*pptr;
	int		i;

	progname = *argv++; argc--;		/* get options if any */
	while (argc > 1 && argv[0][0] == '-')
		switch (argv[0][1]) {
		case 'c':			/* elide comments */
			incl_comments = -1;		/* header only */
			argv++; argc--;
			continue;
		case 's':			/* samples */
			nsamps = 1000*atoi(*++argv);
			argv++; argc -= 2;
			continue;
		default:
			fputs(progname, stderr);
			fputs(": unknown option '", stderr);
			fputs(argv[0], stderr);
			fputs("'\n", stderr);
			goto userr;
		}
	if ((argc < 1) | (argc > 2))
		goto userr;
	origIDF = argv[0];
	revIDF = (argc == 1) ? argv[0] : argv[1];
						/* load Input Data File */
	our_idf = idf_load(origIDF);
	if (our_idf == NULL) {
		fputs(progname, stderr);
		fputs(": cannot load IDF '", stderr);
		fputs(origIDF, stderr);
		fputs("'\n", stderr);
		return(1);
	}
						/* remove existing UVFs */
	if ((pptr = idf_getobject(our_idf, UVF_PNAME)) != NULL) {
		IDF_OBJECT	*pnext;
		fputs(progname, stderr);
		fputs(": removing previous User View Factors\n", stderr);
		do {
			pnext = pptr->pnext;
			idf_delobject(our_idf, pptr);
		} while (pnext != NULL);
	}
						/* add to header */
	if (our_idf->hrem == NULL ||
		(i = strlen(our_idf->hrem)-strlen(ADD_HEADER)) < 0 ||
			strcmp(our_idf->hrem+i, ADD_HEADER))
		idf_add2hdr(our_idf, ADD_HEADER);
						/* gather zone surfaces */
	for (i = 0; surf_type[i].pname != NULL; i++)
		for (pptr = idf_getobject(our_idf, surf_type[i].pname);
				pptr != NULL; pptr = pptr->pnext) {
			IDF_FIELD	*fptr = idf_getfield(pptr,
							surf_type[i].zone_fld);
			if (fptr == NULL) {
				fputs(progname, stderr);
				fputs(": warning - missing zone field\n", stderr);
				continue;
			}
			if (add2zone(pptr, fptr->val) == NULL)
				return(1);
		}
						/* add subsurfaces */
	for (pptr = idf_getobject(our_idf, SUBSURF_PNAME);
			pptr != NULL; pptr = pptr->pnext)
		if (add_subsurf(pptr) == NULL)
			return(1);
						/* run rcontrib on each zone */
	if (!compute_zones())
		return(1);
						/* write out modified IDF */
	if (!idf_write(our_idf, revIDF, incl_comments)) {
		fputs(progname, stderr);
		fputs(": error writing IDF '", stderr);
		fputs(revIDF, stderr);
		fputs("'\n", stderr);
		return(1);
	}
	return(0);				/* finito! */
userr:
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-c][-s Ksamps] Model.idf [Revised.idf]\n", stderr);
	return(1);
}
