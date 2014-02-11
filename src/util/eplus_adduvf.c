#ifndef lint
static const char RCSid[] = "$Id: eplus_adduvf.c,v 2.7 2014/02/11 23:44:11 greg Exp $";
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

char		*progname;			/* global argv[0] */

int		nsamps = NSAMPLES;		/* number of samples to use */

char		temp_octree[128];		/* temporary octree */

const char	UVF_PNAME[] =
			"ZoneProperty:UserViewFactor:bySurfaceName";

const char	ADD_HEADER[] =
			"\n!+++ User View Factors computed by Radiance +++!\n\n";

#define NAME_FLD	1			/* name field always first? */

typedef struct {
	const char	*pname;			/* parameter type name */
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
	IDF_PARAMETER	*pfirst;		/* first matching parameter */
	IDF_PARAMETER	*plast;			/* last matching parameter */
} ZONE;			/* a list of collected zone surfaces */

ZONE		*zone_list = NULL;	/* our list of zones */

IDF_LOADED	*our_idf = NULL;	/* loaded/modified IDF */

typedef struct {
	FVECT		sdir[3];	/* UVW unit sampling vectors */
	double		poff;		/* W-offset for plane of polygon */
	double		area_left;	/* area left to sample */
	int		samp_left;	/* remaining samples */
	int		wd;		/* output file descriptor */
} POLYSAMP;		/* structure for polygon sampling */

/* Create a new zone and push to top of our list */
static ZONE *
new_zone(const char *zname, IDF_PARAMETER *param)
{
	ZONE	*znew = (ZONE *)malloc(sizeof(ZONE));

	if (znew == NULL)
		return(NULL);
	znew->zname = zname;			/* assumes static */
	znew->next = zone_list;
	znew->pfirst = znew->plast = param;
	znew->nsurf = 1;
	return(zone_list = znew);
}

/* Add the detailed surface (polygon) to the named zone */
static ZONE *
add2zone(IDF_PARAMETER *param, const char *zname)
{
	ZONE	*zptr;

	for (zptr = zone_list; zptr != NULL; zptr = zptr->next)
		if (!strcmp(zptr->zname, zname))
			break;
	if (zptr == NULL)
		return(new_zone(zname, param));
						/* keep surfaces together */
	if (!idf_movparam(our_idf, param, zptr->plast))
		return(NULL);
	zptr->plast = param;
	zptr->nsurf++;
	return(zptr);
}

/* Return field for vertices in the given parameter */
static IDF_FIELD *
get_vlist(IDF_PARAMETER *param, const char *zname)
{
	int		i = 0;
	IDF_FIELD	*fptr;
						/* check for surface type */
	while (strcmp(surf_type[i].pname, param->pname))
		if (surf_type[++i].pname == NULL)
			return(NULL);

	if (zname != NULL) {			/* matches specified zone? */
		fptr = idf_getfield(param, surf_type[i].zone_fld);
		if (fptr == NULL || strcmp(fptr->val, zname))
			return(NULL);
	}
						/* return field for #verts */
	return(idf_getfield(param, surf_type[i].vert_fld));
}

/* Convert surface to Radiance with modifier based on unique name */
static int
rad_surface(IDF_PARAMETER *param, FILE *ofp)
{
	const char	*sname = idf_getfield(param,NAME_FLD)->val;
	IDF_FIELD	*fptr = get_vlist(param, NULL);
	int		nvert, i;

	if (fptr == NULL || (nvert = atoi(fptr->val)) < 3) {
		fprintf(stderr, "%s: bad surface '%s'\n", progname, sname);
		return(0);
	}
	fprintf(ofp, "\nvoid glow '%s'\n0\n0\n4 1 1 1 0\n", sname);
	fprintf(ofp, "\n'%s' polygon 's_%s'\n0\n0\n%d\n", sname, sname, 3*nvert);
	while (nvert--) {
		for (i = 3; i--; ) {
			fptr = fptr->next;
			if (fptr == NULL || !isflt(fptr->val)) {
				fprintf(stderr,
					"%s: missing/bad vertex for %s '%s'\n",
						progname, param->pname, sname);
				return(0);
			}
			fputc('\t', ofp);
			fputs(fptr->val, ofp);
		}
		fputc('\n', ofp);
	}
	return(!ferror(ofp));
}

/* Start rcontrib process */
static int
start_rcontrib(SUBPROC *pd, ZONE *zp)
{
#define	BASE_AC		5
	static char	*base_av[BASE_AC] = {
				"rcontrib", "-ff", "-h", "-x", "1"
			};
	char		cbuf[300];
	char		**av;
	FILE		*ofp;
	IDF_PARAMETER	*pptr;
	int		i, n;
						/* start oconv command */
	sprintf(cbuf, "oconv - > '%s'", temp_octree);
	if ((ofp = popen(cbuf, "w")) == NULL) {
		fputs(progname, stderr);
		fputs(": cannot open oconv process\n", stderr);
		return(0);
	}
						/* allocate argument list */
	av = (char **)malloc(sizeof(char *)*(BASE_AC+4+2*zp->nsurf));
	if (av == NULL)
		return(0);
	for (i = 0; i < BASE_AC; i++)
		av[i] = base_av[i];
	sprintf(cbuf, "%d", nsamps);
	av[i++] = "-c";
	av[i++] = cbuf;				/* add modifier arguments */
	for (n = zp->nsurf, pptr = zp->pfirst; n--; pptr = pptr->dnext) {
		IDF_FIELD	*fptr = idf_getfield(pptr,NAME_FLD);
		if (fptr == NULL || !fptr->val[0]) {
			fputs(progname, stderr);
			fputs(": missing name for surface parameter\n", stderr);
			return(0);
		}
		if (!rad_surface(pptr, ofp))	/* add surface to octree */
			return(0);
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
	if (!open_process(pd, av)) {		/* start process */
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
init_poly(POLYSAMP *ps, IDF_FIELD *f0, int nv)
{
	IDF_FIELD	*fptr = f0;
	int		i, j;
	FVECT		*vl3, e1, e2, vc;
	Vert2_list	*vl2 = polyAlloc(nv);

	if (vl2 == NULL)
		return(NULL);
	vl2->p = ps;
					/* get 3-D vertices */
	vl3 = (FVECT *)malloc(sizeof(FVECT)*nv);
	if (vl3 == NULL)
		return(NULL);
	for (i = nv; i--; )		/* reverse vertex ordering */
		for (j = 0; j < 3; j++) {
			if (fptr == NULL) {
				fputs(progname, stderr);
				fputs(": missing vertex in init_poly()\n", stderr);
				return(NULL);
			}
			vl3[i][j] = atof(fptr->val);
			fptr = fptr->next;
		}
					/* compute area and normal */
	ps->sdir[2][0] = ps->sdir[2][1] = ps->sdir[2][2] = 0;
	VSUB(e1, vl3[1], vl3[0]);
	for (i = 2; i < nv; i++) {
		VSUB(e2, vl3[i], vl3[0]); 
		fcross(vc, e1, e2);
		ps->sdir[2][0] += vc[0];
		ps->sdir[2][1] += vc[1];
		ps->sdir[2][2] += vc[2];
		VCOPY(e1, e2);
	}
	ps->area_left = .5 * normalize(ps->sdir[2]);
	if (ps->area_left == .0) {
		fputs(progname, stderr);
		fputs(": degenerate polygon in init_poly()\n", stderr);
		return(0);
	}
					/* create X & Y axes */
	VCOPY(ps->sdir[0], e1);
	normalize(ps->sdir[0]);
	fcross(ps->sdir[1], ps->sdir[2], ps->sdir[0]);
					/* compute plane offset */
	ps->poff = DOT(vl3[0], ps->sdir[2]);
					/* assign 2-D vertices */
	for (i = 0; i < nv; i++) {
		vl2->v[i].mX = DOT(vl3[i], ps->sdir[0]);
		vl2->v[i].mY = DOT(vl3[i], ps->sdir[1]);
	}
	free(vl3);			/* it's ready! */
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
				(ps->poff+.001)*ps->sdir[2][i];
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
static int
sample_surface(IDF_PARAMETER *param, int wd)
{
	IDF_FIELD	*fptr = get_vlist(param, NULL);
	POLYSAMP	psamp;
	int		nv;
	Vert2_list	*vlist2;
					/* set up our polygon sampler */
	if (fptr == NULL || (nv = atoi(fptr->val)) < 3 ||
			(vlist2 = init_poly(&psamp, fptr->next, nv)) == NULL) {
		fprintf(stderr, "%s: bad polygon %s '%s'\n",
				progname, param->pname,
				idf_getfield(param,NAME_FLD)->val);
		return(0);
	}
	psamp.samp_left = nsamps;	/* assign samples & destination */
	psamp.wd = wd;
					/* sample each subtriangle */
	if (!polyTriangulate(vlist2, &sample_triangle))
		return(0);
	polyFree(vlist2);		/* clean up and return */
	return(1);
}

/* Compute User View Factors using open rcontrib process */
static int
compute_uvfs(SUBPROC *pd, ZONE *zp)
{
	IDF_PARAMETER	*pptr, *pout, *pptr1;
	float		*uvfa;
	char		uvfbuf[24];
	int		n, m;
						/* create output parameter */
	pout = idf_newparam(our_idf, UVF_PNAME,
			"    ! computed by Radiance\n        ", zp->plast);
	if (pout == NULL) {
		fputs(progname, stderr);
		fputs(": cannot create new IDF parameter\n", stderr);
		return(0);
	}
	if (!idf_addfield(pout, zp->zname,
			"    !- Zone Name\n        ")) {
		fputs(progname, stderr);
		fputs(": cannot add zone name field\n", stderr);
		return(0);
	}
						/* allocate read buffer */
	uvfa = (float *)malloc(sizeof(float)*3*zp->nsurf);
	if (uvfa == NULL)
		return(0);
						/* UVFs from each surface */
	for (n = zp->nsurf, pptr = zp->pfirst; n--; pptr = pptr->dnext) {
		double	vfsum = 0;
						/* send samples to rcontrib */
		if (!sample_surface(pptr, pd->w))
			return(0);
						/* read results */
		if (readbuf(pd->r, (char *)uvfa, sizeof(float)*3*zp->nsurf) !=
				sizeof(float)*3*zp->nsurf) {
			fputs(progname, stderr);
			fputs(": read error from rcontrib process\n", stderr);
			return(0);
		}
						/* append UVF fields */
		for (m = 0, pptr1 = zp->pfirst;
				m < zp->nsurf; m++, pptr1 = pptr1->dnext) {
			vfsum += uvfa[3*m + 1];
			if (pptr1 == pptr) {
				if (uvfa[3*m + 1] > .001)
					fprintf(stderr,
		"%s: warning - non-zero self-VF (%.1f%%) for surface '%s'\n",
						progname, 100.*uvfa[3*m + 1],
						idf_getfield(pptr,NAME_FLD)->val);
				continue;	/* don't record self-factor */
			}
			sprintf(uvfbuf, "%.4f", uvfa[3*m + 1]);
			if (!idf_addfield(pout,
					idf_getfield(pptr,NAME_FLD)->val, NULL) ||
				!idf_addfield(pout,
					idf_getfield(pptr1,NAME_FLD)->val, NULL) ||
				!idf_addfield(pout, uvfbuf,
						(n || m < zp->nsurf-2) ?
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
	IDF_PARAMETER	*pptr;
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
	if ((pptr = idf_getparam(our_idf, UVF_PNAME)) != NULL) {
		IDF_PARAMETER	*pnext;
		fputs(progname, stderr);
		fputs(": removing previous User View Factors\n", stderr);
		do {
			pnext = pptr->pnext;
			idf_delparam(our_idf, pptr);
		} while (pnext != NULL);
	}
						/* add to header */
	if (our_idf->hrem == NULL ||
		(i = strlen(our_idf->hrem)-strlen(ADD_HEADER)) < 0 ||
			strcmp(our_idf->hrem+i, ADD_HEADER))
		idf_add2hdr(our_idf, ADD_HEADER);
						/* gather zone surfaces */
	for (i = 0; surf_type[i].pname != NULL; i++)
		for (pptr = idf_getparam(our_idf, surf_type[i].pname);
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
