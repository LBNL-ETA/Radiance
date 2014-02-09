#ifndef lint
static const char RCSid[] = "$Id: eplus_adduvf.c,v 2.2 2014/02/09 22:19:30 greg Exp $";
#endif
/*
 * Add User View Factors to EnergyPlus Input Data File
 *
 *	G.Ward for LBNL
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eplus_idf.h"
#include "triangulate.h"
#include "rtprocess.h"

#ifndef NSAMPLES
#define NSAMPLES	100000			/* number of samples to use */
#endif

char		*progname;			/* global argv[0] */

char		temp_octree[128];			/* temporary octree */

const char	UVF_PNAME[] =
			"ZoneProperty:UserViewFactor:bySurfaceName";

const char	ADD_HEADER[] =
			"!+++ User View Factors computed by Radiance +++!\n";

#define NAME_FLD	1			/* name field always first? */

typedef struct {
	const char	*pname;			/* parameter type name */
	short		zone_fld;		/* zone field index */
	short		vert_fld;		/* vertex field index */
} SURF_PTYPE;		/* surface type we're interested in */

const SURF_PTYPE	surf_type[] = {
		{"BuildingSurface:Detailed", 4, 10},
		{NULL}
	};

typedef struct s_zone {
	const char	*zname;			/* zone name */
	struct s_zone	*next;			/* next zone in list */
	int		nsurf;			/* surface count */
	IDF_PARAMETER	*pfirst;		/* first matching parameter */
	IDF_PARAMETER	*plast;			/* last matching parameter */
} ZONE;			/* a list of collected zone surfaces */

ZONE		*zone_list = NULL;	/* our list of zones */

IDF_LOADED	*our_idf = NULL;	/* loaded/modified IDF */

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

/* Determine if a parameter is a surface in the indicated zone */
static int
in_zone(IDF_PARAMETER *param, const char *zname)
{
	int		i = 0;
	IDF_FIELD	*fptr;
						/* check for surface type */
	while (strcmp(surf_type[i].pname, param->pname))
		if (surf_type[++i].pname == NULL)
			return(0);
						/* get zone field */
	fptr = idf_getfield(param, surf_type[i].zone_fld);
	if (fptr == NULL)
		return(0);
						/* check for match */
	if (strcmp(fptr->arg, zname))
		return(0);
						/* return field for #verts */
	return(surf_type[i].vert_fld);
}

/* Convert surface to Radiance with modifier based on unique name */
static int
rad_surface(const char *zname, IDF_PARAMETER *param, FILE *ofp)
{
	const char	*sname = idf_getfield(param, NAME_FLD)->arg;
	IDF_FIELD	*fptr = idf_getfield(param, in_zone(param, zname));
	int		nvert, i;

	if (fptr == NULL || (nvert = atoi(fptr->arg)) < 3) {
		fprintf(stderr, "%s: bad surface '%s'\n", progname, sname);
		return(0);
	}
	fprintf(ofp, "\nvoid glow '%s'\n0\n0\n4 1 1 1 0\n", sname);
	fprintf(ofp, "\n'%s' polygon 's_%s'\n0\n0\n%d\n", sname, sname, 3*nvert);
	while (nvert--) {
		for (i = 3; i--; ) {
			fptr = fptr->next;
			if (fptr == NULL) {
				fprintf(stderr,
				"%s: missing vertex fields in surface '%s'\n",
						progname, sname);
				return(0);
			}
			fputc('\t', ofp);
			fputs(fptr->arg, ofp);
		}
		fputc('\n', ofp);
	}
	return(!ferror(ofp));
}

/* Strat rcontrib process */
static int
start_rcontrib(SUBPROC *pd, ZONE *zp)
{
#define	BASE_AC		3
	static char	*base_av[BASE_AC] = {
				"rcontrib", "-ff", "-h"
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
	sprintf(cbuf, "%d", NSAMPLES);
	av[i++] = "-c";
	av[i++] = cbuf;				/* add modifier arguments */
	for (n = zp->nsurf, pptr = zp->pfirst; n--; pptr = pptr->dnext) {
		IDF_FIELD	*fptr = idf_getfield(pptr, NAME_FLD);
		if (fptr == NULL) {
			fputs(progname, stderr);
			fputs(": missing name for surface parameter\n", stderr);
			return(0);
		}
		av[i++] = "-m";
		av[i++] = fptr->arg;
		if (!rad_surface(zp->zname, pptr, ofp))
			return(0);
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

/* Compute User View Factors using open rcontrib process */
static int
compute_uvfs(SUBPROC *pd, ZONE *sp)
{
	
}

/* Compute zone User View Factors */
static int
compute_zones(void)
{
	ZONE	*zptr;
						/* temporary octree name */
	if (temp_filename(temp_octree, sizeof(temp_octree), TEMPLATE) == NULL) {
		fputs(progname, stderr);
		fputs(": cannot create temporary octree\n", stderr);
		return(0);
	}
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

	progname = argv[0];
	if (argc > 2 && !strcmp(argv[1], "-c")) {
		incl_comments = -1;		/* output header only */
		++argv; --argc;
	}
	if ((argc < 2) | (argc > 3)) {
		fputs("Usage: ", stderr);
		fputs(progname, stderr);
		fputs(" [-c] Model.idf [Revised.idf]\n", stderr);
		return(1);
	}
	origIDF = argv[1];
	revIDF = (argc == 2) ? argv[1] : argv[2];
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
			if (add2zone(pptr, fptr->arg) == NULL)
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
}
