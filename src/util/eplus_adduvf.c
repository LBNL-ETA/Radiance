#ifndef lint
static const char RCSid[] = "$Id: eplus_adduvf.c,v 2.1 2014/02/09 05:49:21 greg Exp $";
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

#ifndef NSAMPLES
#define NSAMPLES	100000			/* number of samples to use */
#endif

#define UVF_PNAME	"ZoneProperty:UserViewFactor:bySurfaceName"

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

struct s_zone {
	const char	*zname;			/* zone name */
	struct s_zone	*next;			/* next zone in list */
	IDF_PARAMETER	*pfirst;		/* first matching parameter */
	IDF_PARAMETER	*plast;			/* last matching parameter */
} *zone_list = NULL;

IDF_LOADED	*our_idf = NULL;	/* loaded/modified IDF */

/* Create a new zone and push to top of our list */
static struct s_zone *
new_zone(const char *zname, IDF_PARAMETER *param)
{
	struct s_zone	*znew = (struct s_zone *)malloc(sizeof(struct s_zone));

	if (znew == NULL)
		return(NULL);
	znew->zname = zname;			/* assumes static */
	znew->next = zone_list;
	znew->pfirst = znew->plast = param;
	return(zone_list = znew);
}

/* Add the detailed surface (polygon) to the named zone */
static struct s_zone *
add2zone(IDF_PARAMETER *param, const char *zname)
{
	struct s_zone	*zptr;

	for (zptr = zone_list; zptr != NULL; zptr = zptr->next)
		if (!strcmp(zptr->zname, zname))
			break;
	if (zptr == NULL)
		return(new_zone(zname, param));
						/* keep surfaces together */
	if (!idf_movparam(our_idf, param, zptr->plast))
		return(NULL);
	zptr->plast = param;
	return(zptr);
}

/* Compute zone User View Factors */
static int
compute_zones(void)
{
	fputs("compute_zones() unimplemented\n", stderr);
	return(0);
}

/* Load IDF and compute User View Factors */
int
main(int argc, char *argv[])
{
	char		*origIDF, *revIDF;
	IDF_PARAMETER	*pptr;
	int		i;

	if ((argc < 2) | (argc > 3)) {
		fputs("Usage: ", stderr);
		fputs(argv[0], stderr);
		fputs(" Model.idf [Revised.idf]\n", stderr);
		return(1);
	}
	origIDF = argv[1];
	revIDF = (argc == 2) ? argv[1] : argv[2];
						/* load Input Data File */
	our_idf = idf_load(origIDF);
	if (our_idf == NULL) {
		fputs(argv[0], stderr);
		fputs(": cannot load IDF '", stderr);
		fputs(origIDF, stderr);
		fputs("'\n", stderr);
		return(1);
	}
						/* remove existing UVFs */
	if ((pptr = idf_getparam(our_idf, UVF_PNAME)) != NULL) {
		IDF_PARAMETER	*pnext;
		fputs(argv[0], stderr);
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
				fputs("Warning: missing zone field!\n", stderr);
				continue;
			}
			if (add2zone(pptr, fptr->arg) == NULL)
				return(1);
		}
						/* run rcontrib on each zone */
	if (!compute_zones())
		return(1);
						/* write out modified IDF */
	if (!idf_write(our_idf, revIDF, 1)) {
		fputs(argv[0], stderr);
		fputs(": error writing IDF '", stderr);
		fputs(revIDF, stderr);
		fputs("'\n", stderr);
		return(1);
	}
	return(0);				/* finito! */
}
