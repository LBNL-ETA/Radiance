#ifndef lint
static const char RCSid[] = "$Id: bsdf2ttree.c,v 2.2 2012/10/20 17:01:26 greg Exp $";
#endif
/*
 * Load measured BSDF interpolant and write out as XML file with tensor tree.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "platform.h"
#include "bsdfrep.h"
				/* global argv[0] */
char			*progname;
				/* percentage to cull (<0 to turn off) */
int			pctcull = 90;
				/* sampling order */
int			samp_order = 6;

/* Interpolate and output isotropic BSDF data */
static void
interp_isotropic()
{
	const int	sqres = 1<<samp_order;
	FILE		*ofp = NULL;
	char		cmd[128];
	int		ix, ox, oy;
	FVECT		ivec, ovec;
	double		bsdf;
#if DEBUG
	fprintf(stderr, "Writing isotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %d%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	if (pctcull >= 0) {			/* begin output */
		sprintf(cmd, "rttree_reduce -h -a -fd -r 3 -t %d -g %d",
				pctcull, samp_order);
		fflush(stdout);
		ofp = popen(cmd, "w");
		if (ofp == NULL) {
			fprintf(stderr, "%s: cannot create pipe to rttree_reduce\n",
					progname);
			exit(1);
		}
		SET_FILE_BINARY(ofp);
	} else
		fputs("{\n", stdout);
						/* run through directions */
	for (ix = 0; ix < sqres/2; ix++) {
		RBFNODE	*rbf;
		SDsquare2disk(ivec, (ix+.5)/sqres, .5);
		ivec[2] = input_orient *
				sqrt(1. - ivec[0]*ivec[0] - ivec[1]*ivec[1]);
		rbf = advect_rbf(ivec);
		for (ox = 0; ox < sqres; ox++)
		    for (oy = 0; oy < sqres; oy++) {
			SDsquare2disk(ovec, (ox+.5)/sqres, (oy+.5)/sqres);
			ovec[2] = output_orient *
				sqrt(1. - ovec[0]*ovec[0] - ovec[1]*ovec[1]);
			bsdf = eval_rbfrep(rbf, ovec) / fabs(ovec[2]);
			if (pctcull >= 0)
				fwrite(&bsdf, sizeof(bsdf), 1, ofp);
			else
				printf("\t%.3e\n", bsdf);
		    }
		free(rbf);
	}
	if (pctcull >= 0) {			/* finish output */
		if (pclose(ofp)) {
			fprintf(stderr, "%s: error running '%s'\n",
					progname, cmd);
			exit(1);
		}
	} else {
		for (ix = sqres*sqres*sqres/2; ix--; )
			fputs("\t0\n", stdout);
		fputs("}\n", stdout);
	}
}

/* Interpolate and output anisotropic BSDF data */
static void
interp_anisotropic()
{
	const int	sqres = 1<<samp_order;
	FILE		*ofp = NULL;
	char		cmd[128];
	int		ix, iy, ox, oy;
	FVECT		ivec, ovec;
	double		bsdf;
#if DEBUG
	fprintf(stderr, "Writing anisotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %d%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	if (pctcull >= 0) {			/* begin output */
		sprintf(cmd, "rttree_reduce -h -a -fd -r 4 -t %d -g %d",
				pctcull, samp_order);
		fflush(stdout);
		ofp = popen(cmd, "w");
		if (ofp == NULL) {
			fprintf(stderr, "%s: cannot create pipe to rttree_reduce\n",
					progname);
			exit(1);
		}
	} else
		fputs("{\n", stdout);
						/* run through directions */
	for (ix = 0; ix < sqres; ix++)
	    for (iy = 0; iy < sqres; iy++) {
		RBFNODE	*rbf;
		SDsquare2disk(ivec, (ix+.5)/sqres, (iy+.5)/sqres);
		ivec[2] = input_orient *
				sqrt(1. - ivec[0]*ivec[0] - ivec[1]*ivec[1]);
		rbf = advect_rbf(ivec);
		for (ox = 0; ox < sqres; ox++)
		    for (oy = 0; oy < sqres; oy++) {
			SDsquare2disk(ovec, (ox+.5)/sqres, (oy+.5)/sqres);
			ovec[2] = output_orient *
				sqrt(1. - ovec[0]*ovec[0] - ovec[1]*ovec[1]);
			bsdf = eval_rbfrep(rbf, ovec) / fabs(ovec[2]);
			if (pctcull >= 0)
				fwrite(&bsdf, sizeof(bsdf), 1, ofp);
			else
				printf("\t%.3e\n", bsdf);
		    }
		free(rbf);
	    }
	if (pctcull >= 0) {			/* finish output */
		if (pclose(ofp)) {
			fprintf(stderr, "%s: error running '%s'\n",
					progname, cmd);
			exit(1);
		}
	} else
		fputs("}\n", stdout);
}

/* Read in BSDF and interpolate as tensor tree representation */
int
main(int argc, char *argv[])
{
	FILE	*fpin = stdin;
	int	i;

	progname = argv[0];			/* get options */
	while (argc > 2 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 't':
			pctcull = atoi(argv[2]);
			break;
		case 'g':
			samp_order = atoi(argv[2]);
			break;
		default:
			goto userr;
		}
		argv += 2; argc -= 2;
	}
	if (argc == 2) {			/* open input if given */
		fpin = fopen(argv[1], "r");
		if (fpin == NULL) {
			fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
					progname, argv[1]);
			return(1);
		}
	} else if (argc != 1)
		goto userr;
	SET_FILE_BINARY(fpin);			/* load BSDF interpolant */
	if (!load_bsdf_rep(fpin))
		return(1);
	draw_edges();
	/* xml_prologue();				/* start XML output */
	if (single_plane_incident)		/* resample dist. */
		interp_isotropic();
	else
		interp_anisotropic();
	/* xml_epilogue();				/* finish XML output */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-t pctcull][-g log2grid] [bsdf.sir] > bsdf.xml\n",
				progname);
	return(1);
}
