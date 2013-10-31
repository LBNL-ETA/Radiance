#ifndef lint
static const char RCSid[] = "$Id: bsdf2rad.c,v 2.3 2013/10/31 18:19:18 greg Exp $";
#endif
/*
 *  Plot 3-D BSDF output based on scattering interpolant representation
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bsdfrep.h"

const float	colarr[6][3] = {
		.7, 1., .7,
		1., .7, .7,
		.7, .7, 1.,
		1., .5, 1.,
		1., 1., .5,
		.5, 1., 1.
	};

char	*progname;

/* Produce a Radiance model plotting the indicated incident direction(s) */
int
main(int argc, char *argv[])
{
	int	showPeaks = 0;
	char	buf[128];
	FILE	*fp;
	RBFNODE	*rbf;
	double	bsdf, min_log;
	FVECT	dir;
	int	i, j, n;

	progname = argv[0];
	if (argc > 1 && !strcmp(argv[1], "-p")) {
		++showPeaks;
		++argv; --argc;
	}
	if (argc < 4) {
		fprintf(stderr, "Usage: %s [-p] bsdf.sir theta1 phi1 .. > output.rad\n", progname);
		return(1);
	}
						/* load input */
	if ((fp = fopen(argv[1], "rb")) == NULL) {
		fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
				progname, argv[1]);
		return(1);
	}
	if (!load_bsdf_rep(fp))
		return(1);
	fclose(fp);
	min_log = log(bsdf_min*.5);
						/* output BSDF rep. */
	for (n = 0; (n < 6) & (2*n+3 < argc); n++) {
		dir[2] = sin((M_PI/180.)*atof(argv[2*n+2]));
		dir[0] = dir[2] * cos((M_PI/180.)*atof(argv[2*n+3]));
		dir[1] = dir[2] * sin((M_PI/180.)*atof(argv[2*n+3]));
		dir[2] = input_orient * sqrt(1. - dir[2]*dir[2]);
#ifdef DEBUG
		fprintf(stderr, "Computing DSF for incident direction (%.1f,%.1f)\n",
				get_theta180(dir), get_phi360(dir));
#endif
		rbf = advect_rbf(dir, 15000);
#ifdef DEBUG
		if (rbf == NULL)
			fputs("NULL RBF\n", stderr);
		else
			fprintf(stderr, "Hemispherical reflectance: %.3f\n", rbf->vtotal);
#endif
		printf("void trans tmat\n0\n0\n7 %f %f %f .04 .04 .9 1\n",
				colarr[n][0], colarr[n][1], colarr[n][2]);
		if (showPeaks && rbf != NULL) {
			printf("void plastic pmat\n0\n0\n5 %f %f %f .04 .08\n",
				1.-colarr[n][0], 1.-colarr[n][1], 1.-colarr[n][2]);
			for (i = 0; i < rbf->nrbf; i++) {
				ovec_from_pos(dir, rbf->rbfa[i].gx, rbf->rbfa[i].gy);
				bsdf = eval_rbfrep(rbf, dir) / (output_orient*dir[2]);
				bsdf = log(bsdf) - min_log;
				printf("pmat sphere p%d\n0\n0\n4 %f %f %f %f\n",
					i+1, dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf,
						.007*bsdf);
			}
		}
		fflush(stdout);
		sprintf(buf, "gensurf tmat bsdf - - - %d %d", GRIDRES-1, GRIDRES-1);
		fp = popen(buf, "w");
		if (fp == NULL) {
			fprintf(stderr, "%s: cannot open '| %s'\n", progname, buf);
			return(1);
		}
		for (i = 0; i < GRIDRES; i++)
		    for (j = 0; j < GRIDRES; j++) {
			ovec_from_pos(dir, i, j);
			bsdf = eval_rbfrep(rbf, dir) / (output_orient*dir[2]);
			bsdf = log(bsdf) - min_log;
			fprintf(fp, "%.8e %.8e %.8e\n",
					dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
		    }
		if (rbf != NULL)
			free(rbf);
		if (pclose(fp))
			return(1);
	}
	return(0);
}
