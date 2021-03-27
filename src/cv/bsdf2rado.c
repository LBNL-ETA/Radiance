#ifndef lint
static const char RCSid[] = "$Id: bsdf2rado.c,v 2.2 2021/03/27 17:50:18 greg Exp $";
#endif
/*
 *  Plot 3-D BSDF output based on scattering interpolant or XML representation
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "rtprocess.h"
#include "bsdfrep.h"

const float	colarr[6][3] = {
		.7, 1., .7,
		1., .7, .7,
		.7, .7, 1.,
		1., .5, 1.,
		1., 1., .5,
		.5, 1., 1.
	};

#if defined(_WIN32) || defined(_WIN64)
char	validf[] = "-e \"valid(s,t)=X`SYS(s,t)^2+Y`SYS(s,t)^2+Z`SYS(s,t)^2-1e-7\"";
#else
char	validf[] = "-e 'valid(s,t)=X`SYS(s,t)^2+Y`SYS(s,t)^2+Z`SYS(s,t)^2-1e-7'";
#endif

char	*progname;

/* Produce a Radiance model plotting the indicated incident direction(s) */
int
main(int argc, char *argv[])
{
	int	showPeaks = 0;
	int	doTrans = 0;
	int	inpXML = -1;
	RBFNODE	*rbf = NULL;
	FILE	*fp;
	char	buf[128];
	SDData	myBSDF;
	double	bsdf, min_log;
	FVECT	idir, odir;
	int	i, j, n;
						/* check arguments */
	progname = argv[0];
	if (argc > 1 && !strcmp(argv[1], "-p")) {
		++showPeaks;
		++argv; --argc;
	}
	if (argc > 1 && !strcmp(argv[1], "-t")) {
		++doTrans;
		++argv; --argc;
	}
	if (argc >= 4 && (n = strlen(argv[1])-4) > 0) {
		if (!strcasecmp(argv[1]+n, ".xml"))
			inpXML = 1;
		else if (!strcasecmp(argv[1]+n, ".sir"))
			inpXML = 0;
	}
	if (inpXML < 0) {
		fprintf(stderr, "Usage: %s [-p] bsdf.sir theta1 phi1 .. > output.rad\n", progname);
		fprintf(stderr, "   Or: %s [-t] bsdf.xml theta1 phi1 .. > output.rad\n", progname);
		return(1);
	}
						/* load input */
	if (inpXML) {
		SDclearBSDF(&myBSDF, argv[1]);
		if (SDreportError(SDloadFile(&myBSDF, argv[1]), stderr))
			return(1);
		bsdf_min = 1./M_PI;
		if (myBSDF.rf != NULL && myBSDF.rLambFront.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.rLambFront.cieY/M_PI;
		if (myBSDF.rb != NULL && myBSDF.rLambBack.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.rLambBack.cieY/M_PI;
		if (myBSDF.tf != NULL && myBSDF.tLambFront.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.tLambFront.cieY/M_PI;
		if (myBSDF.tb != NULL && myBSDF.tLambBack.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.tLambBack.cieY/M_PI;
		if (doTrans && (myBSDF.tf == NULL) & (myBSDF.tb == NULL)) {
			fprintf(stderr, "%s: no transmitted component in '%s'\n",
					progname, argv[1]);
			return(1);
		}
	} else {
		fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
					progname, argv[1]);
			return(1);
		}
		if (!load_bsdf_rep(fp))
			return(1);
		fclose(fp);
	}
#ifdef DEBUG
	fprintf(stderr, "Minimum BSDF set to %.4f\n", bsdf_min);
#endif
	min_log = log(bsdf_min*.5 + 1e-5);
						/* output BSDF rep. */
	for (n = 0; (n < 6) & (2*n+3 < argc); n++) {
		double	theta = (M_PI/180.)*atof(argv[2*n+2]);
		double	phi = (M_PI/180.)*atof(argv[2*n+3]);
		if (theta < -FTINY) {
			fprintf(stderr, "%s: theta values must be positive\n",
					progname);
			return(1);
		}
		if (inpXML) {
			input_orient = (theta <= M_PI/2.) ? 1 : -1;
			output_orient = doTrans ? -input_orient : input_orient;
		}
		idir[2] = sin(theta);
		idir[0] = idir[2] * cos(phi);
		idir[1] = idir[2] * sin(phi);
		idir[2] = input_orient * sqrt(1. - idir[2]*idir[2]);
#ifdef DEBUG
		fprintf(stderr, "Computing BSDF for incident direction (%.1f,%.1f)\n",
				get_theta180(idir), get_phi360(idir));
#endif
		if (!inpXML)
			rbf = advect_rbf(idir, 15000);
#ifdef DEBUG
		if (inpXML)
			fprintf(stderr, "Hemispherical %s: %.3f\n",
				(output_orient > 0 ^ input_orient > 0 ?
					"transmission" : "reflection"),
				SDdirectHemi(idir, SDsampSp|SDsampDf |
					(output_orient > 0 ^ input_orient > 0 ?
						 SDsampT : SDsampR), &myBSDF));
		else if (rbf == NULL)
			fputs("Empty RBF\n", stderr);
		else
			fprintf(stderr, "Hemispherical %s: %.3f\n",
				(output_orient > 0 ^ input_orient > 0 ?
					"transmission" : "reflection"),
				rbf->vtotal);
#endif
		printf("# Incident direction (theta,phi) = (%.2f,%.2f) deg.\n\n",
				(180./M_PI)*theta, (180./M_PI)*phi);
		printf("void trans tmat\n0\n0\n7 %f %f %f .04 .04 .9 1\n",
				colarr[n][0], colarr[n][1], colarr[n][2]);
		if (showPeaks && rbf != NULL) {
			printf("void plastic pmat\n0\n0\n5 %f %f %f .04 .08\n",
				1.-colarr[n][0], 1.-colarr[n][1], 1.-colarr[n][2]);
			for (i = 0; i < rbf->nrbf; i++) {
				ovec_from_pos(odir, rbf->rbfa[i].gx, rbf->rbfa[i].gy);
				bsdf = eval_rbfrep(rbf, odir);
				bsdf = log(bsdf + 1e-5) - min_log;
				printf("pmat sphere p%d\n0\n0\n4 %f %f %f %f\n",
					i+1, odir[0]*bsdf, odir[1]*bsdf, odir[2]*bsdf,
						.007*bsdf);
			}
		}
		fflush(stdout);
		sprintf(buf, "gensurf tmat bsdf%d - - - %d %d %s", n+1,
						GRIDRES-1, GRIDRES-1, validf);
		fp = popen(buf, "w");
		if (fp == NULL) {
			fprintf(stderr, "%s: cannot open '| %s'\n", progname, buf);
			return(1);
		}
		for (i = 0; i < GRIDRES; i++)
		    for (j = 0; j < GRIDRES; j++) {
			ovec_from_pos(odir, i, j);
			if (inpXML) {
				SDValue	sval;
				if (SDreportError(SDevalBSDF(&sval, odir,
							idir, &myBSDF), stderr))
					return(1);
				bsdf = sval.cieY;
			} else
				bsdf = eval_rbfrep(rbf, odir);
			bsdf = log(bsdf + 1e-5) - min_log;
			fprintf(fp, "%.8e %.8e %.8e\n",
					odir[0]*bsdf, odir[1]*bsdf, odir[2]*bsdf);
		    }
		if (rbf != NULL)
			free(rbf);
		if (pclose(fp))
			return(1);
	}
	return(0);
}
