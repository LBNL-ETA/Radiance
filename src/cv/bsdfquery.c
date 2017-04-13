#ifndef lint
static const char RCSid[] = "$Id: bsdfquery.c,v 2.9 2017/04/13 00:14:36 greg Exp $";
#endif
/*
 *  Query values from the given BSDF (scattering interpolant or XML repres.)
 *  Input query is incident and exiting vectors directed away from surface.
 *  We normalize.  Output is a BSDF value for the vector pair.
 *  A zero length in or out vector is ignored, causing output to be flushed.
 *  It is wise to sort the input directions to keep identical ones together
 *  when using a scattering interpolant representation.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rtmath.h"
#include "rtio.h"
#include "bsdfrep.h"

char	*progname;

/* Read in a vector pair */
static int
readIOdir(FVECT idir, FVECT odir, FILE *fp, int fmt)
{
	double	dvec[6];
	float	fvec[6];
tryagain:
	switch (fmt) {
	case 'a':
		if (fscanf(fp, FVFORMAT, &idir[0], &idir[1], &idir[2]) != 3 ||
				fscanf(fp, FVFORMAT, &odir[0], &odir[1], &odir[2]) != 3)
			return(0);
		break;
	case 'd':
		if (getbinary(dvec, sizeof(double), 6, fp) != 6)
			return(0);
		VCOPY(idir, dvec);
		VCOPY(odir, dvec+3);
		break;
	case 'f':
		if (getbinary(fvec, sizeof(float), 6, fp) != 6)
			return(0);
		VCOPY(idir, fvec);
		VCOPY(odir, fvec+3);
		break;
	}
	if ((normalize(idir) == 0) | (normalize(odir) == 0)) {
		fflush(stdout);		/* desired side-effect? */
		goto tryagain;
	}
	return(1);
}

/* Get BSDF values for the input and output directions given on stdin */
int
main(int argc, char *argv[])
{
	int	repXYZ = 0;
	int	inpXML = -1;
	int	inpfmt = 'a';
	int	outfmt = 'a';
	RBFNODE	*rbf = NULL;
	int32	prevInpDir = 0;
	SDData	myBSDF;
	FVECT	idir, odir;
	int	n;
						/* check arguments */
	progname = argv[0];
	while (argc > 2 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 'c':			/* color output */
			repXYZ = 1;
			break;
		case 'f':			/* i/o format */
			if (!argv[1][2] || strchr("afd", argv[1][2]) == NULL)
				goto userr;
			inpfmt = outfmt = argv[1][2];
			if (argv[1][3] && strchr("afd", argv[1][3]) != NULL)
				outfmt = argv[1][3];
			break;
		default:
			goto userr;
		}
		++argv; --argc;
	}
	if (argc > 1 && (n = strlen(argv[1])-4) > 0) {
		if (!strcasecmp(argv[1]+n, ".xml"))
			inpXML = 1;
		else if (!strcasecmp(argv[1]+n, ".sir"))
			inpXML = 0;
	}
	if ((argc != 2) | (inpXML < 0))
		goto userr;
						/* load BSDF representation */
	if (inpXML) {
		SDclearBSDF(&myBSDF, argv[1]);
		if (SDreportError(SDloadFile(&myBSDF, argv[1]), stderr))
			return(1);
	} else {
		FILE	*fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
					progname, argv[1]);
			return(1);
		}
		if (!load_bsdf_rep(fp))
			return(1);
		fclose(fp);
	}
						/* query BSDF values */
	while (readIOdir(idir, odir, stdin, inpfmt)) {
		SDValue	sval;
		if (inpXML) {
			if (SDreportError(SDevalBSDF(&sval, odir,
						idir, &myBSDF), stderr))
				return(1);
		} else {
			int32	inpDir = encodedir(idir);
			if (inpDir != prevInpDir) {
				if (idir[2] > 0 ^ input_orient > 0) {
					fprintf(stderr, "%s: input hemisphere error\n",
							progname);
					return(1);
				}
				if (rbf != NULL) free(rbf);
				rbf = advect_rbf(idir, 15000);
				prevInpDir = inpDir;
			}
			if (odir[2] > 0 ^ output_orient > 0) {
				fprintf(stderr, "%s: output hemisphere error\n",
						progname);
				return(1);
			}
			if (SDreportError(eval_rbfcol(&sval, rbf, odir), stderr))
				return(1);
		}

		switch (outfmt) {		/* write to stdout */
		case 'a':
			if (repXYZ) {
				double	cieX = sval.spec.cx/sval.spec.cy*sval.cieY;
				double	cieZ = (1. - sval.spec.cx - sval.spec.cy) /
						sval.spec.cy * sval.cieY;
				printf("%.6e %.6e %.6e\n", cieX, sval.cieY, cieZ);
			} else
				printf("%.6e\n", sval.cieY);
			break;
		case 'd':
			if (repXYZ) {
				double	cieXYZ[3];
				cieXYZ[0] = sval.spec.cx/sval.spec.cy*sval.cieY;
				cieXYZ[1] = sval.cieY;
				cieXYZ[2] = (1. - sval.spec.cx - sval.spec.cy) /
						sval.spec.cy * sval.cieY;
				putbinary(cieXYZ, sizeof(double), 3, stdout);
			} else
				putbinary(&sval.cieY, sizeof(double), 1, stdout);
			break;
		case 'f':
			if (repXYZ) {
				float	cieXYZ[3];
				cieXYZ[0] = sval.spec.cx/sval.spec.cy*sval.cieY;
				cieXYZ[1] = sval.cieY;
				cieXYZ[2] = (1. - sval.spec.cx - sval.spec.cy) /
						sval.spec.cy * sval.cieY;
				putbinary(cieXYZ, sizeof(float), 3, stdout);
			} else {
				float	cieY = sval.cieY;
				putbinary(&cieY, sizeof(float), 1, stdout);
			}
			break;
		}
	}
	/* if (rbf != NULL) free(rbf); */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-c][-fio] bsdf.{sir|xml}\n", progname);
	return(1);
}
