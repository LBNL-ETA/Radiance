#ifndef lint
static const char RCSid[] = "$Id: bsdfquery.c,v 2.2 2013/11/22 00:00:27 greg Exp $";
#endif
/*
 *  Query values from the given BSDF (scattering interpolant or XML repres.)
 *  Input query is incident and exiting vectors directed away from surface.
 *  We normalize.  Output is a BSDF value for the vector pair.
 *  It is wise to sort the input directions to keep identical ones together
 *  when using a scattering interpolant representation.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rtmath.h"
#include "bsdfrep.h"

char	*progname;

/* Read in a vector pair */
static int
readIOdir(FVECT idir, FVECT odir, FILE *fp, int fmt)
{
	double	dvec[6];
	float	fvec[6];

	switch (fmt) {
	case 'a':
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf",
			dvec, dvec+1, dvec+2, dvec+3, dvec+4, dvec+5) != 6)
			return(0);
		VCOPY(idir, dvec);
		VCOPY(odir, dvec+3);
		break;
	case 'd':
		if (fread(dvec, sizeof(double), 6, fp) != 6)
			return(0);
		VCOPY(idir, dvec);
		VCOPY(odir, dvec+3);
		break;
	case 'f':
		if (fread(fvec, sizeof(float), 6, fp) != 6)
			return(0);
		VCOPY(idir, fvec);
		VCOPY(odir, fvec+3);
		break;
	}
	if ((normalize(idir) == 0) | (normalize(odir) == 0)) {
		fprintf(stderr, "%s: zero input vector!\n", progname);
		return(0);
	}
	return(1);
}

/* Get BSDF values for the input and output directions given on stdin */
int
main(int argc, char *argv[])
{
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
	if (argc > 2 && argv[1][0] == '-' && argv[1][1] == 'f' &&
			argv[1][2] && strchr("afd", argv[1][2]) != NULL) {
		inpfmt = outfmt = argv[1][2];
		if (argv[1][3] && strchr("afd", argv[1][3]) != NULL)
			outfmt = argv[1][3];
		++argv; --argc;
	}
	if (argc > 1 && (n = strlen(argv[1])-4) > 0) {
		if (!strcasecmp(argv[1]+n, ".xml"))
			inpXML = 1;
		else if (!strcasecmp(argv[1]+n, ".sir"))
			inpXML = 0;
	}
	if ((argc != 2) | (inpXML < 0)) {
		fprintf(stderr, "Usage: %s [-fio] bsdf.{sir|xml}\n", progname);
		return(1);
	}
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
		double	bsdf;
		float	fval;
		if (inpXML) {
			SDValue	sval;
			if (SDreportError(SDevalBSDF(&sval, odir,
						idir, &myBSDF), stderr))
				return(1);
			bsdf = sval.cieY;
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
			bsdf = eval_rbfrep(rbf, odir)/(output_orient*odir[2]);
		}
		switch (outfmt) {		/* write to stdout */
		case 'a':
			printf("%.6e\n", bsdf);
			break;
		case 'd':
			fwrite(&bsdf, sizeof(double), 1, stdout);
			break;
		case 'f':
			fval = bsdf;
			fwrite(&fval, sizeof(float), 1, stdout);
			break;
		}
	}
	return(0);
}
