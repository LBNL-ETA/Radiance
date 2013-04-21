#ifndef lint
static const char RCSid[] = "$Id: pkgBSDF.c,v 2.5 2013/04/21 22:56:14 greg Exp $";
#endif
/*
 * Take BSDF XML file and generate a referencing Radiance object
 */

#include "rtio.h"
#include "paths.h"
#include "bsdf.h"

int	do_stdout = 0;			/* send Radiance object to stdout? */
int	do_instance = 0;		/* produce instance/octree pairs? */


/* Return appropriate suffix for Z offset */
static const char *
get_suffix(double zval)
{
	return((zval > FTINY) ? "_f" : (zval < -FTINY) ? "_b" : "");
}


/* Produce a single rectangle corresponding to the given BSDF dimensions */
static void
faceBSDF(const SDData *bsp, double zval)
{
	const char	*sfx = get_suffix(zval);
							/* output material */
	printf("void BSDF m_%s%s\n", bsp->name, sfx);
	printf("6 %g %s.xml 0 1 0 .\n", zval*1.00004, bsp->name);
	printf("0\n0\n\n");
							/* output surface */
	zval *= (zval < 0) + 0.00002;
	printf("m_%s%s polygon %s%s\n", bsp->name, sfx, bsp->name, sfx);
	printf("0\n0\n12\n");
	printf("\t%.6e\t%.6e\t%.6e\n",
				-.5*bsp->dim[0], -.5*bsp->dim[1], zval);
	printf("\t%.6e\t%.6e\t%.6e\n",
				.5*bsp->dim[0], -.5*bsp->dim[1], zval);
	printf("\t%.6e\t%.6e\t%.6e\n",
				.5*bsp->dim[0], .5*bsp->dim[1], zval);
	printf("\t%.6e\t%.6e\t%.6e\n\n",
				-.5*bsp->dim[0], .5*bsp->dim[1], zval);
}


/* Convert BSDF (MGF) geometry to Radiance description */
static int
geomBSDF(const SDData *bsp)
{
	char	tmpfile[64];
	char	command[SDnameLn+64];
	int	fd;
					/* write MGF to temp file */
	fd = open(mktemp(strcpy(tmpfile,TEMPLATE)), O_WRONLY|O_CREAT|O_EXCL, 0600);
	if (fd < 0) {
		fprintf(stderr, "Cannot open temp file '%s'\n", tmpfile);
		return(0);
	}
	(void)write(fd, bsp->mgf, strlen(bsp->mgf));
	close(fd);
					/* set up command */
	if (do_instance) {
		sprintf(command, "mgf2rad %s | oconv -f - > %s.oct",
					tmpfile, bsp->name);
	} else {
		fflush(stdout);
		sprintf(command, "mgf2rad %s", tmpfile);
	}
	if (system(command)) {
		fprintf(stderr, "Error running: %s\n", command);
		return(0);
	}
	unlink(tmpfile);		/* remove temp file */
	if (do_instance) {		/* need to output instance? */
		printf("void instance %s_g\n", bsp->name);
		printf("1 %s.oct\n0\n0\n\n", bsp->name);
	}
	return(1);
}


/* Load a BSDF XML file and produce a corresponding Radiance object */
static int
cvtBSDF(char *fname)
{
	int	retOK;
	SDData	myBSDF;
	char	*pname;
					/* find and load the XML file */
	retOK = strlen(fname);
	if (retOK < 5 || strcmp(fname+retOK-4, ".xml")) {
		fprintf(stderr, "%s: input does not end in '.xml'\n", fname);
		return(0);
	}
	pname = getpath(fname, getrlibpath(), R_OK);
	if (pname == NULL) {
		fprintf(stderr, "%s: cannot find BSDF file\n", fname);
		return(0);
	}
	SDclearBSDF(&myBSDF, fname);
	if (SDreportError(SDloadFile(&myBSDF, pname), stderr))
		return(0);
	retOK = (myBSDF.dim[0] > FTINY) & (myBSDF.dim[1] > FTINY);
	if (!retOK) {
		fprintf(stderr, "%s: zero width or height\n", fname);
	} else {
		if (!do_stdout) {
			char	rname[SDnameLn+4];
			strcpy(rname, myBSDF.name);
			strcat(rname, ".rad");
			retOK = (freopen(rname, "w", stdout) != NULL);
		}
		if (retOK) {
			if (myBSDF.matn[0] && myBSDF.makr[0])
				printf("# Material '%s' by '%s'\n\n",
						myBSDF.matn, myBSDF.makr);
			if (myBSDF.mgf == NULL) {
				faceBSDF(&myBSDF, .0);
			} else {
				faceBSDF(&myBSDF, myBSDF.dim[2]);
				if (myBSDF.rb != NULL)
					faceBSDF(&myBSDF, -myBSDF.dim[2]);
				retOK = geomBSDF(&myBSDF);
			}
		}
	}
	SDfreeBSDF(&myBSDF);		/* clean up and return */
	return(retOK);
}


/* Generate BSDF-referencing Radiance scenes, one per XML input */
int
main(int argc, char *argv[])
{
	int	status = 0;
	int	i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'i':
			do_instance = 1;
			break;
		case 's':
			do_stdout = 1;
			break;
		default:
			goto userr;
		}
	if (i >= argc) {
		fprintf(stderr, "%s: missing XML input\n", argv[0]);
		goto userr;
	}
	if (i < argc-1 && do_stdout) {
		fprintf(stderr, "%s: cannot send multiple BSDFs to stdout\n",
				argv[0]);
		return(1);
	}
	for ( ; i < argc; i++)
		if (!cvtBSDF(argv[i])) {
			fprintf(stderr, "%s: conversion of '%s' failed\n",
					argv[0], argv[i]);
			status = 1;
		}
	return(status);
userr:
	fprintf(stderr, "Usage: %s [-i][-s] input.xml ..\n", argv[0]);
	return(1);
}
