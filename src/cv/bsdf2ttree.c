#ifndef lint
static const char RCSid[] = "$Id: bsdf2ttree.c,v 2.5 2012/11/09 02:16:29 greg Exp $";
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
	float		bsdf;
#if DEBUG
	fprintf(stderr, "Writing isotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %d%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	if (pctcull >= 0) {			/* begin output */
		sprintf(cmd, "rttree_reduce -h -a -ff -r 3 -t %d -g %d",
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
			bsdf = eval_rbfrep(rbf, ovec) * output_orient/ovec[2];
			if (pctcull >= 0)
				fwrite(&bsdf, sizeof(bsdf), 1, ofp);
			else
				printf("\t%.3e\n", bsdf);
		    }
		if (rbf != NULL)
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
	float		bsdf;
#if DEBUG
	fprintf(stderr, "Writing anisotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %d%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	if (pctcull >= 0) {			/* begin output */
		sprintf(cmd, "rttree_reduce -h -a -ff -r 4 -t %d -g %d",
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
			bsdf = eval_rbfrep(rbf, ovec) * output_orient/ovec[2];
			if (pctcull >= 0)
				fwrite(&bsdf, sizeof(bsdf), 1, ofp);
			else
				printf("\t%.3e\n", bsdf);
		    }
		if (rbf != NULL)
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

/* Output XML prologue to stdout */
static void
xml_prologue(int ac, char *av[])
{
	static const char	*prologue0[] = {
"<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
"<WindowElement xmlns=\"http://windows.lbl.gov\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://windows.lbl.gov/BSDF-v1.4.xsd\">",
				NULL};
	static const char	*prologue1[] = {
"<WindowElementType>System</WindowElementType>",
"<FileType>BSDF</FileType>",
"<Optical>",
"<Layer>",
"\t<Material>",
"\t\t<Name>Name</Name>",
"\t\t<Manufacturer>Manufacturer</Manufacturer>",
"\t\t<DeviceType>Other</DeviceType>",
"\t</Material>",
				NULL};
	static const char	*prologue2[] = {
"\t<WavelengthData>",
"\t\t<LayerNumber>System</LayerNumber>",
"\t\t<Wavelength unit=\"Integral\">Visible</Wavelength>",
"\t\t<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>",
"\t\t<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>",
"\t\t<WavelengthDataBlock>",
"\t\t\t<AngleBasis>LBNL/Shirley-Chiu</AngleBasis>",
"\t\t\t<ScatteringDataType>BTDF</ScatteringDataType>",
				NULL};
	static const char	*bsdf_type[4] = {
					"Reflection Back",
					"Transmission Back",
					"Transmission Front",
					"Reflection Front"
				};
	int			i;

	for (i = 0; prologue0[i] != NULL; i++)
		puts(prologue0[i]);
	fputs("<!-- File produced by:", stdout);
	while (ac-- > 0) {
		fputc(' ', stdout);
		fputs(*av++, stdout);
	}
	puts(" -->");
	for (i = 0; prologue1[i] != NULL; i++)
		puts(prologue1[i]);
	puts("\t<DataDefinition>");
	printf("\t\t<IncidentDataStructure>TensorTree%c</IncidentDataStructure>\n",
			single_plane_incident ? '3' : '4');
	puts("\t</DataDefinition>");
	for (i = 0; prologue2[i] != NULL; i++)
		puts(prologue2[i]);
	printf("\t\t\t<WavelengthDataDirection>%s</WavelengthDataDirection>\n",
		bsdf_type[(input_orient>0)<<1 | (output_orient>0)]);
	puts("\t\t\t<ScatteringData>");
}

/* Output XML epilogue to stdout */
static void
xml_epilogue(void)
{
	static const char	*epilogue[] = {
"\t\t\t</ScatteringData>",
"\t\t</WavelengthDataBlock>",
"\t</WavelengthData>",
"</Layer>",
"</Optical>",
"</WindowElement>",
				NULL};
	int			i;

	for (i = 0; epilogue[i] != NULL; i++)
		puts(epilogue[i]);
}

/* Read in BSDF and interpolate as tensor tree representation */
int
main(int argc, char *argv[])
{
	FILE	*fpin = stdin;
	int	i;

	progname = argv[0];
	for (i = 1; i < argc-1 && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {		/* get option */
		case 't':
			pctcull = atoi(argv[++i]);
			break;
		case 'g':
			samp_order = atoi(argv[++i]);
			break;
		default:
			goto userr;
		}

	if (i == argc-1) {			/* open input if given */
		fpin = fopen(argv[i], "r");
		if (fpin == NULL) {
			fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
					progname, argv[1]);
			return(1);
		}
	} else if (i < argc-1)
		goto userr;
	SET_FILE_BINARY(fpin);			/* load BSDF interpolant */
	if (!load_bsdf_rep(fpin))
		return(1);
	fclose(fpin);
	xml_prologue(argc, argv);		/* start XML output */
	if (single_plane_incident)		/* resample dist. */
		interp_isotropic();
	else
		interp_anisotropic();
	xml_epilogue();				/* finish XML output */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-t pctcull][-g log2grid] [bsdf.sir] > bsdf.xml\n",
				progname);
	return(1);
}
