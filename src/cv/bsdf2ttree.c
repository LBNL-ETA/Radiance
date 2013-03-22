#ifndef lint
static const char RCSid[] = "$Id: bsdf2ttree.c,v 2.11 2013/03/22 02:15:06 greg Exp $";
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
#include "calcomp.h"
#include "bsdfrep.h"
				/* global argv[0] */
char			*progname;
				/* percentage to cull (<0 to turn off) */
double			pctcull = 90.;
				/* sampling order */
int			samp_order = 6;

/* Output XML prologue to stdout */
static void
xml_prologue(int ac, char *av[])
{
	puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
	puts("<WindowElement xmlns=\"http://windows.lbl.gov\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://windows.lbl.gov/BSDF-v1.4.xsd\">");
	fputs("<!-- File produced by:", stdout);
	while (ac-- > 0) {
		fputc(' ', stdout);
		fputs(*av++, stdout);
	}
	puts(" -->");
	puts("<WindowElementType>System</WindowElementType>");
	puts("<FileType>BSDF</FileType>");
	puts("<Optical>");
	puts("<Layer>");
	puts("\t<Material>");
	puts("\t\t<Name>Name</Name>");
	puts("\t\t<Manufacturer>Manufacturer</Manufacturer>");
	puts("\t\t<DeviceType>Other</DeviceType>");
	puts("\t</Material>");
	puts("\t<DataDefinition>");
	printf("\t\t<IncidentDataStructure>TensorTree%c</IncidentDataStructure>\n",
			single_plane_incident ? '3' : '4');
	puts("\t</DataDefinition>");
}

/* Output XML data prologue to stdout */
static void
data_prologue()
{
	static const char	*bsdf_type[4] = {
					"Reflection Front",
					"Transmission Front",
					"Transmission Back",
					"Reflection Back"
				};

	puts("\t<WavelengthData>");
	puts("\t\t<LayerNumber>System</LayerNumber>");
	puts("\t\t<Wavelength unit=\"Integral\">Visible</Wavelength>");
	puts("\t\t<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
	puts("\t\t<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>");
	puts("\t\t<WavelengthDataBlock>");
	printf("\t\t\t<WavelengthDataDirection>%s</WavelengthDataDirection>\n",
			bsdf_type[(input_orient>0)<<1 | (output_orient>0)]);
	puts("\t\t\t<AngleBasis>LBNL/Shirley-Chiu</AngleBasis>");
	puts("\t\t\t<ScatteringDataType>BTDF</ScatteringDataType>");
	puts("\t\t\t<ScatteringData>");
}

/* Output XML data epilogue to stdout */
static void
data_epilogue(void)
{
	puts("\t\t\t</ScatteringData>");
	puts("\t\t</WavelengthDataBlock>");
	puts("\t</WavelengthData>");
}

/* Output XML epilogue to stdout */
static void
xml_epilogue(void)
{
	puts("</Layer>");
	puts("</Optical>");
	puts("</WindowElement>");
}

/* Interpolate and output isotropic BSDF data */
static void
eval_isotropic(char *funame)
{
	const int	sqres = 1<<samp_order;
	FILE		*ofp = NULL;
	char		cmd[128];
	int		ix, ox, oy;
	double		iovec[6];
	float		bsdf;
#if DEBUG
	fprintf(stderr, "Writing isotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %.1f%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	data_prologue();			/* begin output */
	if (pctcull >= 0) {
		sprintf(cmd, "rttree_reduce -h -a -ff -r 3 -t %f -g %d",
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
		RBFNODE	*rbf = NULL;
		iovec[0] = (ix+.5)/sqres - 1.;
		iovec[1] = .0;
		iovec[2] = input_orient * sqrt(1. - iovec[0]*iovec[0]);
		if (funame == NULL)
			rbf = advect_rbf(iovec);
		for (ox = 0; ox < sqres; ox++)
		    for (oy = 0; oy < sqres; oy++) {
			SDsquare2disk(iovec+3, (ox+.5)/sqres, (oy+.5)/sqres);
			iovec[5] = output_orient *
				sqrt(1. - iovec[3]*iovec[3] - iovec[4]*iovec[4]);
			if (funame == NULL)
				bsdf = eval_rbfrep(rbf, iovec+3) *
						output_orient/iovec[5];
			else
				bsdf = funvalue(funame, 6, iovec);
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
	data_epilogue();
}

/* Interpolate and output anisotropic BSDF data */
static void
eval_anisotropic(char *funame)
{
	const int	sqres = 1<<samp_order;
	FILE		*ofp = NULL;
	char		cmd[128];
	int		ix, iy, ox, oy;
	double		iovec[6];
	float		bsdf;
#if DEBUG
	fprintf(stderr, "Writing anisotropic order %d ", samp_order);
	if (pctcull >= 0) fprintf(stderr, "data with %.1f%% culling\n", pctcull);
	else fputs("raw data\n", stderr);
#endif
	data_prologue();			/* begin output */
	if (pctcull >= 0) {
		sprintf(cmd, "rttree_reduce -h -a -ff -r 4 -t %f -g %d",
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
		RBFNODE	*rbf = NULL;		/* Klems reversal */
		SDsquare2disk(iovec, (ix+.5)/sqres, (iy+.5)/sqres);
		iovec[0] = -iovec[0]; iovec[1] = -iovec[1];
		iovec[2] = input_orient *
				sqrt(1. - iovec[0]*iovec[0] - iovec[1]*iovec[1]);
		if (funame == NULL)
			rbf = advect_rbf(iovec);
		for (ox = 0; ox < sqres; ox++)
		    for (oy = 0; oy < sqres; oy++) {
			SDsquare2disk(iovec+3, (ox+.5)/sqres, (oy+.5)/sqres);
			iovec[5] = output_orient *
				sqrt(1. - iovec[3]*iovec[3] - iovec[4]*iovec[4]);
			if (funame == NULL)
				bsdf = eval_rbfrep(rbf, iovec+3) *
						output_orient/iovec[5];
			else
				bsdf = funvalue(funame, 6, iovec);
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
	data_epilogue();
}

/* Read in BSDF and interpolate as tensor tree representation */
int
main(int argc, char *argv[])
{
	int	dofwd = 0, dobwd = 1;
	int	i, na;

	progname = argv[0];
	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_INCHAN|E_OUTCHAN);
	scompile("PI:3.14159265358979323846", NULL, 0);
	biggerlib();
	for (i = 1; i < argc-1 && (argv[i][0] == '-') | (argv[i][0] == '+'); i++)
		switch (argv[i][1]) {		/* get options */
		case 'e':
			scompile(argv[++i], NULL, 0);
			break;
		case 'f':
			if (!argv[i][2])
				fcompile(argv[++i]);
			else
				dofwd = (argv[i][0] == '+');
			break;
		case 'b':
			dobwd = (argv[i][0] == '+');
			break;
		case 't':
			switch (argv[i][2]) {
			case '3':
				single_plane_incident = 0;
				break;
			case '4':
				single_plane_incident = 1;
				break;
			case '\0':
				pctcull = atof(argv[++i]);
				break;
			default:
				goto userr;
			}
			break;
		case 'g':
			samp_order = atoi(argv[++i]);
			break;
		default:
			goto userr;
		}
	if (single_plane_incident >= 0) {	/* function-based BSDF? */
		void	(*evf)(char *s) = single_plane_incident ?
				&eval_isotropic : &eval_anisotropic;
		if (i != argc-1 || fundefined(argv[i]) != 6)
			goto userr;
		xml_prologue(argc, argv);	/* start XML output */
		if (dofwd) {
			input_orient = -1;
			output_orient = -1;
			(*evf)(argv[i]);	/* outside reflectance */
			output_orient = 1;
			(*evf)(argv[i]);	/* outside -> inside */
		}
		if (dobwd) {
			input_orient = 1;
			output_orient = 1;
			(*evf)(argv[i]);	/* inside reflectance */
			output_orient = -1;
			(*evf)(argv[i]);	/* inside -> outside */
		}
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	if (i < argc) {				/* open input files if given */
		xml_prologue(argc, argv);	/* start XML output */
		for ( ; i < argc; i++) {	/* interpolate each component */
			FILE	*fpin = fopen(argv[i], "rb");
			if (fpin == NULL) {
				fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
						progname, argv[i]);
				return(1);
			}
			if (!load_bsdf_rep(fpin))
				return(1);
			fclose(fpin);
			if (single_plane_incident)
				eval_isotropic(NULL);
			else
				eval_anisotropic(NULL);
		}
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	SET_FILE_BINARY(stdin);			/* load from stdin */
	if (!load_bsdf_rep(stdin))
		return(1);
	xml_prologue(argc, argv);		/* start XML output */
	if (single_plane_incident)		/* resample dist. */
		eval_isotropic(NULL);
	else
		eval_anisotropic(NULL);
	xml_epilogue();				/* finish XML output & exit */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-g Nlog2][-t pctcull] [bsdf.sir ..] > bsdf.xml\n",
				progname);
	fprintf(stderr,
	"   or: %s -t{3|4} [-g Nlog2][-t pctcull][{+|-}for[ward]][{+|-}b[ackward]][-e expr][-f file] bsdf_func > bsdf.xml\n",
				progname);
	return(1);
}
