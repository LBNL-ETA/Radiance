#ifndef lint
static const char RCSid[] = "$Id: bsdf2klems.c,v 2.1 2013/04/21 23:01:14 greg Exp $";
#endif
/*
 * Load measured BSDF interpolant and write out as XML file with Klems matrix.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "platform.h"
#include "calcomp.h"
#include "bsdfrep.h"
#include "bsdf_m.h"
				/* global argv[0] */
char			*progname;
				/* selected basis function name */
static const char	*kbasis = "LBNL/Klems Full";
				/* number of BSDF samples per patch */
static int		npsamps = 256;

/* Return angle basis corresponding to the given name */
ANGLE_BASIS *
get_basis(const char *bn)
{
	int	n = nabases;
	
	while (n-- > 0)
		if (!strcasecmp(bn, abase_list[n].name))
			return &abase_list[n];
	return NULL;
}

/* Output XML prologue to stdout */
static void
xml_prologue(int ac, char *av[])
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	int		i;

	if (abp == NULL) {
		fprintf(stderr, "%s: unknown angle basis '%s'\n", progname, kbasis);
		exit(1);
	}
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
	puts("\t\t<IncidentDataStructure>Columns</IncidentDataStructure>");
	puts("\t\t<AngleBasis>");
	printf("\t\t\t<AngleBasisName>%s</AngleBasisName>\n", kbasis);
	for (i = 0; abp->lat[i].nphis; i++) {
		puts("\t\t\t<AngleBasisBlock>");
		printf("\t\t\t<Theta>%g</Theta>\n", i ?
				.5*(abp->lat[i].tmin + abp->lat[i+1].tmin) :
				.0 );
		printf("\t\t\t<nPhis>%d</nPhis>", abp->lat[i].nphis);
		puts("\t\t\t<ThetaBounds>");
		printf("\t\t\t\t<LowerTheta>%g</LowerTheta>\n", abp->lat[i].tmin);
		printf("\t\t\t\t<UpperTheta>%g</UpperTheta>\n", abp->lat[i+1].tmin);
		puts("\t\t\t</ThetaBounds>");
		puts("\t\t\t</AngleBasisBlock>");
	}
	puts("\t\t</AngleBasis>");
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
	printf("\t\t\t<ColumnAngleBasis>%s</ColumnAngleBasis>\n", kbasis);
	printf("\t\t\t<RowAngleBasis>%s</RowAngleBasis>\n", kbasis);
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

/* Load and resample XML BSDF description */
static void
eval_bsdf(const char *fname)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	SDData		bsd;
	SDError		ec;
	FVECT		vin, vout;
	SDValue		sv;
	double		sum;
	int		i, j, n;

	SDclearBSDF(&bsd, fname);		/* load BSDF file */
	if ((ec = SDloadFile(&bsd, fname)) != SDEnone)
		goto err;
						/* front reflection */
	if (bsd.rf != NULL || bsd.rLambFront.cieY > .002) {
	    input_orient = 1; output_orient = 1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+(n+frandom())/npsamps, abp);
			ec = SDevalBSDF(&sv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sv.cieY;
		    }
		    printf("\t%.3e\n", sum/npsamps);
		}
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
						/* back reflection */
	if (bsd.rb != NULL || bsd.rLambBack.cieY > .002) {
	    input_orient = -1; output_orient = -1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);
			bi_getvec(vin, i+(n+frandom())/npsamps, abp);
			ec = SDevalBSDF(&sv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sv.cieY;
		    }
		    printf("\t%.3e\n", sum/npsamps);
		}
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
						/* front transmission */
	if (bsd.tf != NULL || bsd.tLamb.cieY > .002) {
	    input_orient = 1; output_orient = -1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+(n+frandom())/npsamps, abp);
			ec = SDevalBSDF(&sv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sv.cieY;
		    }
		    printf("\t%.3e\n", sum/npsamps);
		}
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
						/* back transmission */
	if (bsd.tb != NULL) {
	    input_orient = -1; output_orient = 1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			bi_getvec(vin, i+(n+frandom())/npsamps, abp);
			ec = SDevalBSDF(&sv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sv.cieY;
		    }
		    printf("\t%.3e\n", sum/npsamps);
		}
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
	SDfreeBSDF(&bsd);			/* all done */
	return;
err:
	SDreportError(ec, stderr);
	exit(1);
}

/* Interpolate and output a BSDF function */
static void
eval_function(char *funame)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	double		iovec[6];
	double		sum;
	int		i, j, n;

	data_prologue();			/* begin output */
	for (j = 0; j < abp->nangles; j++) {	/* run through directions */
	    for (i = 0; i < abp->nangles; i++) {
		sum = 0;
		for (n = npsamps; n--; ) {	/* average over patches */
		    if (output_orient > 0)
			fo_getvec(iovec+3, j+(n+frandom())/npsamps, abp);
		    else
			bo_getvec(iovec+3, j+(n+frandom())/npsamps, abp);

		    if (input_orient > 0)
			fi_getvec(iovec, j+(n+frandom())/npsamps, abp);
		    else
			bi_getvec(iovec, j+(n+frandom())/npsamps, abp);

		    sum += funvalue(funame, 6, iovec);
		}
		printf("\t%.3e\n", sum/npsamps);
	    }
	    putchar('\n');
	}
	data_epilogue();			/* finish output */
}

/* Interpolate and output a radial basis function BSDF representation */
static void
eval_rbf(void)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	double		iovec[6];
	double		sum;
	int		i, j, n;

fprintf(stder, "%s: RBF evaluation currently unimplemented\n", progname);
exit(1);
}

/* Read in BSDF and interpolate as Klems matrix representation */
int
main(int argc, char *argv[])
{
	int	dofwd = 0, dobwd = 1;
	char	*cp;
	int	i, na;

	progname = argv[0];
	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_INCHAN|E_OUTCHAN);
	scompile("PI:3.14159265358979323846", NULL, 0);
	biggerlib();
	for (i = 1; i < argc-1 && (argv[i][0] == '-') | (argv[i][0] == '+'); i++)
		switch (argv[i][1]) {		/* get options */
		case 'n':
			npsamps = atoi(argv[++i]);
			if (npsamps <= 0)
				goto userr;
			break;
		case 'e':
			scompile(argv[++i], NULL, 0);
			single_plane_incident = 0;
			break;
		case 'f':
			if (!argv[i][2]) {
				fcompile(argv[++i]);
				single_plane_incident = 0;
			} else
				dofwd = (argv[i][0] == '+');
			break;
		case 'b':
			dobwd = (argv[i][0] == '+');
			break;
		case 'h':
			kbasis = "LBNL/Klems Half";
			break;
		case 'q':
			kbasis = "LBNL/Klems Quarter";
			break;
		default:
			goto userr;
		}
	if (single_plane_incident >= 0) {	/* function-based BSDF? */
		if (i != argc-1 || fundefined(argv[i]) != 6) {
			fprintf(stderr,
	"%s: need single function with 6 arguments: bsdf(ix,iy,iz,ox,oy,oz)\n",
					progname);
			goto userr;
		}
		xml_prologue(argc, argv);	/* start XML output */
		if (dofwd) {
			input_orient = -1;
			output_orient = -1;
			eval_function(argv[i]);		/* outside reflectance */
			output_orient = 1;
			eval_function(argv[i]);		/* outside -> inside */
		}
		if (dobwd) {
			input_orient = 1;
			output_orient = 1;
			eval_function(argv[i]);		/* inside reflectance */
			output_orient = -1;
			eval_function(argv[i]);		/* inside -> outside */
		}
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	if (i == argc-1 && (cp = strstr(argv[i], ".xml")) != NULL &&
			strlen(cp) == 4) {	/* XML input? */
		xml_prologue(argc, argv);	/* start XML output */
		eval_bsdf(argv[i]);		/* load & resample BSDF */
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	if (i < argc) {				/* open input files if given */
		int	nbsdf = 0;
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
			if (!nbsdf++)		/* start XML on first dist. */
				xml_prologue(argc, argv);
			eval_rbf();
		}
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	SET_FILE_BINARY(stdin);			/* load from stdin */
	if (!load_bsdf_rep(stdin))
		return(1);
	xml_prologue(argc, argv);		/* start XML output */
	eval_rbf();				/* resample dist. */
	xml_epilogue();				/* finish XML output & exit */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-n spp][-h|-q][bsdf.sir ..] > bsdf.xml\n",
				progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q] bsdf_in.xml > bsdf_out.xml\n",
				progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q][{+|-}for[ward]][{+|-}b[ackward]][-e expr][-f file] bsdf_func > bsdf.xml\n",
				progname);
	return(1);
}
