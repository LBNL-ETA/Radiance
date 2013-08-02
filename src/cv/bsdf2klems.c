#ifndef lint
static const char RCSid[] = "$Id: bsdf2klems.c,v 2.8 2013/08/02 20:56:19 greg Exp $";
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
				/* assumed maximum # Klems patches */
#define MAXPATCHES	145
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

/* Output XML header to stdout */
static void
xml_header(int ac, char *av[])
{
	puts("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
	puts("<WindowElement xmlns=\"http://windows.lbl.gov\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://windows.lbl.gov/BSDF-v1.4.xsd\">");
	fputs("<!-- File produced by:", stdout);
	while (ac-- > 0) {
		fputc(' ', stdout);
		fputs(*av++, stdout);
	}
	puts(" -->");
}

/* Output XML prologue to stdout */
static void
xml_prologue(const SDData *sd)
{
	const char	*matn = (sd && sd->matn[0]) ? sd->matn : "Name";
	const char	*makr = (sd && sd->makr[0]) ? sd->makr : "Manufacturer";
	ANGLE_BASIS	*abp = get_basis(kbasis);
	int		i;

	if (abp == NULL) {
		fprintf(stderr, "%s: unknown angle basis '%s'\n", progname, kbasis);
		exit(1);
	}
	puts("<WindowElementType>System</WindowElementType>");
	puts("<FileType>BSDF</FileType>");
	puts("<Optical>");
	puts("<Layer>");
	puts("\t<Material>");
	printf("\t\t<Name>%s</Name>\n", matn);
	printf("\t\t<Manufacturer>%s</Manufacturer>\n", makr);
	if (sd && sd->dim[2] > .001) {
		printf("\t\t<Thickness unit=\"meter\">%.3f</Thickness>\n", sd->dim[2]);
		printf("\t\t<Width unit=\"meter\">%.3f</Width>\n", sd->dim[0]);
		printf("\t\t<Height unit=\"meter\">%.3f</Height>\n", sd->dim[1]);
	}
	puts("\t\t<DeviceType>Other</DeviceType>");
	puts("\t</Material>");
	if (sd && sd->mgf != NULL) {
		puts("\t<Geometry format=\"MGF\">");
		puts("\t\t<MGFblock unit=\"meter\">");
		fputs(sd->mgf, stdout);
		puts("</MGFblock>");
		puts("\t</Geometry>");
	}
	puts("\t<DataDefinition>");
	puts("\t\t<IncidentDataStructure>Columns</IncidentDataStructure>");
	puts("\t\t<AngleBasis>");
	printf("\t\t\t<AngleBasisName>%s</AngleBasisName>\n", kbasis);
	for (i = 0; abp->lat[i].nphis; i++) {
		puts("\t\t\t<AngleBasisBlock>");
		printf("\t\t\t<Theta>%g</Theta>\n", i ?
				.5*(abp->lat[i].tmin + abp->lat[i+1].tmin) :
				.0 );
		printf("\t\t\t<nPhis>%d</nPhis>\n", abp->lat[i].nphis);
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

/* Load and resample XML BSDF description using Klems basis */
static void
eval_bsdf(const char *fname)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	float		*trans_mtx = NULL;
	SDData		bsd;
	SDError		ec;
	FVECT		vin, vout;
	SDValue		sv;
	double		sum;
	int		i, j, n;

	SDclearBSDF(&bsd, fname);		/* load BSDF file */
	if ((ec = SDloadFile(&bsd, fname)) != SDEnone)
		goto err;
	xml_prologue(&bsd);			/* pass geometry */
						/* front reflection */
	if (bsd.rf != NULL || bsd.rLambFront.cieY > .002) {
	    input_orient = 1; output_orient = 1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+urand(n), abp);
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
			bi_getvec(vin, i+urand(n), abp);
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
	    if (bsd.tb == NULL)
		trans_mtx = (float *)malloc(sizeof(float) *
						abp->nangles*abp->nangles);
	    input_orient = 1; output_orient = -1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    for (n = npsamps; n-- > 0; ) {
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+urand(n), abp);
			ec = SDevalBSDF(&sv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sv.cieY;
		    }
		    printf("\t%.3e\n", sum/npsamps);
		    if (trans_mtx != NULL)
			trans_mtx[j*abp->nangles + i] = sum/npsamps;
		}
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
						/* back transmission */
	if (bsd.tb != NULL || trans_mtx != NULL) {
	    if (bsd.tf == NULL)
		trans_mtx = (float *)malloc(sizeof(float) *
						abp->nangles*abp->nangles);
	    input_orient = -1; output_orient = 1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++)
		    if (bsd.tb != NULL) {	/* use tb if we have it */
			sum = 0;		/* average over patches */
			for (n = npsamps; n-- > 0; ) {
			    fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			    bi_getvec(vin, i+urand(n), abp);
			    ec = SDevalBSDF(&sv, vout, vin, &bsd);
			    if (ec != SDEnone)
				goto err;
			    sum += sv.cieY;
			}
			printf("\t%.3e\n", sum/npsamps);
			if (trans_mtx != NULL)
			    trans_mtx[i*abp->nangles + j] = sum/npsamps;
		    } else {			/* else transpose tf */
			printf("\t%.3e\n", trans_mtx[i*abp->nangles + j]);
		    }
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
						/* derived front transmission */
	if (bsd.tf == NULL && trans_mtx != NULL) {
	    input_orient = 1; output_orient = -1;
	    data_prologue();
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++)
		    printf("\t%.3e\n", trans_mtx[j*abp->nangles + i]);
		putchar('\n');			/* extra space between rows */
	    }
	    data_epilogue();
	}
	SDfreeBSDF(&bsd);			/* all done */
	if (trans_mtx != NULL)
		free(trans_mtx);
	return;
err:
	SDreportError(ec, stderr);
	exit(1);
}

/* Interpolate and output a BSDF function using Klems basis */
static void
eval_function(char *funame)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	int		assignD = (fundefined(funame) < 6);
	double		iovec[6];
	double		sum;
	int		i, j, n;

	initurand(npsamps);
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
			fi_getvec(iovec, i+urand(n), abp);
		    else
			bi_getvec(iovec, i+urand(n), abp);

		    if (assignD) {
			varset("Dx", '=', -iovec[3]);
			varset("Dy", '=', -iovec[4]);
			varset("Dz", '=', -iovec[5]);
			++eclock;
		    }
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
	float		bsdfarr[MAXPATCHES*MAXPATCHES];
	FVECT		vin, vout;
	RBFNODE		*rbf;
	double		sum;
	int		i, j, n;
						/* sanity check */
	if (abp->nangles > MAXPATCHES) {
		fprintf(stderr, "%s: too many patches!\n", progname);
		exit(1);
	}
	data_prologue();			/* begin output */
	for (i = 0; i < abp->nangles; i++) {
	    if (input_orient > 0)		/* use incident patch center */
		fi_getvec(vin, i+.5*(i>0), abp);
	    else
		bi_getvec(vin, i+.5*(i>0), abp);

	    rbf = advect_rbf(vin);		/* compute radial basis func */

	    for (j = 0; j < abp->nangles; j++) {
	        sum = 0;			/* sample over exiting patch */
		for (n = npsamps; n--; ) {
		    if (output_orient > 0)
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
		    else
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);

		    sum += eval_rbfrep(rbf, vout) / vout[2];
		}
		bsdfarr[j*abp->nangles + i] = sum*output_orient/npsamps;
	    }
	    if (rbf != NULL)
		free(rbf);
	}
	n = 0;					/* write out our matrix */
	for (j = 0; j < abp->nangles; j++) {
	    for (i = 0; i < abp->nangles; i++)
		printf("\t%.3e\n", bsdfarr[n++]);
	    putchar('\n');
	}
	data_epilogue();			/* finish output */
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
	for (i = 1; i < argc && (argv[i][0] == '-') | (argv[i][0] == '+'); i++)
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
			fprintf(stderr, "\tor 3 arguments using Dx,Dy,Dz: bsdf(ix,iy,iz)\n",
					progname);
			goto userr;
		}
		++eclock;
		xml_header(argc, argv);			/* start XML output */
		xml_prologue(NULL);
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
						/* XML input? */
	if (i == argc-1 && (cp = argv[i]+strlen(argv[i])-4) > argv[i] &&
				!strcasecmp(cp, ".xml")) {
		xml_header(argc, argv);		/* start XML output */
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
			if (!nbsdf++) {		/* start XML on first dist. */
				xml_header(argc, argv);
				xml_prologue(NULL);
			}
			eval_rbf();
		}
		xml_epilogue();			/* finish XML output & exit */
		return(0);
	}
	SET_FILE_BINARY(stdin);			/* load from stdin */
	if (!load_bsdf_rep(stdin))
		return(1);
	xml_header(argc, argv);			/* start XML output */
	xml_prologue(NULL);
	eval_rbf();				/* resample dist. */
	xml_epilogue();				/* finish XML output & exit */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-n spp][-h|-q][bsdf.sir ..] > bsdf.xml\n", progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q] bsdf_in.xml > bsdf_out.xml\n", progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q][{+|-}for[ward]][{+|-}b[ackward]][-e expr][-f file] bsdf_func > bsdf.xml\n",
				progname);
	return(1);
}
