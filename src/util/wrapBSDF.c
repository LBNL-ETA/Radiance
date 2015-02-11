#ifndef lint
static const char RCSid[] = "$Id: wrapBSDF.c,v 2.1 2015/02/11 18:08:10 greg Exp $";
#endif
/*
 * Wrap BSDF data in valid WINDOW XML wrapper
 *
 *	G. Ward		January 2015
 */

#include <ctype.h>
#include "rtio.h"
#include "rtprocess.h"
#include "ezxml.h"
#include "bsdf.h"
#include "bsdf_m.h"
					/* XML template file names */
const char	def_template[] = "minimalBSDFt.xml";
const char	win6_template[] = "window6BSDFt.xml";

const char	stdin_name[] = "<stdin>";

					/* input files (can be stdin_name) */
const char	*xml_input = NULL;
					/* unit for materials & geometry */
const char	*attr_unit = "millimeter";
const char	legal_units[] = "meter|foot|inch|centimeter|millimeter";
					/* system materials & geometry */
const char	*mgf_geometry = NULL;

					/* angle basis */
enum { ABdefault=-1, ABklemsFull=0, ABklemsHalf, ABklemsQuarter,
		ABtensorTree3, ABtensorTree4, ABend };

int		angle_basis = ABdefault;

					/* field IDs and nicknames */
struct s_fieldID {
	char		nickName[4];
	int		has_unit;
	const char	*fullName;
} XMLfieldID[] = {
	{"m", 0, "Manufacturer"},
	{"n", 0, "Name"},
	{"c", 0, "ThermalConductivity"},
	{"ef", 0, "EmissivityFront"},
	{"eb", 0, "EmissivityBack"},
	{"tir", 0, "TIR"},
	{"eo", 0, "EffectiveOpennessFraction"},
	{"t", 1, "Thickness"},
	{"h", 1, "Height"},
	{"w", 1, "Width"},
	{"\0", -1, NULL}	/* terminator */
};
					/* field assignments */
#define MAXASSIGN	12
const char	*field_assignment[MAXASSIGN];
int		nfield_assign = 0;
#define FASEP	';'

					/* data file(s) & spectra */
enum { DTtransForward, DTtransBackward, DTreflForward, DTreflBackward };

enum { DSsolar=-1, DSnir=-2, DSxbar31=-3, DSvisible=-4, DSzbar31=-5 };

#define MAXFILES	20

struct s_dfile {
	const char	*fname;		/* input data file name */
	short		type;		/* BSDF data type */
	short		spectrum;	/* BSDF sensor spectrum */
} data_file[MAXFILES];

int		ndataf = 0;		/* number of data files */

const char	*spectr_file[MAXFILES];	/* custom spectral curve input */

const char	top_level_name[] = "WindowElement";

static char	*basis_definition[] = {

	"\t\t<IncidentDataStructure>Columns</IncidentDataStructure>\n"
	"\t\t<AngleBasis>\n"
	"\t\t\t<AngleBasisName>LBNL/Klems Full</AngleBasisName>\n"
	"\t\t</AngleBasis>\n",

	"\t\t<IncidentDataStructure>Columns</IncidentDataStructure>\n"
	"\t\t<AngleBasis>\n"
	"\t\t\t<AngleBasisName>LBNL/Klems Half</AngleBasisName>\n"
	"\t\t</AngleBasis>\n",

	"\t\t<IncidentDataStructure>Columns</IncidentDataStructure>\n"
	"\t\t<AngleBasis>\n"
	"\t\t\t<AngleBasisName>LBNL/Klems Quarter</AngleBasisName>\n"
	"\t\t</AngleBasis>\n",

	"\t\t<IncidentDataStructure>TensorTree3</IncidentDataStructure>\n",

	"\t\t<IncidentDataStructure>TensorTree4</IncidentDataStructure>\n",
};

/* Copy data from file descriptor to stdout and close */
static int
copy_and_close(int fd)
{
	int	ok = 1;
	char	buf[8192];
	int	n;

	if (fd < 0)
		return 0;
	while ((n = read(fd, buf, sizeof(buf))) > 0)
		if (write(fileno(stdout), buf, n) != n) {
			ok = 0;
			break;
		}
	ok &= (n == 0);
	close(fd);
	return ok;
}

/* Allocate and assign string from file or stream */
static char *
input2str(const char *inpspec)
{
	FILE	*fp = NULL;
	char	*str;
	int	len, pos, n;

	if (inpspec == NULL || !*inpspec)
		return "";
	if (inpspec == stdin_name) {		/* read from stdin */
		fp = stdin;
	} else if (*inpspec == '!') {		/* read from command */
		fp = popen(inpspec+1, "r");
		if (fp == NULL) {
			fprintf(stderr, "Cannot start process '%s'\n",
					inpspec+1);
			return "";
		}
	} else {				/* else load file */
		int	fd = open(inpspec, O_RDONLY);
		if (fd < 0) {
			fprintf(stderr, "%s: cannot open\n", inpspec);
			return "";
		}
		len = lseek(fd, 0L, SEEK_END);
		if (len < 0) {
			fprintf(stderr, "%s: seek error\n", inpspec);
			close(fd);
			return "";
		}
		lseek(fd, 0L, SEEK_SET);
		str = (char *)malloc(len+1);
		if (str == NULL) {
			close(fd);
			goto memerr;
		}
		if (read(fd, str, len) != len) {
			fprintf(stderr, "%s: read error\n", inpspec);
			free(str);
			close(fd);
			return "";
		}
		str[len] = '\0';
		close(fd);
		return str;
	}
						/* reading from stream */
	str = (char *)malloc((len=8192)+1);
	if (str == NULL)
		goto memerr;
	pos = 0;
	while ((n = read(fileno(fp), str+pos, len-pos)) > 0)
		if ((pos += n) >= len) {	/* need more space? */
			str = (char *)realloc(str, (len += len>>2) + 1);
			if (str == NULL)
				goto memerr;
		}
	if (n < 0) {
		fprintf(stderr, "%s: read error\n", inpspec);
		free(str);
		str = "";
	} else {				/* tidy up result */
		str[pos] = '\0';
		str = (char *)realloc(str, (len=pos)+1);
		if (str == NULL)
			goto memerr;
	}
	if (*inpspec != '!')
		fclose(fp);
	else if (pclose(fp))
		fprintf(stderr, "Error running command '%s'\n", inpspec+1);
	return str;
memerr:
	fprintf(stderr, "%s: error allocating memory\n", inpspec);
	if (fp != NULL)
		(*inpspec == '!') ? pclose(fp) : fclose(fp);
	return "";
}

/* Make material assignments in field_assignment to XML fields */
static int
mat_assignments(const char *caller, const char *fn, ezxml_t wtl)
{
	int	i;

	if (!nfield_assign)
		return 1;
	wtl = ezxml_child(wtl, "Material");
	if (wtl == NULL) {
		fprintf(stderr, "%s: missing <Material> tag\n", fn);
		return 0;
	}
	for (i = 0; i < nfield_assign; i++) {
		const char	*fnext = field_assignment[i];
		for ( ; ; ) {
			ezxml_t	fld;
			char	sbuf[512];
			int	j;

			while (isspace(*fnext))
				++fnext;
			if (!*fnext)
				break;
			for (j = 0; (*fnext != '=') & !isspace(*fnext); ) {
				if (!*fnext | (*fnext == FASEP) |
						(j >= sizeof(sbuf)-1)) {
					fprintf(stderr,
					"%s: bad tag name in assignment '%s'\n",
						caller, field_assignment[i]);
					return 0;
				}
				sbuf[j++] = *fnext++;
			}
			sbuf[j] = '\0';		/* check nick-names */
			for (j = 0; XMLfieldID[j].nickName[0]; j++)
				if (!strcasecmp(sbuf, XMLfieldID[j].nickName)) {
					strcpy(sbuf, XMLfieldID[j].fullName);
					break;
				}
						/* check if tag exists */
			fld = ezxml_child(wtl, sbuf);
			if (fld == NULL) {	/* otherwise, create one */
				fprintf(stderr, "%s: warning - adding tag <%s>\n",
					fn, sbuf);
				fld = ezxml_add_child_d(wtl, sbuf, strlen(wtl->txt));
			}
			if (XMLfieldID[j].has_unit)
				ezxml_set_attr(fld, "unit", attr_unit);
			while (isspace(*fnext))
				++fnext;
			if (*fnext++ != '=') {
				fprintf(stderr,
					"%s: missing '=' in assignment '%s'\n",
						caller, field_assignment[i]);
				return 0;
			}
			for (j = 0; *fnext & (*fnext != FASEP); ) {
				if (j >= sizeof(sbuf)-1) {
					fprintf(stderr,
					"%s: field too long in '%s'\n",
						caller, field_assignment[i]);
					return 0;
				}
				sbuf[j++] = *fnext++;
			}
			sbuf[j] = '\0';
			ezxml_set_txt_d(fld, sbuf);
			fnext += (*fnext == FASEP);
		}
	}
	return 1;		/* no errors */
}

/* Complete angle basis specification */
static int
finish_angle_basis(ezxml_t ab)
{
	const char	*bn = ezxml_txt(ezxml_child(ab, "AngleBasisName"));
	int		i, n = nabases;
	char		buf[32];

	if (!*bn) {
		fputs("Internal error - missing <AngleBasisName>!\n", stderr);
		return 0;
	}
	while (n-- > 0)
		if (!strcasecmp(bn, abase_list[n].name))
			break;
	if (n < 0) {
		fprintf(stderr, "Internal error - unknown angle basis '%s'", bn);
		return 0;
	}
	for (i = 0; abase_list[n].lat[i].nphis; i++) {
		ezxml_t	tb, abb = ezxml_add_child(ab, "AngleBasisBlock",
						strlen(ab->txt));
		sprintf(buf, "%g", i ?
		.5*(abase_list[n].lat[i].tmin + abase_list[n].lat[i+1].tmin) :
				.0);
		ezxml_set_txt_d(ezxml_add_child(abb,"Theta",strlen(abb->txt)), buf);
		sprintf(buf, "%d", abase_list[n].lat[i].nphis);
		ezxml_set_txt_d(ezxml_add_child(abb,"nPhis",strlen(abb->txt)), buf);
		tb = ezxml_add_child(abb, "ThetaBounds", strlen(abb->txt));
		sprintf(buf, "%g", abase_list[n].lat[i].tmin);
		ezxml_set_txt(ezxml_add_child(tb,"LowerTheta",strlen(tb->txt)), buf);
		sprintf(buf, "%g", abase_list[n].lat[i+1].tmin);
		ezxml_set_txt(ezxml_add_child(tb,"UpperTheta",strlen(tb->txt)), buf);
	}
	return 1;
}

/* Determine our angle basis from current tags */
static int
determine_angle_basis(const char *fn, ezxml_t wtl)
{
	const char	*ids;
	int		i;

	wtl = ezxml_child(wtl, "DataDefinition");
	if (wtl == NULL)
		return -1;
	ids = ezxml_txt(ezxml_child(wtl, "IncidentDataStructure"));
	if (ids == NULL)
		return -1;
	for (i = 0; i < ABend; i++) {
		ezxml_t	parsed = ezxml_parse_str(basis_definition[i],
					strlen(basis_definition[i]));
		int	match = 0;
		if (!strcmp(ids, ezxml_txt(ezxml_child(parsed,
					"IncidentDataStructure")))) {
			const char	*abn0 = ezxml_txt(
					ezxml_child(ezxml_child(wtl,
					"AngleBasis"), "AngleBasisName"));
			const char	*abn1 = ezxml_txt(
					ezxml_child(ezxml_child(parsed,
					"AngleBasis"), "AngleBasisName"));
			match = !strcmp(abn0, abn1);
		}
		ezxml_free(parsed);
		if (match)
			return i;
	}
	return -1;
}

/* Write out BSDF data block with surrounding tags */
static int
writeBSDFblock(const char *caller, struct s_dfile *df)
{
	char	*cp;

	puts("\t<WavelengthData>");
	puts("\t\t<LayerNumber>System</LayerNumber>");
	switch (df->spectrum) {
	case DSvisible:
		puts("\t\t<Wavelength unit=\"Integral\">Visible</Wavelength>");
		puts("\t\tSourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
		puts("\t\t<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>");
		break;
	case DSxbar31:
		puts("\t\t<Wavelength unit=\"Integral\">CIE-X</Wavelength>");
		puts("\t\tSourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
		puts("\t\t<DetectorSpectrum>ASTM E308 1931 X.dsp</DetectorSpectrum>");
		break;
	case DSzbar31:
		puts("\t\t<Wavelength unit=\"Integral\">CIE-Z</Wavelength>");
		puts("\t\tSourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
		puts("\t\t<DetectorSpectrum>ASTM E308 1931 Z.dsp</DetectorSpectrum>");
		break;
	case DSsolar:
		puts("\t\t<Wavelength unit=\"Integral\">Solar</Wavelength>");
		puts("\t\tSourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
		puts("\t\t<DetectorSpectrum>None</DetectorSpectrum>");
		break;
	case DSnir:
		puts("\t\t<Wavelength unit=\"Integral\">NIR</Wavelength>");
		puts("\t\tSourceSpectrum>PLACE_HOLDER</SourceSpectrum>");
		puts("\t\t<DetectorSpectrum>PLACE_HOLDER</DetectorSpectrum>");
		break;
	default:
		cp = strrchr(spectr_file[df->spectrum], '.');
		if (cp != NULL)
			*cp = '\0';
		printf("\t\t<Wavelength unit=\"Integral\">%s</Wavelength>\n",
				spectr_file[df->spectrum]);
		if (cp != NULL)
			*cp = '.';
		puts("\t\tSourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>");
		printf("\t\t<DetectorSpectrum>%s</DetectorSpectrum>\n",
				spectr_file[df->spectrum]);
		break;
	}
	puts("\t\t<WavelengthDataBlock>");
	fputs("\t\t\t<WavelengthDataDirection>", stdout);
	switch (df->type) {
	case DTtransForward:
		fputs("Transmission Front", stdout);
		break;
	case DTtransBackward:
		fputs("Transmission Back", stdout);
		break;
	case DTreflForward:
		fputs("Reflection Front", stdout);
		break;
	case DTreflBackward:
		fputs("Reflection Back", stdout);
		break;
	default:
		fprintf(stderr, "%s: internal - bad BSDF type (%d)\n", caller, df->type);
		return 0;
	}
	puts("</WavelengthDataDirection>");
	switch (angle_basis) {
	case ABklemsFull:
		puts("\t\t\t<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>");
		break;
	case ABklemsHalf:
		puts("\t\t\t<ColumnAngleBasis>LBNL/Klems Half</ColumnAngleBasis>");
		break;
	case ABklemsQuarter:
		puts("\t\t\t<ColumnAngleBasis>LBNL/Klems Quarter</ColumnAngleBasis>");
		break;
	case ABtensorTree3:
	case ABtensorTree4:
		puts("\t\t\t<AngleBasis>LBNL/Shirley-Chiu</AngleBasis>");
		break;
	default:
		fprintf(stderr, "%s: bad angle basis (%d)\n", caller, angle_basis);
		return 0;
	}
	puts("\t\t\t<ScatteringDataType>BTDF</ScatteringDataType>");
	puts("\t\t\t<ScatteringData>");
	fflush(stdout);
	if (df->fname == stdin_name) {
		copy_and_close(fileno(stdin));
	} else if (*df->fname != '!') {
		if (!copy_and_close(open(df->fname, O_RDONLY))) {
			fprintf(stderr, "%s: error reading from '%s'\n",
					caller, df->fname);
			return 0;
		}
	} else if (system(df->fname+1)) {
		fprintf(stderr, "%s: error running '%s'\n", caller, df->fname+1);
		return 0;
	}
	puts("\t\t\t</ScatteringData>");
	puts("\t\t</WavelengthDataBlock>");
	puts("\t</WavelengthData>");
	return 1;
}

/* Write out XML, interpolating BSDF data block(s) */
static int
writeBSDF(const char *caller, ezxml_t fl)
{
	char	*xml = ezxml_toxml(fl);
	int	ei, i;
						/* locate trailer */
	for (ei = strlen(xml)-strlen("</Layer></Optical></WindowElement>");
			ei >= 0; ei--)
		if (!strncmp(xml+ei, "</Layer>", 8))
			break;
	if (ei < 0) {
		fprintf(stderr, "%s: internal - cannot find trailer\n",
				caller);
		free(xml);
		return 0;
	}
	fflush(stdout);				/* write previous XML info. */
	if (write(fileno(stdout), xml, ei) != ei) {
		free(xml);
		return 0;
	}
	for (i = 0; i < ndataf; i++)		/* interpolate new data */
		if (!writeBSDFblock(caller, &data_file[i])) {
			free(xml);
			return 0;
		}
	fputs(xml+ei, stdout);			/* write trailer */
	free(xml);
	return (fflush(stdout) == 0);
}

/* Insert BSDF data into XML wrapper */
static int
wrapBSDF(const char *caller)
{
	const char	*xml_path = getpath((char *)xml_input, getrlibpath(), R_OK);
	ezxml_t		fl, wtl;
					/* load previous XML/template */
	if (xml_path == NULL) {
		fprintf(stderr, "%s: cannot find XML file named '%s'\n",
				caller, xml_input==NULL ? "NULL" : xml_input);
		return 0;
	}
	fl = ezxml_parse_file(xml_path);
	if (fl == NULL) {
		fprintf(stderr, "%s: cannot open XML path '%s'\n",
				caller, xml_path);
		return 0;
	}
	if (ezxml_error(fl)[0]) {
		fprintf(stderr, "%s: error in XML %s: %s\n", caller, xml_path,
				ezxml_error(fl));
		goto failure;
	}
	if (strcmp(ezxml_name(fl), top_level_name)) {
		fprintf(stderr, "%s: top level in XML '%s' not '%s'\n",
				caller, xml_path, top_level_name);
		goto failure;
	}
	wtl = ezxml_child(fl, "FileType");
	if (wtl != NULL && strcmp(ezxml_txt(wtl), "BSDF")) {
		fprintf(stderr, "%s: wrong FileType in XML '%s' (must be 'BSDF')",
				caller, xml_path);
		goto failure;
	}
	wtl = ezxml_child(ezxml_child(fl, "Optical"), "Layer");
	if (wtl == NULL) {
		fprintf(stderr, "%s: no optical layers in XML '%s'",
				caller, xml_path);
		goto failure;
	}
					/* make material assignments */
	if (!mat_assignments(caller, xml_path, wtl)) {
		goto failure;
	}
	if (mgf_geometry != NULL) {	/* add geometry if specified */
		ezxml_t	geom = ezxml_child(wtl, "Geometry");
		if (geom == NULL)
			geom = ezxml_add_child(wtl, "Geometry", 0);
		ezxml_set_attr(geom, "format", "MGF");
		geom = ezxml_child(geom, "MGFblock");
		if (geom == NULL) {
			geom = ezxml_child(wtl, "Geometry");
			geom = ezxml_add_child(geom, "MGFblock", strlen(geom->txt));
		}
		ezxml_set_attr(geom, "unit", attr_unit);
		ezxml_set_txt(geom, input2str(mgf_geometry));
		if (geom->txt[0])
			ezxml_set_flag(geom, EZXML_TXTM);
	}
					/* check basis */
	if (angle_basis != ABdefault) {
		ezxml_t	ab, dd = ezxml_child(wtl, "DataDefinition");
		if (dd == NULL) {
			dd = ezxml_add_child(wtl, "DataDefinition", 0);
		} else {
			if (ezxml_child(dd, "DataDefinition") != NULL)
				fprintf(stderr,
		"%s: warning - replacing existing <DataDefinition> in '%s'\n",
						caller, xml_path);
			while (dd->child != NULL)
				ezxml_remove(dd->child);
		}
		ezxml_insert(ezxml_parse_str(basis_definition[angle_basis],
					strlen(basis_definition[angle_basis])),
				dd, 0);
		if ((ab = ezxml_child(dd, "AngleBasis")) != NULL &&
				!finish_angle_basis(ab))
			goto failure;
	} else if ((angle_basis = determine_angle_basis(xml_path, wtl)) < 0) {
		fprintf(stderr, "%s: need -a option to set angle basis\n",
				caller);
		goto failure;
	}
					/* write & add BSDF data blocks */
	if (!writeBSDF(caller, fl))
		goto failure;
	ezxml_free(fl);			/* all done */
	return 1;
failure:
	ezxml_free(fl);
	return 0;
}

/* Report usage and exit */
static void
UsageExit(const char *pname)
{
	fputs("Usage: ", stderr);
	fputs(pname, stderr);
	fputs(" [options] [input.xml]\n", stderr);
	exit(1);
}

/* Load XML file and use to wrap BSDF data (or modify fields) */
int
main(int argc, char *argv[])
{
	int	cur_spectrum = DSvisible;
	int	ncust_spec = 0;
	int	used_stdin = 0;
	int	units_set = 0;
	int	i;
					/* get/check arguments */
	for (i = 1; i < argc && argv[i][0] == '-'; i++) {
		switch (argv[i][1]) {
		case 'W':		/* customize for WINDOW 6 output */
			xml_input = win6_template;
			continue;
		case 'f':		/* field assignment(s) */
			if (++i >= argc)
				UsageExit(argv[0]);
			if (nfield_assign >= MAXASSIGN) {
				fprintf(stderr, "%s: too many -f options",
						argv[0]);
				return 1;
			}
			field_assignment[nfield_assign++] = argv[i];
			continue;
		case 'u':		/* unit */
			if (++i >= argc)
				UsageExit(argv[0]);
			if (units_set++) {
				fprintf(stderr, "%s: only one -u option allowed\n",
						argv[0]);
				return 1;
			}
			if (strstr(argv[i], legal_units) == NULL) {
				fprintf(stderr, "%s: unit must be one of (%s)\n",
						argv[0], legal_units);
				return 1;
			}
			attr_unit = argv[i];
			continue;
		case 'a':		/* angle basis */
			if (++i >= argc)
				UsageExit(argv[0]);
			if (angle_basis != ABdefault) {
				fprintf(stderr, "%s: only one -a option allowed\n",
						argv[0]);
				return 1;
			}
			if (!strcasecmp(argv[i], "kf"))
				angle_basis = ABklemsFull;
			else if (!strcasecmp(argv[i], "kh"))
				angle_basis = ABklemsHalf;
			else if (!strcasecmp(argv[i], "kq"))
				angle_basis = ABklemsQuarter;
			else if (!strcasecmp(argv[i], "t3"))
				angle_basis = ABtensorTree3;
			else if (!strcasecmp(argv[i], "t4"))
				angle_basis = ABtensorTree4;
			else
				UsageExit(argv[0]);
			continue;
		case 't':		/* transmission */
			if (i >= argc-1)
				UsageExit(argv[0]);
			if (ndataf >= MAXFILES) {
				fprintf(stderr, "%s: too many data files\n",
						argv[0]);
				return 1;
			}
			if (!strcmp(argv[i], "-tf"))
				data_file[ndataf].type = DTtransForward;
			else if (!strcmp(argv[i], "-tb"))
				data_file[ndataf].type = DTtransBackward;
			else
				UsageExit(argv[0]);
			if (!strcmp(argv[++i], "-")) {
				if (used_stdin++) UsageExit(argv[i]);
				argv[i] = (char *)stdin_name;
			}
			data_file[ndataf].fname = argv[i];
			data_file[ndataf].spectrum = cur_spectrum;
			ndataf++;
			continue;
		case 'r':		/* reflection */
			if (i >= argc-1)
				UsageExit(argv[0]);
			if (ndataf >= MAXFILES) {
				fprintf(stderr, "%s: too many data files\n",
						argv[0]);
				return 1;
			}
			if (!strcmp(argv[i], "-rf"))
				data_file[ndataf].type = DTreflForward;
			else if (!strcmp(argv[i], "-rb"))
				data_file[ndataf].type = DTreflBackward;
			else
				UsageExit(argv[0]);
			if (!strcmp(argv[++i], "-")) {
				if (used_stdin++) UsageExit(argv[i]);
				argv[i] = (char *)stdin_name;
			}
			data_file[ndataf].fname = argv[i];
			data_file[ndataf++].spectrum = cur_spectrum;
			continue;
		case 's':		/* spectrum name or input file */
			if (++i >= argc)
				UsageExit(argv[0]);
			if (!strcasecmp(argv[i], "Solar"))
				cur_spectrum = DSsolar;
			else if (!strcasecmp(argv[i], "Visible") ||
					!strcasecmp(argv[i], "CIE-Y"))
				cur_spectrum = DSvisible;
			else if (!strcasecmp(argv[i], "CIE-X"))
				cur_spectrum = DSxbar31;
			else if (!strcasecmp(argv[i], "CIE-Z"))
				cur_spectrum = DSzbar31;
			else if (!strcasecmp(argv[i], "NIR"))
				cur_spectrum = DSnir;
			else {
				if (!strcmp(argv[i], "-")) {
					fprintf(stderr,
					  "%s: cannot read spectra from stdin",
							argv[0]);
					return 1;
				}
				cur_spectrum = ncust_spec;
				spectr_file[ncust_spec++] = argv[i];
			}
			continue;
		case 'g':		/* MGF geometry file */
			if (i >= argc-1)
				UsageExit(argv[0]);
			if (mgf_geometry != NULL) {
				fprintf(stderr, "%s: only one -g option allowed\n",
						argv[0]);
				return 1;
			}
			if (!strcmp(argv[++i], "-")) {
				if (used_stdin++) UsageExit(argv[i]);
				argv[i] = (char *)stdin_name;
			}
			mgf_geometry = argv[i];
			continue;
		case '\0':		/* input XML from stdin */
			break;
		default:
			UsageExit(argv[0]);
			break;
		}
		break;
	}
doneOptions:					/* get XML input */
	if (i >= argc) {
		if (xml_input == NULL)
			xml_input = def_template;
	} else if ((i < argc-1) | (xml_input != NULL)) {
		fprintf(stderr, "%s: only one XML input allowed\n", argv[0]);
		UsageExit(argv[0]);
	} else if (!strcmp(argv[i], "-")) {
		if (used_stdin++) UsageExit(argv[0]);
		xml_input = stdin_name;
	} else {
		xml_input = argv[i];
	}
	if ((xml_input == win6_template) & (angle_basis == ABdefault))
		angle_basis = ABklemsFull;
						/* wrap it! */
	return !wrapBSDF(argv[0]);
}
