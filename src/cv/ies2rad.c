#ifndef lint
static const char	RCSid[] = "$Id: ies2rad.c,v 2.35 2021/11/29 16:07:36 greg Exp $";
#endif
/*
 * ies2rad -- Convert IES luminaire data to Radiance description
 *
 * ies2rad converts an IES LM-63 luminare description to a Radiance
 * luminaire description.  In addition, ies2rad manages a local
 * database of Radiance luminaire files.
 *
 * Ies2rad generates two or three files for each luminaire. For a
 * luminaire named LUM, ies2rad will generate LUM.rad, a Radiance
 * scene description file which describes the light source, LUM.dat,
 * which contains the photometric data from the IES LM-63 file, and
 * (if tilt data is provided) LUM%.dat, which contains the tilt data
 * from the IES file.
 *
 * Ies2rad is supported by the Radiance function files source.cal and
 * tilt.cal, which transform the coordinates in the IES data into
 * Radiance (θ,φ) luminaire coordinates and then apply photometric and
 * tilt data to generate Radiance light. θ is altitude from the
 * negative z-axis and φ is azimuth from the positive x-axis,
 * increasing towards the positive y-axis. This system matches none of
 * the usual goniophotometric conventions, but it is closest to IES
 * type C; V in type C photometry is θ in Radiance and L is -φ.
 *
 * The ies2rad scene description for a luminaire LUM, with tilt data,
 * uses the following Radiance scene description primitives:
 *
 *     void brightdata LUM_tilt
 *     …
 *     LUM_tilt brightdata LUM_dist
 *     …
 *     LUM_dist light LUM_light
 *     …
 *     LUM_light surface1 name1
 *     …
 *     LUM_light surface2 name2
 *     …
 *     LUM_light surface_n name_n
 *
 * Without tilt data, the primitives are:
 *
 *     void brightdata LUM_dist
 *     …
 *     LUM_dist light LUM_light
 *     …
 *     LUM_light surface1 name1
 *     …
 *     LUM_light surface2 name2
 *     …
 *     LUM_light surface_n name_n
 *
 * As many surfaces are given as required to describe the light
 * source. Illum may be used rather than light so that a visible form
 * (impostor) may be given to the luminaire, rather than a simple
 * glowing shape. If an impostor is provided, it must be wholly
 * contained within the illum and if it provides impostor light
 * sources, those must be given with glow, so that they do not
 * themselves illuminate the scene, providing incorrect results.
 *
 * Overview of the LM-63 file format
 * =================================
 * Here we offer a summary of the IESNA LM-63 photometry file format
 * for the perplexed reader.  Dear reader, do remember that this is
 * our interpretation of the five different versions of the standard.
 * When our interpretation of the standard conflicts with the official
 * standard, the official document is to be respected.  In conflicts
 * with practice, do take into account the robustness principle and be
 * permissive, accepting reasonable deviations from the standard.
 *
 * LM-63 files are organized as a version tag, followed by a series of
 * luminaire data sets.  The luminaire data sets, in turn, are
 * organized into a label, a tilt data section, and a photometric data
 * section.  Finally, the data sections are organized into records,
 * which are made up of lines of numeric values delimited by spaces or
 * commas.  Lines are delimited by CR LF sequences.  Records are made
 * up of one or more lines, and every record must be made up of some
 * number of complete lines, but there is no delimiter which makes the
 * end of a record.  The first records of the tilt and photometric
 * data sections have fixed numbers of numeric values; the initial
 * records contain counts that describe the remaining records.
 *
 * Ies2rad allows only one luminaire data set per file.
 *
 * The tilt section is made up of exactly four records; the second gives
 * the number of values in the third and fourth records.
 *
 * The photometric section begins with two records, which give both the
 * number of records following and the number of values in each of the
 * following records.
 *
 * The original 1986 version of LM-63 does not have a version tag.
 * 
 * The 1986, 1991, and 1995 versions allow 80 characters for the label
 * lines and the "TILT=" line which begins the tilt data section, and
 * 132 characters thereafter.  (Those counts do not include the CR LF
 * line terminator.)  The 2002 version dispenses with those limits,
 * allowing 256 characters per line, including the CR LF line
 * terminator.  The 2019 version does not specify a line length at
 * all.  Ies2rad allows lines of up to 256 characters and will accept
 * CR LF or LF alone as line terminators.
 *
 * In the 1986 version, the label is a series of free-form lines of up
 * to 80 characters.  In later versions, the label is a series of
 * lines of beginning with keywords in brackets with interpretation
 * rules which differ between versions.
 *
 * The tilt data section begins with a line beginning with "TILT=",
 * optionally followed by either a file name or four records of
 * numerical data.  The 2019 version no longer allows a file name to
 * be given.
 *
 * The main photometric data section contains two header records
 * followed by a record of vertical angles, a record of horizontal
 * angles, and one record of candela values for each horizontal angle.
 * Each record of candela values contains exactly one value for each
 * vertical angle. Data values in records are separated by spaces or
 * commas.  In keeping with the robustness principle, commas
 * surrounded by spaces will also be accepted as separators.
 *
 * The first header record of the photometric data section contains
 * exactly 10 values.  The second contains exactly 3 values.  Most of
 * the data values are floating point numbers; the exceptions are
 * various counts and enumerators, which are integers: the number of
 * lamps, the numbers of vertical and horizontal angles, the
 * photometric type identifier, and the units type identifier.  In the
 * 2019 version, a field with information about how the file was
 * generated has replaced a field unused since 1995; it is a textual
 * representation of a bit string, but may - we hope! - safely be
 * interpreted as a floating point number and decoded later.
 *
 * Style Note
 * ==========
 * The ies2rad code uses the "bsd" style. For emacs, this is set up
 * automatically in the "Local Variables" section at the end of the
 * file. For vim, use ":set tabstop=8 shiftwidth=8".
 *
 * History
 * =======
 * 
 *	07Apr90		Greg Ward
 *
 *  Fixed correction factor for flat sources 29Oct2001 GW
 *  Extensive comments added by Randolph Fritz May2018
 */

#include <math.h>
#include <ctype.h>

#include "rtio.h"
#include "color.h"
#include "paths.h"

#define PI		3.14159265358979323846

#define FAIL            (-1)
#define SUCCESS         0

/* floating point comparisons -- floating point numbers within FTINY
 * of each other are considered equal */
#define FTINY		1e-6
#define FEQ(a,b)	((a)<=(b)+FTINY&&(a)>=(b)-FTINY)

#define IESFIRSTVER 1986
#define IESLASTVER 2019

/* tilt specs
 *
 * This next series of definitions address metal-halide lamps, which
 * change their brightness depending on the angle at which they are
 * mounted. The section begins with "TILT=".  The constants in this
 * section are all defined in LM-63.
 *
 */

#define TLTSTR		"TILT="
#define TLTSTRLEN	5
#define TLTNONE		"NONE"
#define TLTINCL		"INCLUDE"
#define TLT_VERT	1
#define TLT_H0		2
#define TLT_H90		3

/* Constants from LM-63 files */

/* photometric types
 *
 * This enumeration reflects three different methods of measuring the
 * distribution of light from a luminaire -- "goniophotometry" -- and
 * the different coordinate systems related to these
 * goniophotometers.  All are described in IES standard LM-75-01.
 * Earlier and shorter descriptions may be found the LM-63 standards
 * from 1986, 1991, and 1995.
 *
 * ies2rad does not support type A photometry.
 *
 * In the 1986 file format, LM-63-86, 1 is used for type C and type A
 * photometric data.
 *
 */
#define PM_C		1
#define PM_B		2
#define PM_A		3

/* unit types */
#define U_FEET		1
#define U_METERS	2

/* string lengths */
/* Maximum length of a keyword, including brackets and NUL */
#define MAXKW           21
/* Maximum input line is 256 characters including CR LF and NUL at end. */
#define MAXLINE		257
#define MAXUNITNAME     64
#define RMAXWORD	76

/* Shapes defined in the IES LM-63 standards
 *
 * PH stands for photometric horizontal
 * PPH stands for perpendicular to photometric horizontal
 * Cylinders are vertical and circular unless otherwise stated
 *
 * The numbers assigned here are not part of any LM-63 standard; they
 * are for programming convenience.
 */
/* Error and not-yet-assigned constants */
#define IESERROR        -2
#define IESNONE         -1
/* Shapes */
#define IESPT            0
#define IESRECT          1
#define IESBOX           2
#define IESDISK          3
#define IESELLIPSE       4
#define IESVCYL          5
#define IESVECYL         6
#define IESSPHERE        7
#define IESELLIPSOID     8
#define IESHCYL_PH       9
#define IESHECYL_PH     10
#define IESHCYL_PPH     11
#define IESHECYL_PPH    12
#define IESVDISK_PH     13
#define IESVEL_PH       14

/* End of LM-63 related #defines */

/* file extensions */
#define T_RAD		".rad"
#define T_DST		".dat"
#define T_TLT		"%.dat"
#define T_OCT		".oct"

/* Radiance shape types
 * These #defines enumerate the shapes of the Radiance objects which
 * emit the light.
 */
#define RECT		1
#define DISK		2
#define SPHERE		3

/* 1mm.  The diameter of a point source luminaire model. Also the minimum
 * size (in meters) that the luminous opening of a luminaire must have
 * to be treated as other than a point source. */
#define MINDIM		.001

/* feet to meters */
/* length_in_meters = length_in_feet * F_M */
#define F_M		.3048

/* abspath - return true if a path begins with a directory separator
 * or a '.' (current directory) */
#define abspath(p)	(ISDIRSEP((p)[0]) || (p)[0] == '.')

/* LM-63 related constants */
typedef struct {
	char *tag;
	int yr; } IESversions;

IESversions IESFILEVERSIONS[] = {
	{ "IESNA91", 1991 },
	{ "IESNA:LM-63-1995", 1995 },
	{ "IESNA:LM-63-2002", 2002 },
	{ "IES:LM-63-2019", 2019 },
	{ NULL, 1986 }
};

char *IESHAPENAMES[] = {
	"point", "rectangle", "box", "disk", "ellipse", "vertical cylinder",
	"vertical elliptical cylinder", "sphere", "ellipsoid",
	"horizontal cylinder along photometric horizontal",
	"horizontal elliptical cylinder along photometric horizontal",
	"horizontal cylinder perpendicular to photometric horizontal",
	"horizontal elliptical cylinder perpendicular to photometric horizontal",
	"vertical disk facing photometric horizontal",
	"vertical ellipse facing photometric horizontal" };

/* end of LM-63 related constants */

/* Radiance shape names */
char *RADSHAPENAMES[] = { "rectangle or box", "disk or cylinder", "sphere" };

/* Global variables.
 *
 * Mostly, these are a way of communicating command line parameters to
 * the rest of the program.
 */
static char	default_name[] = "default";

char	*libdir = NULL;			/* library directory location */
char	*prefdir = NULL;		/* subdirectory */
char	*lampdat = "lamp.tab";		/* lamp data file */

double	meters2out = 1.0;		/* conversion from meters to output */
char	*lamptype = NULL;		/* selected lamp type */
char	*deflamp = NULL;		/* default lamp type */
float	defcolor[3] = {1.,1.,1.};	/* default lamp color */
float	*lampcolor = defcolor;		/* pointer to current lamp color */
double	multiplier = 1.0;		/* multiplier for all light sources */
char	units[MAXUNITNAME] = "meters";	/* output units */
int	out2stdout = 0;			/* put out to stdout r.t. file */
int	instantiate = 0;		/* instantiate geometry */
double	illumrad = 0.0;			/* radius for illum sphere */

/* This struct describes the Radiance source object */
typedef struct {
	int	isillum;			/* do as illum */
	int	type;				/* RECT, DISK, SPHERE */
	double	mult;				/* candela multiplier */
	double	w, l, h;			/* width, length, height */
	double	area;				/* max. projected area */
	int	filerev;			/* IES file version */
	int     havelamppos;	                /* Lamp position was given */
	float   lamppos[2];	                /* Lamp position */
	int     iesshape;			/* Shape number */
	char   *warn;				/* Warning message */
} SRCINFO;				/* a source shape (units=meters) */

/* A count and pointer to the list of input file names */
int	gargc;				/* global argc */
char	**gargv;			/* global argv */

/* macros to scan numbers out of IES files
 *
 * fp is a file pointer.  scnint() places the number in the integer
 * indicated by ip; scnflt() places the number in the double indicated
 * by rp. The macros return 1 if successful, 0 if not.
 *
 */
#define scnint(fp,ip)	cvtint(ip,getword(fp))
#define scnflt(fp,rp)	cvtflt(rp,getword(fp))

/* The original (1986) version of LM-63 allows decimals points in
 * integers, so that, for instance, the number of lamps may be written
 * 3.0 (the number, obviously, must still be an integer.) This
 * confusing define accommodates that.  */
#define isint		isflt

/* IES file conversion functions */
static int ies2rad(char *inpname, char *outname);
static void initlamps(void);
static int dosource(SRCINFO *sinf, FILE *in, FILE *out, char *mod, char *name);
static int dotilt(FILE *in, FILE *out, char *dir, char *tltspec,
		char *dfltname, char *tltid);
static int cvgeometry(char *inpname, SRCINFO *sinf, char *outname, FILE *outfp);
static int cvtint(int *ip, char *wrd);
static int cvdata(FILE *in, FILE *out, int ndim, int npts[], double mult,
		double lim[][2]);
static int cvtflt(double *rp, char *wrd);
static int makeiesshape(SRCINFO *shp, double length, double width, double height);
static int makeillumsphere(SRCINFO *shp);
static int makeshape(SRCINFO *shp, double width, double length, double height);
static void makecylshape(SRCINFO *shp, double diam, double height);
static void makeelshape(SRCINFO *shp, double width, double length, double height);
static void makeecylshape(SRCINFO *shp, double width, double length, double height);
static void makeelshape(SRCINFO *shp, double width, double length, double height);
static void makeboxshape(SRCINFO *shp, double length, double width, double height);
static int makepointshape(SRCINFO *shp);
static int putsource(SRCINFO *shp, FILE *fp, char *mod, char *name,
		int dolower, int doupper, int dosides);
static void putrectsrc(SRCINFO *shp, FILE *fp, char *mod, char *name, int up);
static void putsides(SRCINFO *shp, FILE *fp, char *mod, char *name);
static void putdisksrc(SRCINFO *shp, FILE *fp, char *mod, char *name, int up);
static void putspheresrc(SRCINFO *shp, FILE *fp, char *mod, char *name);
static void putrect(SRCINFO *shp, FILE *fp, char *mod, char *name, char *suffix,
		int a, int b, int c, int d);
static void putpoint(SRCINFO *shp, FILE *fp, int p);
static void putcyl(SRCINFO *shp, FILE *fp, char *mod, char *name);
static void shapearea(SRCINFO *shp);

/* string and filename functions */
static int isprefix(char *p, char *s);
static char * matchprefix(char *p, char *s);
static char * tailtrunc(char *name);
static char * filename(char *path);
static char * libname(char *path, char *fname, char *suffix);
static char * getword(FILE *fp);
static char * fullnam(char *path, char *fname, char *suffix);

/* output function */
static void fpcomment(FILE *fp, char *prefix, char *s);

/* main - process arguments and run the conversion
 *
 * Refer to the man page for details of the arguments.
 *
 * Following Unix environment conventions, main() exits with 0 on
 * success and 1 on failure.
 *
 * ies2rad outputs either two or three files for a given IES
 * file. There is always a .rad file containing Radiance scene
 * description primitives and a .dat file for the photometric data. If
 * tilt data is given, that is placed in a separate .dat file.  So
 * ies2rad must have a filename to operate. Sometimes this name is the
 * input file name, shorn of its extension; sometimes it is given in
 * the -o option. But an output file name is required for ies2rad to
 * do its work.
 *
 * Older versions of the LM-63 standard allowed inclusion of multiple
 * luminaires in one IES file; this is not supported by ies2rad.
 *
 * This code sometimes does not check to make sure it has not run out
 * of arguments; this can lead to segmentation faults and perhaps
 * other errors.
 *
 */
int
main(
	int	argc,
	char	*argv[]
)
{
	char	*outfile = NULL;
	int	status;
	char	outname[RMAXWORD];
	double	d1;
	int	i;

	/* Scan the options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'd':		/* dimensions */
			if (argv[i][2] == '\0')
				goto badopt;
			if (argv[i][3] == '\0')
				d1 = 1.0;
			else if (argv[i][3] == '/') {
				d1 = atof(argv[i]+4);
				if (d1 <= FTINY)
					goto badopt;
			} else
				goto badopt;
			switch (argv[i][2]) {
			case 'c':		/* centimeters */
				if (FEQ(d1,10.))
					strcpy(units,"millimeters");
				else {
					strcpy(units,"centimeters");
					strcat(units,argv[i]+3);
				}
				meters2out = 100.*d1;
				break;
			case 'm':		/* meters */
				if (FEQ(d1,1000.))
					strcpy(units,"millimeters");
				else if (FEQ(d1,100.))
					strcpy(units,"centimeters");
				else {
					strcpy(units,"meters");
					strcat(units,argv[i]+3);
				}
				meters2out = d1;
				break;
			case 'i':		/* inches */
				strcpy(units,"inches");
				strcat(units,argv[i]+3);
				meters2out = d1*(12./F_M);
				break;
			case 'f':		/* feet */
				if (FEQ(d1,12.))
					strcpy(units,"inches");
				else {
					strcpy(units,"feet");
					strcat(units,argv[i]+3);
				}
				meters2out = d1/F_M;
				break;
			default:
				goto badopt;
			}
			break;
		case 'l':		/* library directory */
			libdir = argv[++i];
			break;
		case 'p':		/* prefix subdirectory */
			prefdir = argv[++i];
			break;
		case 'f':		/* lamp data file */
			lampdat = argv[++i];
			break;
		case 'o':		/* output file root name */
			outfile = argv[++i];
			break;
		case 's':		/* output to stdout */
			out2stdout = !out2stdout;
			break;
		case 'i':		/* illum */
			illumrad = atof(argv[++i]);
			break;
		case 'g':		/* instantiate geometry? */
			instantiate = !instantiate;
			break;
		case 't':		/* override lamp type */
			lamptype = argv[++i];
			break;
		case 'u':		/* default lamp type */
			deflamp = argv[++i];
			break;
		case 'c':		/* default lamp color */
			defcolor[0] = atof(argv[++i]);
			defcolor[1] = atof(argv[++i]);
			defcolor[2] = atof(argv[++i]);
			break;
		case 'm':		/* multiplier */
			multiplier = atof(argv[++i]);
			break;
		default:
		badopt:
			fprintf(stderr, "%s: bad option: %s\n",
					argv[0], argv[i]);
			exit(1);
		}
	/* Save pointers to the list of input file names */
	gargc = i;
	gargv = argv;

	/* get lamp data (if needed) */
	initlamps();

        /* convert ies file(s) */
	/* If an output file name is specified */
	if (outfile != NULL) {
		if (i == argc)
			/* If no input filename is given, use stdin as
			 * the source for the IES file */
			exit(ies2rad(NULL, outfile) == 0 ? 0 : 1);
		else if (i == argc-1)
			/* If exactly one input file name is given, use it. */
			exit(ies2rad(argv[i], outfile) == 0 ? 0 : 1);
		else
			goto needsingle; /* Otherwise, error. */
	} else if (i >= argc) {
		/* If an output file and an input file are not give, error. */
		fprintf(stderr, "%s: missing output file specification\n",
				argv[0]);
		exit(1);
	}
	/* If no input or output file is given, error. */
	if (out2stdout && i != argc-1)
		goto needsingle;
	/* Otherwise, process each input file in turn. */
	status = 0;
	for ( ; i < argc; i++) {
		tailtrunc(strcpy(outname,filename(argv[i])));
		if (ies2rad(argv[i], outname) != 0)
			status = 1;
	}
	exit(status);
needsingle:
	fprintf(stderr, "%s: single input file required\n", argv[0]);
	exit(1);
}

/* Initlamps -- If necessary, read lamp data table */
void
initlamps(void)				/* set up lamps */
{
	float	*lcol;
	int	status;

	/* If the lamp name is set to default, don't bother to read
	 * the lamp data table. */
	if (lamptype != NULL && !strcmp(lamptype, default_name) &&
			deflamp == NULL)
		return;

	if ((status = loadlamps(lampdat)) < 0) /* Load the lamp data table */
		exit(1);		       /* Exit if problems
						* with the file. */
	if (status == 0) {
                /* If can't open the file, just use the standard default lamp */
		fprintf(stderr, "%s: warning - no lamp data\n", lampdat);
		lamptype = default_name;
		return;
	}
	if (deflamp != NULL) {
                /* Look up the specified default lamp type */
		if ((lcol = matchlamp(deflamp)) == NULL)
			/* If it can't be found, use the default */
			fprintf(stderr,
				"%s: warning - unknown default lamp type\n",
					deflamp);
		else
			/* Use the selected default lamp color */
			copycolor(defcolor, lcol);
	}
	/* If a lamp type is specified and can be found, use it, and
	 * release the lamp data table memory; it won't be needed any more. */
	if (lamptype != NULL) {
		if (strcmp(lamptype, default_name)) {
			if ((lcol = matchlamp(lamptype)) == NULL) {
				fprintf(stderr,
					"%s: warning - unknown lamp type\n",
						lamptype);
				lamptype = default_name;
			} else
				copycolor(defcolor, lcol);
		}
		freelamps();			/* all done with data */
	}
	/* else keep lamp data */
}

/*
 * String functions
 */

/*
 * isprefix - return 1 (true) if p is a prefix of s, 0 otherwise
 * 
 * For this to work properly, s must be as long or longer than p.
 */
int
isprefix(char *p, char *s) {
	return matchprefix(p,s) != NULL;
}

/*
 * matchprefix - match p against s
 *
 * If p is a prefix of s, return a pointer to the character of s just
 * past p.
 *
 * For this to work properly, s must be as long or longer than p.
 */ 
char *
matchprefix(char *p, char *s) {
	int c;
	
	while ((c = *p++)) {
		if (c != *s++)
			return NULL;
	}
	return s;
}

/*
 * skipws - skip whitespace
 */
char *
skipws(char *s) {
	while (isspace(*s))
		s++;
	return s;
}

/*
 * streq - test strings for equality
 */
int
streq(char *s1, char *s2) {
	return strcmp(s1,s2) == 0;
}

/*
 * strneq - test strings for equality, with a length limit
 */
int
strneq(char *s1, char *s2, int n) {
	return strncmp(s1,s2,n) == 0;
}

/*
 * IES (LM-63) file functions
 */

/*
 * prockwd - process keywords on a label line
 * 
 * We're looking for four keywords: LAMP, LAMPCAT, LAMPPOSITION, and
 * LUMINOUSGEOMETRY.  Any other keywords are ignored.
 * 
 * LAMP and LAMPCAT are searched for a known lamp type name.
 * LAMPPOSITION is stored.
 * LUMINOUSGEOMETRY contains the name of an MGF file, which is stored.
 */
void
prockwd(char *bp, char *geomfile, char *inpname, SRCINFO *srcinfo) {
	char *kwbegin;
	int kwlen;
	
	bp = skipws(bp);	/* Skip leading whitespace. */
	if (*bp != '[')
		return;		/* If there's no keyword on this line,
				 * do nothing */
	kwbegin = bp;
	while (*bp && *bp != ']')	/* Skip to the end of the keyword or
				 * end of the buffer. */
		bp++;
	if (!(*bp))		/* If the keyword doesn't have a
				 * terminating ']', return. */
		return;
	kwlen = bp - kwbegin + 1;
	bp++;
	if (lampcolor == NULL && strneq("[LAMP]", kwbegin, kwlen)) 
		lampcolor = matchlamp(bp);
	else if (lampcolor == NULL && strneq("[LAMPCAT]", kwbegin, kwlen))
		lampcolor = matchlamp(bp);
	else if (strneq("[LUMINOUSGEOMETRY]", kwbegin, kwlen)) {
		bp = skipws(bp);	/* Skip leading whitespace. */
		strcpy(geomfile, inpname); /* Copy the input file path */
		/* Replace the filename in the input file path with
		 * the name of the MGF file.  Trailing spaces were
		 * trimmed before this routine was called. */
		strcpy(filename(geomfile), bp); 
		srcinfo->isillum = 1;
	}
	else if (strneq("[LAMPPOSITION]", kwbegin, kwlen)) {
		srcinfo->havelamppos = 1;
		sscanf(bp,"%f%f", &(srcinfo->lamppos[0]),
		       &(srcinfo->lamppos[1]));
	}
}

/*
 * iesversion - examine the first line of an IES file and return the version
 *
 * Returns the year of the version.  If the version is unknown,
 * returns 1986, since the first line of a 1986-format IES file can be
 * anything.
 */
int
iesversion(char *buf) {
	IESversions *v;

	for(v = IESFILEVERSIONS; v != NULL; v++)
		if (streq(v->tag,buf))
			return v->yr;
	return v->yr;
}


/*
 * File path operations
 *
 * These provide file path operations that operate on both MS-Windows
 * and *nix. They will ignore and pass, but will not necessarily
 * process correctly, Windows drive letters. Paths including Windows
 * UNC network names (\\server\folder\file) may also cause problems.
 *
 */

/*
 * stradd()
 *
 * Add a string to the end of a string, optionally concatenating a
 * file path separator character.  If the path already ends with a
 * path separator, no additional separator is appended.
 *
 */
char *
stradd(			/* add a string at dst */
	char	*dst,
	char	*src,
	int	sep
)
{
	if (src && *src) {
		do
			*dst++ = *src++;
		while (*src);
		if (sep && dst[-1] != sep)
			*dst++ = sep;
	}
	*dst = '\0';
	return(dst);
}

/*
 * fullnam () - return a usable path name for an output file
 */
char *
fullnam(
	char	*path,		/* The base directory path */
	char	*fname,		/* The file name */
	char	*suffix		/* A suffix, which usually contains
				 * a file name extension. */
)
{
	extern char *prefdir;
	extern char *libdir;

	if (prefdir != NULL && abspath(prefdir))
		/* If the subdirectory path is absolute or '.', just
		 * concatenate the names together */
		libname(path, fname, suffix);
	else if (abspath(fname))
		/* If there is no subdirectory, and the file name is
		 * an absolute path or '.', concatenate the path,
		 * filename, and suffix. */
		strcpy(stradd(path, fname, 0), suffix);
	else
		/* If the file name is relative, concatenate path,
		 * library directory, directory separator, file name,
		 * and suffix.  */
		libname(stradd(path, libdir, DIRSEP), fname, suffix);

	return(path);
}


/*
 * libname - convert a file name to a path
 */
char *
libname(
	char	*path,		/* The base directory path */
	char	*fname,		/* The file name */
	char	*suffix		/* A suffix, which usually contains
				 * a file name extension. */
)
{
	extern char *prefdir;	/* The subdirectory where the file
				 * name is stored. */

	if (abspath(fname))
		/* If the file name begins with '/' or '.', combine
		 * it with the path and attach the suffix */
		strcpy(stradd(path, fname, 0), suffix);
	else
		/* If the file name is relative, attach it to the
		 * path, include the subdirectory, and append the suffix. */
		strcpy(stradd(stradd(path, prefdir, DIRSEP), fname, 0), suffix);

	return(path);
}

/* filename - pointer to filename in buffer containing path 
 *
 * Scan the path, recording directory separators.  Return the location
 * of the character past the last one.  If no directory separators are
 * found, returns a pointer to beginning of the path.
 */
char *
filename(
	char	*path
)
{
	char	*cp = path;

	for (; *path; path++)
		if (ISDIRSEP(*path))
			cp = path+1;
	return(cp);
}


/* filetrunc() - return the directory portion of a path
 *
 * The path is passed in in a pointer to a buffer; a null character is
 * inserted in the buffer after the last directory separator
 *
 */
char *
filetrunc(
	char	*path
)
{
	char	*p1, *p2;

	for (p1 = p2 = path; *p2; p2++)
		if (ISDIRSEP(*p2))
			p1 = p2;
	if (p1 == path && ISDIRSEP(*p1))
		p1++;
	*p1 = '\0';
	return(path);
}

/* tailtrunc() - trim a file name extension, if any.
 *
 * The file name is passed in in a buffer indicated by *name; the
 * period which begins the extension is replaced with a 0 byte.
 */
char *
tailtrunc(
	char	*name
)
{
	char	*p1, *p2;

	/* Skip leading periods */
	for (p1 = filename(name); *p1 == '.'; p1++)
		;
	/* Find the last period in a file name */
	p2 = NULL;
	for ( ; *p1; p1++)
		if (*p1 == '.')
			p2 = p1;
	/* If present, trim the filename at that period */
	if (p2 != NULL)
		*p2 = '\0';
	return(name);
}

/* blanktrunc() - trim spaces at the end of a string
 *
 * the string is passed in a character array, which is modified
 */
void
blanktrunc(
	char	*s
)
{
	char	*cp;

	for (cp = s; *cp; cp++)
		;
	while (cp-- > s && isspace(*cp))
		;
	*++cp = '\0';
}

/* fpcomment - output a multi-line comment 
 *
 * The comment may be multiple lines, with each line separated by a
 * newline.  Each line is prefixed by prefix.  If the last line isn't
 * terminated by a newline, no newline will be output.
 */
void
fpcomment(FILE *fp, char *prefix, char *s) {
  while (*s) {			/* While there are characters left to output */
    fprintf(fp, "%s", prefix);	/* Output the prefix */
    for (; *s && *s != '\n'; s++) /* Output a line */
      putc(*s, fp);
    if (*s == '\n') {		/* Including the newline, if any */
      putc(*s, fp);
      s++;
    }
  }
}

/* putheader - output the header of the .rad file
 *
 * Header is:
 *   # <file> <file> <file> (all files from input line)
 *   # Dimensions in [feet,meters,etc.]
 *
 * ??? Is listing all the input file names correct behavior?
 *
 */
void

putheader(
	FILE	*out
)
{
	int	i;

	putc('#', out);
	for (i = 0; i < gargc; i++) {
		putc(' ', out);
		fputs(gargv[i], out);
	}
	fputs("\n# Dimensions in ", out);
	fputs(units, out);
	putc('\n', out);
}

/* ies2rad - convert an IES LM-63 file to a Radiance light source desc.
 *
 * Return -1 in case of failure, 0 in case of success.
 *
 */
int
ies2rad(		/* convert IES file */
	char	*inpname,
	char	*outname
)
{
	SRCINFO	srcinfo;
	char	buf[MAXLINE], tltid[RMAXWORD];
	char	geomfile[MAXLINE];
	FILE	*inpfp, *outfp;
	int	lineno = 0;


	/* Initialize srcinfo */
	srcinfo.filerev = IESFIRSTVER;
	srcinfo.iesshape = IESNONE;
	srcinfo.warn = NULL;
	srcinfo.isillum = 0;
	srcinfo.havelamppos = 0;
        /* Open input and output files */
	geomfile[0] = '\0';
	if (inpname == NULL) {
		inpname = "<stdin>";
		inpfp = stdin;
	} else if ((inpfp = fopen(inpname, "r")) == NULL) {
		perror(inpname);
		return(-1);
	}
	if (out2stdout)
		outfp = stdout;
	else if ((outfp = fopen(fullnam(buf,outname,T_RAD), "w")) == NULL) {
		perror(buf);
		fclose(inpfp);
		return(-1);
	}

	/* Output the output file header */
	putheader(outfp);

	/* If the lamp type wasn't given on the command line, mark
	 * the lamp color as missing */
	if (lamptype == NULL)
		lampcolor = NULL;

	/* Read the input file header, copying lines to the .rad file
	 * and looking for a lamp type. Stop at EOF or a line
	 * beginning with "TILT=". */
	while (fgets(buf,sizeof(buf),inpfp) != NULL
			&& strncmp(buf,TLTSTR,TLTSTRLEN)) {
		blanktrunc(buf); /* Trim trailing spaces, CR, LF. */
		if (!buf[0])	 /* Skip blank lines */
			continue;
		/* increment the header line count. If we are on the
		 * first line of the file, check for a version tag. If
		 * one is not found, assume the first version of the
		 * file. */
		if (!lineno++)
			srcinfo.filerev = iesversion(buf);
		/* Output the header line as a comment in the .rad file. */
		fputs("#<", outfp);
		fputs(buf, outfp);
		putc('\n', outfp);

		/* For post-1986 version files, process a keyword
		 * line.  Otherwise, just scan the line for a lamp
		 * name */
		if (srcinfo.filerev != 1986)
			prockwd(buf, geomfile, inpname, &srcinfo);
		else if (lampcolor == NULL)
			lampcolor = matchlamp(buf);
	}

	/* Done reading header information. If a lamp color still
	 * hasn't been found, print a warning and use the default
	 * color; if a lamp type hasn't been found, but a color has
	 * been specified, used the specified color. */
	if (lampcolor == NULL) {
		fprintf(stderr, "%s: warning - no lamp type\n", inpname);
		fputs("# Unknown lamp type (used default)\n", outfp);
		lampcolor = defcolor;
	} else if (lamptype == NULL)
		fprintf(outfp,"# CIE(x,y) = (%f,%f)\n# Depreciation = %.1f%%\n",
				lampcolor[3], lampcolor[4], 100.*lampcolor[5]);

	/* If the file ended before a "TILT=" line, that's an error. */
	if (feof(inpfp)) {
		fprintf(stderr, "%s: not in IES format\n", inpname);
		goto readerr;
	}

	/* Process the tilt section of the file. */
	/* Get the tilt file name, or the keyword "INCLUDE". */
	atos(tltid, RMAXWORD, buf+TLTSTRLEN);
	if (inpfp == stdin)
		buf[0] = '\0';
	else
		filetrunc(strcpy(buf, inpname));
	/* Process the tilt data. */
	if (dotilt(inpfp, outfp, buf, tltid, outname, tltid) != 0) {
		fprintf(stderr, "%s: bad tilt data\n", inpname);
		goto readerr;
	}

	/* Process the luminaire data. */
	if (dosource(&srcinfo, inpfp, outfp, tltid, outname) != 0) {
		fprintf(stderr, "%s: bad luminaire data\n", inpname);
		goto readerr;
	}

	/* Close the input file */
	fclose(inpfp);

	/* Process an MGF file, if present. cvgeometry() closes outfp. */
	if (cvgeometry(geomfile, &srcinfo, outname, outfp) != 0) {
		fprintf(stderr, "%s: bad geometry file\n", geomfile);
		return(-1);
	}
	return(0);

readerr:
	/* If there is an error reading the file, close the input and
	 * .rad output files, and delete the .rad file, returning -1. */
	fclose(inpfp);
	fclose(outfp);
	unlink(fullnam(buf,outname,T_RAD));
	return(-1);
}

/* dotilt -- process tilt data
 *
 * Generate a brightdata primitive which describes the effect of
 * luminaire tilt on luminaire output and return its identifier in tltid.
 *
 * Tilt data (if present) is given as a number 1, 2, or 3, which
 * specifies the orientation of the lamp within the luminaire, a
 * number, n, of (angle, multiplier) pairs, followed by n angles and n
 * multipliers.
 *
 * returns 0 for success, -1 for error
 */
int
dotilt(
	FILE	*in,
	FILE	*out,
	char	*dir,
	char	*tltspec,
	char	*dfltname,
	char	*tltid
)
{
	int	nangles, tlt_type;
	double	minmax[1][2];
	char	buf[PATH_MAX], tltname[RMAXWORD];
	FILE	*datin, *datout;

	/* Decide where the tilt data is; if the luminaire description
	 * doesn't have a tilt section, set the identifier to "void". */
	if (!strcmp(tltspec, TLTNONE)) {
		/* If the line is "TILT=NONE", set the input file
		 * pointer to NULL and the identifier to "void". */
		datin = NULL;
		strcpy(tltid, "void");
	} else if (!strcmp(tltspec, TLTINCL)) {
		/* If the line is "TILT=INCLUDE" use the main IES
		 * file as the source of tilt data. */
		datin = in;
		strcpy(tltname, dfltname);
	} else {
		/* If the line is "TILT=<filename>", use that file
		 * name as the source of tilt data. */
		if (ISDIRSEP(tltspec[0]))
			strcpy(buf, tltspec);
		else
			strcpy(stradd(buf, dir, DIRSEP), tltspec);
		if ((datin = fopen(buf, "r")) == NULL) {
			perror(buf);
			return(-1);
		}
		tailtrunc(strcpy(tltname,filename(tltspec)));
	}
	/* If tilt data is present, read, process, and output it. */
	if (datin != NULL) {
		/* Try to open the output file */
		if ((datout = fopen(fullnam(buf,tltname,T_TLT),"w")) == NULL) {
			perror(buf);
			if (datin != in)
				fclose(datin);
			return(-1);
		}
		/* Try to copy the tilt data to the tilt data file */
		if (!scnint(datin,&tlt_type) || !scnint(datin,&nangles)
			|| cvdata(datin,datout,1,&nangles,1.,minmax) != 0) {
			fprintf(stderr, "%s: data format error\n", tltspec);
			fclose(datout);
			if (datin != in)
				fclose(datin);
			unlink(fullnam(buf,tltname,T_TLT));
			return(-1);
		}
		fclose(datout);
		if (datin != in)
			fclose(datin);

		/* Generate the identifier of the brightdata; the filename
		 * with "_tilt" appended. */
		strcat(strcpy(tltid, filename(tltname)), "_tilt");
		/* Write out the brightdata primitive */
		fprintf(out, "\nvoid brightdata %s\n", tltid);
		libname(buf,tltname,T_TLT);
		/* Generate the tilt description */
		switch (tlt_type) {
		case TLT_VERT:
			/* The lamp is mounted vertically; either
			 * base up or base down. */
			fprintf(out, "4 noop %s tilt.cal %s\n", buf,
				minmax[0][1]>90.+FTINY ? "tilt_ang" : "tilt_ang2");
			break;
		case TLT_H0:
			/* The lamp is mounted horizontally and
			 * rotates but does not tilt when the
			 * luminaire is tilted. */
			fprintf(out, "6 noop %s tilt.cal %s -rz 90\n", buf,
			minmax[0][1]>90.+FTINY ? "tilt_xang" : "tilt_xang2");
			break;
		case TLT_H90:
			/* The lamp is mounted horizontally, and
			 * tilts when the luminaire is tilted. */
			fprintf(out, "4 noop %s tilt.cal %s\n", buf,
			minmax[0][1]>90.+FTINY ? "tilt_xang" : "tilt_xang2");
			break;
		default:
			/* otherwise, this is a bad IES file */
			fprintf(stderr,
				"%s: illegal lamp to luminaire geometry (%d)\n",
				tltspec, tlt_type);
			return(-1);
		}
		/* And finally output the numbers of integer and real
		 * arguments, of which there are none. */
		fprintf(out, "0\n0\n");
	}
	return(0);
}

/* dosource -- create the source and distribution primitives */
int
dosource(
	SRCINFO	*sinf,
	FILE	*in,
	FILE	*out,
	char	*mod,
	char	*name
)
{
	char	buf[PATH_MAX], id[RMAXWORD];
	FILE	*datout;
	double	mult, bfactor, pfactor, width, length, height, wattage;
	double	bounds[2][2];
	int	nangles[2], pmtype, unitype;
	double	d1;
	int	doupper, dolower, dosides;

	/* Read in the luminaire description header */
	if (!isint(getword(in)) || !isflt(getword(in)) || !scnflt(in,&mult)
			|| !scnint(in,&nangles[0]) || !scnint(in,&nangles[1])
			|| !scnint(in,&pmtype) || !scnint(in,&unitype)
			|| !scnflt(in,&width) || !scnflt(in,&length)
			|| !scnflt(in,&height) || !scnflt(in,&bfactor)
			|| !scnflt(in,&pfactor) || !scnflt(in,&wattage)) {
		fprintf(stderr, "dosource: bad lamp specification\n");
		return(-1);
	}
	
	/* pfactor is only provided in 1986 and 1991 format files, and
	 * is something completely different in 2019 files.  If the
	 * file version is 1995 or later, set it to 1.0 to avoid
	 * error. */
	if (sinf->filerev >= 1995)
		pfactor = 1.0;

	/* Type A photometry is not supported */
	if (pmtype != PM_C && pmtype != PM_B) {
		fprintf(stderr, "dosource: unsupported photometric type (%d)\n",
				pmtype);
		return(-1);
	}

	/* Multiplier = the multiplier from the -m option, times the
	 * multiplier from the IES file, times the ballast factor,
	 * times the "ballast lamp photometric factor," (pfactor)
	 * which was part of the 1986 and 1991 standards. In the 1995
	 * and 2002 standards, it is always supposed to be 1 and in
	 * the 2019 standard it encodes information about the source
	 * of the file.  For those files, pfactor is set to 1.0,
	 * above.  */
	sinf->mult = multiplier*mult*bfactor*pfactor;

	/* If the count of angles is wrong, raise an error and quit. */
	if (nangles[0] < 2 || nangles[1] < 1) {
		fprintf(stderr, "dosource: too few measured angles\n");
		return(-1);
	}

	/* For internal computation, convert units to meters. */
	if (unitype == U_FEET) {
		width *= F_M;
		length *= F_M;
		height *= F_M;
	}

	/* Make decisions about the shape of the light source
	 * geometry, and store them in sinf. */
	if (makeshape(sinf, width, length, height) != 0) {
		fprintf(stderr, "dosource: illegal source dimensions\n");
		return(-1);
	}
	/* If any warning messages were generated by makeshape(), output them */
	if ((sinf->warn) != NULL)
		fputs(sinf->warn, stderr);

	/* Copy the candela values into a Radiance data file. */
	if ((datout = fopen(fullnam(buf,name,T_DST), "w")) == NULL) {
		perror(buf);
		return(-1);
	}
	if (cvdata(in, datout, 2, nangles, 1./WHTEFFICACY, bounds) != 0) {
		fprintf(stderr, "dosource: bad distribution data\n");
		fclose(datout);
		unlink(fullnam(buf,name,T_DST));
		return(-1);
	}
	fclose(datout);

	/* Output explanatory comment */
	fprintf(out, "\n# %g watt luminaire, lamp*ballast factor = %g\n",
			wattage, bfactor*pfactor);
	if (sinf->iesshape >= 0)
		fprintf(out, "# IES file shape = %s\n",
			IESHAPENAMES[sinf->iesshape]);
	else
		fprintf(out, "# IES file shape overridden\n");
	fprintf(out, "# Radiance geometry shape = %s\n",
		RADSHAPENAMES[sinf->type - 1]);
	if (sinf->warn != NULL)
		fpcomment(out, "# ", sinf->warn);

	/* Output distribution "brightdata" primitive. Start handling
	   the various cases of symmetry of the distribution.  This
	   code reflects the complexity of the LM-63 format, as
	   described under "<horizontal angles>" in the various
	   versions of the standard. */
	strcat(strcpy(id, filename(name)), "_dist");
	fprintf(out, "\n'%s' brightdata '%s'\n", mod, id);
	if (nangles[1] < 2)
		 /* if it's a radially-symmetric type C distribution */
		fprintf(out, "4 ");
	else if (pmtype == PM_B)
		/* Photometry type B */
		fprintf(out, "5 ");
	else if (FEQ(bounds[1][0],90.) && FEQ(bounds[1][1],270.))
		/* Symmetric around the 90-270 degree plane */
		fprintf(out, "7 ");
	else
		/* Just regular type C photometry */
		fprintf(out, "5 ");

	/* If the generated source geometry will be a box, a flat
	 * rectangle, or a disk figure out if it needs a top, a
	 * bottom, and/or sides. */
	dolower = (bounds[0][0] < 90.-FTINY); /* Smallest vertical angle */
	doupper = (bounds[0][1] > 90.+FTINY); /* Largest vertical angle */
	dosides = (doupper & dolower && sinf->h > MINDIM); /* Sides */

	/* Select the appropriate function and parameters from source.cal */
	fprintf(out, "%s '%s' source.cal ",
			sinf->type==SPHERE ? "corr" :
			!dosides ? "flatcorr" :
			sinf->type==DISK ? "cylcorr" : "boxcorr",
			libname(buf,name,T_DST));
	if (pmtype == PM_B) {
		/* Type B photometry */
		if (FEQ(bounds[1][0],0.))
			/* laterally symmetric around a vertical plane */
			fprintf(out, "srcB_horiz2 ");
		else
			fprintf(out, "srcB_horiz ");
		fprintf(out, "srcB_vert ");
	} else /* pmtype == PM_C */ {
		if (nangles[1] >= 2) {
			/* Not radially symmetric */
			d1 = bounds[1][1] - bounds[1][0];
			if (d1 <= 90.+FTINY)
				/* Data for a quadrant */
				fprintf(out, "src_phi4 ");
			else if (d1 <= 180.+FTINY) {
				/* Data for a hemisphere */
				if (FEQ(bounds[1][0],90.))
					fprintf(out, "src_phi2+90 ");
				else
					fprintf(out, "src_phi2 ");
			} else	/* Data for a whole sphere */
				fprintf(out, "src_phi ");
			fprintf(out, "src_theta ");
			/* For the hemisphere around the 90-270 degree plane */
			if (FEQ(bounds[1][0],90.) && FEQ(bounds[1][1],270.))
				fprintf(out, "-rz -90 ");
		} else		/* Radially symmetric */
			fprintf(out, "src_theta ");
	}
	/* finish the brightdata primitive with appropriate data */
	if (!dosides || sinf->type == SPHERE)
		fprintf(out, "\n0\n1 %g\n", sinf->mult/sinf->area);
	else if (sinf->type == DISK)
		fprintf(out, "\n0\n3 %g %g %g\n", sinf->mult,
				sinf->w, sinf->h);
	else
		fprintf(out, "\n0\n4 %g %g %g %g\n", sinf->mult,
				sinf->l, sinf->w, sinf->h);
	/* Brightdata primitive written out. */

	/* Finally, output the descriptions of the actual radiant
	 * surfaces. */
	if (putsource(sinf, out, id, filename(name),
			dolower, doupper, dosides) != 0)
		return(-1);
	return(0);
}

/* putsource - output the actual light emitting geometry
 *
 * Three kinds of geometry are produced: rectangles and boxes, disks
 * ("ring" primitive, but the radius of the hole is always zero) and
 * cylinders, and spheres.
 */
int
putsource(
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name,
	int	dolower,
	int	doupper,
	int	dosides
)
{
	char	lname[RMAXWORD];

	/* First, describe the light. If a materials and geometry
	 * file is given, generate an illum instead. */
	strcat(strcpy(lname, name), "_light");
	fprintf(fp, "\n'%s' %s '%s'\n", mod,
			shp->isillum ? "illum" : "light", lname);
	fprintf(fp, "0\n0\n3 %g %g %g\n",
			lampcolor[0], lampcolor[1], lampcolor[2]);
	switch (shp->type) {
	case RECT:
		/* Output at least one rectangle. If light is radiated
		 * from the sides of the luminaire, output rectangular
		 * sides as well. */
		if (dolower)
			putrectsrc(shp, fp, lname, name, 0);
		if (doupper)
			putrectsrc(shp, fp, lname, name, 1);
		if (dosides)
			putsides(shp, fp, lname, name);
		break;
	case DISK:
		/* Output at least one disk. If light is radiated from
		 * the sides of luminaire, output a cylinder as well. */
		if (dolower)
			putdisksrc(shp, fp, lname, name, 0);
		if (doupper)
			putdisksrc(shp, fp, lname, name, 1);
		if (dosides)
			putcyl(shp, fp, lname, name);
		break;
	case SPHERE:
		/* Output a sphere. */
		putspheresrc(shp, fp, lname, name);
		break;
	}
	return(0);
}

/* makeshape -- decide what shape will be used
 *
 * Makeshape decides what Radiance geometry will be used to represent
 * the light source and stores information about it in shp.
 * 
 * The height, width, and length parameters are values from the
 * IES file, given in meters.
 *
 * The various versions of the IES LM-63 standard give a "luminous
 * opening" (really a crude shape) a width, a length (or depth), and a
 * height.  If all three values are positive, they describe a box.  If
 * they are all zero, they describe a point.  Various combinations of
 * negative values are used to denote disks, circular or elliptical
 * cylinders, spheres, and ellipsoids.  This encoding differs from
 * version to version of LM-63.
 *
 * Ies2rad simplifies this, reducing the geometry of LM-63 files to
 * three forms which can be easily represented by Radiance primitives:
 * boxes (RECT), cylinders or disks (DISK), and spheres (SPHERE.)  A
 * point is necessarily represented by a small sphere, since a point
 * is not a Radiance object.
 *
 * Makeshape() returns 0 if it succeeds in choosing a shape, and -1 if
 * it fails.
 *
 */
int
makeshape(
	SRCINFO	*shp,
	double	width,
	double	length,
	double	height
)
{
	int rc;
	
	if (illumrad != 0.0)
		rc = makeillumsphere(shp);
	else
		rc = makeiesshape(shp, length, width, height);
	if (rc == SUCCESS)
		shapearea(shp);
	return rc;
}

/* 
 * Return 1 if d < 0, 2 if d == 0, 3 if d > 0.  This is used to encode
 * the signs of IES file dimensions for quick lookup.  As usual with
 * macros, don't use an expression with side effects as an argument.
 */
#define CONVSGN(d) ((d) < 0 ? 1 : ((d) == 0 ? 2 : 3))

/* makeiesshape - convert IES shape to Radiance shape
 * 
 * Some 34 cases in the various versions of the IES LM-63 standard are
 * handled, though some only by approximation.  For each case which is
 * processed a Radiance box, cylinder, or sphere is selected.
 *
 * Shapes are categorized by version year of the standard and the
 * signs of the LM-63 length, width (depth), and height fields.  These
 * are combined and converted to an integer, which is then used as the
 * argument to switch().  The last two digits of the IES file version
 * year are used and the signs of length, width, and height are
 * encoded, in that order, as 1 for negative, 2 for zero, and 3 for
 * positive.  These are then combined into a numeric key by the
 * following formula: 
 *
 *   version * 1000 + sgn(length) * 100 + sgn(width) * 10 + sgn(height).
 *
 * Since the 1991 version uses the same encoding as the 1986 version,
 * and the 2019 version uses the same encoding as the 2002 version,
 * these are collapsed into the earlier years.
 *
 * In the cases of the switch() statement, further processing takes
 * place. Circles and ellipses are distinguished by comparisons.  Then
 * routines are called to fill out the fields of the shp structure.
 *
 * As per the conventions of the rest of ies2rad, makeiesshape()
 * returns 0 on success and -1 on failure.  -1 reflects an error in
 * the IES file and is unusual.
 *
 * By convention, the shape generating routines are always given
 * positive values for dimensions and always succeed; all errors are
 * caught before they are called.  The absolute values of all three
 * dimensions are calculated at the beginning of makeiesshape() and
 * used throughout the function, this has a low cost and eliminates
 * the chance of sign errors.
 * 
 * There is one extension to the ies standard here, devised to
 * accomdate wall-mounted fixtures; vertical rectangles, not formally
 * supported by any version of LM-63, are treated as boxes.
 *
 * The code is complicated by the way that earlier versions of the
 * standard (1986 and 1991) prioritize width in their discussions, and
 * later versions prioritize length.  It is not always clear which to
 * write first and there is hesitation between the older code which
 * invokes makeiesshape() and makeiesshape() itself.
 */
int
makeiesshape(SRCINFO *shp, double l, double w, double h) {
	int rc = SUCCESS;
	int shape = IESNONE;
        /* Get the last two digits of the standard year  */
	int ver = shp->filerev % 100;
	/* Make positive versions of all dimensions, for clarity in
	 * function calls.  If you like, read this as l', w', and h'. */
	double lp = fabs(l), wp = fabs(w), hp = fabs(h);
	int thumbprint;

	/* Change 1991 into 1986 and 2019 in 2002 */
	switch (ver) {
	case 91:
		ver = 86;
		break;
	case 19:
		ver = 02;
		break;
	}
	
	thumbprint =
		ver * 1000 + CONVSGN(l) * 100 + CONVSGN(w) * 10 + CONVSGN(h);
	switch(thumbprint) {
	case 86222: case 95222: case 2222:
		shp->iesshape = IESPT;
		shp->type = SPHERE;
		shp->w = shp->l = shp->h = MINDIM;
		break;
	case 86332: case 95332: case 2332:
		shp->iesshape = IESRECT;
		makeboxshape(shp, lp, wp, hp);
		break;
	case 86333: case 86233: case 86323:
	case 95333: case 95233: case 95323:
	case 2333: case 2233: case 2323:
		shp->iesshape = IESBOX;
		makeboxshape(shp, lp, wp, hp);
		break;
	case 86212: case 95212:
		shp->iesshape = IESDISK;
		makecylshape(shp, wp, hp);
		break;
	case 86213:
		shp->iesshape = IESVCYL;
		makecylshape(shp, wp, hp);
		break;
	case 86312:
		shp->iesshape = IESELLIPSE;
		makeecylshape(shp, lp, wp, 0);
		break;
	case 86313:
		shp->iesshape = IESELLIPSOID;
		makeelshape(shp, wp, lp, hp);
		break;
	case 95211:
		shp->iesshape = FEQ(lp,hp) ? IESSPHERE : IESNONE;
		if (shp->iesshape == IESNONE) {
			shp->warn = "makeshape: cannot determine shape\n";
			rc = FAIL;
			break;
		}
		shp->type = SPHERE;
		shp->w = shp->l = shp->h = wp;
		break;
	case 95213:
		shp->iesshape = IESVCYL;
		makecylshape(shp, wp, hp);
		break;
	case 95321:
		shp->iesshape = IESHCYL_PH;
		shp->warn = "makeshape: shape is a horizontal cylinder, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, lp, wp, hp);
		break;
	case 95231:
		shp->iesshape = IESHCYL_PPH;
		shp->warn = "makeshape: shape is a horizontal cylinder, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, lp, wp, hp);
		break;
	case 95133: case 95313:
		shp->iesshape = IESVECYL;
		makeecylshape(shp, lp, wp, hp);
		break;
	case 95131: case 95311:
		shp->iesshape = IESELLIPSOID;
		makeelshape(shp, lp, wp, hp);
		break;
	case 2112:
		shp->iesshape = FEQ(l,w) ? IESDISK : IESELLIPSE;
		if (shp->iesshape == IESDISK)
			makecylshape(shp, wp, hp);
		else
			makeecylshape(shp, wp, lp, hp);
		break;
	case 2113:
		shp->iesshape = FEQ(l,w) ? IESVCYL : IESVECYL;
		if (shp->iesshape == IESVCYL)
			makecylshape(shp, wp, hp);
		else
			makeecylshape(shp, wp, lp, hp);
		break;
	case 2111:
		shp->iesshape = FEQ(l,w) && FEQ(l,h) ? IESSPHERE : IESELLIPSOID;
		if (shp->iesshape == IESSPHERE) {
			shp->type = SPHERE;
			shp->w = shp->l = shp->h = wp;
		}
		else
			makeelshape(shp, lp, wp, hp);
		break;
	case 2311:
		shp->iesshape = FEQ(w,h) ? IESHCYL_PH : IESHECYL_PH;
		shp->warn = "makeshape: shape is a horizontal cylinder, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, lp, wp, hp);
		break;
	case 2131:
		shp->iesshape = FEQ(l,h) ? IESHCYL_PPH : IESHECYL_PPH;
		shp->warn = "makeshape: shape is a horizontal cylinder, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, lp, wp, hp);
		break;
	case 2121:
		shp->iesshape = FEQ(w,h) ? IESVDISK_PH : IESVEL_PH;
		shp->warn = "makeshape: shape is a vertical ellipse, which is not supported.\nmakeshape: replaced with rectangle\n";
		makeboxshape(shp, lp, wp, hp);
		break;
	default:
		/* We don't recognize the shape - report an error. */
		rc = FAIL;
	}
	return rc;
}

/* makeillumsphere - create an illum sphere */
int
makeillumsphere(SRCINFO *shp) {
	/* If the size is too small or negative, error. */
	if (illumrad/meters2out < MINDIM/2.) {
		fprintf(stderr, "makeillumsphere: -i argument is too small or negative\n");
		return FAIL;
	}
	shp->isillum = 1;
	shp->type = SPHERE;
	shp->w = shp->l = shp->h = 2.*illumrad / meters2out;
	return SUCCESS;
}

/* makeboxshape - create a box */
void
makeboxshape(SRCINFO *shp, double l, double w, double h) {
	shp->type = RECT;
	shp->l = fmax(l, MINDIM);
	shp->w = fmax(w, MINDIM);
	shp->h = fmax(h, .5*MINDIM);
}

/* makecylshape - output a vertical cylinder or disk 
 *
 * If the shape has no height, make it a half-millimeter.
 */
void
makecylshape(SRCINFO *shp, double diam, double height) {
	shp->type = DISK;
	shp->w = shp->l = diam;
	shp->h = fmax(height, .5*MINDIM);
}

/* makeelshape - create a substitute for an ellipsoid
 * 
 * Because we don't actually support ellipsoids, and they don't seem
 * to be common in actual IES files.
 */
void
makeelshape(SRCINFO *shp, double w, double l, double h) {
	float avg = (w + l + h) / 3;
	float bot = .5 * avg;
	float top = 1.5 * avg;

	if (bot < w && w < top
	    && bot < l && l < top
	    && bot < h && h > top) {
		/* it's sort of spherical, replace it with a sphere */
		shp->warn = "makeshape: shape is an ellipsoid, which is not supported.\nmakeshape: replaced with sphere\n";
		shp->type = SPHERE;
		shp->w = shp->l = shp->h = avg;
	} else if (bot < w && w < top
		   && bot < l && l < top
		   && h <= .5*MINDIM) {
		/* It's flat and sort of circular, replace it
		 * with a disk. */
		shp->warn = "makeshape: shape is an ellipse, which is not supported.\nmakeshape: replaced with disk\n";
		makecylshape(shp, w, 0);
	} else {
		shp->warn = "makeshape: shape is an ellipsoid, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, w, l, h);
	}
}

/* makeecylshape - create a substitute for an elliptical cylinder or disk */
void
makeecylshape(SRCINFO *shp, double l, double w, double h) {
	float avg = (w + l) / 2;
	float bot = .5 * avg;
	float top = 1.5 * avg;

	if (bot < w && w < top
		   && bot < l && l < top) {
		/* It's sort of circular, replace it
		 * with a circular cylinder. */
		shp->warn = "makeshape: shape is a vertical elliptical cylinder, which is not supported.\nmakeshape: replaced with circular cylinder\n";
		makecylshape(shp, w, h);
	} else {
		shp->warn = "makeshape: shape is a  vertical elliptical cylinder, which is not supported.\nmakeshape: replaced with box\n";
		makeboxshape(shp, w, l, h);
	}
}

void
shapearea(SRCINFO *shp) {
	switch (shp->type) {
	case RECT:
		shp->area = shp->w * shp->l;
		break;
	case DISK:
	case SPHERE:
		shp->area = PI/4. * shp->w * shp->w;
		break;
	}
}	

/* Rectangular or box-shaped light source.
 *
 * putrectsrc, putsides, putrect, and putpoint are used to output the
 * Radiance description of a box.  The box is centered on the origin
 * and has the dimensions given in the IES file.  The coordinates
 * range from [-1/2*length, -1/2*width, -1/2*height] to [1/2*length,
 * 1/2*width, 1/2*height].
 *
 * The location of the point is encoded in the low-order three bits of
 * an integer. If the integer is p, then: bit 0 is (p & 1),
 * representing length (x), bit 1 is (p & 2) representing width (y),
 * and bit 2 is (p & 4), representing height (z).
 *
 * Looking down from above (towards -z), the vertices of the box or
 * rectangle are numbered so:
 *
 *     2,6					  3,7
 *        +--------------------------------------+
 *        |					 |
 *        |					 |
 *        |					 |
 *        |					 |
 *        +--------------------------------------+
 *     0,4					  1,5
 *
 * The higher number of each pair is above the x-y plane (positive z),
 * the lower number is below the x-y plane (negative z.)
 *
 */

/* putrecsrc - output a rectangle parallel to the x-y plane 
 *
 * Putrecsrc calls out the vertices of a rectangle parallel to the x-y
 * plane.  The order of the vertices is different for the upper and
 * lower rectangles of a box, since a right-hand rule based on the
 * order of the vertices is used to determine the surface normal of
 * the rectangle, and the surface normal determines the direction the
 * light radiated by the rectangle.
 *
 */
void
putrectsrc(
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name,
	int	up
)
{
	if (up)
		putrect(shp, fp, mod, name, ".u", 4, 5, 7, 6);
	else
		putrect(shp, fp, mod, name, ".d", 0, 2, 3, 1);
}

/* putsides - put out sides of box */
void
putsides(
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name
)
{
	putrect(shp, fp, mod, name, ".1", 0, 1, 5, 4);
	putrect(shp, fp, mod, name, ".2", 1, 3, 7, 5);
	putrect(shp, fp, mod, name, ".3", 3, 2, 6, 7);
	putrect(shp, fp, mod, name, ".4", 2, 0, 4, 6);
}

/* putrect - put out a rectangle 
 *
 * putrect generates the "polygon" primitive which describes a
 * rectangle.
 */
void
putrect(
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name,
	char	*suffix,
	int	a,
	int b,
	int c,
	int d
)
{
	fprintf(fp, "\n'%s' polygon '%s%s'\n0\n0\n12\n", mod, name, suffix);
	putpoint(shp, fp, a);
	putpoint(shp, fp, b);
	putpoint(shp, fp, c);
	putpoint(shp, fp, d);
}

/* putpoint -- output a the coordinates of a vertex
 *
 * putpoint maps vertex numbers to coordinates and outputs the
 * coordinates.
 */
void
putpoint(
	SRCINFO	*shp,
	FILE	*fp,
	int	p
)
{
	static double	mult[2] = {-.5, .5};

	fprintf(fp, "\t%g\t%g\t%g\n",
			mult[p&1]*shp->l*meters2out,
			mult[p>>1&1]*shp->w*meters2out,
			mult[p>>2]*shp->h*meters2out);
}

/* End of routines to output a box-shaped light source */

/* Routines to output a cylindrical or disk shaped light source 
 * 
 * As with other shapes, the light source is centered on the origin.
 * The "ring" and "cylinder" primitives are used.
 *
 */
void
putdisksrc(		/* put out a disk source */
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name,
	int	up
)
{
	if (up) {
		fprintf(fp, "\n'%s' ring '%s.u'\n", mod, name);
		fprintf(fp, "0\n0\n8\n");
		fprintf(fp, "\t0 0 %g\n", .5*shp->h*meters2out);
		fprintf(fp, "\t0 0 1\n");
		fprintf(fp, "\t0 %g\n", .5*shp->w*meters2out);
	} else {
		fprintf(fp, "\n'%s' ring '%s.d'\n", mod, name);
		fprintf(fp, "0\n0\n8\n");
		fprintf(fp, "\t0 0 %g\n", -.5*shp->h*meters2out);
		fprintf(fp, "\t0 0 -1\n");
		fprintf(fp, "\t0 %g\n", .5*shp->w*meters2out);
	}
}


void
putcyl(			/* put out a cylinder */
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name
)
{
	fprintf(fp, "\n'%s' cylinder '%s.c'\n", mod, name);
	fprintf(fp, "0\n0\n7\n");
	fprintf(fp, "\t0 0 %g\n", .5*shp->h*meters2out);
	fprintf(fp, "\t0 0 %g\n", -.5*shp->h*meters2out);
	fprintf(fp, "\t%g\n", .5*shp->w*meters2out);
}

/* end of of routines to output cylinders and disks */

void
putspheresrc(		/* put out a sphere source */
	SRCINFO	*shp,
	FILE	*fp,
	char	*mod,
	char	*name
)
{
	fprintf(fp, "\n'%s' sphere '%s.s'\n", mod, name);
	fprintf(fp, "0\n0\n4 0 0 0 %g\n", .5*shp->w*meters2out);
}

/* cvdata - convert LM-63 tilt and candela data to Radiance brightdata format
 *
 * The files created by this routine are intended for use with the Radiance
 * "brightdata" material type.
 *
 * Two types of data are converted; one-dimensional tilt data, which
 * is given in polar coordinates, and two-dimensional candela data,
 * which is given in spherical co-ordinates.
 *
 * Return 0 for success, -1 for failure.
 *
 */
int
cvdata(
	FILE	*in,		/* Input file */
	FILE	*out,		/* Output file */
	int	ndim,		/* Number of dimensions; 1 for
				 * tilt data, 2 for photometric data. */
	int	npts[],		/* Number of points in each dimension */
	double	mult,		/* Multiple each value by this
				 * number. For tilt data, always
				 * 1. For candela values, the
				 * efficacy of white Radiance light.  */
	double	lim[][2]	/* The range of angles in each dimension. */
)
{
	double	*pt[4];		/* Four is the expected maximum of ndim. */
	int	i, j;
	double	val;
	int	total;

	/* Calculate and output the number of data values */
	total = 1; j = 0;
	for (i = 0; i < ndim; i++)
		if (npts[i] > 1) {
			total *= npts[i];
			j++;
		}
	fprintf(out, "%d\n", j);

	/* Read in the angle values, and note the first and last in
	 * each dimension, if there is a place to store them. In the
	 * case of tilt data, there is only one list of angles. In the
	 * case of candela values, vertical angles appear first, and
	 * horizontal angles occur second. */
	for (i = 0; i < ndim; i++) {
		/* Allocate space for the angle values. */
		pt[i] = (double *)malloc(npts[i]*sizeof(double));
		for (j = 0; j < npts[i]; j++)
			if (!scnflt(in, &pt[i][j]))
				return(-1);
		if (lim != NULL) {
			lim[i][0] = pt[i][0];
			lim[i][1] = pt[i][npts[i]-1];
		}
	}

	/* Output the angles. If this is candela data, horizontal
	 * angles output first. There are two cases: the first where
	 * the angles are evenly spaced, the second where they are
	 * not.
	 *
	 * When the angles are evenly spaced, three numbers are
	 * output: the first angle, the last angle, and the number of
	 * angles.  When the angles are not evenly spaced, instead
	 * zero, zero, and the count of angles is given, followed by a
	 * list of angles.  In this case, angles are output four to a line.
	 */
	for (i = ndim-1; i >= 0; i--) {
		if (npts[i] > 1) {
			/* Determine if the angles are evenly spaces */
			for (j = 1; j < npts[i]-1; j++)
				if (!FEQ(pt[i][j]-pt[i][j-1],
						pt[i][j+1]-pt[i][j]))
					break;
			/* If they are, output the first angle, the
			 * last angle, and a count */
			if (j == npts[i]-1)
				fprintf(out, "%g %g %d\n", pt[i][0], pt[i][j],
						npts[i]);
			else {
				/* otherwise, output 0, 0, and a
				 * count, followed by the list of
				 * angles, one to a line. */
				fprintf(out, "0 0 %d", npts[i]);
				for (j = 0; j < npts[i]; j++) {
					if (j%4 == 0)
						putc('\n', out);
					fprintf(out, "\t%g", pt[i][j]);
				}
				putc('\n', out);
			}
		}
		/* Free the storage containing the angle values. */
		free((void *)pt[i]);
	}

	/* Finally, read in the data values (candela or multiplier values,
	 * depending on the part of the file) and output them four to
	 * a line. */
	for (i = 0; i < total; i++) {
		if (i%4 == 0)
			putc('\n', out);
		if (!scnflt(in, &val))
			return(-1);
		fprintf(out, "\t%g", val*mult);
	}
	putc('\n', out);
	return(0);
}

/* getword - get an LM-63 delimited word from fp
 *
 * Getword gets a word from an IES file delimited by either white
 * space or a comma surrounded by white space. A pointer to the word
 * is returned, which will persist only until getword is called again.
 * At EOF, return NULL instead.
 *
 */
char *
getword(			/* scan a word from fp */
	FILE	*fp
)
{
	static char	wrd[RMAXWORD];
	char	*cp;
	int	c;

	/* Skip initial spaces */
	while (isspace(c=getc(fp)))
		;
	/* Get characters to a delimiter or until wrd is full */
	for (cp = wrd; c != EOF && cp < wrd+RMAXWORD-1;
			*cp++ = c, c = getc(fp))
		if (isspace(c) || c == ',') {
			/* If we find a delimiter */
			/* Gobble up whitespace */
			while (isspace(c))
				c = getc(fp);
			/* If it's not a comma, put the first
			 * character of the next data item back */
			if ((c != EOF) & (c != ','))
				ungetc(c, fp);
			/* Close out the strimg */
			*cp = '\0';
			/* return it */
			return(wrd);
		}
	/* If we ran out of space or are at the end of the file,
	 * return either the word or NULL, as appropriate. */
	*cp = '\0';
	return(cp > wrd ? wrd : NULL);
}

/* cvtint - convert an IES word to an integer
 *
 * A pointer to the word is passed in wrd; ip is expected to point to
 * an integer.  cvtint() will silently truncate a floating point value
 * to an integer; "1", "1.0", and "1.5" will all return 1.
 *
 * cvtint() returns 0 if it fails, 1 if it succeeds.
 */
int
cvtint(
	int	*ip,
	char	*wrd
)
{
	if (wrd == NULL || !isint(wrd))
		return(0);
	*ip = atoi(wrd);
	return(1);
}


/* cvtflt - convert an IES word to a double precision floating-point number
 *
 * A pointer to the word is passed in wrd; rp is expected to point to
 * a double.
 *
 * cvtflt returns 0 if it fails, 1 if it succeeds.
 */
int
cvtflt(
	double	*rp,
	char	*wrd
)
{
	if (wrd == NULL || !isflt(wrd))
		return(0);
	*rp = atof(wrd);
	return(1);
}

/* cvgeometry - process materials and geometry format luminaire data
 *
 * The materials and geometry format (MGF) for describing luminaires
 * was a part of Radiance that was first adopted and then retracted by
 * the IES as part of LM-63.  It provides a way of describing
 * luminaire geometry similar to the Radiance scene description
 * format.
 *
 * cvgeometry() generates an mgf2rad command and then, if "-g" is given
 * on the command line, an oconv command, both of which are then
 * executed with the system() function.
 *
 * The generated commands are:
 *   mgf2rad -e <multiplier> -g <size> <mgf_filename> \
 *     | xform -s <scale_factor> \
 *     >> <luminare_scene_description_file
 * or:
 *   mgf2rad -e <multiplier> -g <size> <mgf_filename> \
 *     oconv - > <instance_filename>
 */
int
cvgeometry(
	char	*inpname,
	SRCINFO	*sinf,
	char	*outname,
	FILE	*outfp			/* close output file upon return */
)
{
	char	buf[256];
	char	*cp;

	if (inpname == NULL || !inpname[0]) {	/* no geometry file */
		fclose(outfp);
		return(0);
	}
	putc('\n', outfp);
	strcpy(buf, "mgf2rad ");		/* build mgf2rad command */
	cp = buf+8;
	if (!FEQ(sinf->mult, 1.0)) {
                /* if there's an output multiplier, include in the
		 * mgf2rad command */
		sprintf(cp, "-e %f ", sinf->mult);
		cp += strlen(cp);
	}
	/* Include the glow distance for the geometry */
	sprintf(cp, "-g %f %s ",
		sqrt(sinf->w*sinf->w + sinf->h*sinf->h + sinf->l*sinf->l),
			inpname);
	cp += strlen(cp);
	if (instantiate) {		/* instantiate octree */
		/* If "-g" is given on the command line, include an
		 * "oconv" command in the pipe. */
		strcpy(cp, "| oconv - > ");
		cp += 12;
		fullnam(cp,outname,T_OCT);
		/* Only update if the input file is newer than the
		 * output file */
		if (fdate(inpname) > fdate(outname) &&
				system(buf)) {		/* create octree */
			fclose(outfp);
			return(-1);
		}
		/* Reference the instance file in the scene description */
		fprintf(outfp, "void instance %s_inst\n", outname);
		/* If the geometry isn't in meters, scale it appropriately. */
		if (!FEQ(meters2out, 1.0))
			fprintf(outfp, "3 %s -s %f\n",
					libname(buf,outname,T_OCT),
					meters2out);
		else
			fprintf(outfp, "1 %s\n", libname(buf,outname,T_OCT));
		/* Close off the "instance" primitive. */
		fprintf(outfp, "0\n0\n");
		/* And the Radiance scene description. */
		fclose(outfp);
	} else {			/* else append to luminaire file */
		if (!FEQ(meters2out, 1.0)) {	/* apply scalefactor */
			sprintf(cp, "| xform -s %f ", meters2out);
			cp += strlen(cp);
		}
		if (!out2stdout) {
			fclose(outfp);
			strcpy(cp, ">> ");	/* append works for DOS? */
			cp += 3;
			fullnam(cp,outname,T_RAD);
		}
		if (system(buf))
			return(-1);
	}
	return(0);
}

/* Set up emacs indentation */
/* Local Variables: */
/*   c-file-style: "bsd" */
/* End: */

/* For vim, use ":set tabstop=8 shiftwidth=8" */
