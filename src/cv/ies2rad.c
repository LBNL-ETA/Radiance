#ifndef lint
static const char	RCSid[] = "$Id: ies2rad.c,v 2.21 2003/07/03 22:41:44 schorsch Exp $";
#endif
/*
 * Convert IES luminaire data to Radiance description
 *
 *	07Apr90		Greg Ward
 *
 *  Fixed correction factor for flat sources 29Oct2001 GW
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <ctype.h>
#include "color.h"
#include "paths.h"

#define PI		3.14159265358979323846
					/* floating comparisons */
#define FTINY		1e-6
#define FEQ(a,b)	((a)<=(b)+FTINY&&(a)>=(b)-FTINY)
					/* keywords */
#define MAGICID		"IESNA"
#define LMAGICID	5
#define FIRSTREV	86
#define LASTREV		95

#define D86		0		/* keywords defined in LM-63-1986 */

#define K_TST		0
#define K_MAN		1
#define K_LMC		2
#define K_LMN		3
#define K_LPC		4
#define K_LMP		5
#define K_BAL		6
#define K_MTC		7
#define K_OTH		8
#define K_SCH		9
#define K_MOR		10
#define K_BLK		11
#define K_EBK		12

#define D91		((1L<<13)-1)	/* keywords defined in LM-63-1991 */

#define K_LMG		13

#define D95		((1L<<14)-1)	/* keywords defined in LM-63-1995 */

char	k_kwd[][20] = {"TEST", "MANUFAC", "LUMCAT", "LUMINAIRE", "LAMPCAT",
			"LAMP", "BALLAST", "MAINTCAT", "OTHER", "SEARCH",
			"MORE", "BLOCK", "ENDBLOCK", "LUMINOUSGEOMETRY"};

long k_defined[] = {D86, D86, D86, D86, D86, D91, D91, D91, D91, D95};

int	filerev = FIRSTREV;

#define keymatch(i,s)	(k_defined[filerev-FIRSTREV]&1L<<(i) &&\
				k_match(k_kwd[i],s))

#define checklamp(s)	(!(k_defined[filerev-FIRSTREV]&(1<<K_LMP|1<<K_LPC)) ||\
				keymatch(K_LMP,s) || keymatch(K_LPC,s))

					/* tilt specs */
#define TLTSTR		"TILT="
#define TLTSTRLEN	5
#define TLTNONE		"NONE"
#define TLTINCL		"INCLUDE"
#define TLT_VERT	1
#define TLT_H0		2
#define TLT_H90		3
					/* photometric types */
#define PM_C		1
#define PM_B		2
					/* unit types */
#define U_FEET		1
#define U_METERS	2
					/* string lengths */
#define MAXLINE		132
#define RMAXWORD	76
					/* file types */
#define T_RAD		".rad"
#define T_DST		".dat"
#define T_TLT		"%.dat"
#define T_OCT		".oct"
					/* shape types */
#define RECT		1
#define DISK		2
#define SPHERE		3

#define MINDIM		.001		/* minimum dimension (point source) */

#define F_M		.3048		/* feet to meters */

#define abspath(p)	(ISDIRSEP((p)[0]) || (p)[0] == '.')

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
char	units[64] = "meters";		/* output units */
int	out2stdout = 0;			/* put out to stdout r.t. file */
int	instantiate = 0;		/* instantiate geometry */
double	illumrad = 0.0;			/* radius for illum sphere */

typedef struct {
	int	isillum;			/* do as illum */
	int	type;				/* RECT, DISK, SPHERE */
	double	mult;				/* candela multiplier */
	double	w, l, h;			/* width, length, height */
	double	area;				/* max. projected area */
} SRCINFO;				/* a source shape (units=meters) */

int	gargc;				/* global argc (minus filenames) */
char	**gargv;			/* global argv */

extern char	*stradd(), *tailtrunc(), *filetrunc(),
		*filename(), *libname(), *fullnam(), *getword(), *atos();
extern float	*matchlamp();
extern time_t	fdate();

#define scnint(fp,ip)	cvtint(ip,getword(fp))
#define scnflt(fp,rp)	cvtflt(rp,getword(fp))
#define isint		isflt			/* IES allows real as integer */


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*outfile = NULL;
	int	status;
	char	outname[RMAXWORD];
	double	d1;
	int	i;
	
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
		case 'g':		/* instatiate geometry? */
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
	gargc = i;
	gargv = argv;
	initlamps();			/* get lamp data (if needed) */
					/* convert ies file(s) */
	if (outfile != NULL) {
		if (i == argc)
			exit(ies2rad(NULL, outfile) == 0 ? 0 : 1);
		else if (i == argc-1)
			exit(ies2rad(argv[i], outfile) == 0 ? 0 : 1);
		else
			goto needsingle;
	} else if (i >= argc) {
		fprintf(stderr, "%s: missing output file specification\n",
				argv[0]);
		exit(1);
	}
	if (out2stdout && i != argc-1)
		goto needsingle;
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


initlamps()				/* set up lamps */
{
	float	*lcol;
	int	status;

	if (lamptype != NULL && !strcmp(lamptype, default_name) &&
			deflamp == NULL)
		return;				/* no need for data */
						/* else load file */
	if ((status = loadlamps(lampdat)) < 0)
		exit(1);
	if (status == 0) {
		fprintf(stderr, "%s: warning - no lamp data\n", lampdat);
		lamptype = default_name;
		return;
	}
	if (deflamp != NULL) {			/* match default type */
		if ((lcol = matchlamp(deflamp)) == NULL)
			fprintf(stderr,
				"%s: warning - unknown default lamp type\n",
					deflamp);
		else
			copycolor(defcolor, lcol);
	}
	if (lamptype != NULL) {			/* match selected type */
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


char *
stradd(dst, src, sep)			/* add a string at dst */
register char	*dst, *src;
int	sep;
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


char *
fullnam(path, fname, suffix)		/* return full path name */
char	*path, *fname, *suffix;
{
	if (prefdir != NULL && abspath(prefdir))
		libname(path, fname, suffix);
	else if (abspath(fname))
		strcpy(stradd(path, fname, 0), suffix);
	else
		libname(stradd(path, libdir, DIRSEP), fname, suffix);

	return(path);
}


char *
libname(path, fname, suffix)		/* return library relative name */
char	*path, *fname, *suffix;
{
	if (abspath(fname))
		strcpy(stradd(path, fname, 0), suffix);
	else
		strcpy(stradd(stradd(path, prefdir, DIRSEP), fname, 0), suffix);

	return(path);
}


char *
filename(path)			/* get final component of pathname */
register char	*path;
{
	register char	*cp;

	for (cp = path; *path; path++)
		if (ISDIRSEP(*path))
			cp = path+1;
	return(cp);
}


char *
filetrunc(path)				/* truncate filename at end of path */
char	*path;
{
	register char	*p1, *p2;

	for (p1 = p2 = path; *p2; p2++)
		if (ISDIRSEP(*p2))
			p1 = p2;
	if (p1 == path && ISDIRSEP(*p1))
		p1++;
	*p1 = '\0';
	return(path);
}


char *
tailtrunc(name)				/* truncate tail of filename */
char	*name;
{
	register char	*p1, *p2;

	for (p1 = filename(name); *p1 == '.'; p1++)
		;
	p2 = NULL;
	for ( ; *p1; p1++)
		if (*p1 == '.')
			p2 = p1;
	if (p2 != NULL)
		*p2 = '\0';
	return(name);
}


blanktrunc(s)				/* truncate spaces at end of line */
char	*s;
{
	register char	*cp;

	for (cp = s; *cp; cp++)
		;
	while (cp-- > s && isspace(*cp))
		;
	*++cp = '\0';
}


k_match(kwd, hdl)			/* header line matches keyword? */
register char	*kwd, *hdl;
{
	if (!*hdl++ == '[')
		return(0);
	while (islower(*hdl) ? toupper(*hdl) == *kwd++ : *hdl == *kwd++)
		if (!*hdl++)
			return(0);
	return(!*kwd & *hdl == ']');
}


char *
keyargs(hdl)				/* return keyword arguments */
register char	*hdl;
{
	while (*hdl && *hdl++ != ']')
		;
	while (isspace(*hdl))
		hdl++;
	return(hdl);
}


putheader(out)				/* print header */
FILE	*out;
{
	register int	i;
	
	putc('#', out);
	for (i = 0; i < gargc; i++) {
		putc(' ', out);
		fputs(gargv[i], out);
	}
	fputs("\n# Dimensions in ", out);
	fputs(units, out);
	putc('\n', out);
}


ies2rad(inpname, outname)		/* convert IES file */
char	*inpname, *outname;
{
	SRCINFO	srcinfo;
	char	buf[MAXLINE], tltid[RMAXWORD];
	char	geomfile[128];
	FILE	*inpfp, *outfp;
	int	lineno = 0;

	geomfile[0] = '\0';
	srcinfo.isillum = 0;
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
	putheader(outfp);
	if (lamptype == NULL)
		lampcolor = NULL;
	while (fgets(buf,sizeof(buf),inpfp) != NULL
			&& strncmp(buf,TLTSTR,TLTSTRLEN)) {
		blanktrunc(buf);
		if (!buf[0])
			continue;
		if (!lineno++ && !strncmp(buf, MAGICID, LMAGICID)) {
			filerev = atoi(buf+LMAGICID);
			if (filerev < FIRSTREV)
				filerev = FIRSTREV;
			else if (filerev > LASTREV)
				filerev = LASTREV;
		}
		fputs("#<", outfp);
		fputs(buf, outfp);
		putc('\n', outfp);
		if (lampcolor == NULL && checklamp(buf))
			lampcolor = matchlamp( buf[0] == '[' ?
						keyargs(buf) : buf );
		if (keymatch(K_LMG, buf)) {		/* geometry file */
			strcpy(geomfile, inpname);
			strcpy(filename(geomfile), keyargs(buf));
			srcinfo.isillum = 1;
		}
	}
	if (lampcolor == NULL) {
		fprintf(stderr, "%s: warning - no lamp type\n", inpname);
		fputs("# Unknown lamp type (used default)\n", outfp);
		lampcolor = defcolor;
	} else if (lamptype == NULL)
		fprintf(outfp,"# CIE(x,y) = (%f,%f)\n# Depreciation = %.1f%%\n",
				lampcolor[3], lampcolor[4], 100.*lampcolor[5]);
	if (feof(inpfp)) {
		fprintf(stderr, "%s: not in IES format\n", inpname);
		goto readerr;
	}
	atos(tltid, RMAXWORD, buf+TLTSTRLEN);
	if (inpfp == stdin)
		buf[0] = '\0';
	else
		filetrunc(strcpy(buf, inpname));
	if (dotilt(inpfp, outfp, buf, tltid, outname, tltid) != 0) {
		fprintf(stderr, "%s: bad tilt data\n", inpname);
		goto readerr;
	}
	if (dosource(&srcinfo, inpfp, outfp, tltid, outname) != 0) {
		fprintf(stderr, "%s: bad luminaire data\n", inpname);
		goto readerr;
	}
	fclose(inpfp);
					/* cvgeometry closes outfp */
	if (cvgeometry(geomfile, &srcinfo, outname, outfp) != 0) {
		fprintf(stderr, "%s: bad geometry file\n", geomfile);
		return(-1);
	}
	return(0);
readerr:
	fclose(inpfp);
	fclose(outfp);
	unlink(fullnam(buf,outname,T_RAD));
	return(-1);
}


dotilt(in, out, dir, tltspec, dfltname, tltid)	/* convert tilt data */
FILE	*in, *out;
char	*dir, *tltspec, *dfltname, *tltid;
{
	int	nangles, tlt_type;
	double	minmax[2];
	char	buf[PATH_MAX], tltname[RMAXWORD];
	FILE	*datin, *datout;

	if (!strcmp(tltspec, TLTNONE)) {
		datin = NULL;
		strcpy(tltid, "void");
	} else if (!strcmp(tltspec, TLTINCL)) {
		datin = in;
		strcpy(tltname, dfltname);
	} else {
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
	if (datin != NULL) {
		if ((datout = fopen(fullnam(buf,tltname,T_TLT),"w")) == NULL) {
			perror(buf);
			if (datin != in)
				fclose(datin);
			return(-1);
		}
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
		strcat(strcpy(tltid, filename(tltname)), "_tilt");
		fprintf(out, "\nvoid brightdata %s\n", tltid);
		libname(buf,tltname,T_TLT);
		switch (tlt_type) {
		case TLT_VERT:			/* vertical */
			fprintf(out, "4 noop %s tilt.cal %s\n", buf,
				minmax[1]>90.+FTINY ? "tilt_ang" : "tilt_ang2");
			break;
		case TLT_H0:			/* horiz. in 0 deg. plane */
			fprintf(out, "6 noop %s tilt.cal %s -rz 90\n", buf,
			minmax[1]>90.+FTINY ? "tilt_xang" : "tilt_xang2");
			break;
		case TLT_H90:
			fprintf(out, "4 noop %s tilt.cal %s\n", buf,
			minmax[1]>90.+FTINY ? "tilt_xang" : "tilt_xang2");
			break;
		default:
			fprintf(stderr,
				"%s: illegal lamp to luminaire geometry (%d)\n",
				tltspec, tlt_type);
			return(-1);
		}
		fprintf(out, "0\n0\n");
	}
	return(0);
}


dosource(sinf, in, out, mod, name)	/* create source and distribution */
SRCINFO	*sinf;
FILE	*in, *out;
char	*mod, *name;
{
	char	buf[PATH_MAX], id[RMAXWORD];
	FILE	*datout;
	double	mult, bfactor, pfactor, width, length, height, wattage;
	double	bounds[2][2];
	int	nangles[2], pmtype, unitype;
	double	d1;
	int	doupper, dolower, dosides; 

	if (!isint(getword(in)) || !isflt(getword(in)) || !scnflt(in,&mult)
			|| !scnint(in,&nangles[0]) || !scnint(in,&nangles[1])
			|| !scnint(in,&pmtype) || !scnint(in,&unitype)
			|| !scnflt(in,&width) || !scnflt(in,&length)
			|| !scnflt(in,&height) || !scnflt(in,&bfactor)
			|| !scnflt(in,&pfactor) || !scnflt(in,&wattage)) {
		fprintf(stderr, "dosource: bad lamp specification\n");
		return(-1);
	}
	sinf->mult = multiplier*mult*bfactor*pfactor;
	if (nangles[0] < 2 || nangles[1] < 1) {
		fprintf(stderr, "dosource: too few measured angles\n");
		return(-1);
	}
	if (unitype == U_FEET) {
		width *= F_M;
		length *= F_M;
		height *= F_M;
	}
	if (makeshape(sinf, width, length, height) != 0) {
		fprintf(stderr, "dosource: illegal source dimensions");
		return(-1);
	}
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
	fprintf(out, "# %g watt luminaire, lamp*ballast factor = %g\n",
			wattage, bfactor*pfactor);
	strcat(strcpy(id, filename(name)), "_dist");
	fprintf(out, "\n%s brightdata %s\n", mod, id);
	if (nangles[1] < 2)
		fprintf(out, "4 ");
	else if (pmtype == PM_B)
		fprintf(out, "5 ");
	else if (FEQ(bounds[1][0],90.) && FEQ(bounds[1][1],270.))
		fprintf(out, "7 ");
	else
		fprintf(out, "5 ");
	dolower = (bounds[0][0] < 90.);
	doupper = (bounds[0][1] > 90.);
	dosides = (doupper & dolower && sinf->h > MINDIM);
	fprintf(out, "%s %s source.cal ",
			sinf->type==SPHERE ? "corr" :
			!dosides ? "flatcorr" :
			sinf->type==DISK ? "cylcorr" : "boxcorr",
			libname(buf,name,T_DST));
	if (pmtype == PM_B) {
		if (FEQ(bounds[1][0],0.))
			fprintf(out, "srcB_horiz2 ");
		else
			fprintf(out, "srcB_horiz ");
		fprintf(out, "srcB_vert ");
	} else /* pmtype == PM_A */ {
		if (nangles[1] >= 2) {
			d1 = bounds[1][1] - bounds[1][0];
			if (d1 <= 90.+FTINY)
				fprintf(out, "src_phi4 ");
			else if (d1 <= 180.+FTINY) {
				if (FEQ(bounds[1][0],90.))
					fprintf(out, "src_phi2+90 ");
				else
					fprintf(out, "src_phi2 ");
			} else
				fprintf(out, "src_phi ");
			fprintf(out, "src_theta ");
			if (FEQ(bounds[1][0],90.) && FEQ(bounds[1][1],270.))
				fprintf(out, "-rz -90 ");
		} else
			fprintf(out, "src_theta ");
	}
	if (!dosides || sinf->type == SPHERE)
		fprintf(out, "\n0\n1 %g\n", sinf->mult/sinf->area);
	else if (sinf->type == DISK)
		fprintf(out, "\n0\n3 %g %g %g\n", sinf->mult,
				sinf->w, sinf->h);
	else
		fprintf(out, "\n0\n4 %g %g %g %g\n", sinf->mult,
				sinf->l, sinf->w, sinf->h);
	if (putsource(sinf, out, id, filename(name),
			dolower, doupper, dosides) != 0)
		return(-1);
	return(0);
}


putsource(shp, fp, mod, name, dolower, doupper, dosides) /* put out source */
SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
int	dolower, doupper;
{
	char	lname[RMAXWORD];
	
	strcat(strcpy(lname, name), "_light");
	fprintf(fp, "\n%s %s %s\n", mod,
			shp->isillum ? "illum" : "light", lname);
	fprintf(fp, "0\n0\n3 %g %g %g\n",
			lampcolor[0], lampcolor[1], lampcolor[2]);
	switch (shp->type) {
	case RECT:
		if (dolower)
			putrectsrc(shp, fp, lname, name, 0);
		if (doupper)
			putrectsrc(shp, fp, lname, name, 1);
		if (dosides)
			putsides(shp, fp, lname, name);
		break;
	case DISK:
		if (dolower)
			putdisksrc(shp, fp, lname, name, 0);
		if (doupper)
			putdisksrc(shp, fp, lname, name, 1);
		if (dosides)
			putcyl(shp, fp, lname, name);
		break;
	case SPHERE:
		putspheresrc(shp, fp, lname, name);
		break;
	}
	return(0);
}


makeshape(shp, width, length, height)		/* make source shape */
register SRCINFO	*shp;
double	width, length, height;
{
	if (illumrad/meters2out >= MINDIM/2.) {
		shp->isillum = 1;
		shp->type = SPHERE;
		shp->w = shp->l = shp->h = 2.*illumrad / meters2out;
	} else if (width < MINDIM) {
		width = -width;
		if (width < MINDIM) {
			shp->type = SPHERE;
			shp->w = shp->l = shp->h = MINDIM;
		} else if (height < .5*width) {
			shp->type = DISK;
			shp->w = shp->l = width;
			if (height >= MINDIM)
				shp->h = height;
			else
				shp->h = .5*MINDIM;
		} else {
			shp->type = SPHERE;
			shp->w = shp->l = shp->h = width;
		}
	} else {
		shp->type = RECT;
		shp->w = width;
		if (length >= MINDIM)
			shp->l = length;
		else
			shp->l = MINDIM;
		if (height >= MINDIM)
			shp->h = height;
		else
			shp->h = .5*MINDIM;
	}
	switch (shp->type) {
	case RECT:
		shp->area = shp->w * shp->l;
		break;
	case DISK:
	case SPHERE:
		shp->area = PI/4. * shp->w * shp->w;
		break;
	}
	return(0);
}


putrectsrc(shp, fp, mod, name, up)		/* rectangular source */
SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
int	up;
{
	if (up)
		putrect(shp, fp, mod, name, ".u", 4, 5, 7, 6);
	else
		putrect(shp, fp, mod, name, ".d", 0, 2, 3, 1);
}


putsides(shp, fp, mod, name)			/* put out sides of box */
register SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
{
	putrect(shp, fp, mod, name, ".1", 0, 1, 5, 4);
	putrect(shp, fp, mod, name, ".2", 1, 3, 7, 5);
	putrect(shp, fp, mod, name, ".3", 3, 2, 6, 7);
	putrect(shp, fp, mod, name, ".4", 2, 0, 4, 6);
}
	

putrect(shp, fp, mod, name, suffix, a, b, c, d)	/* put out a rectangle */
SRCINFO	*shp;
FILE	*fp;
char	*mod, *name, *suffix;
int	a, b, c, d;
{
	fprintf(fp, "\n%s polygon %s%s\n0\n0\n12\n", mod, name, suffix);
	putpoint(shp, fp, a);
	putpoint(shp, fp, b);
	putpoint(shp, fp, c);
	putpoint(shp, fp, d);
}


putpoint(shp, fp, p)				/* put out a point */
register SRCINFO	*shp;
FILE	*fp;
int	p;
{
	static double	mult[2] = {-.5, .5};

	fprintf(fp, "\t%g\t%g\t%g\n",
			mult[p&1]*shp->l*meters2out,
			mult[p>>1&1]*shp->w*meters2out,
			mult[p>>2]*shp->h*meters2out);
}


putdisksrc(shp, fp, mod, name, up)		/* put out a disk source */
register SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
int	up;
{
	if (up) {
		fprintf(fp, "\n%s ring %s.u\n", mod, name);
		fprintf(fp, "0\n0\n8\n");
		fprintf(fp, "\t0 0 %g\n", .5*shp->h*meters2out);
		fprintf(fp, "\t0 0 1\n");
		fprintf(fp, "\t0 %g\n", .5*shp->w*meters2out);
	} else {
		fprintf(fp, "\n%s ring %s.d\n", mod, name);
		fprintf(fp, "0\n0\n8\n");
		fprintf(fp, "\t0 0 %g\n", -.5*shp->h*meters2out);
		fprintf(fp, "\t0 0 -1\n");
		fprintf(fp, "\t0 %g\n", .5*shp->w*meters2out);
	}
}


putcyl(shp, fp, mod, name)			/* put out a cylinder */
register SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
{
	fprintf(fp, "\n%s cylinder %s.c\n", mod, name);
	fprintf(fp, "0\n0\n7\n");
	fprintf(fp, "\t0 0 %g\n", .5*shp->h*meters2out);
	fprintf(fp, "\t0 0 %g\n", -.5*shp->h*meters2out);
	fprintf(fp, "\t%g\n", .5*shp->w*meters2out);
}


putspheresrc(shp, fp, mod, name)		/* put out a sphere source */
SRCINFO	*shp;
FILE	*fp;
char	*mod, *name;
{
	fprintf(fp, "\n%s sphere %s.s\n", mod, name);
	fprintf(fp, "0\n0\n4 0 0 0 %g\n", .5*shp->w*meters2out);
}


cvdata(in, out, ndim, npts, mult, lim)		/* convert data */
FILE	*in, *out;
int	ndim, npts[];
double	mult, lim[][2];
{
	double	*pt[4];
	register int	i, j;
	double	val;
	int	total;

	total = 1; j = 0;
	for (i = 0; i < ndim; i++)
		if (npts[i] > 1) {
			total *= npts[i];
			j++;
		}
	fprintf(out, "%d\n", j);
					/* get coordinates */
	for (i = 0; i < ndim; i++) {
		pt[i] = (double *)malloc(npts[i]*sizeof(double));
		for (j = 0; j < npts[i]; j++)
			if (!scnflt(in, &pt[i][j]))
				return(-1);
		if (lim != NULL) {
			lim[i][0] = pt[i][0];
			lim[i][1] = pt[i][npts[i]-1];
		}
	}
					/* write out in reverse */
	for (i = ndim-1; i >= 0; i--) {
		if (npts[i] > 1) {
			for (j = 1; j < npts[i]-1; j++)
				if (!FEQ(pt[i][j]-pt[i][j-1],
						pt[i][j+1]-pt[i][j]))
					break;
			if (j == npts[i]-1)
				fprintf(out, "%g %g %d\n", pt[i][0], pt[i][j],
						npts[i]);
			else {
				fprintf(out, "0 0 %d", npts[i]);
				for (j = 0; j < npts[i]; j++) {
					if (j%4 == 0)
						putc('\n', out);
					fprintf(out, "\t%g", pt[i][j]);
				}
				putc('\n', out);
			}
		}
		free((void *)pt[i]);
	}
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


char *
getword(fp)			/* scan a word from fp */
register FILE	*fp;
{
	static char	wrd[RMAXWORD];
	register char	*cp;
	register int	c;

	while (isspace(c=getc(fp)))
		;
	for (cp = wrd; c != EOF && cp < wrd+RMAXWORD-1;
			*cp++ = c, c = getc(fp))
		if (isspace(c) || c == ',') {
			while (isspace(c))
				c = getc(fp);
			if (c != EOF & c != ',')
				ungetc(c, fp);
			*cp = '\0';
			return(wrd);
		}
	*cp = '\0';
	return(cp > wrd ? wrd : NULL);
}


cvtint(ip, wrd)			/* convert a word to an integer */
int	*ip;
char	*wrd;
{
	if (wrd == NULL || !isint(wrd))
		return(0);
	*ip = atoi(wrd);
	return(1);
}


cvtflt(rp, wrd)			/* convert a word to a double */
double	*rp;
char	*wrd;
{
	if (wrd == NULL || !isflt(wrd))
		return(0);
	*rp = atof(wrd);
	return(1);
}


cvgeometry(inpname, sinf, outname, outfp)
char	*inpname;
register SRCINFO	*sinf;
char	*outname;
FILE	*outfp;			/* close output file upon return */
{
	char	buf[256];
	register char	*cp;

	if (inpname == NULL || !inpname[0]) {	/* no geometry file */
		fclose(outfp);
		return(0);
	}
	putc('\n', outfp);
	strcpy(buf, "mgf2rad ");		/* build mgf2rad command */
	cp = buf+8;
	if (!FEQ(sinf->mult, 1.0)) {
		sprintf(cp, "-m %f ", sinf->mult);
		cp += strlen(cp);
	}
	sprintf(cp, "-g %f %s ",
		sqrt(sinf->w*sinf->w + sinf->h*sinf->h + sinf->l*sinf->l),
			inpname);
	cp += strlen(cp);
	if (instantiate) {		/* instantiate octree */
		strcpy(cp, "| oconv - > ");
		cp += 12;
		fullnam(cp,outname,T_OCT);
		if (fdate(inpname) > fdate(outname) &&
				system(buf)) {		/* create octree */
			fclose(outfp);
			return(-1);
		}
		fprintf(outfp, "void instance %s_inst\n", outname);
		if (!FEQ(meters2out, 1.0))
			fprintf(outfp, "3 %s -s %f\n",
					libname(buf,outname,T_OCT),
					meters2out);
		else
			fprintf(outfp, "1 %s\n", libname(buf,outname,T_OCT));
		fprintf(outfp, "0\n0\n");
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
