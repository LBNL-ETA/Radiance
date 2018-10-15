#ifndef lint
static const char RCSid[] = "$Id: radcompare.c,v 2.4 2018/10/15 21:47:52 greg Exp $";
#endif
/*
 * Compare Radiance files for significant differences
 *
 *	G. Ward
 */

#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "platform.h"
#include "rtio.h"
#include "resolu.h"
#include "color.h"
#include "lookup.h"
			/* Reporting levels */
#define REP_QUIET	0	/* no reporting */
#define REP_ERROR	1	/* report errors only */
#define REP_WARN	2	/* report warnings as well */
#define REP_VERBOSE	3	/* verbose reporting */

int	report = REP_WARN;	/* reporting level */

int	ign_header = 0;		/* ignore header differences? */

double	rel_min = 1e-5;		/* positive for relative comparisons */

double	rms_lim = 0.01;		/* RMS difference limit */

double	max_lim = 0.25;		/* difference limit if non-negative */

int	lin1cnt=0, lin2cnt=0;	/* file line position */

				/* file types */
const char	*file_type[] = {
			"Unrecognized",
			"TEXT_generic",
			"ascii",
			COLRFMT,
			CIEFMT,
			"float",
			"double",
			"BSDF_RBFmesh",
			"Radiance_octree",
			"Radiance_tmesh",
			"BINARY_unknown",
			NULL	/* terminator */
		};

enum {TYP_UNKNOWN, TYP_TEXT, TYP_ASCII, TYP_RGBE, TYP_XYZE, TYP_FLOAT,
		TYP_DOUBLE, TYP_RBFMESH, TYP_OCTREE, TYP_TMESH, TYP_BINARY};
		
#define has_header(t)	(!( 1L<<(t) & (1L<<TYP_TEXT | 1L<<TYP_BINARY) ))

				/* header variables to always ignore */
const char	*hdr_ignkey[] = {
			"SOFTWARE",
			"CAPDATE",
			"GMT",
			NULL	/* terminator */
		};
				/* header variable settings */
LUTAB	hdr1 = LU_SINIT(free,free);
LUTAB	hdr2 = LU_SINIT(free,free);

				/* advance appropriate file line count */
#define adv_linecnt(htp)	(lin1cnt += (htp == &hdr1), \
					lin2cnt += (htp == &hdr2))

				/* input files */
char		*progname = NULL;
const char	stdin_name[] = "<stdin>";
const char	*f1name=NULL, *f2name=NULL;
FILE		*f1in=NULL, *f2in=NULL;

				/* running real differences */
double	diff2sum = 0;
int	nsum = 0;

/* Report usage and exit */
static void
usage()
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-h][-s|-w|-v][-rel min_test][-rms epsilon][-max epsilon] reference test\n",
			stderr);
	exit(1);
}

/* Get type ID from name (or 0 if not found) */
static int
xlate_type(const char *nm)
{
	int	i;

	if (!nm || !*nm)
		return(TYP_UNKNOWN);
	for (i = 1; file_type[i]; i++)
		if (!strcmp(nm, file_type[i]))
			return(i);
	return(TYP_UNKNOWN);
}

/* Compare real values and keep track of differences */
static int
real_check(double r1, double r2)
{
	double	diff2 = (r1 - r2)*(r1 - r2);

	if (rel_min > 0) {	/* doing relative differences? */
		double	av2 = .25*(r1*r1 + 2.*fabs(r1*r2) + r2*r2);
		if (av2 > rel_min*rel_min)
			diff2 /= av2;
	}
	if (max_lim >= 0 && diff2 > max_lim*max_lim) {
		if (report != REP_QUIET)
			printf(
			"%s: %sdifference between %.8g and %.8g exceeds epsilon\n",
					progname,
					(rel_min > 0) ? "relative " : "",
					r1, r2);
		return(0);
	}
	diff2sum += diff2;
	nsum++;
	return(1);
}

/* Compare two strings for equivalence */
static int
equiv_string(char *s1, char *s2)
{
#define CLS_STR		0
#define CLS_INT		1
#define CLS_FLT		2
				/* skip whitespace at beginning */
	while (isspace(*s1)) s1++;
	while (isspace(*s2)) s2++;
	while (*s1) {		/* check each word */
		int	inquote;
		if (!*s2)	/* unexpected EOL in s2? */
			return(0);
		inquote = *s1;
		if ((inquote != '\'') & (inquote != '"'))
			inquote = 0;
		if (inquote) {	/* quoted text must match exactly */
			if (*s1++ != *s2++)
				return(0);
			while (*s1 != inquote) {
				if (!*s1)
					return(0);
				if (*s1++ != *s2++)
					return(0);
			}
			s1++;
			if (*s2++ != inquote)
				return(0);
		} else {	/* else classify word type */
			char	*s1s = s1;
			char	*s2s = s2;
			int	cls = CLS_STR;
			s1 = sskip(s1);
			s2 = sskip(s2);
			if (iskip(s1s) == s1) {
				if (iskip(s2s) == s2)
					cls = CLS_INT;
				else if (fskip(s2s) == s2)
					cls = CLS_FLT;
			} else if (fskip(s1s) == s1) {
				if (fskip(s2s) != s2)
					return(0);
				cls = CLS_FLT;
			}
			switch (cls) {
			case CLS_INT:		/* strncmp() faster */
			case CLS_STR:
				if (s1 - s1s != s2 - s2s)
					return(0);
				if (strncmp(s1s, s2s, s1 - s1s))
					return(0);
				break;
			case CLS_FLT:
				if (!real_check(atof(s1s), atof(s2s)))
					return(0);
				break;
			}
		}
		while (isspace(*s1)) s1++;
		while (isspace(*s2)) s2++;
	}
	return(!*s2);		/* match if we reached the end of s2, too */
#undef CLS_STR
#undef CLS_INT
#undef CLS_FLT
}

/* Check if string is var=value pair and set if not in ignore list */
static int
setheadvar(char *val, void *p)
{
	LUTAB	*htp = (LUTAB *)p;
	LUENT	*tep;
	char	*key;
	int	kln, vln;
	int	n;

	adv_linecnt(htp);	/* side-effect is to count lines */
	if (!isalpha(*val))	/* key must start line */
		return(0);
	key = val++;
	while (*val && !isspace(*val) & (*val != '='))
		val++;
	kln = val - key;
	while (isspace(*val))	/* check for value */
		*val++ = '\0';
	if (*val != '=')
		return(0);
	*val++ = '\0';
	while (isspace(*val))
		val++;
	if (!*val)		/* nothing set? */
		return(0);
	vln = strlen(val);	/* eliminate space and newline at end */
	while (--vln > 0 && isspace(val[vln]))
		;
	val[++vln] = '\0';
				/* check if key to ignore */
	for (n = 0; hdr_ignkey[n]; n++)
		if (!strcmp(key, hdr_ignkey[n]))
			return(0);
	if (!(tep = lu_find(htp, key)))
		return(-1);	/* memory allocation error */
	if (!tep->key)
		tep->key = strcpy(malloc(kln+1), key);
	if (tep->data)
		free(tep->data);
	tep->data = strcpy(malloc(vln+1), val);
	return(1);
}

/* Lookup correspondent in other header */
static int
match_val(const LUENT *ep1, void *p2)
{
	const LUENT	*ep2 = lu_find((LUTAB *)p2, ep1->key);
	if (!ep2 || !ep2->data) {
		if (report != REP_QUIET)
			printf("%s: Variable '%s' missing in '%s'\n",
					progname, ep1->key, f2name);
		return(-1);
	}
	if (!equiv_string((char *)ep1->data, (char *)ep2->data)) {
		if (report != REP_QUIET) {
			printf("%s: Header variable '%s' has different values\n",
					progname, ep1->key);
			if (report >= REP_VERBOSE) {
				printf("%s: %s=%s\n", f1name,
						ep1->key, (char *)ep1->data);
				printf("%s: %s=%s\n", f2name,
						ep2->key, (char *)ep2->data);
			}
		}
		return(-1);
	}
	return(1);		/* good match */
}

/* Compare two sets of header variables */
static int
headers_match(LUTAB *hp1, LUTAB *hp2)
{
	int	ne = lu_doall(hp1, match_val, hp2);
	if (ne < 0)
		return(0);	/* something didn't match! */
				/* non-fatal if second header has extra */
	if (report >= REP_WARN && (ne = lu_doall(hp2, NULL, NULL) - ne))
		printf("%s: '%s' has %d extra header setting(s)\n",
					progname, f2name, ne);
	return(1);		/* good match */
}

/* Check generic input to determine if it is binary, -1 on error */
static int
input_is_binary(FILE *fin)
{
	int	n = 0;
	int	c = 0;

	while ((c = getc(fin)) != EOF) {
		if (!c | (c > 127))
			break;			/* non-ascii character */
		if (++n >= 10240)
			break;			/* enough to be confident */
	}
	if (!n)
		return(-1);			/* first read failed */
	if (fseek(fin, 0L, 0) < 0)
		return(-1);			/* rewind failed */
	return(!c | (c > 127));
}

/* Identify and return data type from header (if one) */
static int
identify_type(const char *name, FILE *fin, LUTAB *htp)
{
	extern const char	HDRSTR[];
	int			c;
						/* check magic header start */
	if ((c = getc(fin)) != HDRSTR[0]) {
		if (c == EOF) goto badeof;
		ungetc(c, fin);
		c = 0;
	} else if ((c = getc(fin)) != HDRSTR[1]) {
		if (c == EOF) goto badeof;
		ungetc(c, fin); ungetc(HDRSTR[0], fin);
		c = 0;
	}
	if (c) {				/* appears to have a header */
		char	sbuf[32];
		if (!fgets(sbuf, sizeof(sbuf), fin))
			goto badeof;
		adv_linecnt(htp);		/* for #?ID string */
		if (report >= REP_WARN && strncmp(sbuf, "RADIANCE", 8)) {
			fputs(name, stdout);
			fputs(": warning - unexpected header ID: ", stdout);
			fputs(sbuf, stdout);
		}
		if (getheader(fin, setheadvar, htp) < 0) {
			fputs(name, stderr);
			fputs(": Unknown error reading header\n", stderr);
			return(-1);
		}
		adv_linecnt(htp);		/* for trailing emtpy line */
		return(xlate_type((const char *)lu_find(htp,"FORMAT")->data));
	}
	c = input_is_binary(fin);		/* else peek to see if binary */
	if (c < 0) {
		fputs(name, stderr);
		fputs(": read/seek error\n", stderr);
		return(-1);
	}
	if (c)
		return(TYP_BINARY);
	SET_FILE_TEXT(fin);			/* originally set to binary */
	return(TYP_TEXT);
badeof:
	if (report != REP_QUIET) {
		fputs(name, stdout);
		fputs(": Unexpected end-of-file\n", stdout);
	}
	return(-1);
}

/* Check that overall RMS error is below threshold */
static int
good_RMS()
{
	if (!nsum)
		return(1);
	if (diff2sum/(double)nsum > rms_lim*rms_lim) {
		if (report != REP_QUIET)
			printf(
	"%s: %sRMS difference between '%s' and '%s' of %.5g exceeds limit\n",
					progname,
					(rel_min > 0) ? "relative " : "",
					f1name, f2name,
					sqrt(diff2sum/(double)nsum));
		return(0);
	}
	if (report >= REP_VERBOSE)
		printf("%s: %sRMS difference of reals in '%s' and '%s' is %.5g\n",
				progname, (rel_min > 0) ? "relative " : "",
				f1name, f2name, sqrt(diff2sum/(double)nsum));
	return(1);
}

/* Compare two inputs as generic binary files */
static int
compare_binary()
{
	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as binary\n", stdout);
	}
	for ( ; ; ) {				/* exact byte matching */
		int	c1 = getc(f1in);
		int	c2 = getc(f2in);
		if (c1 == EOF) {
			if (c2 == EOF)
				return(1);	/* success! */
			if (report != REP_QUIET) {
				fputs(f1name, stdout);
				fputs(": Unexpected end-of-file\n", stdout);
			}
			return(0);
		}
		if (c2 == EOF) {
			if (report != REP_QUIET) {
				fputs(f2name, stdout);
				fputs(": Unexpected end-of-file\n", stdout);
			}
			return(0);
		}
		if (c1 != c2)
			break;			/* quit and report difference */
	}
	if (report == REP_QUIET)
		return(0);
	printf("%s: binary files '%s' and '%s' differ at offset %ld|%ld\n",
			progname, f1name, f2name, ftell(f1in), ftell(f2in));
	return(0);
}

/* Compare two inputs as generic text files */
static int
compare_text()
{
	char	l1buf[4096], l2buf[4096];

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as ASCII text\n", stdout);
	}
						/* compare a line at a time */
	while (fgets(l1buf, sizeof(l1buf), f1in)) {
		lin1cnt++;
		if (!*sskip2(l1buf,0))
			continue;		/* ignore empty lines */
		while (fgets(l2buf, sizeof(l2buf), f2in)) {
			lin2cnt++;
			if (*sskip2(l2buf,0))
				break;		/* found other non-empty line */
		}
		if (feof(f2in)) {
			if (report != REP_QUIET) {
				fputs(f2name, stdout);
				fputs(": Unexpected end-of-file\n", stdout);
			}
			return(0);
		}
						/* compare non-empty lines */
		if (!equiv_string(l1buf, l2buf)) {
			if (report != REP_QUIET) {
				printf("%s: inputs '%s' and '%s' differ at line %d|%d\n",
						progname, f1name, f2name,
						lin1cnt, lin2cnt);
				if (report >= REP_VERBOSE) {
					fputs("------------- Mismatch -------------\n", stdout);
					printf("%s(%d):\t%s", f1name,
							lin1cnt, l1buf);
					printf("%s(%d):\t%s", f2name,
							lin2cnt, l2buf);
				}
			}
			return(0);
		}
	}
						/* check for EOF on input 2 */
	while (fgets(l2buf, sizeof(l2buf), f2in)) {
		if (!*sskip2(l2buf,0))
			continue;
		if (report != REP_QUIET) {
			fputs(f1name, stdout);
			fputs(": Unexpected end-of-file\n", stdout);
		}
		return(0);
	}
	return(good_RMS());			/* final check for reals */
}

/* Compare two inputs that are known to be RGBE or XYZE images */
static int
compare_hdr()
{
	RESOLU	rs1, rs2;
	COLOR	*scan1, *scan2;
	int	x, y;

	fgetsresolu(&rs1, f1in);
	fgetsresolu(&rs2, f2in);
	if (rs1.rt != rs2.rt) {
		if (report != REP_QUIET)
			printf(
			"%s: Images '%s' and '%s' have different pixel ordering\n",
					progname, f1name, f2name);
		return(0);
	}
	if ((rs1.xr != rs2.xr) | (rs1.yr != rs2.yr)) {
		if (report != REP_QUIET)
			printf(
			"%s: Images '%s' and '%s' are different sizes\n",
					progname, f1name, f2name);
		return(0);
	}
	scan1 = (COLOR *)malloc(scanlen(&rs1)*sizeof(COLOR));
	scan2 = (COLOR *)malloc(scanlen(&rs2)*sizeof(COLOR));
	if (!scan1 | !scan2) {
		fprintf(stderr, "%s: out of memory in compare_hdr()\n", progname);
		exit(2);
	}
	for (y = 0; y < numscans(&rs1); y++) {
		if ((freadscan(scan1, scanlen(&rs1), f1in) < 0) |
				(freadscan(scan2, scanlen(&rs2), f2in) < 0)) {
			if (report != REP_QUIET)
				printf("%s: Unexpected end-of-file\n",
						progname);
			free(scan1);
			free(scan2);
			return(0);
		}
		for (x = scanlen(&rs1); x--; ) {
			if (real_check(colval(scan1[x],RED),
						colval(scan2[x],RED)) &
					real_check(colval(scan1[x],GRN),
						colval(scan2[x],GRN)) &
					real_check(colval(scan1[x],BLU),
						colval(scan2[x],BLU)))
				continue;
			if (report != REP_QUIET) {
				printf(
				"%s: pixels at scanline %d offset %d differ\n",
					progname, y, x);
				if (report >= REP_VERBOSE) {
					printf("%s: (R,G,B)=(%g,%g,%g)\n",
						f1name, colval(scan1[x],RED),
						colval(scan1[x],GRN),
						colval(scan1[x],BLU));
					printf("%s: (R,G,B)=(%g,%g,%g)\n",
						f2name, colval(scan2[x],RED),
						colval(scan2[x],GRN),
						colval(scan2[x],BLU));
				}
			}
			free(scan1);
			free(scan2);
			return(0);
		}
	}
	free(scan1);
	free(scan2);
	return(good_RMS());			/* final check of RMS */
}

const char	nsuffix[10][3] = {		/* 1st, 2nd, 3rd, etc. */
			"th","st","nd","rd","th","th","th","th","th","th"
		};

/* Compare two inputs that are known to be 32-bit floating-point data */
static int
compare_float()
{
	long	nread = 0;
	float	f1, f2;

	while (getbinary(&f1, sizeof(f1), 1, f1in)) {
		if (!getbinary(&f2, sizeof(f2), 1, f2in))
			goto badeof;
		++nread;
		if (real_check(f1, f2))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s float values differ\n",
					progname, nread, nsuffix[nread%10]);
		return(0);
	}
	if (!getbinary(&f2, sizeof(f2), 1, f2in))
		return(good_RMS());		/* final check of RMS */
badeof:
	if (report != REP_QUIET)
		printf("%s: Unexpected end-of-file\n", progname);
	return(0);
}

/* Compare two inputs that are known to be 64-bit floating-point data */
static int
compare_double()
{
	long	nread = 0;
	double	f1, f2;

	while (getbinary(&f1, sizeof(f1), 1, f1in)) {
		if (!getbinary(&f2, sizeof(f2), 1, f2in))
			goto badeof;
		++nread;
		if (real_check(f1, f2))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s float values differ\n",
					progname, nread, nsuffix[nread%10]);
		return(0);
	}
	if (!getbinary(&f2, sizeof(f2), 1, f2in))
		return(good_RMS());		/* final check of RMS */
badeof:
	if (report != REP_QUIET)
		printf("%s: Unexpected end-of-file\n", progname);
	return(0);
}

/* Compare two Radiance files for equivalence */
int
main(int argc, char *argv[])
{
	int	typ1, typ2;
	int	a;

	progname = argv[0];
	for (a = 1; a < argc && argv[a][0] == '-'; a++) {
		switch (argv[a][1]) {
		case 'h':			/* ignore header info. */
			ign_header = !ign_header;
			continue;
		case 's':			/* silent operation */
			report = REP_QUIET;
			continue;
		case 'w':			/* turn off warnings */
			if (report > REP_ERROR)
				report = REP_ERROR;
			continue;
		case 'v':			/* turn on verbose mode */
			report = REP_VERBOSE;
			continue;
		case 'm':			/* maximum epsilon */
			max_lim = atof(argv[++a]);
			continue;
		case 'r':
			if (argv[a][2] == 'e')		/* relative difference */
				rel_min = atof(argv[++a]);
			else if (argv[a][2] == 'm')	/* RMS limit */
				rms_lim = atof(argv[++a]);
			else
				usage();
			continue;
		case '\0':			/* first file from stdin */
			f1in = stdin;
			f1name = stdin_name;
			break;
		default:
			usage();
		}
		break;
	}
	if (a != argc-2)
		usage();
	if (!f1name) f1name = argv[a];
	if (!f2name) f2name = argv[a+1];

	if (!strcmp(f1name, f2name)) {		/* inputs are same? */
		if (report >= REP_WARN)
			printf("%s: warning - identical inputs given\n",
					progname);
		return(0);
	}
						/* open inputs */
	SET_FILE_BINARY(stdin);			/* in case we're using it */
	if (!f1in && !(f1in = fopen(f1name, "rb"))) {
		fprintf(stderr, "%s: cannot open for reading\n", f1name);
		return(1);
	}
	if (!strcmp(f2name, "-")) {
		f2in = stdin;
		f2name = stdin_name;
	} else if (!(f2in = fopen(f2name, "rb"))) {
		fprintf(stderr, "%s: cannot open for reading\n", f2name);
		return(1);
	}
						/* load headers */
	if ((typ1 = identify_type(f1name, f1in, &hdr1)) < 0)
		return(1);
	if ((typ2 = identify_type(f2name, f2in, &hdr2)) < 0)
		return(1);
	if (typ1 != typ2) {
		if (report != REP_QUIET)
			printf("%s: '%s' is %s and '%s' is %s\n",
					progname, f1name, file_type[typ1],
					f2name, file_type[typ2]);
		return(1);
	}
	ign_header |= !has_header(typ1);	/* check headers if indicated */
	if (!ign_header && !headers_match(&hdr1, &hdr2))
		return(1);
	lu_done(&hdr1); lu_done(&hdr2);
	if (!ign_header & (report >= REP_WARN)) {
		if (typ1 == TYP_UNKNOWN)
			printf("%s: warning - unrecognized format, comparing as binary\n",
					progname);
		if (lin1cnt != lin2cnt)
			printf("%s: warning - headers are different lengths\n",
					progname);
	}
	if (report >= REP_VERBOSE)
		printf("%s: input file type is %s\n",
				progname, file_type[typ1]);

	switch (typ1) {				/* compare based on type */
	case TYP_BINARY:
	case TYP_TMESH:
	case TYP_OCTREE:
	case TYP_RBFMESH:
	case TYP_UNKNOWN:
		return( !compare_binary() );
	case TYP_TEXT:
	case TYP_ASCII:
		return( !compare_text() );
	case TYP_RGBE:
	case TYP_XYZE:
		return( !compare_hdr() );
	case TYP_FLOAT:
		return( !compare_float() );
	case TYP_DOUBLE:
		return( !compare_double() );
	}
	return(1);
}
