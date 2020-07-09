#ifndef lint
static const char RCSid[] = "$Id: radcompare.c,v 2.25 2020/07/09 17:29:04 greg Exp $";
#endif
/*
 * Compare Radiance files for significant differences
 *
 *	G. Ward
 */

#include <stdlib.h>
#include <ctype.h>
#include "rtmath.h"
#include "platform.h"
#include "rtio.h"
#include "resolu.h"
#include "color.h"
#include "depthcodec.h"
#include "normcodec.h"
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

const char	nsuffix[10][3] = {		/* 1st, 2nd, 3rd, etc. */
			"th","st","nd","rd","th","th","th","th","th","th"
		};
#define	num_sfx(n)	nsuffix[(n)%10]

				/* file types */
const char	*file_type[] = {
			"Unrecognized",
			"TEXT_generic",
			"ascii",
			COLRFMT,
			CIEFMT,
			DEPTH16FMT,
			NORMAL32FMT,
			"float",
			"double",
			"BSDF_RBFmesh",
			"Radiance_octree",
			"Radiance_tmesh",
			"8-bit_indexed_name",
			"16-bit_indexed_name",
			"24-bit_indexed_name",
			"BINARY_unknown",
			NULL	/* terminator */
		};
				/* keep consistent with above */
enum {TYP_UNKNOWN, TYP_TEXT, TYP_ASCII, TYP_RGBE, TYP_XYZE,
		TYP_DEPTH, TYP_NORM, TYP_FLOAT, TYP_DOUBLE,
		TYP_RBFMESH, TYP_OCTREE, TYP_TMESH,
		TYP_ID8, TYP_ID16, TYP_ID24, TYP_BINARY};
		
#define has_header(t)	(!( 1L<<(t) & (1L<<TYP_TEXT | 1L<<TYP_BINARY) ))

				/* header variables to always ignore */
const char	*hdr_ignkey[] = {
			"SOFTWARE",
			"CAPDATE",
			"GMT",
			"FRAME",
			NULL	/* terminator */
		};
				/* header variable settings */
LUTAB	hdr1 = LU_SINIT(free,free);
LUTAB	hdr2 = LU_SINIT(free,free);

				/* advance appropriate file line count */
#define adv_linecnt(htp)	(lin1cnt += (htp == &hdr1), \
					lin2cnt += (htp == &hdr2))

typedef struct {		/* dynamic line buffer */
	char	*str;
	int	len;
	int	siz;
} LINEBUF;

#define init_line(bp)	((bp)->str = NULL, (bp)->siz = 0)
				/* 100 MByte limit on line buffer */
#define MAXBUF		(100L<<20)

				/* input files */
char		*progname = NULL;
const char	stdin_name[] = "<stdin>";
const char	*f1name=NULL, *f2name=NULL;
FILE		*f1in=NULL, *f2in=NULL;
int		f1swap=0, f2swap=0;

				/* running real differences */
double	diff2sum = 0;
long	nsum = 0;

/* Report usage and exit */
static void
usage()
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-h][-s|-w|-v][-rel min_test][-rms epsilon][-max epsilon] reference test\n",
			stderr);
	exit(2);
}

/* Read a text line, increasing buffer size as necessary */
static int
read_line(LINEBUF *bp, FILE *fp)
{
	static int	doneWarn = 0;

	bp->len = 0;
	if (!bp->str) {
		bp->str = (char *)malloc(bp->siz = 512);
		if (!bp->str)
			goto memerr;
	}
	while (fgets(bp->str + bp->len, bp->siz - bp->len, fp)) {
		bp->len += strlen(bp->str + bp->len);
		if (bp->str[bp->len-1] == '\n')
			break;		/* found EOL */
		if (bp->len < bp->siz - 4)
			continue;	/* at EOF? */
		if (bp->siz >= MAXBUF) {
			if ((report >= REP_WARN) & !doneWarn) {
				fprintf(stderr,
			"%s: warning - input line(s) past %ld MByte limit\n",
					progname, MAXBUF>>20);
				doneWarn++;
			}
			break;		/* return MAXBUF partial line */
		}
		if ((bp->siz += bp->siz/2) > MAXBUF)
			bp->siz = MAXBUF;
		bp->str = (char *)realloc(bp->str, bp->siz);
		if (!bp->str)
			goto memerr;
	}
	return(bp->len);
memerr:
	fprintf(stderr,
		"%s: out of memory in read_line() allocating %d byte buffer\n",
			progname, bp->siz);
	exit(2);
}

/* Free line buffer */
static void
free_line(LINEBUF *bp)
{
	if (bp->str) free(bp->str);
	init_line(bp);
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
			"%s: %sdifference between %.8g and %.8g exceeds epsilon of %.8g\n",
					progname,
					(rel_min > 0) ? "relative " : "",
					r1, r2, max_lim);
		return(0);
	}
	diff2sum += diff2;
	nsum++;
	return(1);
}

/* Compare two color values for equivalence */
static int
color_check(COLOR c1, COLOR c2)
{
	int	p;

	if (!real_check((colval(c1,RED)+colval(c1,GRN)+colval(c1,BLU))*(1./3.),
			(colval(c2,RED)+colval(c2,GRN)+colval(c2,BLU))*(1./3.)))
		return(0);

	p = (colval(c1,GRN) > colval(c1,RED)) ? GRN : RED;
	if (colval(c1,BLU) > colval(c1,p)) p = BLU;
	
	return(real_check(colval(c1,p), colval(c2,p)));
}

/* Compare two normal directions for equivalence */
static int
norm_check(FVECT nv1, FVECT nv2)
{
	double	max2 = nv1[2]*nv2[2];
	int	imax = 2;
	int	i = 2;
					/* identify largest component */
	while (i--) {
		double	tm2 = nv1[i]*nv2[i];
		if (tm2 > max2) {
			imax = i;
			max2 = tm2;
		}
	}
	i = 3;				/* compare smaller components */
	while (i--) {
		if (i == imax)
			continue;
		if (!real_check(nv1[i], nv2[i]))
			return(0);
	}
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
	char	newval[128];
	LUTAB	*htp = (LUTAB *)p;
	LUENT	*tep;
	char	*key;
	int	kln, vln;
	int	n;

	adv_linecnt(htp);	/* side-effect is to count lines */
	if (!isalpha(*val))	/* key must start line */
		return(0);
				/* check if we need to swap binary data */
	if ((n = isbigendian(val)) >= 0) {
		if (nativebigendian() == n)
			return(0);
		f1swap += (htp == &hdr1);
		f2swap += (htp == &hdr2);
		return(0);
	}
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
				/* check if key to ignore */
	for (n = 0; hdr_ignkey[n]; n++)
		if (!strcmp(key, hdr_ignkey[n]))
			return(0);
	vln = strlen(val);	/* eliminate space and newline at end */
	while (isspace(val[--vln]))
		;
	val[++vln] = '\0';
	if (!(tep = lu_find(htp, key)))
		return(-1);	/* memory allocation error */
	if (!tep->key)
		tep->key = strcpy(malloc(kln+1), key);
	if (tep->data) {	/* check for special cases */
		if (!strcmp(key, "EXPOSURE")) {
			sprintf(newval, "%f", atof(tep->data)*atof(val));
			vln = strlen(val = newval);
		}
		free(tep->data);
	}
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
			printf("%s: variable '%s' missing in '%s'\n",
					progname, ep1->key, f2name);
		return(-1);
	}
	if (!equiv_string((char *)ep1->data, (char *)ep2->data)) {
		if (report != REP_QUIET) {
			printf("%s: header variable '%s' has different values\n",
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
headers_match()
{
	int	ne = lu_doall(&hdr1, match_val, &hdr2);
	if (ne < 0)
		return(0);	/* something didn't match! */
				/* non-fatal if second header has extra */
	if (report >= REP_WARN && (ne = lu_doall(&hdr2, NULL, NULL) - ne))
		printf("%s: warning - '%s' has %d extra header setting(s)\n",
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
		++n;
		if (!c | (c > 127))
			break;			/* non-ascii character */
		if (n >= 10240)
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
			fputs(": unknown error reading header\n", stderr);
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
	return(TYP_TEXT);
badeof:
	if (report != REP_QUIET) {
		fputs(name, stdout);
		fputs(": unexpected end-of-file\n", stdout);
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
	"%s: %sRMS difference between '%s' and '%s' of %.5g exceeds limit of %.5g\n",
					progname,
					(rel_min > 0) ? "relative " : "",
					f1name, f2name,
					sqrt(diff2sum/(double)nsum), rms_lim);
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
	int	c1=0, c2=0;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as binary\n", stdout);
	}
	for ( ; ; ) {				/* exact byte matching */
		c1 = getc(f1in);
		c2 = getc(f2in);
		if (c1 == EOF) {
			if (c2 == EOF)
				return(1);	/* success! */
			if (report != REP_QUIET) {
				fputs(f1name, stdout);
				fputs(": unexpected end-of-file\n", stdout);
			}
			return(0);
		}
		if (c2 == EOF) {
			if (report != REP_QUIET) {
				fputs(f2name, stdout);
				fputs(": unexpected end-of-file\n", stdout);
			}
			return(0);
		}
		if (c1 != c2)
			break;			/* quit and report difference */
	}
	if (report == REP_QUIET)
		return(0);
	printf("%s: binary files '%s' and '%s' differ at byte offset %ld|%ld\n",
			progname, f1name, f2name, ftell(f1in), ftell(f2in));
	if (report >= REP_VERBOSE)
		printf("%s: byte in '%s' is 0x%X, byte in '%s' is 0x%X\n",
				progname, f1name, c1, f2name, c2);
	return(0);
}

/* Compare two inputs as generic text files */
static int
compare_text()
{
	LINEBUF	l1buf, l2buf;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as ASCII text\n", stdout);
	}
	SET_FILE_TEXT(f1in);			/* originally set to binary */
	SET_FILE_TEXT(f2in);
	init_line(&l1buf); init_line(&l2buf);	/* compare a line at a time */
	while (read_line(&l1buf, f1in)) {
		lin1cnt++;
		if (!*sskip2(l1buf.str,0))
			continue;		/* ignore empty lines */

		while (read_line(&l2buf, f2in)) {
			lin2cnt++;
			if (*sskip2(l2buf.str,0))
				break;		/* found other non-empty line */
		}
		if (!l2buf.len) {		/* input 2 EOF? */
			if (report != REP_QUIET) {
				fputs(f2name, stdout);
				fputs(": unexpected end-of-file\n", stdout);
			}
			free_line(&l1buf); free_line(&l2buf);
			return(0);
		}
						/* compare non-empty lines */
		if (!equiv_string(l1buf.str, l2buf.str)) {
			if (report != REP_QUIET) {
				printf("%s: inputs '%s' and '%s' differ at line %d|%d\n",
						progname, f1name, f2name,
						lin1cnt, lin2cnt);
				if ( report >= REP_VERBOSE &&
						(l1buf.len < 256) &
						(l2buf.len < 256) ) {
					fputs("------------- Mismatch -------------\n", stdout);
					printf("%s@%d:\t%s", f1name,
							lin1cnt, l1buf.str);
					printf("%s@%d:\t%s", f2name,
							lin2cnt, l2buf.str);
				}
			}
			free_line(&l1buf); free_line(&l2buf);
			return(0);
		}
	}
	free_line(&l1buf);			/* check for EOF on input 2 */
	while (read_line(&l2buf, f2in)) {
		if (!*sskip2(l2buf.str,0))
			continue;
		if (report != REP_QUIET) {
			fputs(f1name, stdout);
			fputs(": unexpected end-of-file\n", stdout);
		}
		free_line(&l2buf);
		return(0);
	}
	free_line(&l2buf);
	return(good_RMS());			/* final check for reals */
}

/* Check image/map resolutions */
static int
check_resolu(const char *class, RESOLU *r1p, RESOLU *r2p)
{
	if (r1p->rt != r2p->rt) {
		if (report != REP_QUIET)
			printf(
			"%s: %ss '%s' and '%s' have different pixel ordering\n",
					progname, class, f1name, f2name);
		return(0);
	}
	if ((r1p->xr != r2p->xr) | (r1p->yr != r2p->yr)) {
		if (report != REP_QUIET)
			printf(
			"%s: %ss '%s' and '%s' are different sizes\n",
					progname, class, f1name, f2name);
		return(0);
	}
	return(1);
}

/* Compare two inputs that are known to be RGBE or XYZE images */
static int
compare_hdr()
{
	RESOLU	rs1, rs2;
	COLOR	*scan1, *scan2;
	int	x, y;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as HDR images\n", stdout);
	}
	fgetsresolu(&rs1, f1in);
	fgetsresolu(&rs2, f2in);
	if (!check_resolu("HDR image", &rs1, &rs2))
		return(0);
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
				printf("%s: unexpected end-of-file\n",
						progname);
			free(scan1);
			free(scan2);
			return(0);
		}
		for (x = 0; x < scanlen(&rs1); x++) {
			if (color_check(scan1[x], scan2[x]))
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

/* Set reference depth based on header variable */
static int
set_refdepth(DEPTHCODEC *dcp, LUTAB *htp)
{
	static char	depthvar[] = DEPTHSTR;
	const char	*drval;

	depthvar[LDEPTHSTR-1] = '\0';
	drval = (const char *)lu_find(htp, depthvar)->data;
	if (!drval)
		return(0);
	dcp->refdepth = atof(drval);
	if (dcp->refdepth <= 0) {
		if (report != REP_QUIET) {
			fputs(dcp->inpname, stderr);
			fputs(": bad reference depth '", stderr);
			fputs(drval, stderr);
			fputs("'\n", stderr);
		}
		return(-1);
	}
	return(1);
}

/* Compare two encoded depth maps */
static int
compare_depth()
{
	long		nread = 0;
	DEPTHCODEC	dc1, dc2;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as depth maps\n", stdout);
	}
	set_dc_defaults(&dc1);
	dc1.hdrflags = HF_RESIN;
	dc1.finp = f1in;
	dc1.inpname = f1name;
	set_dc_defaults(&dc2);
	dc2.hdrflags = HF_RESIN;
	dc2.finp = f2in;
	dc2.inpname = f2name;
	if (report != REP_QUIET) {
		dc1.hdrflags |= HF_STDERR;
		dc2.hdrflags |= HF_STDERR;
	}
	if (!process_dc_header(&dc1, 0, NULL))
		return(0);
	if (!process_dc_header(&dc2, 0, NULL))
		return(0);
	if (!check_resolu("Depth map", &dc1.res, &dc2.res))
		return(0);
	if (set_refdepth(&dc1, &hdr1) < 0)
		return(0);
	if (set_refdepth(&dc2, &hdr2) < 0)
		return(0);
	while (nread < dc1.res.xr*dc1.res.yr) {
		double	d1 = decode_depth_next(&dc1);
		double	d2 = decode_depth_next(&dc2);
		if ((d1 < 0) | (d2 < 0)) {
			if (report != REP_QUIET)
				printf("%s: unexpected end-of-file\n",
						progname);
			return(0);
		}
		++nread;
		if (real_check(d1, d2))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s depth values differ\n",
					progname, nread, num_sfx(nread));
		return(0);
	}
	return(good_RMS());		/* final check of RMS */
}

/* Compare two encoded normal maps */
static int
compare_norm()
{
	long		nread = 0;
	NORMCODEC	nc1, nc2;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as normal maps\n", stdout);
	}
	set_nc_defaults(&nc1);
	nc1.hdrflags = HF_RESIN;
	nc1.finp = f1in;
	nc1.inpname = f1name;
	set_nc_defaults(&nc2);
	nc2.hdrflags = HF_RESIN;
	nc2.finp = f2in;
	nc2.inpname = f2name;
	if (report != REP_QUIET) {
		nc1.hdrflags |= HF_STDERR;
		nc2.hdrflags |= HF_STDERR;
	}
	if (!process_nc_header(&nc1, 0, NULL))
		return(0);
	if (!process_nc_header(&nc2, 0, NULL))
		return(0);
	if (!check_resolu("Normal map", &nc1.res, &nc2.res))
		return(0);
	while (nread < nc1.res.xr*nc1.res.yr) {
		FVECT	nv1, nv2;
		int	rv1 = decode_normal_next(nv1, &nc1);
		int	rv2 = decode_normal_next(nv2, &nc2);
		if ((rv1 < 0) | (rv2 < 0)) {
			if (report != REP_QUIET)
				printf("%s: unexpected end-of-file\n",
						progname);
			return(0);
		}
		++nread;
		if (rv1 == rv2 && (!rv1 || norm_check(nv1, nv2)))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s normal vectors differ\n",
					progname, nread, num_sfx(nread));
		return(0);
	}
	return(good_RMS());		/* final check of RMS */
}

/* Compare two inputs that are known to be 32-bit floating-point data */
static int
compare_float()
{
	long	nread = 0;
	float	f1, f2;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as 32-bit IEEE floats\n", stdout);
	}
	while (getbinary(&f1, sizeof(f1), 1, f1in)) {
		if (!getbinary(&f2, sizeof(f2), 1, f2in))
			goto badeof;
		++nread;
		if (f1swap) swap32((char *)&f1, 1);
		if (f2swap) swap32((char *)&f2, 1);
		if (real_check(f1, f2))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s float values differ\n",
					progname, nread, num_sfx(nread));
		return(0);
	}
	if (!getbinary(&f2, sizeof(f2), 1, f2in))
		return(good_RMS());		/* final check of RMS */
badeof:
	if (report != REP_QUIET)
		printf("%s: unexpected end-of-file\n", progname);
	return(0);
}

/* Compare two inputs that are known to be 64-bit floating-point data */
static int
compare_double()
{
	long	nread = 0;
	double	f1, f2;

	if (report >= REP_VERBOSE) {
		fputs(progname, stdout);
		fputs(": comparing inputs as 64-bit IEEE doubles\n", stdout);
	}
	while (getbinary(&f1, sizeof(f1), 1, f1in)) {
		if (!getbinary(&f2, sizeof(f2), 1, f2in))
			goto badeof;
		++nread;
		if (f1swap) swap64((char *)&f1, 1);
		if (f2swap) swap64((char *)&f2, 1);
		if (real_check(f1, f2))
			continue;
		if (report != REP_QUIET)
			printf("%s: %ld%s double values differ\n",
					progname, nread, num_sfx(nread));
		return(0);
	}
	if (!getbinary(&f2, sizeof(f2), 1, f2in))
		return(good_RMS());		/* final check of RMS */
badeof:
	if (report != REP_QUIET)
		printf("%s: unexpected end-of-file\n", progname);
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
	if (a != argc-2)			/* make sure of two inputs */
		usage();
	if (!strcmp(argv[a], argv[a+1])) {	/* inputs are same? */
		if (report >= REP_WARN)
			printf("%s: warning - identical inputs given\n",
					progname);
		return(0);
	}
	if (!f1name) f1name = argv[a];
	if (!f2name) f2name = argv[a+1];
						/* open inputs */
	SET_FILE_BINARY(stdin);			/* in case we're using it */
	if (!f1in && !(f1in = fopen(f1name, "rb"))) {
		fprintf(stderr, "%s: cannot open for reading\n", f1name);
		return(2);
	}
	if (!strcmp(f2name, "-")) {
		f2in = stdin;
		f2name = stdin_name;
	} else if (!(f2in = fopen(f2name, "rb"))) {
		fprintf(stderr, "%s: cannot open for reading\n", f2name);
		return(2);
	}
						/* load headers */
	if ((typ1 = identify_type(f1name, f1in, &hdr1)) < 0)
		return(2);
	if ((typ2 = identify_type(f2name, f2in, &hdr2)) < 0)
		return(2);
	if (typ1 != typ2) {
		if (report != REP_QUIET)
			printf("%s: '%s' format is %s and '%s' is %s\n",
					progname, f1name, file_type[typ1],
					f2name, file_type[typ2]);
		return(1);
	}
	ign_header |= !has_header(typ1);	/* check headers if indicated */
	if (!ign_header && !headers_match())
		return(1);
	if (!ign_header & (report >= REP_WARN)) {
		if (lin1cnt != lin2cnt)
			printf("%s: warning - headers are different lengths\n",
					progname);
		if (typ1 == TYP_UNKNOWN)
			printf("%s: warning - unrecognized format\n",
					progname);
	}
	if (report >= REP_VERBOSE) {
		printf("%s: data format is %s\n", progname, file_type[typ1]);
		if ((typ1 == TYP_FLOAT) | (typ1 == TYP_DOUBLE)) {
			if (f1swap)
				printf("%s: input '%s' is byte-swapped\n",
						progname, f1name);
			if (f2swap)
				printf("%s: input '%s' is byte-swapped\n",
						progname, f2name);
		}
	}
	switch (typ1) {				/* compare based on type */
	case TYP_BINARY:
	case TYP_TMESH:
	case TYP_OCTREE:
	case TYP_RBFMESH:
	case TYP_ID8:
	case TYP_ID16:
	case TYP_ID24:
	case TYP_UNKNOWN:
		return( !compare_binary() );
	case TYP_TEXT:
	case TYP_ASCII:
		return( !compare_text() );
	case TYP_RGBE:
	case TYP_XYZE:
		return( !compare_hdr() );
	case TYP_DEPTH:
		return( !compare_depth() );
	case TYP_NORM:
		return( !compare_norm() );
	case TYP_FLOAT:
		return( !compare_float() );
	case TYP_DOUBLE:
		return( !compare_double() );
	}
	return(1);
}
