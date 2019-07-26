#ifndef lint
static const char RCSid[] = "$Id: rcode_ident.c,v 2.9 2019/07/26 18:38:39 greg Exp $";
#endif
/*
 * Create or read identifier index map
 *
 */

#include <stdlib.h>
#include <ctype.h>
#include "rtio.h"
#include "platform.h"
#include "lookup.h"
#include "idmap.h"

#ifndef MAXIDLEN
#define	MAXIDLEN	256		/* longest ID length supported */
#endif

#define HF_TABLOUT	0x100		/* output ID table after header */

char		*progname;		/* global argv[0] */

static int	unbuffered = 0;
static int	numeric = 0;
static int	sepc = '\n';


/* Report usage error and exit */
static void
usage_exit(int code)
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-8|-16|-24][-h][-tS][-x xr -y yr] [input.txt [output.idx]]\n",
			stderr);
	fputs("   Or: ", stderr);
	fputs(progname, stderr);
	fputs(" -r [-i][-u][-h][-H][-l][-n][-tS] input.idx [output.txt]\n",
			stderr);
	exit(code);
}


/* read identifier from input, trimming white space on either end */
static int
scan_ident(char ident[MAXIDLEN], FILE *fp)
{
	char	*cp = ident;
	int	c;

	while ((c = getc(fp)) != EOF) {
		if (cp == ident && isspace(c))
			continue;
		if (c == sepc) {
			while (cp > ident && isspace(cp[-1]))
				cp--;
			*cp = '\0';
			return cp - ident;	/* return string length */
		}
		if (cp - ident >= MAXIDLEN-1) {
			fputs(progname, stderr);
			fputs(": identifier too long in scan_ident()\n", stderr);
			return -1;
		}
		*cp++ = c;
	}
	return -1;
}


/* elide format from header lines */
static int
headline(char *s, void *p)
{
	if (isformat(s))
		return 0;

	fputs(s, stdout);
	return 1;
}


/* Create an identifier index from the given input file (stdin if NULL) */
int
create_index(const char *fname, int hdrflags, int ndxbytes, int xres, int yres)
{
	FILE	*fp = stdin;
	char	**idmap;
	int	idmlen;
	int	nextID = 0;
	LUTAB	hashtab;
	RESOLU	rs;
	long	n;
	int	ndx;
					/* open file if not stdin */
	if (!fname)
		fname = "<stdin>";
	else if (!(fp = fopen(fname, "r"))) {
		fputs(fname, stderr);
		fputs(": cannot open for input\n", stderr);
		return 0;
	}
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(fp);
#endif
					/* load/copy header */
	if (hdrflags & HF_HEADOUT && getheader(fp, headline, NULL) < 0) {
		fputs(fname, stderr);
		fputs(": bad header or wrong input format\n", stderr);
		fclose(fp);
		return 0;
	}
	if (!(hdrflags & HF_HEADOUT))
		newheader("RADIANCE", stdout);
					/* finish & allocate ID table */
	switch (ndxbytes) {
	case 1:
		fputformat(IDMAP8FMT, stdout);
		idmap = (char **)malloc((idmlen=256)*sizeof(char *));
		break;
	case 2:
		fputformat(IDMAP16FMT, stdout);
		idmap = (char **)malloc((idmlen=4096)*sizeof(char *));
		break;
	case 3:
		fputformat(IDMAP24FMT, stdout);
		idmap = (char **)malloc((idmlen=1L<<17)*sizeof(char *));
		break;
	default:
		fputs(progname, stderr);
		fputs(": unsupported bits/pixel\n", stderr);
		return 0;
	}
	memset(&hashtab, 0, sizeof(hashtab));
	hashtab.hashf = lu_shash;
	hashtab.keycmp = strcmp;
	hashtab.freek = free;
	if (!idmap || !lu_init(&hashtab, idmlen))
		goto memerr;
	fputc('\n', stdout);		/* end of info header */

	rs.rt = PIXSTANDARD;		/* get/write map resolution */
	if (hdrflags & HF_RESOUT) {
		if (!fgetsresolu(&rs, fp)) {
			fputs(fname, stderr);
			fputs(": bad resolution string\n", stderr);
			return 0;
		}
	} else if (((rs.xr = xres) <= 0) | ((rs.yr = yres) <= 0)) {
		fputs(fname, stderr);
		fputs(": missing/illegal resolution\n", stderr);
		return 0;
	}
	fputsresolu(&rs, stdout);
					/* read/assign/write each ID */
	for (n = (long)rs.xr*rs.yr; n-- > 0; ) {
		LUENT	*lep;
		char	ident[MAXIDLEN];
		int	len = scan_ident(ident, fp);
		if (len < 0) {
			if (feof(fp)) {
				fputs(fname, stderr);
				fputs(": unexpected EOF\n", stderr);
			}
			return 0;
		}
		lep = lu_find(&hashtab, ident);
		if (!lep->key) {			/* new entry */
			if (nextID >= 1<<(ndxbytes<<3)) {
				fputs(progname, stderr);
				fputs(
			": too many unique identifiers for index size\n",
						stderr);
				return 0;
			}
			if (nextID >= idmlen) {		/* reallocate? */
				idmlen += idmlen>>2;	/* +25% growth */
				idmap = (char **)realloc(idmap, idmlen);
				if (!idmap)
					goto memerr;
			}
			lep->key = (char *)malloc(len+1);
			if (!lep->key)
				goto memerr;
			idmap[nextID] = strcpy(lep->key, ident);
			lep->data = (char *)nextID++;
		}
						/* write this id index */
		putint((int)lep->data, ndxbytes, stdout);
	}
	while ((ndx = getc(fp)) != EOF)		/* verify we got them all */
		if (!isspace(ndx)) {
			fputs(fname, stderr);
			fputs(": warning -  non-white characters past end\n",
					stderr);
			break;
		}
	if (fp != stdin)			/* done with input */
		fclose(fp);
	for (ndx = 0; ndx < nextID; ndx++) {	/* append string table */
		fputs(idmap[ndx], stdout);
		putchar('\0');			/* nul-terminated IDs */
	}
	free(idmap);				/* clean up */
	lu_done(&hashtab);
	return 1;
memerr:
	fputs(progname, stderr);
	fputs(": out of memory!\n", stderr);
	return 0;
}


/* print out ID table */
void
print_IDs(IDMAP *idmp)
{
	int	i;

	printf("============ %d IDs ============\n", idmp->nids);

	for (i = 0; i < idmp->nids; i++) {
		fputs(mapID(idmp, i), stdout);
		putchar(sepc);
	}
	if (sepc != '\n')
		fputc('\n', stdout);
	puts("============= END =============");
}


/* Load selected pixels from identifier index file */
int
decode_select(const char *fname, int hdrflags)
{
	IDMAP	*idmp = idmap_ropen(fname, hdrflags);
	int	i, j;

	if (!idmp)
		return 0;

	if (hdrflags & HF_TABLOUT)
		print_IDs(idmp);

	if (idmp->res.rt != PIXSTANDARD) {
		fputs(progname, stderr);
		fputs(": can only handle standard pixel ordering\n", stderr);
		idmap_close(idmp);
		return 0;
	}
	while (scanf("%d %d", &i, &j) == 2) {
		i = idmap_seek(idmp, i, idmp->res.yr-1 - j);
		if (!i) {
			fputs(progname, stderr);
			fputs(": warning - pixel index is off map\n", stderr);
			continue;
		}
		if (i > 0)
			i = idmap_next_i(idmp);
		if (i < 0) {
			fputs(fname, stderr);
			fputs(": read/seek error in decode_select()\n", stderr);
			idmap_close(idmp);
			return 0;
		}
		if (numeric) {
			printf("%d", i);
		} else {
			const char	*id = mapID(idmp, i);
			if (!id) {
				fputs(fname, stderr);
				fputs(": bad ID index in file\n", stderr);
				idmap_close(idmp);
				return 0;
			}
			fputs(id, stdout);
		}
		putchar(sepc);
		if (unbuffered && fflush(stdout) == EOF) {
			fputs(progname, stderr);
			fputs(": write error on output\n", stderr);
			idmap_close(idmp);
			return 0;
		}
	}
	idmap_close(idmp);

	if (!feof(stdin)) {
		fputs(progname, stderr);
		fputs(": non-numeric input for pixel index\n", stderr);
		return 0;
	}
	return 1;
}


/* Load and convert identifier index file to ASCII */
int
decode_all(const char *fname, int hdrflags)
{
	IDMAP	*idmp = idmap_ropen(fname, hdrflags);
	long	n;

	if (!idmp)
		return 0;

	if (hdrflags & HF_TABLOUT)
		print_IDs(idmp);

	for (n = idmp->res.xr*idmp->res.yr; n-- > 0; ) {
		const int	ndx = idmap_next_i(idmp);
		if (ndx < 0) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			idmap_close(idmp);
			return 0;
		}
		if (numeric) {
			printf("%d", ndx);
		} else {
			const char	*id = mapID(idmp, ndx);
			if (!id) {
				fputs(fname, stderr);
				fputs(": bad ID index in file\n", stderr);
				idmap_close(idmp);
				return 0;
			}
			fputs(id, stdout);
		}
		putchar(sepc);
	}
	idmap_close(idmp);
	return 1;
}


int
main(int argc, char *argv[])
{
	char	*inpfname = NULL;
	int	reverse = 0;
	int	bypixel = 0;
	int	hdrflags = HF_HEADOUT|HF_RESOUT|HF_STDERR;
	int	xres=0, yres=0;
	int	ndxbytes = 2;
	int	a;
	
	progname = argv[0];
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'r':
			reverse++;
			break;
		case 'i':
			bypixel++;
			break;
		case 'u':
			unbuffered++;
			break;
		case 'n':
			numeric++;
			break;
		case 'l':
			hdrflags |= HF_TABLOUT;
			break;
		case 'h':
			hdrflags &= ~HF_HEADOUT;
			break;
		case 'H':
			hdrflags &= ~HF_RESOUT;
			break;
		case 't':
			if (!(sepc = argv[a][2]))
				sepc = '\n';
			break;
		case 'x':
			xres = atoi(argv[++a]);
			break;
		case 'y':
			yres = atoi(argv[++a]);
			break;
		case '8':
		case '1':
		case '2':
			switch (atoi(argv[a]+1)) {
			case 8:
				ndxbytes = 1;
				break;
			case 16:
				ndxbytes = 2;
				break;
			case 24:
				ndxbytes = 3;
				break;
			default:
				usage_exit(1);
			}
			break;
		default:
			usage_exit(1);
		}
	if ((xres > 0) & (yres > 0))	/* user specified resolution? */
		hdrflags &= ~HF_RESOUT;
	if (reverse && a >= argc) {
		fputs(progname, stderr);
		fputs(": -r option requires named input file\n", stderr);
		return 1;
	}
	if (a < argc-2) {
		fputs(progname, stderr);
		fputs(": too many file arguments\n", stderr);
		return 1;
	}
	if (a < argc)
		inpfname = argv[a];
	if (a+1 < argc && !freopen(argv[a+1], "w", stdout)) {
		fputs(argv[a+1], stderr);
		fputs(": cannot open output\n", stderr);
		return 1;
	}
	if (!reverse)
		SET_FILE_BINARY(stdout);
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(stdout);
#endif
	if (reverse) {
		if (!(bypixel ? decode_select(inpfname, hdrflags) :
				decode_all(inpfname, hdrflags)))
			return 1;
	} else if (!create_index(inpfname, hdrflags, ndxbytes, xres, yres))
		return 1;

	if (fflush(stdout) == EOF) {
		fputs(progname, stderr);
		fputs(": error writing output\n", stderr);
		return 1;
	}
	return 0;
}
