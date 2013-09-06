#ifndef lint
static const char RCSid[] = "$Id: rcollate.c,v 2.4 2013/09/06 21:34:39 greg Exp $";
#endif
/*
 * Utility to re-order records in a binary or ASCII data file (matrix)
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "platform.h"
#include "rtio.h"
#include "resolu.h"
#ifndef _WIN32
#include <sys/mman.h>
#endif

#ifdef getc_unlocked            /* avoid horrendous overhead of flockfile */
#undef getc
#undef putc
#define getc    getc_unlocked
#define putc    putc_unlocked
#endif

typedef struct {
	void	*base;		/* pointer to base memory */
	size_t	len;		/* allocated memory length */
	int	mapped;		/* memory-mapped file? */
} MEMLOAD;		/* file loaded/mapped into memory */

typedef struct {
	int	nw_rec;		/* number of words per record */
	int	nrecs;		/* number of records we found */
	char	*rec[1];	/* record array (extends struct) */
} RECINDEX;

int		warnings = 1;	/* report warnings? */

/* free loaded file */
static void
free_load(MEMLOAD *mp)
{
	if (mp == NULL || (mp->base == NULL) | (mp->len <= 0))
		return;
#ifdef MAP_FILE
	if (mp->mapped)
		munmap(mp->base, mp->len);
	else
#endif
		free(mp->base);
	mp->base = NULL;
	mp->len = 0;
}

/* load a file into memory */
static int
load_file(MEMLOAD *mp, FILE *fp)
{
	int	fd;
	off_t	skip, flen;

	if (mp == NULL)
		return(-1);
	mp->base = NULL;
	mp->len = 0;
	mp->mapped = 0;
	if (fp == NULL)
		return(-1);
	fd = fileno(fp);
	skip = ftello(fp);
	flen = lseek(fd, 0, SEEK_END);
	if (flen <= skip)
		return((int)(flen - skip));
	mp->len = (size_t)(flen - skip);
#ifdef MAP_FILE
	if (mp->len > 1L<<20) {		/* map file if > 1 MByte */
		mp->base = mmap(NULL, mp->len, PROT_READ, MAP_PRIVATE, fd, skip);
		if (mp->base != MAP_FAILED) {
			mp->mapped = 1;
			return(1);	/* mmap() success */
		}
		mp->base = NULL;	/* fall back to reading it in... */
	}
#endif
	if (lseek(fd, skip, SEEK_SET) != skip ||
			(mp->base = malloc(mp->len)) == NULL) {
		mp->len = 0;
		return(-1);
	}
	if (read(fd, (char *)mp->base, mp->len) != mp->len) {
		free_load(mp);
		return(-1);
	}
	return(1);
}

/* load memory from an input stream, starting from current position */
static int
load_stream(MEMLOAD *mp, FILE *fp)
{
	char	buf[8192];
	size_t	nr;

	if (mp == NULL)
		return(-1);
	mp->base = NULL;
	mp->len = 0;
	mp->mapped = 0;
	if (fp == NULL)
		return(-1);
	while ((nr = fread(buf, 1, sizeof(buf), fp)) > 0) {
		if (!mp->len)
			mp->base = malloc(nr);
		else
			mp->base = realloc(mp->base, mp->len+nr);
		if (mp->base == NULL)
			return(-1);
		memcpy((char *)mp->base + mp->len, buf, nr);
		mp->len += nr;
	}
	if (ferror(fp)) {
		free_load(mp);
		return(-1);
	}
	return(mp->len > 0);
}

/* free a record index */
#define free_records(rp)	free(rp)

/* compute record index */
static RECINDEX *
index_records(const MEMLOAD *mp, int nw_rec)
{
	RECINDEX	*rp;
	char		*cp, *mend;
	int		n;

	if (mp == NULL || (mp->base == NULL) | (mp->len <= 0))
		return(NULL);
	if (nw_rec <= 0)
		return(NULL);
	rp = (RECINDEX *)malloc(sizeof(RECINDEX) + mp->len/(2*nw_rec)*sizeof(char *));
	if (rp == NULL)
		return(NULL);
	rp->nw_rec = nw_rec;
	rp->nrecs = 0;
	cp = (char *)mp->base;
	mend = cp + mp->len;
	for ( ; ; ) {				/* whitespace-separated words */
		while (cp < mend && !*cp | isspace(*cp))
			++cp;
		if (cp >= mend)
			break;
		rp->rec[rp->nrecs++] = cp;	/* point to first non-white */
		n = rp->nw_rec;
		while (++cp < mend)		/* find end of record */
			if (!*cp | isspace(*cp)) {
				if (--n <= 0)
					break;	/* got requisite # words */
				do {		/* else find next word */
					if (*cp == '\n') {
						fprintf(stderr,
						"Unexpected EOL in record!\n");
						free_records(rp);
						return(NULL);
					}
					if (++cp >= mend)
						break;
				} while (!*cp | isspace(*cp));
			}
	}
	rp->rec[rp->nrecs] = mend;		/* reallocate to save space */
	rp = (RECINDEX *)realloc(rp,
			sizeof(RECINDEX) + rp->nrecs*sizeof(char *));
	return(rp);
}

/* count number of columns based on first EOL */
static int
count_columns(const RECINDEX *rp)
{
	char	*cp = rp->rec[0];
	char	*mend = rp->rec[rp->nrecs];
	int	i;

	while (*cp != '\n')
		if (++cp >= mend)
			return(0);
	for (i = 0; i < rp->nrecs; i++)
		if (rp->rec[i] >= cp)
			break;
	return(i);
}

/* copy nth record from index to stdout */
static int
print_record(const RECINDEX *rp, int n)
{
	int	words2go = rp->nw_rec;
	char	*scp;

	if ((n < 0) | (n >= rp->nrecs))
		return(0);
	scp = rp->rec[n];
	do {
		putc(*scp++, stdout);
		if (!*scp | isspace(*scp)) {
			if (--words2go <= 0)
				break;
			putc(' ', stdout);	/* single space btwn. words */
			do
				if (++scp >= rp->rec[n+1])
					break;
			while (!*scp | isspace(*scp));
		}
	} while (scp < rp->rec[n+1]);
						/* caller adds record sep. */
	return(1);
}

/* copy a stream to stdout */
static int
output_stream(FILE *fp)
{
	char	buf[8192];
	ssize_t	n;

	if (fp == NULL)
		return(0);
	fflush(stdout);			/* assumes nothing in input buffer */
	while ((n = read(fileno(fp), buf, sizeof(buf))) > 0)
		if (write(fileno(stdout), buf, n) != n)
			return(0);
	return(n >= 0);
}

/* get next word from stream, leaving stream on EOL or start of next word */
static char *
fget_word(char buf[256], FILE *fp)
{
	int	c;
	char	*cp;
					/* skip nul's and white space */
	while (!(c = getc(fp)) || isspace(c))
		;
	if (c == EOF)
		return(NULL);
	cp = buf;
	do
		*cp++ = c;
	while ((c = getc(fp)) != EOF && !isspace(c) && cp < buf+255);
	*cp = '\0';
	while (isspace(c) & (c != '\n'))
		c = getc(fp);
	if (c != EOF)
		ungetc(c, fp);
	return(buf);
}

char		*fmtid = "ascii";		/* format id */
int		record_width = 3;		/* words/record (<0 binary) */
int		ni_columns = 0;			/* number of input columns */
int		ni_rows = 0;			/* number of input rows */
int		no_columns = 0;			/* number of output columns */
int		no_rows = 0;			/* number of output rows */

/* output transposed ASCII or binary data from memory */
static int
do_transpose(const MEMLOAD *mp)
{
	static const char	tabEOL[2] = {'\t','\n'};
	RECINDEX		*rp = NULL;
	long			nrecords;
	int			i, j;
						/* get # records (& index) */
	if (record_width > 0) {
		if ((rp = index_records(mp, record_width)) == NULL)
			return(0);
		if (ni_columns <= 0)
			ni_columns = count_columns(rp);
		nrecords = rp->nrecs;
	} else if ((ni_rows > 0) & (ni_columns > 0)) {
		nrecords = ni_rows*ni_columns;
		if (nrecords > mp->len / -record_width) {
			fprintf(stderr,
			    "Input too small for specified size and type\n");
			return(0);
		}
	} else
		nrecords = mp->len / -record_width;
						/* check sizes */
	if (ni_rows <= 0)
		ni_rows = no_columns;
	if (ni_columns <= 0)
		ni_columns = no_rows;
	if ((ni_rows <= 0) & (ni_columns > 0))
		ni_rows = nrecords/ni_columns;
	if ((ni_columns <= 0) & (ni_rows > 0))
		ni_columns = nrecords/ni_rows;
	if (nrecords != ni_rows*ni_columns)
		goto badspec;
	if (no_columns <= 0)
		no_columns = ni_rows;
	if (no_rows <= 0)
		no_rows = ni_columns;
	if ((no_rows != ni_columns) | (no_columns != ni_rows))
		goto badspec;
						/* transpose records */
	for (i = 0; i < no_rows; i++) {
	    for (j = 0; j < no_columns; j++)
		if (rp != NULL) {		/* ASCII output */
			print_record(rp, j*ni_columns + i);
			putc(tabEOL[j >= no_columns-1], stdout);
		} else {			/* binary output */
			fwrite((char *)mp->base +
					-record_width*(j*ni_columns + i),
					-record_width, 1, stdout);
		}
	    if (ferror(stdout)) {
		fprintf(stderr, "Error writing to stdout\n");
		return(0);
	    }
	}
	if (rp != NULL)
		free_records(rp);
	return(1);
badspec:
	fprintf(stderr, "Bad transpose specification -- check dimension(s)\n");
	return(0);
}

/* resize ASCII stream input by ignoring EOLs between records */
static int
do_resize(FILE *fp)
{
	long	records2go = ni_rows*ni_columns;
	int	columns2go = no_columns;
	char	word[256];
						/* sanity checks */
	if (record_width <= 0) {
		fprintf(stderr, "Bad call to do_resize (record_width = %d)\n",
				record_width);
		return(0);
	}
	if (no_columns <= 0) {
		fprintf(stderr, "Missing -oc specification\n");
		return(0);
	}
	if ((records2go <= 0) & (no_rows > 0))
		records2go = no_rows*no_columns;
	else if (no_rows*no_columns != records2go) {
		fprintf(stderr,
			"Input and output data sizes disagree (%dx%d != %dx%d)\n",
				ni_rows, ni_columns, no_rows, no_columns);
		return(0);
	}
	do {					/* reshape records */
		int	n;

		for (n = record_width; n--; ) {
			if (fget_word(word, fp) == NULL) {
				if (records2go > 0 || n < record_width-1)
					break;
				goto done;	/* normal EOD */
			}
			fputs(word, stdout);
			if (n) {		/* mid-record? */
				int	c = getc(fp);
				if ((c == '\n') | (c == EOF))
					break;
				ungetc(c, fp);
				putc(' ', stdout);
			}
		}
		if (n >= 0) {
			fprintf(stderr, "Incomplete record / unexpected EOF\n");
			return(0);
		}
		if (--columns2go <= 0) {	/* time to end output row? */
			putc('\n', stdout);
			columns2go = no_columns;
		} else				/* else separate records */
			putc('\t', stdout);
	} while (--records2go);			/* expected EOD? */
done:
	if (warnings && columns2go != no_columns)
		fprintf(stderr, "Warning -- incomplete final row\n");
	if (warnings && fget_word(word, fp) != NULL)
		fprintf(stderr, "Warning -- characters beyond expected EOD\n");
	return(1);
}

/* process a header line and copy to stdout */
static int
headline(char *s, void *p)
{
	char	fmt[32];

	if (formatval(fmt, s)) {
		if (!strcmp(fmt, fmtid))
			return(0);
		fprintf(stderr, "Input format '%s' != '%s'\n", fmt, fmtid);
		return(-1);
	}
	fputs(s, stdout);			/* copy header info. */
	return(0);
}

/* main routine for converting rows/columns in data file */
int
main(int argc, char *argv[])
{
	int	do_header = 1;			/* header i/o? */
	int	transpose = 0;			/* transpose rows & cols? */
	int	i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'i':			/* input */
			if (argv[i][2] == 'c')	/* columns */
				ni_columns = atoi(argv[++i]);
			else if (argv[i][2] == 'r')
				ni_rows = atoi(argv[++i]);
			else
				goto userr;
			break;
		case 'o':			/* output */
			if (argv[i][2] == 'c')	/* columns */
				no_columns = atoi(argv[++i]);
			else if (argv[i][2] == 'r')
				no_rows = atoi(argv[++i]);
			else
				goto userr;
			break;
		case 'h':			/* header on/off */
			do_header = !do_header;
			break;
		case 't':			/* transpose on/off */
			transpose = !transpose;
			break;
		case 'f':			/* format */
			switch (argv[i][2]) {
			case 'a':		/* ASCII */
			case 'A':
				fmtid = "ascii";
				record_width = 1;
				break;
			case 'f':		/* float */
			case 'F':
				fmtid = "float";
				record_width = -(int)sizeof(float);
				break;
			case 'd':		/* double */
			case 'D':
				fmtid = "double";
				record_width = -(int)sizeof(double);
				break;
			case 'b':		/* binary (bytes) */
			case 'B':
				fmtid = "byte";
				record_width = -1;
				break;
			default:
				goto userr;
			}
			if (argv[i][3]) {
				if (!isdigit(argv[i][3]))
					goto userr;
				record_width *= atoi(argv[i]+3);
			}
			break;
		case 'w':			/* warnings on/off */
			warnings = !warnings;
			break;
		default:
			goto userr;
		}
	if (!record_width)
		goto userr;
	if (i < argc-1)				/* arg count OK? */
		goto userr;
						/* open input file? */
	if (i == argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: cannot open for reading\n", argv[i]);
		return(1);
	}
	if (record_width < 0) {
		SET_FILE_BINARY(stdin);
		SET_FILE_BINARY(stdout);
	}
						/* check for no-op */
	if (!transpose && (record_width < 0 ||
			(no_columns == ni_columns) & (no_rows == ni_rows))) {
		if (warnings)
			fprintf(stderr, "%s: no-op -- copying input verbatim\n",
				argv[0]);
		if (!output_stream(stdin))
			return(1);
		return(0);
	}
	if (do_header) {			/* read/write header */
		if (getheader(stdin, &headline, NULL) < 0)
			return(1);
		printargs(argc, argv, stdout);
		fputformat(fmtid, stdout);
		fputc('\n', stdout);		/* finish new header */
	}
	if (transpose) {			/* transposing rows & columns? */
		MEMLOAD	myMem;			/* need to load into memory */
		if (i == argc-1) {
			if (load_file(&myMem, stdin) <= 0) {
				fprintf(stderr, "%s: error loading file into memory\n",
						argv[i]);
				return(1);
			}
		} else if (load_stream(&myMem, stdin) <= 0) {
			fprintf(stderr, "%s: error loading stdin into memory\n",
						argv[0]);
			return(1);
		}
		if (!do_transpose(&myMem))
			return(1);
		/* free_load(&myMem); */
	} else if (!do_resize(stdin))		/* just reshaping input */
		return(1);
	return(0);
userr:
	fprintf(stderr,
"Usage: %s [-h][-w][-f[afdb][N]][-t][-ic in_col][-ir in_row][-oc out_col][-or out_row] [input.dat]\n",
			argv[0]);
	return(1);
}
