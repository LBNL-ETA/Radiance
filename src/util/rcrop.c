#ifndef lint
static const char RCSid[] = "$Id: rcrop.c,v 1.12 2022/03/21 20:19:19 greg Exp $";
#endif
/*
 * rcrop.c - crop a Radiance picture or matrix data
 */

#include <ctype.h>
#include "rtio.h"
#include "platform.h"
#include "color.h"
#include "fvect.h"
#include "view.h"

char	*progname;		/* global argv[0] */

VIEW	vw = STDVIEW;
int	gotvw = 0;
char	fmt[MAXFMTLEN] = "ascii";	/* assumed when unspecified */
int	ncomp = 0;
RESOLU	res;
int	rmin, cmin, nrows, ncols;

/* Process header line, copying to stdout when appropriate */
static int
headline(char *s, void *p)
{
	if (formatval(fmt, s))
		return(0);
	if (!strncmp(s, "NCOMP=", 6)) {
		ncomp = atoi(s+6);
		return(-(ncomp <= 0));
	}
	if (!strncmp(s, "NROWS=", 6)) {
		res.rt = PIXSTANDARD;
		res.yr = atoi(s+6);
		return(-(res.yr <= 0));
	}
	if (!strncmp(s, "NCOLS=", 6)) {
		res.rt = PIXSTANDARD;
		res.xr = atoi(s+6);
		return(-(res.xr <= 0));
	}
	if (isview(s)) {
		gotvw += sscanview(&vw, s);
		return(0);
	}
	fputs(s, stdout);		/* copy other header info. */
	return(0);
}

/* Copy routine for COLR data */
static int
colr_copyf(FILE *fp)
{
	const int	width = scanlen(&res);
	COLR		*scan = (COLR *)malloc(sizeof(COLR)*width);
	int		y;

	if (!scan) {
		fputs(progname, stderr);
		fputs(": out of memory!\n", stderr);
		return(0);
	}
	for (y = 0; y < rmin; y++)	/* initial skip */
		if (freadcolrs(scan, width, fp) < 0)
			goto readerr;
					/* scanlines to copy */
	for (y = 0; y < nrows; y++) {
		if (freadcolrs(scan, width, fp) < 0)
			goto readerr;
		if (fwritecolrs(scan+cmin, ncols, stdout) < 0)
			goto writerr;
	}
	free(scan);
	if (fflush(stdout) == 0)
		return(1);
writerr:
	fputs(progname, stderr);
	fputs(": error writing picture\n", stderr);
	return(0);
readerr:
	fputs(progname, stderr);
	fputs(": error reading picture\n", stderr);
	return(0);
}

/* Copy routine for binary data (asize = sizeof(type)) */
static int
binary_copyf(FILE *fp, int asize)
{
	const int	skip_thresh = 8192;
	const size_t	elsiz = asize*ncomp;
	const int	width = scanlen(&res);
	const long	skip_len = (width-ncols)*elsiz;
	char		*buf;
	int		y;
					/* check if fseek() useful */
	if (skip_len > skip_thresh &&
			fseek(fp, (rmin*width + cmin)*elsiz, SEEK_CUR) == 0) {
		buf = (char *)malloc(ncols*elsiz);
		if (!buf)
			goto memerr;
		for (y = nrows; y-- > 0; ) {
			if (getbinary(buf, elsiz, ncols, fp) != ncols)
				goto readerr;
			if (putbinary(buf, elsiz, ncols, stdout) != ncols)
				goto writerr;
			if (y && fseek(fp, skip_len, SEEK_CUR) < 0) {
				fputs(progname, stderr);
				fputs(": unexpected seek error on input\n", stderr);
				return(0);
			}
		}
		free(buf);
		if (fflush(stdout) == EOF)
			goto writerr;
		return(1);		/* success! */
	}				/* else need to read it all... */
	buf = (char *)malloc(width*elsiz);
	if (!buf)
		goto memerr;
					/* skip rows as requested */
	if (skip_len > skip_thresh ||
			(rmin && fseek(fp, rmin*width*elsiz, SEEK_CUR) < 0))
		for (y = 0; y < rmin; y++)
			if (getbinary(buf, elsiz, width, fp) != width)
				goto readerr;
	for (y = 0; y < nrows; y++) {	/* copy portion */
		if (getbinary(buf, elsiz, width, fp) != width)
			goto readerr;
		if (putbinary(buf+cmin*elsiz, elsiz, ncols, stdout) != ncols)
			goto writerr;
	}
	free(buf);			/* we're done */
	if (fflush(stdout) == 0)
		return(1);
writerr:
	fputs(progname, stderr);
	fputs(": error writing binary data\n", stderr);
	return(0);
readerr:
	fputs(progname, stderr);
	fputs(": error reading binary data\n", stderr);
	return(0);
memerr:
	fputs(progname, stderr);
	fputs(": out of memory!\n", stderr);
	return(0);
}

/* Read (and copy) specified number of white-space-separated words */
static int
readwords(FILE *finp, int nwords, FILE *fout)
{
	while (nwords-- > 0) {
		int	c;
		do {
			c = getc(finp);
		} while (isspace(c));
		if (c == EOF)
			return(-1);
		if (fout && fputc(' ', fout) == EOF)
			return(-1);
		do {
			if (fout)
				putc(c, fout);
		} while ((c = getc(finp)) != EOF && !isspace(c));
	}
	return(0);
}

/* Copy routine for ascii data */
static int
ascii_copyf(FILE *fp)
{
	const int	width = scanlen(&res);
	int		x, y;

	SET_FILE_TEXT(fp);		/* started as binary */
	SET_FILE_TEXT(stdout);
					/* skip rows as requested */
	if (readwords(fp, rmin*width*ncomp, NULL) < 0)
		goto io_err;
	for (y = 0; y < nrows; y++) {	/* copy part */
		if (readwords(fp, cmin*ncomp, NULL) < 0)
			goto io_err;
		if (readwords(fp, ncols*ncomp, stdout) < 0)
			goto io_err;
		fputc('\n', stdout);	/* newline per row */
		if (readwords(fp, (width-ncols-cmin)*ncomp, NULL) < 0)
			goto io_err;
	}
	if (fflush(stdout) == 0)
		return(1);
io_err:
	fputs(progname, stderr);
	fputs(": error copying ascii data\n", stderr);
	return(0);
}

/* Adjust (crop) our view */
static int
adjust_view(void)
{
	double		p0[2], p1[2];
	const char	*err;

	if (res.rt & YMAJOR) {
		p0[0] = cmin/(double)res.xr;
		p0[1] = rmin/(double)res.yr;
		p1[0] = (cmin+ncols)/(double)res.xr;
		p1[1] = (rmin+nrows)/(double)res.yr;
	} else {
		p0[0] = rmin/(double)res.xr;
		p0[1] = cmin/(double)res.yr;
		p1[0] = (rmin+nrows)/(double)res.xr;
		p1[1] = (cmin+ncols)/(double)res.yr;
	}
	if (res.rt & XDECR) {
		p0[0] = 1. - p0[0];
		p1[0] = 1. - p1[0];
	}
	if (res.rt & YDECR) {
		p0[1] = 1. - p0[1];
		p1[1] = 1. - p1[1];
	}
	err = cropview(&vw, p0[0], p0[1], p1[0], p1[1]);

	if (!err)
		return(1);	/* success! */

	fputs(progname, stderr);
	fputs(": view error - ", stderr);
	fputs(err, stderr);
	fputc('\n', stderr);
	return(0);		/* something went wrong */
}


/* Main routine -- load header and call processor */
int
main(int argc, char *argv[])
{
	FILE	*fp = stdin;
	int	asiz = 0;
	int	gotdims;

	progname = argv[0];
				/* get input and output */
	if ((argc < 5) | (argc > 7))
		goto usage;
	if (!isint(argv[1]) | !isint(argv[2]) |
			!isint(argv[3]) | !isint(argv[4]))
		goto usage;
	rmin = atoi(argv[1]);
	cmin = atoi(argv[2]);
	nrows = atoi(argv[3]);
	ncols = atoi(argv[4]);
	if ((rmin < 0) | (cmin < 0) | (nrows < 0) | (ncols < 0))
		goto usage;
	if (argc <= 5)
		SET_FILE_BINARY(fp);
	else if (!(fp = fopen(argv[5], "rb"))) {
		fputs(argv[5], stderr);
		fputs(": cannot open for reading\n", stderr);
		return(1);
	}
	if (argc <= 6)
		SET_FILE_BINARY(stdout);
	else if (!freopen(argv[6], "wb", stdout)) {
		fputs(argv[6], stderr);
		fputs(": cannot open for writing\n", stderr);
		return(1);
	}
#ifdef getc_unlocked		/* avoid stupid semaphores */
	flockfile(fp);
	flockfile(stdout);
#endif
				/* process information header */
	if (getheader(fp, headline, NULL) < 0) {
		fputs(progname, stderr);
		fputs(": bad input header\n", stderr);
		return(1);
	}
	gotdims = (res.rt == PIXSTANDARD) & (res.xr > 0) & (res.yr > 0);
	if (!gotdims && !fgetsresolu(&res, fp)) {
		fputs(progname, stderr);
		fputs(": missing input dimensions\n", stderr);
		return(1);
	}
	if (!nrows)
		nrows = numscans(&res) - rmin;
	if (!ncols)
		ncols = scanlen(&res) - cmin;
	if ((nrows <= 0) | (ncols <= 0) |
			(rmin+nrows > numscans(&res)) |
			(cmin+ncols > scanlen(&res))) {
		fputs(progname, stderr);
		fputs(": illegal crop\n", stderr);
		return(1);
	}
	printargs(5, argv, stdout);	/* add to header */
	if (gotvw && adjust_view()) {
		fputs(VIEWSTR, stdout);	/* write adjusted view */
		fprintview(&vw, stdout);
		fputc('\n', stdout);
	}
	if (gotdims)			/* dimensions + format */
		printf("NROWS=%d\nNCOLS=%d\n", nrows, ncols);
	if (ncomp)
		printf("NCOMP=%d\n", ncomp);
	fputformat(fmt, stdout);	/* will align bytes if it can */
	fputc('\n', stdout);		/* end of new header */
	if (!gotdims) {			/* add resolution string? */
		RESOLU	newres;
		if (res.rt & YMAJOR) {
			newres.xr = ncols;
			newres.yr = nrows;
		} else {
			newres.xr = nrows;
			newres.yr = ncols;
		}
		newres.rt = res.rt;
		fputsresolu(&newres, stdout);
	}
					/* call appropriate processor */
	if (!strcmp(fmt, "float")) {
		asiz = sizeof(float);
	} else if (!strcmp(fmt, "double")) {
		asiz = sizeof(double);
	} else if (!strcmp(fmt, "32-bit_encoded_normal")) {
		asiz = 4;
		ncomp = 1;
	} else if (!strcmp(fmt, "16-bit_encoded_depth")) {
		asiz = 2;
		ncomp = 1;
	} else if (globmatch(PICFMT, fmt)) {
		asiz = -1;
		if (!ncomp) ncomp = 3;
		else ncomp *= (ncomp == 3);
	} else if (strcasecmp(fmt, "ascii")) {
		fputs(progname, stderr);
		fputs(": unsupported format - ", stderr);
		fputs(fmt, stderr);
		fputc('\n', stderr);
		return(1);
	}
	if (ncomp <= 0) {
		fputs(progname, stderr);
		fputs(": illegal number of components\n", stderr);
		return(1);
	}
	if (!(asiz < 0 ? colr_copyf(fp) :
			asiz ? binary_copyf(fp, asiz) : ascii_copyf(fp)))
		return(1);
					/* need to consume the rest? */
	if (fp == stdin && rmin+nrows < numscans(&res) &&
			fseek(fp, 0L, SEEK_END) < 0)
		while (getc(fp) != EOF)
			;
	return(0);
usage:
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" row0 col0 nrows ncols [input [output]]\n", stderr);
	return(1);
}
