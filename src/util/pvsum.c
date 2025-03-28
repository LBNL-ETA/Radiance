#ifndef lint
static const char RCSid[] = "$Id: pvsum.c,v 2.3 2025/03/28 21:36:31 greg Exp $";
#endif
/*
 *	pvsum.c - add together spectral and/or float pictures
 *			based on vector or matrix, similar to dctimestep
 */

#include <math.h>
#include "rtio.h"
#include "resolu.h"
#include "platform.h"
#include "paths.h"
#include "random.h"
#include "rmatrix.h"
#if !defined(_WIN32) && !defined(_WIN64)
#include <sys/mman.h>
#endif

int	nprocs = 1;			/* # of calculation processes (Unix) */
int	in_type = DTfromHeader;		/* input data type */
int	out_type = DTfromHeader;	/* output data type */
char	*in_spec = NULL;		/* input specification */
char	*out_spec = NULL;		/* output file specification */

int	iswapped = 0;			/* input data is byte-swapped? */
int	ncomp = 3;			/* # input components */
int	xres=0, yres=0;			/* input image dimensions */

RMATRIX	*cmtx = NULL;			/* coefficient matrix */

/* does the given spec contain integer format? */
int
hasFormat(const char *s)
{
restart:
	if (s) s = strchr(s, '%');
	if (!s) return(0);
	if (s[1] == '%') {		/* "%%" ? */
		s += 2;
		goto restart;
	}
	while (*++s) {
		if (strchr("diouxX", *s))
			return(1);	/* integer format */
		if (strchr("%fFeEgGaAcsb", *s))
			break;		/* non-integer format */
	}
	return(0);
}

/* get first input header data we'll need */
int
iheadline(char *s, void *p)
{
	char	fmt[MAXFMTLEN];
	int	i;

	if (!strncmp(s, "NCOLS=", 6)) {
		xres = atoi(s+6);
		return(1);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		yres = atoi(s+6);
		return(1);
	}
	if (isncomp(s)) {
		ncomp = ncompval(s);
		return(1);
	}
	if (formatval(fmt, s)) {
		for (in_type = DTend; --in_type > DTfromHeader; )
			if (!strcmp(fmt, cm_fmt_id[in_type]))
				return(1);
		return(-1);
	}
	i = isbigendian(s);
	if (i >= 0) {
		iswapped = (i != nativebigendian());
		return(1);
	}
	return(0);
}

/* open initial file and get relevant dimensions and type */
int
get_iotypes(void)
{
	char	fbuf[256];
	FILE	*fp;

	sprintf(fbuf, in_spec, 0);
	fp = fopen(fbuf, "rb");
	if (!fp) {
		fprintf(stderr, "%s: cannot open for reading\n", fbuf);
		return(0);
	}
	if (getheader(fp, iheadline, NULL) < 0) {
		fprintf(stderr, "%s: bad header - wrong format?\n", fbuf);
		fclose(fp);
		return(0);
	}
	if ((xres <= 0) | (yres <= 0) && !fscnresolu(&xres, &yres, fp)) {
		fprintf(stderr, "%s: missing input resolution\n", fbuf);
		fclose(fp);
		return(0);
	}
	if (nprocs > 1 && (in_type==DTrgbe) | (in_type==DTxyze)) {
		long	data_start = ftell(fp);		/* make sure input flat */
		off_t	dend = lseek(fileno(fp), 0, SEEK_END);
		if (dend < data_start + 4L*xres*yres) {
			fputs("Warning - multi-processing requires flat input files\n",
					stderr);
			nprocs = 1;
		}
	}
	fclose(fp);
	if ((cmtx->ncomp == 1) & (ncomp != 1)) {
		double	xfm[MAXCSAMP];
		RMATRIX	*nmtx;
		int	i;
		for (i = ncomp; i--; )
			xfm[i] = 1.;
		nmtx = rmx_transform(cmtx, ncomp, xfm);
		if (!nmtx)
			return(0);
		rmx_free(cmtx);
		cmtx = nmtx;
	} else if (cmtx->ncomp != ncomp) {
		fprintf(stderr, "operation %s needs %d components, has %d\n",
				cmtx->ncols == 1 ? "vector" : "matrix",
				ncomp, cmtx->ncomp);
		return(0);
	}
	if ((in_type != DTrgbe) & (in_type != DTxyze) & (in_type != DTspec) &
			(in_type != DTfloat)) {
		fprintf(stderr, "%s unsupported input data type: %s\n",
				fbuf, cm_fmt_id[in_type]);
		return(0);
	}
	if ((out_type == DTrgbe) & (ncomp > 3))
		out_type = DTspec;
	else if (out_type == DTfromHeader ||
			(out_type == DTrgbe) & (in_type != DTfloat))
		out_type = in_type;
	return(1);
}

/* check subsequent headers match initial file */
int
checkline(char *s, void *p)
{
	static int	exposWarned = 0;
	int		*xyres = (int *)p;
	char		fmt[MAXFMTLEN];

	if (!strncmp(s, "NCOLS=", 6)) {
		xyres[0] = atoi(s+6);
		if (xyres[0] <= 0)
			return(-1);
		return(1);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		xyres[1] = atoi(s+6);
		if (xyres[1] <= 0)
			return(-1);
		return(1);
	}
	if (isncomp(s)) {
		if (ncompval(s) != ncomp)
			return(-1);
		return(1);
	}
	if (isexpos(s)) {
		if (!exposWarned && fabs(1. - exposval(s)) > 0.04) {
			fputs("Warning - ignoring EXPOSURE setting(s)\n",
					stderr);
			exposWarned++;
		}
		return(1);
	}
	if (formatval(fmt, s)) {
		if (strcmp(fmt, cm_fmt_id[in_type]))
			return(-1);
		return(1);
	}
	return(0);
}

/* open and check input file */
FILE *
open_input(char *fname)
{
	int	xyres[2];
	FILE	*fp = fopen(fname, "rb");

	if (!fp) {
		fprintf(stderr, "%s: cannot open for reading\n", fname);
		return(NULL);
	}
	xyres[0] = xyres[1] = 0;
	if (getheader(fp, checkline, xyres) < 0) {
		fprintf(stderr, "%s: bad/inconsistent header\n", fname);
		fclose(fp);
		return(NULL);
	}
	if ((xyres[0] <= 0) | (xyres[1] <= 0) &&
			!fscnresolu(&xyres[0], &xyres[1], fp)) {
		fprintf(stderr, "%s: missing resolution\n", fname);
		fclose(fp);
		return(NULL);
	}
	if ((xyres[0] != xres) | (xyres[1] != yres)) {
		fprintf(stderr, "%s: mismatched resolution\n", fname);
		fclose(fp);
		return(NULL);
	}
	return(fp);
}

/* open output file or command (if !NULL) and write info header */
FILE *
open_output(char *ospec, int fno)
{
	FILE	*fp;

	if (!ospec) {
		ospec = "<stdout>";
		fp = stdout;
	} else if (ospec[0] == '!') {
		if (!(fp = popen(ospec+1, "w"))) {
			fprintf(stderr, "Cannot start: %s\n", ospec);
			return(NULL);
		}
	} else if (!(fp = fopen(ospec, "w"))) {
		fprintf(stderr, "%s: cannot open for writing\n", ospec);
		return(NULL);
	}
	SET_FILE_BINARY(fp);
	newheader("RADIANCE", fp);
	if (cmtx->info)			/* prepend matrix metadata */
		fputs(cmtx->info, fp);
	else
		fputnow(fp);
	if (fno >= 0)
		fprintf(fp, "FRAME=%d\n", fno);
	switch (out_type) {
	case DTfloat:
	case DTdouble:
		fprintf(fp, "NCOLS=%d\nNROWS=%d\n", xres, yres);
		fputncomp(ncomp, fp);
		fputendian(fp);
		fputformat(cm_fmt_id[out_type], fp);
		fputc('\n', fp);
		break;
	case DTrgbe:
		fputformat(COLRFMT, fp);
		fputc('\n', fp);
		fprtresolu(xres, yres, fp);
		break;
	case DTxyze:
		fputformat(CIEFMT, fp);
		fputc('\n', fp);
		fprtresolu(xres, yres, fp);
		break;
	case DTspec:
		fputncomp(ncomp, fp);
		fputwlsplit(cmtx->wlpart, fp);
		fputformat(SPECFMT, fp);
		fputc('\n', fp);
		fprtresolu(xres, yres, fp);
		break;
	default:
		fputs("Unsupported output type!\n", stderr);
		return(NULL);
	}
	if (fflush(fp) < 0) {
		fprintf(stderr, "%s: write error\n", ospec);
		fclose(fp);
		return(NULL);
	}
	return(fp);
}

/* run calculation from a single process */
int
solo_process(void)
{
	float	*osum = (float *)calloc((size_t)xres*yres, sizeof(float)*ncomp);
	COLORV	*iscan = (COLORV *)malloc(sizeof(COLORV)*ncomp*xres);
	char	fbuf[512];
	int	c;

	if (!osum | !iscan) {
		fprintf(stderr, "Cannot allocate %dx%d %d-component accumulator\n",
				xres, yres, ncomp);
		return(0);
	}
	if (sizeof(float) != sizeof(COLORV)) {
		fputs("Code Error 1 in solo_process()\n", stderr);
		return(0);
	}
	for (c = 0; c < cmtx->ncols; c++) {	/* run through each column/output */
		FILE	*fout;
		int	y;
		int	rc = cmtx->nrows;
		if (c > 0)		/* clear accumulator? */
			memset(osum, 0, sizeof(float)*ncomp*xres*yres);
		while (rc-- > 0) {	/* run through each input file */
			const int	r = c&1 ? rc : cmtx->nrows-1 - rc;
			const rmx_dtype	*cval = rmx_val(cmtx, r, c);
			FILE		*finp;
			int		i, x;
			for (i = ncomp; i--; )
				if (cval[i] != 0) break;
			if (i < 0)	/* this coefficient is zero, skip */
				continue;
			sprintf(fbuf, in_spec, r);
			finp = open_input(fbuf);
			if (!finp)
				return(0);
			for (y = 0; y < yres; y++) {
				float	*dst = osum + y*xres*ncomp;
				if (in_type == DTfloat	? getbinary(iscan, sizeof(float)*ncomp,
									xres, finp) != xres
							: freadsscan(iscan, ncomp, xres, finp) < 0)
					goto readerr;
				if ((in_type == DTfloat) & iswapped)
					swap32((char *)iscan, ncomp*xres);
				for (x = 0; x < xres; x++, dst += ncomp)
					for (i = ncomp; i--; )
						dst[i] += cval[i]*iscan[x*ncomp + i];
			}
			fclose(finp);
		}			/* write out accumulated column result... */
		if (out_spec) {		/* ...to file or command */
			if (cmtx->ncols > 1 && !hasFormat(out_spec)) {
				fputs("Sequential result must go to stdout\n", stderr);
				return(0);
			}
			sprintf(fbuf, out_spec, c);
			fout = open_output(fbuf, c-(cmtx->ncols==1));
		} else {		/* ...to stdout */
			if ((out_type == DTfloat) & (cmtx->ncols > 1)) {
				fputs("Float outputs must have separate destinations\n",
						stderr);
				return(0);
			}
			strcpy(fbuf, "<stdout>");
			fout = open_output(NULL, c-(cmtx->ncols==1));
		}
		if (!fout)
			return(0);	/* assume error was reported */
		if (out_type != DTfloat) {
			for (y = 0; y < yres; y++)
				if (fwritesscan(osum + (size_t)y*xres*ncomp,
							ncomp, xres, fout) < 0)
					goto writerr;
		} else if (fwrite(osum, sizeof(float)*ncomp, (size_t)xres*yres, fout) !=
					(size_t)xres*yres)
			goto writerr;

		if (fbuf[0] == '!') {
			if (pclose(fout) != 0) {
				fprintf(stderr, "Bad status from: %s\n", fbuf);
				return(0);
			}
		} else if (fout != stdout && fclose(fout) == EOF)
			goto writerr;
	}
	free(osum);			/* clean up on success */
	free(iscan);
	return(1);
readerr:
	fprintf(stderr, "%s: read error\n", fbuf);
	return(0);
writerr:
	fprintf(stderr, "%s: write error\n", fbuf);
	return(0);
}

/* allocate a scrambled index array of the specified length */
int *
scramble(int n)
{
	int	*scarr = (int *)malloc(sizeof(int)*n);
	int	i;

	if (!scarr) {
		fprintf(stderr, "Out of memory in scramble(%d)\n", n);
		exit(1);
	}
	for (i = n; i--; )
		scarr[i] = i;
					/* perform Fisher-Yates shuffle */
	for (i = 0; i < n-1; i++) {
		int	ix = irandom(n-i) + i;
		int	ndx = scarr[i];
		scarr[i] = scarr[ix];
		scarr[ix] = ndx;
	}
	return(scarr);
}

/* run calculation on multiple processes, using memory maps and fork() */
int
multi_process(void)
{
#if defined(_WIN32) || defined(_WIN64)
	fputs("Bad call to multi_process()\n", stderr);
	return(0);
#else
	int	coff = nprocs;
	int	odd = 0;
	char	fbuf[512];
	float	*osum;
	int	*syarr;
	int	c;
					/* sanity check */
	if (sizeof(float) != sizeof(COLORV)) {
		fputs("Code Error 1 in multi_process()\n", stderr);
		return(0);
	}
	while (--coff > 0) {		/* parent births children */
		int	pid = fork();
		if (pid < 0) {
			fputs("fork() call failed!\n", stderr);
			return(0);
		}
		if (pid == 0) break;	/* child gets to work */
	}
	osum = (float *)calloc((size_t)xres*yres, sizeof(float)*ncomp);
	if (!osum) {
		fprintf(stderr, "Cannot allocate %dx%d %d-component accumulator\n",
				xres, yres, ncomp);
		return(0);
	}
	srandom(113*coff + 5669);	/* randomize row access for this process */
	syarr = scramble(yres);
					/* run through our unique set of columns */
	for (c = coff; c < cmtx->ncols; c += nprocs) {
		FILE	*fout;
		int	y;
		int	rc = cmtx->nrows;
		if (c > coff)		/* clear accumulator? */
			memset(osum, 0, sizeof(float)*ncomp*xres*yres);
		while (rc-- > 0) {	/* map & sum each input file */
			const int	r = odd ? rc : cmtx->nrows-1 - rc;
			const rmx_dtype	*cval = rmx_val(cmtx, r, c);
			long		dstart;
			size_t		maplen;
			void		*imap;
			FILE		*finp;
			float		*dst;
			int		i, x;
			for (i = ncomp; i--; )
				if (cval[i] != 0) break;
			if (i < 0)	/* this coefficient is zero, skip */
				continue;
			sprintf(fbuf, in_spec, r);
			finp = open_input(fbuf);
			if (!finp)
				return(0);
			dstart = ftell(finp);
			if (dstart < 0) {
				fprintf(stderr, "%s: ftell() failed!\n", fbuf);
				return(0);
			}
			if (in_type == DTfloat && dstart%sizeof(float)) {
				fprintf(stderr, "%s: float header misalignment\n", fbuf);
				return(0);
			}
			i = in_type==DTfloat ? ncomp*(int)sizeof(float) : ncomp+1;
			maplen = dstart + yres*xres*i;
			imap = mmap(NULL, maplen, PROT_READ,
					MAP_FILE|MAP_SHARED, fileno(finp), 0);
			fclose(finp);		/* will read from map (randomly) */
			if (imap == MAP_FAILED) {
				fprintf(stderr, "%s: unable to map input file\n", fbuf);
				return(0);
			}
			if (in_type == DTfloat)
			    for (y = yres; y-- > 0; ) {
			    	const float	*fvp = (float *)((char *)imap + dstart) +
			    				(size_t)ncomp*xres*syarr[y];
				dst = osum + (size_t)ncomp*xres*syarr[y];
				for (x = xres; x-- > 0; dst += ncomp, fvp += ncomp)
				    for (i = ncomp; i--; )
					dst[i] += cval[i]*fvp[i];
			    }
			else
			    for (y = yres; y-- > 0; ) {
				const COLRV	*cvp = (COLRV *)((char *)imap + dstart) +
							(ncomp+1L)*xres*syarr[y];
				dst = osum + (size_t)ncomp*xres*syarr[y];
				for (x = xres; x-- > 0; dst += ncomp, cvp += ncomp+1) {
			    	    const rmx_dtype	fe = cxponent[cvp[ncomp]];
				    for (i = ncomp; i--; )
					dst[i] += cval[i]*(cvp[i]+(rmx_dtype).5)*fe;
				}
			    }
			munmap(imap, maplen);
		}			/* write accumulated column picture/matrix */
		sprintf(fbuf, out_spec, c);
		fout = open_output(fbuf, c);
		if (!fout)
			return(0);	/* assume error was reported */
		if (out_type != DTfloat) {
			for (y = 0; y < yres; y++)
				if (fwritesscan(osum + (size_t)y*xres*ncomp,
							ncomp, xres, fout) < 0)
					goto writerr;
		} else if (fwrite(osum, sizeof(float)*ncomp, (size_t)xres*yres, fout) !=
					(size_t)xres*yres)
			goto writerr;

		if (fbuf[0] == '!') {
			if (pclose(fout) != 0) {
				fprintf(stderr, "Bad status from: %s\n", fbuf);
				return(0);
			}
		} else if (fclose(fout) == EOF)
			goto writerr;
		odd = !odd;		/* go back & forth to milk page cache */
	}
	free(osum);
	free(syarr);
	if (coff)			/* child processes return here... */
		return(1);
	c = 0;				/* ...but parent waits for children */
	while (++coff < nprocs) {
		int	st;
		if (wait(&st) < 0)
			break;
		if (st) c = st;
	}
	return(c == 0);
writerr:
	fprintf(stderr, "%s: write error\n", fbuf);
	return(0);
#endif
}

int
main(int argc, char *argv[])
{
	int	a;

	for (a = 1; a < argc-1 && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'o':		/* output spec/format */
			switch (argv[a][2]) {
			case '\0':
				out_spec = argv[++a];
				break;
			case 'f':
				out_type = DTfloat;
				break;
			case 'c':
				out_type = DTrgbe;
				break;
			default:
				goto badopt;
			}
			break;
		case 'N':		/* number of desired processes */
		case 'n':		/* quietly supported alternate */
			nprocs = atoi(argv[++a]);
			if (nprocs <= 0)
				goto userr;
			break;
		default:;
badopt:			fprintf(stderr, "%s: bad option: %s\n", argv[0], argv[a]);
			goto userr;
			return(1);
		}
	if ((argc-a < 1) | (argc-a > 2) || argv[a][0] == '-')
		goto userr;
	in_spec = argv[a];
	cmtx = rmx_load(argv[a+1], RMPnone);	/* loads from stdin if a+1==argc */
	if (cmtx == NULL)
		return(1);		/* error reported */
	if (nprocs > cmtx->ncols)
		nprocs = cmtx->ncols;
#if defined(_WIN32) || defined(_WIN64)
	if (nprocs > 1) {
		fprintf(stderr, "%s: warning - Windows only allows -N 1\n", argv[0]);
		nprocs = 1;
	}
#else
	if ((nprocs > 1) & !out_spec) {
		fprintf(stderr, "%s: multi-processing result cannot go to stdout\n",
				argv[0]);
		nprocs = 1;
	}
	if ((nprocs > 1)  & iswapped && (in_type==DTfloat) | (in_type==DTdouble)) {
		fprintf(stderr, "%s: multi-processing unsupported on swapped input\n",
				argv[0]);
		nprocs = 1;
	}
#endif
	if (cmtx->nrows > 1 && !hasFormat(in_spec)) {
		fprintf(stderr, "%s: input specification '%s' needs %%d format\n",
				argv[0], in_spec);
		goto userr;
	}
	if (!get_iotypes())
		return(1);
	if (!(nprocs == 1 ? solo_process() : multi_process()))
		return(1);
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-oc | -of][-o ospec][-N nproc] inpspec [mtx]\n",
				argv[0]);
	return(1);
}
