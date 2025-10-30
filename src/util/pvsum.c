#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 *	pvsum.c - add together spectral and/or float pictures
 *			based on vector or matrix, similar to dctimestep
 */

#include <math.h>
#include "rtio.h"
#include "resolu.h"
#include "platform.h"
#include "random.h"
#include "rmatrix.h"
#if !defined(_WIN32) && !defined(_WIN64)
#include <sys/mman.h>
#include <sys/wait.h>
#endif

#define  VIEWSTR	"VIEW="		/* borrowed from common/view.h */
#define  VIEWSTRL	5

int	nprocs = 1;			/* # of calculation processes (Unix) */
int	in_type = DTfromHeader;		/* input data type */
int	out_type = DTfromHeader;	/* output data type */
char	*in_spec = NULL;		/* input specification */
char	*out_spec = NULL;		/* output file specification */

int	iswapped = 0;			/* input data is byte-swapped? */
int	ncomp = 3;			/* # input components */
int	xres=0, yres=0;			/* input image dimensions */
char	viewspec[128] = "";		/* VIEW= line from first header */
char	pixasp[48] = "";		/* PIXASPECT= line from header */

int	gargc;				/* global argc */
char	**gargv;			/* global argv */

RMATRIX	*cmtx;				/* coefficient matrix */
int	row0, rowN;			/* rows for current pass */

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
	if (!strncmp(s, VIEWSTR, VIEWSTRL)) {
		strcpy(viewspec, s);
		return(1);
	}
	if (isaspect(s)) {
		strcpy(pixasp, s);
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
			fprintf(stderr, "%s: warning - multi-processing requires flat input files\n",
					gargv[0]);
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
		fprintf(stderr, "%s: operation %s needs %d components, has %d\n",
				gargv[0], cmtx->ncols == 1 ? "vector" : "matrix",
				ncomp, cmtx->ncomp);
		return(0);
	}
	if ((in_type != DTrgbe) & (in_type != DTxyze) & (in_type != DTspec) &
			(in_type != DTfloat)) {
		fprintf(stderr, "%s: unsupported input data type '%s'\n",
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

struct hdata {
	int	xr, yr;		/* resolution */
	int	fno;		/* frame # */
	char	fmt[MAXFMTLEN];	/* format */
};

/* check subsequent headers match initial file */
int
checkline(char *s, void *p)
{
	static int	exposWarned = 0;
	struct hdata	*hp = (struct hdata *)p;

	if (!strncmp(s, "NCOLS=", 6)) {
		hp->xr = atoi(s+6);
		if (hp->xr <= 0)
			return(-1);
		return(1);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		hp->yr = atoi(s+6);
		if (hp->yr <= 0)
			return(-1);
		return(1);
	}
	if (!strncmp(s, "FRAME=", 6)) {
		hp->fno = atoi(s+6);
		return(1);
	}
	if (isncomp(s)) {
		if (ncompval(s) != ncomp)
			return(-1);
		return(1);
	}
	if (isexpos(s)) {
		if (!exposWarned && fabs(1. - exposval(s)) > 0.04) {
			fprintf(stderr, "%s: warning - ignoring EXPOSURE setting(s)\n",
					gargv[0]);
			exposWarned++;
		}
		return(1);
	}
	if (formatval(hp->fmt, s))
		return(1);

	return(0);
}

/* open and check input/output file, read/write mode if fno >= 0 */
FILE *
open_iofile(char *fname, int fno)
{
	struct hdata	hd;
	FILE		*fp = fopen(fname, fno>=0 ? "r+b" : "rb");

	if (!fp) {
		fprintf(stderr, "%s: cannot open for reading%s\n",
				fname, fno>=0 ? "/writing" : "");
		return(NULL);
	}
	hd.xr = hd.yr = 0;
	hd.fno = -1;
	hd.fmt[0] = '\0';
	if (getheader(fp, checkline, &hd) < 0) {
		fprintf(stderr, "%s: bad/inconsistent header\n", fname);
		fclose(fp);
		return(NULL);
	}
	if ((hd.fno >= 0) & (fno >= 0) & (hd.fno != fno)) {
		fprintf(stderr, "%s: unexpected frame number (%d != %d)\n",
				fname, hd.fno, fno);
		fclose(fp);
		return(NULL);
	}
	if (strcmp(hd.fmt, cm_fmt_id[fno>=0 ? out_type : in_type])) {
		fprintf(stderr, "%s: wrong format\n", fname);
		fclose(fp);
		return(NULL);
	}
	if ((hd.xr <= 0) | (hd.yr <= 0) &&
			!fscnresolu(&hd.xr, &hd.yr, fp)) {
		fprintf(stderr, "%s: missing resolution\n", fname);
		fclose(fp);
		return(NULL);
	}
	if ((hd.xr != xres) | (hd.yr != yres)) {
		fprintf(stderr, "%s: mismatched resolution\n", fname);
		fclose(fp);
		return(NULL);
	}
	return(fp);
}

/* read in previous pixel data from output and rewind to data start */
int
reload_data(float *osum, FILE *fp)
{
	long	dstart;

	if (!osum | !fp)
		return(0);
	if ((dstart = ftell(fp)) < 0) {
		fprintf(stderr, "%s: ftell() error in reload_data()\n",
				gargv[0]);
		return(0);
	}
	if (out_type == DTfloat) {
		if (fread(osum, sizeof(float)*ncomp, (size_t)xres*yres, fp) !=
				(size_t)xres*yres) {
			fprintf(stderr, "%s: fread() error\n", gargv[0]);
			return(0);
		}
	} else {
		int	y;
		for (y = 0; y < yres; y++, osum += ncomp*xres)
			if (freadsscan(osum, ncomp, xres, fp) < 0) {
				fprintf(stderr, "%s: freadsscan() error\n", gargv[0]);
				return(0);
			}
	}
	if (fseek(fp, dstart, SEEK_SET) < 0) {
		fprintf(stderr, "%s: fseek() error in reload_data()\n",
				gargv[0]);
		return(0);
	}
	return(1);
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
			fprintf(stderr, "%s: cannot start: %s\n", gargv[0], ospec);
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
	printargs(gargc, gargv, fp);	/* this command */
	if (fno >= 0)
		fprintf(fp, "FRAME=%d\n", fno);
	if (viewspec[0])
		fputs(viewspec, fp);
	if (pixasp[0])
		fputs(pixasp, fp);
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
		fprintf(stderr, "%s: unsupported output type!\n", gargv[0]);
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
		fprintf(stderr, "%s: annot allocate %dx%d %d-component accumulator\n",
				gargv[0], xres, yres, ncomp);
		return(0);
	}
	if (sizeof(float) != sizeof(COLORV)) {
		fprintf(stderr, "%s: Code Error 1 in solo_process()\n", gargv[0]);
		return(0);
	}
					/* run through each column/output */
	for (c = 0; c < cmtx->ncols; c++) {
		int	rc = rowN - row0;
		FILE	*fout;
		int	y;
					/* open output (load if multipass) */
		if (out_spec) {		/* file or command */
			if (cmtx->ncols > 1 && !hasFormat(out_spec)) {
				fprintf(stderr, "%s: sequential result must go to stdout\n",
						gargv[0]);
				return(0);
			}
			sprintf(fbuf, out_spec, c);
			if (row0) {	/* another pass -- get prev. data */
				fout = open_iofile(fbuf, c);
				if (!reload_data(osum, fout))
					return(0);
			} else		/* else new output (clobbers prev. file) */
				fout = open_output(fbuf, c-(cmtx->ncols==1));
		} else {			/* else stdout */
			if ((out_type == DTfloat) & (cmtx->ncols > 1)) {
				fprintf(stderr, "%s: float outputs must have separate destinations\n",
						gargv[0]);
				return(0);
			}
			strcpy(fbuf, "<stdout>");
			fout = open_output(NULL, c-(cmtx->ncols==1));
		}
		if (!fout)
			return(0);	/* assume error was reported */
		if (!row0 & (c > 0))	/* clear accumulator? */
			memset(osum, 0, sizeof(float)*ncomp*xres*yres);
		while (rc-- > 0) {	/* run through each input file */
			const int	r = c&1 ? row0 + rc : rowN-1 - rc;
			const rmx_dtype	*cval = rmx_val(cmtx, r, c);
			FILE		*finp;
			int		i, x;
			for (i = ncomp; i--; )
				if (cval[i] != 0) break;
			if (i < 0)	/* this coefficient is zero, skip */
				continue;
			sprintf(fbuf, in_spec, r);
			finp = open_iofile(fbuf, -1);
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
		}			/* write accumulated picture */
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
				fprintf(stderr, "%s: bad status from: %s\n", gargv[0], fbuf);
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

#if defined(_WIN32) || defined(_WIN64)
#define	multi_process	solo_process
#else

/* allocate a scrambled index array of the specified length */
int *
scramble(int n)
{
	int	*scarr = (int *)malloc(sizeof(int)*n);
	int	i;

	if (!scarr) {
		fprintf(stderr, "%s: out of memory in scramble(%d)\n", gargv[0], n);
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
	int	coff = nprocs;
	int	odd = 0;
	float	*osum = NULL;
	int	*syarr = NULL;
	char	fbuf[512];
	int	c;
					/* sanity check */
	if (sizeof(float) != sizeof(COLORV)) {
		fprintf(stderr, "%s: code Error 1 in multi_process()\n", gargv[0]);
		return(0);
	}
	fflush(NULL);			/* parent births helper subprocs */
	while (--coff > 0) {
		int	pid = fork();
		if (pid < 0) {
			fprintf(stderr, "%s: fork() call failed!\n", gargv[0]);
			return(0);
		}
		if (!pid) break;	/* new child gets to work */
	}
	osum = (float *)calloc((size_t)xres*yres, sizeof(float)*ncomp);
	if (!osum) {
		fprintf(stderr, "%s: cannot allocate %dx%d %d-component accumulator\n",
				gargv[0], xres, yres, ncomp);
		return(0);
	}
	srandom(113*coff + 5669);	/* randomize row access for this process */
	syarr = scramble(yres);
					/* run through our unique set of columns */
	for (c = coff; c < cmtx->ncols; c += nprocs) {
		int	rc = rowN - row0;
		FILE	*fout;
		int	y;
					/* create/load output */
		sprintf(fbuf, out_spec, c);
		if (row0) {		/* making another pass? */
			fout = open_iofile(fbuf, c);
			if (!reload_data(osum, fout))
				return(0);
		} else {		/* else new output (clobbers prev. file) */
			fout = open_output(fbuf, c);
			if (!fout) return(0);
			if (c > coff)	/* clear accumulator? */
				memset(osum, 0, sizeof(float)*ncomp*xres*yres);
		}
		while (rc-- > 0) {	/* map & sum each input file */
			const int	r = odd ? row0 + rc : rowN-1 - rc;
			const rmx_dtype	*cval = rmx_val(cmtx, r, c);
			long		dstart;
			size_t		imaplen;
			void		*imap;
			FILE		*finp;
			float		*dst;
			int		i, x;
			for (i = ncomp; i--; )
				if (cval[i] != 0) break;
			if (i < 0)	/* this coefficient is zero, skip */
				continue;
			sprintf(fbuf, in_spec, r);
			finp = open_iofile(fbuf, -1);
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
			imaplen = dstart + (size_t)yres*xres*i;
			imap = mmap(NULL, imaplen, PROT_READ,
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
			munmap(imap, imaplen);
		}			/* write accumulated column picture */
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
				fprintf(stderr, "%s: bad status from: %s\n", gargv[0], fbuf);
				return(0);
			}
		} else if (fclose(fout) == EOF)
			goto writerr;
		odd = !odd;		/* go back & forth to milk page cache */
	}
	if (coff) _exit(0);		/* child exits here */
					/* but parent waits for children */
	free(osum);
	free(syarr);
	c = 0;
	while (++coff < nprocs) {
		int	st;
		if (wait(&st) < 0) {
			fprintf(stderr, "%s: warning - child disappeared\n", gargv[0]);
			break;
		}
		if (st) {
			fprintf(stderr, "%s: bad exit status from child\n", gargv[0]);
			c = st;
		}
	}
	return(c == 0);
writerr:
	fprintf(stderr, "%s: write error\n", fbuf);
	return(0);
}

#endif		/* ! Windows */

int
main(int argc, char *argv[])
{
	double	cacheGB = 0;
	int	rintvl;
	int	a;

	gargc = argc;			/* for header output */
	gargv = argv;

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
		case 'm':		/* cache size in GigaBytes */
			cacheGB = atof(argv[++a]);
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
	cmtx = rmx_load(argv[a+1]);	/* loads from stdin if a+1==argc */
	if (cmtx == NULL)
		return(1);		/* error reported */
	cacheGB *= (cmtx->ncols > 1);
	if (cacheGB > 0 && (!out_spec || *out_spec == '!')) {
		fprintf(stderr, "%s: -m option incompatible with output to %s\n",
				argv[0], out_spec ? "command" : "stdout");
		return(1);
	}
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
	if (cacheGB > 1e-4) {		/* figure out # of passes => rintvl */
		size_t	inp_bytes = (in_type==DTfloat ? sizeof(float)*ncomp
						: (size_t)(ncomp+1)) * xres*yres;
		size_t	over_bytes = sizeof(float)*ncomp*xres*yres +
					2*(out_type==DTfloat ? sizeof(float)*ncomp
						: (size_t)(ncomp+1)) * xres*yres;
		int	npasses = (double)inp_bytes*cmtx->nrows /
				(cacheGB*(1L<<30) - (double)over_bytes*nprocs) + 1;
		if ((npasses <= 0) | (npasses*6 >= cmtx->nrows)) {
			fprintf(stderr,
			    "%s: warning - insufficient cache space for multi-pass\n",
					argv[0]);
			npasses = 1;
		}
		rintvl = cmtx->nrows / npasses;
		rintvl += (rintvl*npasses < cmtx->nrows);
	} else
		rintvl = cmtx->nrows;
					/* make our output accumulation passes */
	for (row0 = 0; row0 < cmtx->nrows; row0 += rintvl) {
		if ((rowN = row0 + rintvl) > cmtx->nrows)
			rowN = cmtx->nrows;
		if (nprocs==1 ? !solo_process() : !multi_process())
			return(1);
	}
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-oc | -of][-o ospec][-N nproc][-m cacheGB] inpspec [mtx]\n",
				argv[0]);
	return(1);
}
