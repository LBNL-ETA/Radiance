#ifndef lint
static const char RCSid[] = "$Id: normcodec.c,v 2.10 2025/06/07 05:09:45 greg Exp $";
#endif
/*
 * Routines to encode/decode 32-bit normals
 */

#include "copyright.h"

#include "rtio.h"
#include "paths.h"
#include "rtmath.h"
#include "normcodec.h"


/* Set codec defaults */
void
set_nc_defaults(NORMCODEC *ncp)
{
	memset(ncp, 0, sizeof(NORMCODEC));
	ncp->finp = stdin;
	ncp->inpname = "<stdin>";
	ncp->format = 'a';
	ncp->res.rt = PIXSTANDARD;
	if (!progname) progname = "norm_codec";
}


/* process header line */
static int
headline(char *s, void *p)
{
	NORMCODEC	*ncp = (NORMCODEC *)p;
	int		rv;

	if (formatval(ncp->inpfmt, s))	/* don't pass format */
		return 0;

	if ((rv = isbigendian(s)) >= 0) {
		ncp->swapped = (nativebigendian() != rv);
		return 0;
	}
	if (!strncmp(s, "NCOMP=", 6)) {
		if (atoi(s+6) != 3) {
			if (ncp->hdrflags & HF_STDERR) {
				fputs(ncp->inpname, stderr);
				fputs(": NCOMP must equal 3\n", stderr);
			}
			return -1;
		}
		return 0;
	}
	if (!strncmp(s, "NROWS=", 6)) {
		ncp->res.yr = atoi(s+6);
		return 0;
	}
	if (!strncmp(s, "NCOLS=", 6)) {
		ncp->res.xr = atoi(s+6);
		return 0;
	}
	if (ncp->hdrflags & HF_HEADOUT)
		fputs(s, stdout);	/* copy to standard output */
	return 1;
}


/* Load/copy header */
int
process_nc_header(NORMCODEC *ncp, int ac, char *av[])
{
	if (ncp->hdrflags & HF_HEADIN &&
			getheader(ncp->finp, headline, ncp) < 0) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(ncp->inpname, stderr);
			fputs(": bad header\n", stderr);
		}
		return 0;
	}
					/* get resolution string? */
	if (ncp->hdrflags & HF_RESIN &&
			(ncp->res.xr <= 0) | (ncp->res.yr <= 0) &&
			!fgetsresolu(&ncp->res, ncp->finp)) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(ncp->inpname, stderr);
			fputs(": bad resolution string\n", stderr);
		}
		return 0;
	}
	if (ncp->hdrflags & HF_HEADOUT) {	/* finish header */
		if (!(ncp->hdrflags & HF_HEADIN))
			newheader("RADIANCE", stdout);
		if (ac > 0)
			printargs(ac, av, stdout);
		if (ncp->hdrflags & HF_ENCODE) {
			fputformat(NORMAL32FMT, stdout);
		} else {
			fputs("NCOMP=3\n", stdout);
			if ((ncp->hdrflags & (HF_RESIN|HF_RESOUT)) == HF_RESIN)
				printf("NCOLS=%d\nNROWS=%d\n",
						scanlen(&ncp->res),
						numscans(&ncp->res));
			switch (ncp->format) {
			case 'a':
				fputformat("ascii", stdout);
				break;
			case 'f':
				fputendian(stdout);
				fputformat("float", stdout);
				break;
			case 'd':
				fputendian(stdout);
				fputformat("double", stdout);
				break;
			}
		}
		fputc('\n', stdout);
	}
	if (ncp->hdrflags & HF_RESOUT)	/* put resolution string? */
		fputsresolu(&ncp->res, stdout);

	ncp->dstart = ftell(ncp->finp);
	return 1;
}


/* Check that we have what we need to decode normals */
int
check_decode_normals(NORMCODEC *ncp)
{
	if (ncp->hdrflags & HF_ENCODE) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(progname, stderr);
			fputs(": wrong header mode for decode\n", stderr);
		}
		return 0;
	}
	if (ncp->inpfmt[0] && strcmp(ncp->inpfmt, NORMAL32FMT)) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(ncp->inpname, stderr);
			fputs(": unexpected input format: ", stderr);
			fputs(ncp->inpfmt, stderr);
			fputc('\n', stderr);
		}
		return 0;
	}
	return 1;
}


/* Decode next normal from input */
int
decode_normal_next(FVECT nrm, NORMCODEC *ncp)
{
	static int32	lastc;
	static FVECT	lastv;
	int32		c = getint(4, ncp->finp);

	if (c == EOF && feof(ncp->finp))
		return -1;

	if (c == lastc) {			/* optimization */
		VCOPY(nrm, lastv);
	} else {
		decodedir(nrm, c);
		if (c) {
			lastc = c;
			VCOPY(lastv, nrm);
		}
	}
	return (c != 0);
}


/* Seek to the indicated pixel position */
int
seek_nc_pix(NORMCODEC *ncp, int x, int y)
{
	if ((ncp->res.xr <= 0) | (ncp->res.yr <= 0)) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(progname, stderr);
			fputs(": need map resolution to seek\n", stderr);
		}
		return -1;
	}
	if ((x < 0) | (y < 0) ||
			(x >= scanlen(&ncp->res)) | (y >= numscans(&ncp->res))) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(ncp->inpname, stderr);
			fputs(": warning - pixel index off map\n", stderr);
		}
		return 0;
	}
	if (fseek(ncp->finp, ncp->dstart + 4*((long)y*scanlen(&ncp->res) + x),
			SEEK_SET) == EOF) {
		if (ncp->hdrflags & HF_STDERR) {
			fputs(ncp->inpname, stderr);
			fputs(": seek error\n", stderr);
		}
		return -1;
	}
	return 1;
}


/* Read and decode normal for the given pixel */
int
decode_normal_pix(FVECT nrm, NORMCODEC *ncp, int x, int y)
{
	int	rval = seek_nc_pix(ncp, x, y);

	if (rval < 0)
		return -1;

	if (!rval) {
		nrm[0] = nrm[1] = nrm[2] = .0;
		return 0;
	}
	return decode_normal_next(nrm, ncp);
}
