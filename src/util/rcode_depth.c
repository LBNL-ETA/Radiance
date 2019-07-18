#ifndef lint
static const char RCSid[] = "$Id: rcode_depth.c,v 2.1 2019/07/18 18:51:56 greg Exp $";
#endif
/*
 * Encode and decode depth values using 16-bit integers
 */

#include "copyright.h"

#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "platform.h"
#include "rtio.h"
#include "fvect.h"
#include "depthcodec.h"

enum {CV_FWD, CV_REV, CV_PTS};


/* Report usage error and exit */
static void
usage_exit(int code)
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(
	" [-d ref_depth/unit][-h[io]][-H[io]][-f[afd]] [input [output]]\n",
			stderr);
	fputs("   Or: ", stderr);
	fputs(progname, stderr);
	fputs(
	" {-r|-p} [-i][-u][-h[io]][-H[io]][-f[afd]] [input [output]]\n",
			stderr);
	exit(code);
}


/* Encode depths from input stream */
static int
encode_depths(DEPTHCODEC *dcp)
{
	long	nexpected = 0;

	if (dcp->inpfmt[0]) {
		if (strcasestr(dcp->inpfmt, "ascii") != NULL)
			dcp->format = 'a';
		else if (strcasestr(dcp->inpfmt, "float") != NULL)
			dcp->format = 'f';
		else if (strcasestr(dcp->inpfmt, "double") != NULL)
			dcp->format = 'd';
		else {
			fputs(dcp->inpname, stderr);
			fputs(": unsupported input format: ", stderr);
			fputs(dcp->inpfmt, stderr);
			fputc('\n', stderr);
			return 0;
		}
	}
	if (dcp->hdrflags & HF_RESIN)
		nexpected = (long)dcp->res.xr * dcp->res.yr;
	do {
		int	ok = 0;
		float	f;
		double	d;

		switch (dcp->format) {
		case 'a':
			ok = (fscanf(dcp->finp, "%lf", &d) == 1);
			break;
		case 'f':
			ok = (getbinary(&f, sizeof(f), 1, dcp->finp) == 1);
			d = f;
			break;
		case 'd':
			ok = (getbinary(&d, sizeof(d), 1, dcp->finp) == 1);
			break;
		}
		if (!ok)
			break;
		putint(depth2code(d, dcp->refdepth), 2, stdout);
	} while (--nexpected);

	if (nexpected > 0) {
		fputs(dcp->inpname, stderr);
		fputs(": fewer depths than expected\n", stderr);
		return 0;
	}
	return 1;
}


/* Convert and output the given depth code to stdout */
static void
output_depth(DEPTHCODEC *dcp, double d)
{
	float	f;

	switch (dcp->format) {
	case 'a':
		printf("%.5e\n", d);
		break;
	case 'f':
		f = d;
		putbinary(&f, sizeof(f), 1, stdout);
		break;
	case 'd':
		putbinary(&d, sizeof(d), 1, stdout);
		break;
	}
}


/* Decode depths from input stream */
static int
decode_depths(DEPTHCODEC *dcp)
{
	long	nexpected = 0;

	if (!check_decode_depths(dcp))
		return 0;

	if (dcp->hdrflags & HF_RESIN)
		nexpected = (long)dcp->res.xr * dcp->res.yr;
	do {
		double	d = decode_depth_next(dcp);		
		if (d < -FTINY)
			break;
		output_depth(dcp, d);

	} while (--nexpected);

	if (nexpected > 0) {
		fputs(dcp->inpname, stderr);
		fputs(": unexpected EOF\n", stderr);
		return 0;
	}
	return 1;
}


/* Decode depths at particular pixels given on standard input */
static int
pixel_depths(DEPTHCODEC *dcp, int unbuf)
{
	int	xy[2];
	double	d;

	if (dcp->finp == stdin) {
		fputs(progname, stderr);
		fputs(": <stdin> used for encoded depths\n", stderr);
		return 0;
	}
	if (!check_decode_depths(dcp))
		return 0;

	while (scanf("%d %d", &xy[0], &xy[1]) == 2) {
		loc2pix(xy, &dcp->res,
			(xy[0]+.5)/dcp->res.xr, (xy[1]+.5)/dcp->res.yr);
		d = decode_depth_pix(dcp, xy[0], xy[1]);
		if (d < -FTINY)
			return 0;
		output_depth(dcp, d);
		if (unbuf && fflush(stdout) == EOF) {
			fputs(progname, stderr);
			fputs(": write error on output\n", stderr);
			return 0;
		}
	}
	if (!feof(stdin)) {
		fputs(progname, stderr);
		fputs(": non-integer as pixel index\n", stderr);
		return 0;
	}
	return 1;
}


/* Output the given world position */
static void
output_worldpos(DEPTHCODEC *dcp, FVECT wpos)
{
	switch (dcp->format) {
	case 'a':
		fprintf(stdout, "%.5e %.5e %.5e\n",
				wpos[0], wpos[1], wpos[2]);
		break;
#ifdef SMLFLT
	case 'f':
		putbinary(wpos, sizeof(*wpos), 3, stdout);
		break;
	case 'd': {
			double	wpd[3];
			VCOPY(wpd, wpos);
			putbinary(wpd, sizeof(*wpd), 3, stdout);
		}
		break;
#else
	case 'f': {
			float	wpf[3];
			VCOPY(wpf, wpos);
			putbinary(wpf, sizeof(*wpf), 3, stdout);
		}
		break;
	case 'd':
		putbinary(wpos, sizeof(*wpos), 3, stdout);
		break;
#endif
	}
}


/* Decode world points from input stream */
static int
decode_points(DEPTHCODEC *dcp)
{
	long	n = (long)dcp->res.xr * dcp->res.yr;

	if (!check_decode_worldpos(dcp))
		return 0;

	while (n-- > 0) {
		FVECT	wpos;
		if (decode_worldpos_next(wpos, dcp) < 0) {
			if (feof(dcp->finp)) {
				fputs(dcp->inpname, stderr);
				fputs(": unexpected EOF\n", stderr);
			}
			return 0;
		}
		output_worldpos(dcp, wpos);
	}
	return 1;	
}


/* Decode world points at particular pixels given on standard input */
static int
pixel_points(DEPTHCODEC *dcp, int unbuf)
{
	int	xy[2];
	FVECT	wpos;

	if (dcp->finp == stdin) {
		fputs(progname, stderr);
		fputs(": <stdin> used for encoded depths\n", stderr);
		return 0;
	}
	if (!check_decode_worldpos(dcp))
		return 0;
	
	while (scanf("%d %d", &xy[0], &xy[1]) == 2) {
		loc2pix(xy, &dcp->res,
			(xy[0]+.5)/dcp->res.xr, (xy[1]+.5)/dcp->res.yr);
		if (get_worldpos_pix(wpos, dcp, xy[0], xy[1]) < 0)
			return 0;
		output_worldpos(dcp, wpos);
		if (unbuf && fflush(stdout) == EOF) {
			fputs(progname, stderr);
			fputs(": write error on output\n", stderr);
			return 0;
		}
	}
	if (!feof(stdin)) {
		fputs(progname, stderr);
		fputs(": non-integer as pixel index\n", stderr);
		return 0;
	}
	return 1;
}


int
main(int argc, char *argv[])
{
	int		conversion = CV_FWD;
	int		bypixel = 0;
	int		unbuffered = 0;
	DEPTHCODEC	dc;
	int		a;

	progname = argv[0];
	set_dc_defaults(&dc);
	dc.hdrflags = HF_ALL;
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'd':
			strlcpy(dc.depth_unit, argv[++a], sizeof(dc.depth_unit));
			dc.refdepth = atof(dc.depth_unit);
			if (dc.refdepth <= .0) {
				fputs(progname, stderr);
				fputs(": illegal -d reference depth\n", stderr);
				return 1;
			}
			break;
		case 'r':
			conversion = CV_REV;
			break;
		case 'p':
			conversion = CV_PTS;
			break;
		case 'f':
			switch (argv[a][2]) {
			case 'f':
			case 'd':
			case 'a':
				dc.format = argv[a][2];
				break;
			default:
				fputs(progname, stderr);
				fputs(": unsupported -f? format\n", stderr);
				usage_exit(1);
			}
			break;
		case 'h':
			switch (argv[a][2]) {
			case '\0':
				dc.hdrflags &= ~(HF_HEADIN|HF_HEADOUT);
				break;
			case 'i':
			case 'I':
				dc.hdrflags &= ~HF_HEADIN;
				break;
			case 'o':
			case 'O':
				dc.hdrflags &= ~HF_HEADOUT;
				break;
			default:
				usage_exit(1);
			}
			break;
		case 'H':
			switch (argv[a][2]) {
			case '\0':
				dc.hdrflags &= ~(HF_RESIN|HF_RESOUT);
				break;
			case 'i':
			case 'I':
				dc.hdrflags &= ~HF_RESIN;
				break;
			case 'o':
			case 'O':
				dc.hdrflags &= ~HF_RESOUT;
				break;
			default:
				usage_exit(1);
			}
			break;
		case 'i':
			bypixel++;
			break;
		case 'u':
			unbuffered++;
			break;
		default:
			usage_exit(1);
		}
	dc.hdrflags |= (conversion == CV_FWD) * HF_ENCODE;

	if ((dc.hdrflags & (HF_RESIN|HF_RESOUT)) == HF_RESOUT) {
		fputs(progname, stderr);
		fputs(": unknown resolution for output\n", stderr);
		return 1;
	}
	if (bypixel) {
		if (conversion == CV_FWD) {
			fputs(progname, stderr);
			fputs(": -i only supported with -r or -p option\n",
					stderr);
			usage_exit(1);
		}
		if (a >= argc) {
			fputs(progname, stderr);
			fputs(": -i option requires input file\n", stderr);
			usage_exit(1);
		}
		if (!(dc.hdrflags & HF_RESIN)) {
			fputs(progname, stderr);
			fputs(": -i option requires input resolution\n", stderr);
			usage_exit(1);
		}
		dc.hdrflags &= ~(HF_HEADOUT|HF_RESOUT);
	}
	if (a < argc-2) {
		fputs(progname, stderr);
		fputs(": too many file arguments\n", stderr);
		usage_exit(1);
	}
	if (a < argc && !(dc.finp = fopen(dc.inpname=argv[a], "r"))) {
		fputs(dc.inpname, stderr);
		fputs(": cannot open for reading\n", stderr);
		return 1;
	}
	if (a+1 < argc && !freopen(argv[a+1], "w", stdout)) {
		fputs(argv[a+1], stderr);
		fputs(": cannot open for writing\n", stderr);
		return 1;
	}
	SET_FILE_BINARY(dc.finp);
	SET_FILE_BINARY(stdout);
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(dc.finp);
	flockfile(stdout);
#endif
					/* read/copy header */
	if (!process_dc_header(&dc, argc, argv))
		return 1;
					/* process data */
	switch (conversion) {
	case CV_FWD:			/* distance -> depth code */
		if (!encode_depths(&dc))
			return 1;
		break;
	case CV_REV:			/* depth code -> distance */
		if (!(bypixel ? pixel_depths(&dc, unbuffered) :
				decode_depths(&dc)))
			return 1;
		break;
	case CV_PTS:			/* depth code -> world points */
		if (!(bypixel ? pixel_points(&dc, unbuffered) :
				decode_points(&dc)))
			return 1;
		break;
	}
	if (fflush(stdout) == EOF) {
		fputs(progname, stderr);
		fputs(": error writing output\n", stderr);
		return 1;
	}
	return 0;
}
