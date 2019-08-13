#ifndef lint
static const char RCSid[] = "$Id: rcode_norm.c,v 2.3 2019/08/13 16:31:35 greg Exp $";
#endif
/*
 * Encode and decode surface normal map using 32-bit integers
 */

#include "copyright.h"

#include <stdlib.h>
#include "platform.h"
#include "rtio.h"
#include "rtmath.h"
#include "normcodec.h"

char		*progname;		/* set in main() */


/* Report usage error and exit */
static void
usage_exit(int code)
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-h[io]][-H[io]][-f[afd]] [input [output.nrm]]\n", stderr);
	fputs("   Or: ", stderr);
	fputs(progname, stderr);
	fputs(" -r [-i][-u][-h[io]][-H[io]][-f[afd]] [input.nrm [output]]\n",
			stderr);
	exit(code);
}


/* Encode normals from input stream */
static int
encode_normals(NORMCODEC *ncp)
{
	long	nexpected = (long)ncp->res.xr * ncp->res.yr;

	if (ncp->inpfmt[0]) {
		if (strstr(ncp->inpfmt, "ascii") != NULL)
			ncp->format = 'a';
		else if (strstr(ncp->inpfmt, "float") != NULL)
			ncp->format = 'f';
		else if (strstr(ncp->inpfmt, "double") != NULL)
			ncp->format = 'd';
		else {
			fputs(ncp->inpname, stderr);
			fputs(": unsupported input format: ", stderr);
			fputs(ncp->inpfmt, stderr);
			fputc('\n', stderr);
			return 0;
		}
	}

	do {
		int	ok = 0;
		FVECT	nrm;

		switch (ncp->format) {
		case 'a':
			ok = (fscanf(ncp->finp, FVFORMAT,
					&nrm[0], &nrm[1], &nrm[2]) == 3);
			break;
#ifdef SMLFLT
		case 'f':
			ok = (getbinary(nrm, sizeof(*nrm), 3, ncp->finp) == 3);
			break;
		case 'd': {
				double	nrmd[3];
				ok = (getbinary(nrmd, sizeof(*nrmd),
						3, ncp->finp) == 3);
				if (ok) VCOPY(nrm, nrmd);
			}
			break;
#else
		case 'f': {
				float	nrmf[3];
				ok = (getbinary(nrmf, sizeof(*nrmf),
						3, ncp->finp) == 3);
				if (ok) VCOPY(nrm, nrmf);
			}
			break;
		case 'd':
			ok = (getbinary(nrm, sizeof(*nrm), 3, ncp->finp) == 3);
			break;
#endif
		}
		if (!ok)
			break;

		normalize(nrm);
		putint(encodedir(nrm), 4, stdout);

	} while (--nexpected);

	if (nexpected > 0) {
		fputs(ncp->inpname, stderr);
		fputs(": fewer normals than expected\n", stderr);
		return 0;
	}
	return 1;
}


/* Output the given normal to stdout */
static void
output_normal(NORMCODEC *ncp, FVECT nrm)
{
	switch (ncp->format) {
	case 'a':
		printf("%.9f %.9f %.9f\n", nrm[0], nrm[1], nrm[2]);
		break;
#ifdef SMLFLT
	case 'f':
		putbinary(nrm, sizeof(*nrm), 3, stdout);
		break;
	case 'd': {
			double	nrmd[3];
			VCOPY(nrmd, nrm);
			putbinary(nrmd, sizeof(*nrmd), 3, stdout);
		}
		break;
#else
	case 'f': {
			float	nrmf[3];
			VCOPY(nrmf, nrm);
			putbinary(nrmf, sizeof(*nrmf), 3, stdout);
		}
		break;
	case 'd':
		putbinary(nrm, sizeof(*nrm), 3, stdout);
		break;
#endif
	}
}


/* Decode normals from input stream */
static int
decode_normals(NORMCODEC *ncp)
{
	long	nexpected = (long)ncp->res.xr * ncp->res.yr;

	if (!check_decode_normals(ncp))
		return 0;

	do {
		FVECT	nrm;

		if (decode_normal_next(nrm, ncp) < 0)
			break;

		output_normal(ncp, nrm);

	} while (--nexpected);

	if (nexpected > 0) {
		fputs(ncp->inpname, stderr);
		fputs(": unexpected EOF\n", stderr);
		return 0;
	}
	return 1;
}


/* Decode normals at particular pixels given on standard input */
static int
pixel_normals(NORMCODEC *ncp, int unbuf)
{
	int	x, y;
	FVECT	nrm;

	if (ncp->finp == stdin) {
		fputs(progname, stderr);
		fputs(": <stdin> used for encoded normals\n", stderr);
		return 0;
	}
	if (!check_decode_normals(ncp))
		return 0;
	if (ncp->res.rt != PIXSTANDARD) {
		fputs(progname, stderr);
		fputs(": can only handle standard pixel ordering\n", stderr);
		return 0;
	}
	while (scanf("%d %d", &x, &y) == 2) {

		y = ncp->res.yr-1 - y;

		if (decode_normal_pix(nrm, ncp, x, y) < 0)
			return 0;

		output_normal(ncp, nrm);

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
	int		reverse = 0;
	int		bypixel = 0;
	int		unbuffered = 0;
	NORMCODEC	nc;
	int		a;

	progname = argv[0];
	set_nc_defaults(&nc);
	nc.hdrflags = HF_ALL;
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'r':
			reverse++;
			break;
		case 'f':
			switch (argv[a][2]) {
			case 'f':
			case 'd':
			case 'a':
				nc.format = argv[a][2];
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
				nc.hdrflags &= ~(HF_HEADIN|HF_HEADOUT);
				break;
			case 'i':
			case 'I':
				nc.hdrflags &= ~HF_HEADIN;
				break;
			case 'o':
			case 'O':
				nc.hdrflags &= ~HF_HEADOUT;
				break;
			default:
				usage_exit(1);
			}
			break;
		case 'H':
			switch (argv[a][2]) {
			case '\0':
				nc.hdrflags &= ~(HF_RESIN|HF_RESOUT);
				break;
			case 'i':
			case 'I':
				nc.hdrflags &= ~HF_RESIN;
				break;
			case 'o':
			case 'O':
				nc.hdrflags &= ~HF_RESOUT;
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
	nc.hdrflags |= !reverse * HF_ENCODE;

	if ((nc.hdrflags & (HF_RESIN|HF_RESOUT)) == HF_RESOUT) {
		fputs(progname, stderr);
		fputs(": unknown resolution for output\n", stderr);
		return 1;
	}
	if (bypixel) {
		if (!reverse) {
			fputs(progname, stderr);
			fputs(": -i only supported with -r option\n",
					stderr);
			usage_exit(1);
		}
		if (a >= argc) {
			fputs(progname, stderr);
			fputs(": -i option requires input file\n", stderr);
			usage_exit(1);
		}
		if (!(nc.hdrflags & HF_RESIN)) {
			fputs(progname, stderr);
			fputs(": -i option requires input resolution\n", stderr);
			usage_exit(1);
		}
		nc.hdrflags &= ~HF_RESOUT;
	}
	if (a < argc-2) {
		fputs(progname, stderr);
		fputs(": too many file arguments\n", stderr);
		usage_exit(1);
	}
	if (a < argc && !(nc.finp = fopen(nc.inpname=argv[a], "r"))) {
		fputs(nc.inpname, stderr);
		fputs(": cannot open for reading\n", stderr);
		return 1;
	}
	if (a+1 < argc && !freopen(argv[a+1], "w", stdout)) {
		fputs(argv[a+1], stderr);
		fputs(": cannot open for writing\n", stderr);
		return 1;
	}
	SET_FILE_BINARY(nc.finp);
	if (reverse || nc.format != 'a')
		SET_FILE_BINARY(stdout);
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(nc.finp);
	flockfile(stdout);
#endif
					/* read/copy header */
	if (!process_nc_header(&nc, a, argv))
		return 1;
					/* process data */
	if (!reverse) {
		if (!encode_normals(&nc))
			return 1;
	} else if (bypixel) {
		if (!pixel_normals(&nc, unbuffered))
			return 1;
	} else if (!decode_normals(&nc))
		return 1;

	if (fflush(stdout) == EOF) {
		fputs(progname, stderr);
		fputs(": error writing output\n", stderr);
		return 1;
	}
	return 0;
}
