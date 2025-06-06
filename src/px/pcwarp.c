#ifndef lint
static const char	RCSid[] = "$Id: pcwarp.c,v 3.8 2025/06/06 19:11:21 greg Exp $";
#endif
/*
 * Warp colors in Radiance picture to correct for input/output changes.
 */

#include "paths.h"
#include "rtio.h"
#include "resolu.h"
#include "color.h"
#include "warp3d.h"

FILE	*infp = NULL;			/* input stream */
int	xres, yres;			/* input picture resolution */

WARP3D	*cwarp;				/* our warp map */
int	iclip = CGAMUT_UPPER;		/* input value gamut clipping */
int	oclip = CGAMUT_LOWER;		/* output value gamut clipping */

static void syserror(char *s);
static void picwarp(void);


int
main(
	int	argc,
	char	*argv[]
)
{
	static char	picfmt[MAXFMTLEN] = PICFMT;
	int	cwflags = 0;
	int	rval;
	int	i;

	fixargv0(argv[0]);
	infp = stdin;
					/* get options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'e':
			cwflags = W3EXACT;
			break;
		case 'f':
			cwflags = W3FAST;
			break;
		case 'i':
			iclip = 0;
			break;
		case 'o':
			oclip = CGAMUT;
			break;
		default:
			goto userr;
		}
					/* load warp map */
	if (i >= argc)
		goto userr;
	if ((cwarp = load3dw(argv[i++], NULL)) == NULL)
		syserror("load3dw");
	set3dwfl(cwarp, cwflags);
					/* open input and output pictures */
	if (i < argc && (infp = fopen(argv[i], "r")) == NULL)
		syserror(argv[i]);
	if (i < argc-1 && freopen(argv[i+1], "w", stdout) == NULL)
		syserror(argv[i+1]);
					/* transfer header */
	if ((rval = checkheader(infp, picfmt, stdout)) < 0) {
		fprintf(stderr, "%s: input not a Radiance picture\n",
				progname);
		exit(1);
	}
	if (rval)
		fputformat(picfmt, stdout);
					/* add new header info. */
	printargs(i, argv, stdout);
	putchar('\n');
					/* get picture size */
	if ((rval = fgetresolu(&xres, &yres, infp)) < 0) {
		fprintf(stderr, "%s: bad picture size\n", progname);
		exit(1);
	}
					/* new picture size the same */
	fputresolu(rval, xres, yres, stdout);
					/* warp those colors! */
	picwarp();
	exit(0);
userr:
	fprintf(stderr,
		"Usage: %s [-i][-o][-e|-f] map.cwp [input.hdr [output.hdr]]\n",
			progname);
	exit(1);
}


static void
syserror(			/* print system error and exit */
	char	*s
)
{
	fprintf(stderr, "%s: ", progname);
	perror(s);
	exit(2);
}


static void
picwarp(void)			/* warp our picture scanlines */
{
	register COLOR	*scan;
	long	ngamut = 0;
	int	rval;
	int	y;
	register int	x;

	scan = (COLOR *)malloc(xres*sizeof(COLOR));
	if (scan == NULL)
		syserror("picwarp");
	for (y = 0; y < yres; y++) {
		if (freadscan(scan, xres, infp) < 0) {
			fprintf(stderr, "%s: error reading input picture\n",
					progname);
			exit(1);
		}
		for (x = 0; x < xres; x++) {
			if (iclip)
				clipgamut(scan[x], bright(scan[x]), iclip,
						cblack, cwhite);
			rval = warp3d(scan[x], scan[x], cwarp);
			if (rval & W3ERROR)
				syserror("warp3d");
			if (rval & W3BADMAP) {
				fprintf(stderr, "%s: singular color mapping\n",
						progname);
				exit(1);
			}
			if (rval & W3GAMUT)
				ngamut++;
			if (oclip)
				clipgamut(scan[x], bright(scan[x]), oclip,
						cblack, cwhite);
		}
		if (fwritescan(scan, xres, stdout) < 0) {
			fprintf(stderr, "%s: error writing output picture\n",
					progname);
			exit(1);
		}
	}
	if (ngamut >= (long)xres*yres/100)
		fprintf(stderr, "%s: warning - %ld%% of pixels out of gamut\n",
				progname, 100*ngamut/((long)xres*yres));
	free((void *)scan);
}
