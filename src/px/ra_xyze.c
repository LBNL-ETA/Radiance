#ifndef lint
static const char	RCSid[] = "$Id: ra_xyze.c,v 2.17 2025/06/07 05:09:46 greg Exp $";
#endif
/*
 *  Program to convert between RADIANCE RGBE and XYZE formats
 *  Added white-balance adjustment 10/01 (GW).
 */

#include  <math.h>
#include  "platform.h"
#include  "color.h"
#include  "resolu.h"
#include  "rtio.h"

enum {InpUNK, InpRGB, InpXYZ, InpSPEC};	/* input format */
int  infmt = InpUNK;
int  rgbout = 0;			/* output should be RGBE? */
RGBPRIMS  inprims = STDPRIMS;		/* input primaries */
RGBPRIMS  outprims = STDPRIMS;		/* output primaries */
double	expcomp = 1.0;			/* exposure compensation */
int  doflat = -1;			/* produce flat file? */
double  origexp = -1.0;			/* original exposure */

static gethfunc headline;
static void quiterr(char *err);
static void myreadscan(COLOR *scn, int len);
static void convert(void);



static int
headline(				/* process header line */
	char	*s,
	void	*p
)
{
	char	fmt[MAXFMTLEN];

	if (formatval(fmt, s)) {	/* check if format string */
		if (!strcmp(fmt,COLRFMT))
			infmt = InpRGB;
		else if (!strcmp(fmt,CIEFMT))
			infmt = InpXYZ;
		else if (!strcmp(fmt,SPECFMT))
			infmt = InpSPEC;
		else
			infmt = InpUNK;
		return(0);		/* don't echo */
	}
	if (origexp > 0.0 && isexpos(s)) {
		origexp *= exposval(s);
		return(0);		/* don't echo */
	}
	if (isprims(s)) {		/* get input primaries */
		primsval(inprims, s);
		return(0);		/* don't echo */
	}
	if (iswlsplit(s)) {		/* get wavelength limits */
		wlsplitval(WLPART, s);
		return(0);		/* don't echo */
	}
	if (isncomp(s)) {		/* number of color samples */
		NCSAMP = ncompval(s);
		return(0);		/* don't echo */
	}
					/* should I grok colcorr also? */
	return(fputs(s, stdout));
}


int
main(int  argc, char  *argv[])
{
	int  i;

	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	fixargv0(argv[0]);

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'c':		/* rle-compressed output */
				doflat = 0;
				break;
			case 'u':		/* flat output */
				doflat = 1;
				break;
			case 'r':		/* RGBE output */
				rgbout = 1;
				break;
			case 'p':		/* RGB primaries */
				if (i+8 >= argc)
					goto userr;
				outprims[RED][CIEX] = atof(argv[++i]);
				outprims[RED][CIEY] = atof(argv[++i]);
				outprims[GRN][CIEX] = atof(argv[++i]);
				outprims[GRN][CIEY] = atof(argv[++i]);
				outprims[BLU][CIEX] = atof(argv[++i]);
				outprims[BLU][CIEY] = atof(argv[++i]);
				outprims[WHT][CIEX] = atof(argv[++i]);
				outprims[WHT][CIEY] = atof(argv[++i]);
				break;
			case 'o':		/* original exposure */
				origexp = 1.0;
				break;
			case 'e':		/* exposure compensation */
				expcomp = atof(argv[++i]);
				if (argv[i][0] == '+' || argv[i][0] == '-')
					expcomp = pow(2., expcomp);
				break;
			default:
				goto userr;
			}
		else
			break;

	if (doflat < 0)
		doflat = !rgbout;
	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				progname, argv[i]);
		exit(1);
	}
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		fprintf(stderr, "%s: can't open output \"%s\"\n",
				progname, argv[i+1]);
		exit(1);
	}
	getheader(stdin, headline, NULL);
	if (infmt == InpUNK)
		quiterr("unrecognized/missing input file format");
	if (infmt == InpSPEC ? (NCSAMP <= 3) | (NCSAMP > MAXCSAMP) :
			NCSAMP != 3)
		quiterr("bad number of color components");
	printargs(argc, argv, stdout);		/* add to header */
	convert();				/* convert picture */
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-r][-o][-e exp][-c|-u]", progname);
	fprintf(stderr, "[-p rx ry gx gy bx by wx wy] [input [output]]\n");
	exit(1);
}


static void
quiterr(		/* print message and exit */
	char  *err
)
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


static void
myreadscan(COLOR *scn, int len)
{
	if (infmt == InpSPEC) {		/* read & convert to XYZ */
		static COLOR	*scomp = NULL;
		SCOLR		sclr;
		SCOLOR		scol;
		COLOR		xyz;
		int		n;
		if (scomp == NULL) {	/* initialize conversion */
			scomp = (COLOR *)malloc(sizeof(COLOR)*NCSAMP);
			if (scomp == NULL)
				quiterr("out of memory in myreadscan");
			for (n = NCSAMP; n--; )
				spec_cie(scomp[n],
					WLPART[0] + (WLPART[3] - WLPART[0])*(n+1)/NCSAMP,
					WLPART[0] + (WLPART[3] - WLPART[0])*n/NCSAMP);
		}
		while (len-- > 0) {
			if (getbinary(sclr, LSCOLR, 1, stdin) != 1)
				goto readerr;
			scolr_scolor(scol, sclr);
			setcolor(*scn, 0, 0, 0);
			for (n = NCSAMP; n--; ) {
				copycolor(xyz, scomp[n]);
				scalecolor(xyz, scol[n]);
				addcolor(*scn, xyz);
			}
			scn++;
		}
		return;
	}				/* else read as RGBE/XYZE */
	if (freadscan(scn, len, stdin) >= 0)
		return;
readerr:
	quiterr("error reading input picture");
}


static void
convert(void)				/* convert to XYZE or RGBE picture */
{
	int	order;
	int	xmax, ymax;
	COLORMAT	xfm;
	COLOR	*scanin;
	COLR	*scanout;
	double	exp2do = expcomp;
	double	exp2report = expcomp;
	int	y;
	int	x;
						/* recover original? */
	if (origexp > 0.0)
		exp2do /= origexp;
						/* compute transform */
	if (rgbout) {
		if (infmt == InpRGB) {		/* RGBE -> RGBE */
			comprgb2rgbWBmat(xfm, inprims, outprims);
		} else {			/* XYZE/Spectral -> RGBE */
			compxyz2rgbWBmat(xfm, outprims);
		}
		if (infmt == InpXYZ) {
			if (origexp > 0.0)
				exp2do /= WHTEFFICACY;
			else
				exp2report *= WHTEFFICACY;
		}
	} else {
		if (infmt == InpRGB) {		/* RGBE -> XYZE */
			comprgb2xyzWBmat(xfm, inprims);
		} else {			/* XYZE/Spectral -> XYZE */
			memset(xfm, 0, sizeof(xfm));
			for (x = 3; x--; )
				xfm[x][x] = 1.;
		}
		if (infmt != InpXYZ) {
			if (origexp > 0)
				exp2do *= WHTEFFICACY;
			else
				exp2report /= WHTEFFICACY;
		}
	}
	for (y = 0; y < 3; y++)
		for (x = 0; x < 3; x++)
			xfm[y][x] *= exp2do;
						/* get input resolution */
	if ((order = fgetresolu(&xmax, &ymax, stdin)) < 0)
		quiterr("bad picture format");
						/* complete output header */
	if ((exp2report < 0.99) | (exp2report > 1.01))
		fputexpos(exp2report, stdout);
	if (rgbout) {
		fputprims(outprims, stdout);
		fputformat(COLRFMT, stdout);
	} else
		fputformat(CIEFMT, stdout);
	putc('\n', stdout);
	fputresolu(order, xmax, ymax, stdout);
						/* allocate scanline */
	scanin = (COLOR *)malloc(xmax*sizeof(COLOR));
	if (scanin == NULL)
		quiterr("out of memory in convert");
	scanout = doflat ? (COLR *)malloc(xmax*sizeof(COLR)) : (COLR *)NULL;
						/* convert image */
	for (y = 0; y < ymax; y++) {
		myreadscan(scanin, xmax);
		for (x = 0; x < xmax; x++) {
			colortrans(scanin[x], xfm, scanin[x]);
			if (rgbout)
				clipgamut(scanin[x], bright(scanin[x]),
						CGAMUT_LOWER, cblack, cwhite);
		}
		if (scanout != NULL) {
			for (x = 0; x < xmax; x++)
				setcolr(scanout[x], colval(scanin[x],RED),
						colval(scanin[x],GRN),
						colval(scanin[x],BLU));
			putbinary(scanout, sizeof(COLR), xmax, stdout);
		} else
			fwritescan(scanin, xmax, stdout);
		if (ferror(stdout))
			quiterr("error writing output picture");
	}
						/* free scanline */
	free(scanin);
	if (scanout != NULL)
		free(scanout);
}
