#ifndef lint
static const char	RCSid[] = "$Id: pcond.c,v 3.24 2006/06/07 17:52:04 schorsch Exp $";
#endif
/*
 * Condition Radiance picture for display/output
 *  Added white-balance adjustment 10/01 (GW).
 */

#include "platform.h"
#include "paths.h"
#include "rtprocess.h"
#include "pcond.h"


#define	LDMAX		100		/* default max. display luminance */
#define LDDYN		32		/* default dynamic range */

int	what2do = 0;			/* desired adjustments */

double	ldmax = LDMAX;			/* maximum output luminance */
double	lddyn = LDDYN;			/* display dynamic range */
double	Bldmin, Bldmax;			/* Bl(ldmax/lddyn) and Bl(ldmax) */

char	*progname;			/* global argv[0] */

char	*infn;				/* input file name */
FILE	*infp;				/* input stream */
FILE	*mapfp = NULL;			/* tone-mapping function stream */
VIEW	ourview = STDVIEW;		/* picture view */
int	gotview = 0;			/* picture has view */
double	pixaspect = 1.0;		/* pixel aspect ratio */
double	fixfrac = 0.;			/* histogram share due to fixations */
RESOLU	inpres;				/* input picture resolution */

COLOR	*fovimg;			/* foveal (1 degree) averaged image */
int	fvxr, fvyr;			/* foveal image resolution */
float	*crfimg;			/* contrast reduction factors */
short	(*fixlst)[2];			/* fixation history list */
int	nfixations;			/* number of fixation points */
double	bwhist[HISTRES];		/* luminance histogram */
double	histot;				/* total count of histogram */
double	bwmin, bwmax;			/* histogram limits */
double	bwavg;				/* mean brightness */

double	scalef = 0.;			/* linear scaling factor */

static gethfunc headline;
static void getahead(void);
static void mapimage(void);
static void getfovimg(void);
static void check2do(void);



int
main(
	int	argc,
	char	*argv[]
)
{
	static RGBPRIMS	outprimS;
	int	i;
#define	bool(flg)		switch (argv[i][2]) { \
				case '\0': what2do ^= flg; break; \
				case 'y': case 'Y': case 't': case 'T': \
				case '+': case '1': what2do |= flg; break; \
				case 'n': case 'N': case 'f': case 'F': \
				case '-': case '0': what2do &= ~(flg); break; \
				default: goto userr; }

	progname = argv[0];

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'h':
			bool(DO_HUMAN);
			break;
		case 'a':
			bool(DO_ACUITY);
			break;
		case 'v':
			bool(DO_VEIL);
			break;
		case 's':
			bool(DO_HSENS);
			break;
		case 'c':
			bool(DO_COLOR);
			break;
		case 'w':
			bool(DO_CWEIGHT);
			break;
		case 'i':
			if (i+1 >= argc) goto userr;
			fixfrac = atof(argv[++i]);
			if (fixfrac > FTINY) what2do |= DO_FIXHIST;
			else what2do &= ~DO_FIXHIST;
			break;
		case 'I':
			bool(DO_PREHIST);
			break;
		case 'l':
			bool(DO_LINEAR);
			break;
		case 'p':
			if (i+8 >= argc) goto userr;
			outprimS[RED][CIEX] = atof(argv[++i]);
			outprimS[RED][CIEY] = atof(argv[++i]);
			outprimS[GRN][CIEX] = atof(argv[++i]);
			outprimS[GRN][CIEY] = atof(argv[++i]);
			outprimS[BLU][CIEX] = atof(argv[++i]);
			outprimS[BLU][CIEY] = atof(argv[++i]);
			outprimS[WHT][CIEX] = atof(argv[++i]);
			outprimS[WHT][CIEY] = atof(argv[++i]);
			outprims = outprimS;
			break;
		case 'e':
			if (i+1 >= argc) goto userr;
			scalef = atof(argv[++i]);
			if ((argv[i][0] == '+') | (argv[i][0] == '-'))
				scalef = pow(2.0, scalef);
			what2do |= DO_LINEAR;
			break;
		case 'f':
			if (i+1 >= argc) goto userr;
			mbcalfile = argv[++i];
			break;
		case 'm':
			if (i+1 >= argc) goto userr;
			cwarpfile = argv[++i];
			break;
		case 'u':
			if (i+1 >= argc) goto userr;
			ldmax = atof(argv[++i]);
			if (ldmax <= FTINY)
				goto userr;
			break;
		case 'd':
			if (i+1 >= argc) goto userr;
			lddyn = atof(argv[++i]);
			break;
		case 'x':
			if (i+1 >= argc) goto userr;
			if ((mapfp = fopen(argv[++i], "w")) == NULL) {
				fprintf(stderr,
					"%s: cannot open for writing\n",
						argv[i]);
				exit(1);
			}
			break;
		default:
			goto userr;
		}
	if ((what2do & (DO_FIXHIST|DO_PREHIST)) == (DO_FIXHIST|DO_PREHIST)) {
		fprintf(stderr, "%s: only one of -i or -I option\n", progname);
		exit(1);
	}
	if ((mbcalfile != NULL) + (cwarpfile != NULL) +
			(outprims != stdprims) > 1) {
		fprintf(stderr,
			"%s: only one of -p, -m or -f option supported\n",
				progname);
		exit(1);
	}
	if ((outprims == stdprims) & (inprims != stdprims))
		outprims = inprims;
	Bldmin = Bl(ldmax/lddyn);
	Bldmax = Bl(ldmax);
	if (i >= argc || i+2 < argc)
		goto userr;
					/* open input file */
	if ((infp = fopen(infn=argv[i], "r")) == NULL)
		syserror(infn);
					/* open output file */
	if (i+2 == argc && freopen(argv[i+1], "w", stdout) == NULL)
		syserror(argv[i+1]);
	SET_FILE_BINARY(infp);
	SET_FILE_BINARY(stdout);
	getahead();			/* load input header */
	printargs(argc, argv, stdout);	/* add to output header */
	if ((mbcalfile == NULL) & (outprims != stdprims))
		fputprims(outprims, stdout);
	if ((what2do & (DO_PREHIST|DO_VEIL|DO_ACUITY)) != DO_PREHIST)
		getfovimg();		/* get foveal sample image? */
	if (what2do&DO_PREHIST)		/* get histogram? */
		gethisto(stdin);
	else if (what2do&DO_FIXHIST)	/* get fixation history? */
		getfixations(stdin);
	mapimage();			/* map the picture */
	if (mapfp != NULL)		/* write out basic mapping */
		putmapping(mapfp);
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-{h|a|v|s|c|l|w}[+-]][-I|-i ffrac][-e ev][-p xr yr xg yg xb yb xw yw|-f mbf.cal|-m rgb.cwp][-u Ldmax][-d Lddyn][-x mapfile] inpic [outpic]\n",
			progname);
	exit(1);
	return 1; /* pro forma return */
#undef bool
}


extern void
syserror(				/* report system error and exit */
	char	*s
)
{
	fprintf(stderr, "%s: ", progname);
	perror(s);
	exit(2);
}


static int
headline(				/* process header line */
	char	*s,
	void	*p
)
{
	static RGBPRIMS	inprimS;
	char	fmt[32];

	if (formatval(fmt, s)) {	/* check if format string */
		if (!strcmp(fmt,COLRFMT)) lumf = rgblum;
		else if (!strcmp(fmt,CIEFMT)) lumf = cielum;
		else lumf = NULL;
		return(0);		/* don't echo */
	}
	if (isprims(s)) {		/* get input primaries */
		primsval(inprimS, s);
		inprims= inprimS;
		return(0);		/* don't echo */
	}
	if (isexpos(s)) {		/* picture exposure */
		inpexp *= exposval(s);
		return(0);		/* don't echo */
	}
	if (isaspect(s))		/* pixel aspect ratio */
		pixaspect *= aspectval(s);
	if (isview(s))			/* image view */
		gotview += sscanview(&ourview, s);
	return(fputs(s, stdout));
}


static void
getahead(void)			/* load picture header */
{
	char	*err;

	getheader(infp, headline, NULL);
	if (lumf == NULL || !fgetsresolu(&inpres, infp)) {
		fprintf(stderr, "%s: %s: not a Radiance picture\n",
			progname, infn);
		exit(1);
	}
	if (lumf == rgblum)
		comprgb2xyzWBmat(inrgb2xyz, inprims);
	else if (mbcalfile != NULL) {
		fprintf(stderr, "%s: macbethcal only works with RGB pictures\n",
				progname);
		exit(1);
	}
	if (!gotview || ourview.type == VT_PAR ||
			(ourview.horiz <= 5.) | (ourview.vert <= 5.)) {
		ourview = stdview;
		ourview.type = VT_PER;
		if (pixaspect*inpres.yr < inpres.xr) {
			ourview.horiz = 40.0;
			ourview.vert = 2.*180./PI *
				atan(.3639702*pixaspect*inpres.yr/inpres.xr);
		} else {
			ourview.vert = 40.0;
			ourview.horiz = 2.*180./PI *
				atan(.3639702*inpres.xr/pixaspect/inpres.yr);
		}
	}
	if ((err = setview(&ourview)) != NULL) {
		fprintf(stderr, "%s: view error in picture \"%s\": %s\n",
			progname, infn, err);
		exit(1);
	}
}


static void
mapimage(void)				/* map picture and send to stdout */
{
	COLOR	*scan;

	comphist();			/* generate adaptation histogram */
	check2do();			/* modify what2do flags */
	if (what2do&DO_VEIL)
		compveil();		/* compute veil image */
	if (!(what2do&DO_LINEAR))
		if (mkbrmap() < 0)	/* make tone map */
			what2do |= DO_LINEAR;	/* failed! -- use linear */
#if ADJ_VEIL
		else if (what2do&DO_VEIL)
			adjveil();	/* else adjust veil image */
#endif
	if (what2do&DO_LINEAR) {
		if (scalef <= FTINY) {
			if (what2do&DO_HSENS)
				scalef = htcontrs(Lb(0.5*(Bldmax+Bldmin))) /
						htcontrs(Lb(bwavg));
			else
				scalef = Lb(0.5*(Bldmax+Bldmin)) / Lb(bwavg);
			scalef *= WHTEFFICACY/(inpexp*ldmax);
		}
		fputexpos(inpexp*scalef, stdout);	/* record exposure */
		if (lumf == cielum) scalef /= WHTEFFICACY;
	}
	fputformat(COLRFMT, stdout);	/* complete header */
	putchar('\n');
	fputsresolu(&inpres, stdout);	/* resolution doesn't change */
					/* condition our image */
	for (scan = firstscan(); scan != NULL; scan = nextscan())
		if (fwritescan(scan, scanlen(&inpres), stdout) < 0) {
			fprintf(stderr, "%s: scanline write error\n",
					progname);
			exit(1);
		}
}


static void
getfovimg(void)			/* load foveal sampled image */
{
	char	combuf[PATH_MAX];
	FILE	*fp;
	int	x, y;
						/* compute image size */
	fvxr = sqrt(ourview.hn2)/FOVDIA + 0.5;
	if (fvxr < 2) fvxr = 2;
	fvyr = sqrt(ourview.vn2)/FOVDIA + 0.5;
	if (fvyr < 2) fvyr = 2;
	if (!(inpres.rt & YMAJOR)) {		/* picture is rotated? */
		y = fvyr;
		fvyr = fvxr;
		fvxr = y;
	}
	if ((fovimg = (COLOR *)malloc(fvxr*fvyr*sizeof(COLOR))) == NULL)
		syserror("malloc");
	sprintf(combuf, "pfilt -1 -b -pa 0 -x %d -y %d \"%s\"", fvxr, fvyr, infn);
	if ((fp = popen(combuf, "r")) == NULL)
		syserror("popen");
	getheader(fp, NULL, NULL);	/* skip header */
	if (fgetresolu(&x, &y, fp) < 0 || (x != fvxr) | (y != fvyr))
		goto readerr;
	for (y = 0; y < fvyr; y++)
		if (freadscan(fovscan(y), fvxr, fp) < 0)
			goto readerr;
	pclose(fp);
	return;
readerr:
	fprintf(stderr, "%s: error reading from pfilt process in fovimage\n",
			progname);
	exit(1);
}


static void
check2do(void)		/* check histogram to see what isn't worth doing */
{
	double	sum;
	double	b, l;
	register int	i;

					/* check for within display range */
	if (bwmax - bwmin <= Bldmax - Bldmin)
		what2do |= DO_LINEAR;
					/* determine if veiling significant */
	if (bwmax - bwmin < 4.5)		/* heuristic */
		what2do &= ~DO_VEIL;

	if (!(what2do & (DO_ACUITY|DO_COLOR)))
		return;
					/* find 5th percentile */
	sum = histot*0.05;
	for (i = 0; i < HISTRES; i++)
		if ((sum -= bwhist[i]) <= 0)
			break;
	b = (i+.5)*(bwmax-bwmin)/HISTRES + bwmin;
	l = Lb(b);
					/* determine if acuity adj. useful */
	if (what2do&DO_ACUITY &&
			hacuity(l) >= (inpres.xr/sqrt(ourview.hn2) +
			inpres.yr/sqrt(ourview.vn2))/(2.*180./PI))
		what2do &= ~DO_ACUITY;
					/* color sensitivity loss? */
	if (l >= TopMesopic)
		what2do &= ~DO_COLOR;
}
