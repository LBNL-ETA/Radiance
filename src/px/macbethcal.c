#ifndef lint
static const char	RCSid[] = "$Id: macbethcal.c,v 2.25 2008/11/10 19:08:19 greg Exp $";
#endif
/*
 * Calibrate a scanned MacBeth Color Checker Chart
 *
 * Produce a .cal file suitable for use with pcomb,
 * or .cwp file suitable for use with pcwarp.
 *
 * Warping code depends on conformance of COLOR and W3VEC types.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "platform.h"
#include "rtprocess.h"
#include "rtio.h"
#include "color.h"
#include "resolu.h"
#include "pmap.h"
#include "warp3d.h"
#include "mx3.h"

				/* MacBeth colors */
#define	DarkSkin	0
#define	LightSkin	1
#define	BlueSky		2
#define	Foliage		3
#define	BlueFlower	4
#define	BluishGreen	5
#define	Orange		6
#define	PurplishBlue	7
#define	ModerateRed	8
#define	Purple		9
#define	YellowGreen	10
#define	OrangeYellow	11
#define	Blue		12
#define	Green		13
#define	Red		14
#define	Yellow		15
#define	Magenta		16
#define	Cyan		17
#define	White		18
#define	Neutral8	19
#define	Neutral65	20
#define	Neutral5	21
#define	Neutral35	22
#define	Black		23
				/* computed from 5nm spectral measurements */
				/* CIE 1931 2 degree obs, equal-energy white */
float	mbxyY[24][3] = {
		{0.462, 0.3769, 0.0932961},	/* DarkSkin */
		{0.4108, 0.3542, 0.410348},	/* LightSkin */
		{0.2626, 0.267, 0.181554},	/* BlueSky */
		{0.36, 0.4689, 0.108447},	/* Foliage */
		{0.2977, 0.2602, 0.248407},	/* BlueFlower */
		{0.2719, 0.3485, 0.401156},	/* BluishGreen */
		{0.52, 0.4197, 0.357899},	/* Orange */
		{0.229, 0.1866, 0.103911},	/* PurplishBlue */
		{0.4909, 0.3262, 0.242615},	/* ModerateRed */
		{0.3361, 0.2249, 0.0600102},	/* Purple */
		{0.3855, 0.4874, 0.42963},	/* YellowGreen */
		{0.4853, 0.4457, 0.476343},	/* OrangeYellow */
		{0.2026, 0.1369, 0.0529249},	/* Blue */
		{0.3007, 0.4822, 0.221226},	/* Green */
		{0.5805, 0.3238, 0.162167},	/* Red */
		{0.4617, 0.472, 0.64909},	/* Yellow */
		{0.4178, 0.2625, 0.233662},	/* Magenta */
		{0.2038, 0.2508, 0.167275},	/* Cyan */
		{0.3358, 0.337, 0.916877},	/* White */
		{0.3338, 0.3348, 0.604678},	/* Neutral.8 */
		{0.3333, 0.3349, 0.364566},	/* Neutral.65 */
		{0.3353, 0.3359, 0.200238},	/* Neutral.5 */
		{0.3363, 0.336, 0.0878721},	/* Neutral.35 */
		{0.3346, 0.3349, 0.0308383}	/* Black */
	};

COLOR	mbRGB[24];		/* MacBeth RGB values */

#define	NMBNEU		6	/* Number of MacBeth neutral colors */
short	mbneu[NMBNEU] = {Black,Neutral35,Neutral5,Neutral65,Neutral8,White};

#define  NEUFLGS	(1L<<White|1L<<Neutral8|1L<<Neutral65| \
				1L<<Neutral5|1L<<Neutral35|1L<<Black)

#define  SATFLGS	(1L<<Red|1L<<Green|1L<<Blue|1L<<Magenta|1L<<Yellow| \
			1L<<Cyan|1L<<Orange|1L<<Purple|1L<<PurplishBlue| \
			1L<<YellowGreen|1<<OrangeYellow|1L<<BlueFlower)

#define  UNSFLGS	(1L<<DarkSkin|1L<<LightSkin|1L<<BlueSky|1L<<Foliage| \
			1L<<BluishGreen|1L<<ModerateRed)

#define  REQFLGS	NEUFLGS			/* need these colors */
#define  MODFLGS	(NEUFLGS|UNSFLGS)	/* should be in gamut */

#define  RG_BORD	0	/* patch border */
#define  RG_CENT	01	/* central region of patch */
#define  RG_ORIG	02	/* original color region */
#define  RG_CORR	04	/* corrected color region */

#ifndef  DISPCOM
#define  DISPCOM	"ximage -op \"%s\""
#endif

int	scanning = 1;		/* scanned input (or recorded output)? */
double	irrad = 1.0;		/* irradiance multiplication factor */
int	rawmap = 0;		/* put out raw color mapping? */

int	xmax, ymax;		/* input image dimensions */
int	bounds[4][2];		/* image coordinates of chart corners */
double	imgxfm[3][3];		/* coordinate transformation matrix */

COLOR	inpRGB[24];		/* measured or scanned input colors */
long	inpflags = 0;		/* flags of which colors were input */
long	gmtflags = 0;		/* flags of out-of-gamut colors */

COLOR	bramp[NMBNEU][2];	/* brightness ramp (per primary) */
COLORMAT	solmat;		/* color mapping matrix */
COLOR	colmin, colmax;		/* gamut limits */

WARP3D	*wcor = NULL;		/* color space warp */

FILE	*debugfp = NULL;	/* debug output picture */
char	*progname;

static void init(void);
static int chartndx(int	x, int y, int	*np);
static void getpicture(void);
static void getcolors(void);
static void bresp(COLOR	y, COLOR	x);
static void compute(void);
static void putmapping(void);
static void compsoln(COLOR	cin[], COLOR	cout[], int	n);
static void cwarp(void);
static int cvtcolor(COLOR	cout, COLOR	cin);
static int cresp(COLOR	cout, COLOR	cin);
static void xyY2RGB(COLOR	rgbout, float	xyYin[3]);
static void picdebug(void);
static void clrdebug(void);
static void getpos(char	*name, int	bnds[2], FILE	*fp);
static void pickchartpos(char	*pfn);


int
main(
	int	argc,
	char	**argv
)
{
	int	i;

	progname = argv[0];
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'd':				/* debug output */
			i++;
			if (badarg(argc-i, argv+i, "s"))
				goto userr;
			if ((debugfp = fopen(argv[i], "w")) == NULL) {
				perror(argv[i]);
				exit(1);
			}
			SET_FILE_BINARY(debugfp);
			newheader("RADIANCE", debugfp);		/* start */
			printargs(argc, argv, debugfp);		/* header */
			break;
		case 'p':				/* picture position */
			if (badarg(argc-i-1, argv+i+1, "iiiiiiii"))
				goto userr;
			bounds[0][0] = atoi(argv[++i]);
			bounds[0][1] = atoi(argv[++i]);
			bounds[1][0] = atoi(argv[++i]);
			bounds[1][1] = atoi(argv[++i]);
			bounds[2][0] = atoi(argv[++i]);
			bounds[2][1] = atoi(argv[++i]);
			bounds[3][0] = atoi(argv[++i]);
			bounds[3][1] = atoi(argv[++i]);
			scanning = 2;
			break;
		case 'P':				/* pick position */
			scanning = 3;
			break;
		case 'i':				/* irradiance factor */
			i++;
			if (badarg(argc-i, argv+i, "f"))
				goto userr;
			irrad = atof(argv[i]);
			break;
		case 'm':				/* raw map output */
			rawmap = 1;
			break;
		case 'c':				/* color input */
			scanning = 0;
			break;
		default:
			goto userr;
		}
							/* open files */
	if (i < argc && freopen(argv[i], "r", stdin) == NULL) {
		perror(argv[i]);
		exit(1);
	}
	if (i+1 < argc && freopen(argv[i+1], "w", stdout) == NULL) {
		perror(argv[i+1]);
		exit(1);
	}
	if (scanning) {			/* load input picture header */
		SET_FILE_BINARY(stdin);
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				fgetresolu(&xmax, &ymax, stdin) < 0) {
			fprintf(stderr, "%s: bad input picture\n", progname);
			exit(1);
		}
		if (scanning == 3) {
			if (i >= argc)
				goto userr;
			pickchartpos(argv[i]);
			scanning = 2;
		}
	} else {			/* else set default xmax and ymax */
		xmax = 512;
		ymax = 2*512/3;
	}
	if (scanning != 2) {		/* use default boundaries */
		bounds[0][0] = bounds[2][0] = .029*xmax + .5;
		bounds[0][1] = bounds[1][1] = .956*ymax + .5;
		bounds[1][0] = bounds[3][0] = .971*xmax + .5;
		bounds[2][1] = bounds[3][1] = .056*ymax + .5;
	}
	init();				/* initialize */
	if (scanning)			/* get picture colors */
		getpicture();
	else
		getcolors();
	compute();			/* compute color mapping */
	if (rawmap) {			/* print out raw correspondence */
		register int	j;

		printf("# Color correspondence produced by:\n#\t\t");
		printargs(argc, argv, stdout);
		printf("#\tUsage: pcwarp %s uncorrected.hdr > corrected.hdr\n",
				i+1 < argc ? argv[i+1] : "{this_file}");
		printf("#\t   Or: pcond [options] -m %s orig.hdr > output.hdr\n",
				i+1 < argc ? argv[i+1] : "{this_file}");
		for (j = 0; j < 24; j++)
			printf("%f %f %f    %f %f %f\n",
				colval(inpRGB[j],RED), colval(inpRGB[j],GRN),
				colval(inpRGB[j],BLU), colval(mbRGB[j],RED),
				colval(mbRGB[j],GRN), colval(mbRGB[j],BLU));
		if (scanning && debugfp != NULL)
			cwarp();		/* color warp for debugging */
	} else {			/* print color mapping */
						/* print header */
		printf("{\n\tColor correction file computed by:\n\t\t");
		printargs(argc, argv, stdout);
		printf("\n\tUsage: pcomb -f %s uncorrected.hdr > corrected.hdr\n",
				i+1 < argc ? argv[i+1] : "{this_file}");
		if (!scanning)
			printf("\t   Or: pcond [options] -f %s orig.hdr > output.hdr\n",
					i+1 < argc ? argv[i+1] : "{this_file}");
		printf("}\n");
		putmapping();			/* put out color mapping */
	}
	if (debugfp != NULL) {		/* put out debug picture */
		if (scanning)
			picdebug();
		else
			clrdebug();
	}
	exit(0);
userr:
	fprintf(stderr,
"Usage: %s [-d dbg.hdr][-P | -p xul yul xur yur xll yll xlr ylr][-i irrad][-m] input.hdr [output.{cal|cwp}]\n",
			progname);
	fprintf(stderr, "   or: %s [-d dbg.hdr][-i irrad][-m] -c [xyY.dat [output.{cal|cwp}]]\n",
			progname);
	exit(1);
	return 1; /* pro forma return */
}


static void
init(void)				/* initialize */
{
	double	quad[4][2];
	register int	i;
					/* make coordinate transformation */
	quad[0][0] = bounds[0][0];
	quad[0][1] = bounds[0][1];
	quad[1][0] = bounds[1][0];
	quad[1][1] = bounds[1][1];
	quad[2][0] = bounds[3][0];
	quad[2][1] = bounds[3][1];
	quad[3][0] = bounds[2][0];
	quad[3][1] = bounds[2][1];

	if (pmap_quad_rect(0., 0., 6., 4., quad, imgxfm) == PMAP_BAD) {
		fprintf(stderr, "%s: bad chart boundaries\n", progname);
		exit(1);
	}
					/* map MacBeth colors to RGB space */
	for (i = 0; i < 24; i++) {
		xyY2RGB(mbRGB[i], mbxyY[i]);
		scalecolor(mbRGB[i], irrad);
	}
}


static int
chartndx(			/* find color number for position */
	int	x,
	int y,
	int	*np
)
{
	double	ipos[3], cpos[3];
	int	ix, iy;
	double	fx, fy;

	ipos[0] = x;
	ipos[1] = y;
	ipos[2] = 1;
	mx3d_transform(ipos, imgxfm, cpos);
	cpos[0] /= cpos[2];
	cpos[1] /= cpos[2];
	if (cpos[0] < 0. || cpos[0] >= 6. || cpos[1] < 0. || cpos[1] >= 4.)
		return(RG_BORD);
	ix = cpos[0];
	iy = cpos[1];
	fx = cpos[0] - ix;
	fy = cpos[1] - iy;
	*np = iy*6 + ix;
	if (fx >= 0.35 && fx < 0.65 && fy >= 0.35 && fy < 0.65)
		return(RG_CENT);
	if (fx < 0.05 || fx >= 0.95 || fy < 0.05 || fy >= 0.95)
		return(RG_BORD);
	if (fx >= 0.5)			/* right side is corrected */
		return(RG_CORR);
	return(RG_ORIG);		/* left side is original */
}


static void
getpicture(void)				/* load in picture colors */
{
	COLR	*scanln;
	COLOR	pval;
	int	ccount[24];
	double	d;
	int	y, i;
	register int	x;

	scanln = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanln == NULL) {
		perror(progname);
		exit(1);
	}
	for (i = 0; i < 24; i++) {
		setcolor(inpRGB[i], 0., 0., 0.);
		ccount[i] = 0;
	}
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanln, xmax, stdin) < 0) {
			fprintf(stderr, "%s: error reading input picture\n",
					progname);
			exit(1);
		}
		for (x = 0; x < xmax; x++)
			if (chartndx(x, y, &i) == RG_CENT) {
				colr_color(pval, scanln[x]);
				addcolor(inpRGB[i], pval);
				ccount[i]++;
			}
	}
	for (i = 0; i < 24; i++) {		/* compute averages */
		if (ccount[i] == 0)
			continue;
		d = 1./ccount[i];
		scalecolor(inpRGB[i], d);
		inpflags |= 1L<<i;
	}
	free((void *)scanln);
}


static void
getcolors(void)			/* get xyY colors from standard input */
{
	int	gotwhite = 0;
	COLOR	whiteclr;
	int	n;
	float	xyYin[3];

	while (fgetval(stdin, 'i', &n) == 1) {		/* read colors */
		if ((n < 0) | (n > 24) ||
				fgetval(stdin, 'f', &xyYin[0]) != 1 ||
				fgetval(stdin, 'f', &xyYin[1]) != 1 ||
				fgetval(stdin, 'f', &xyYin[2]) != 1 ||
				(xyYin[0] < 0.) | (xyYin[1] < 0.) ||
				xyYin[0] + xyYin[1] > 1.) {
			fprintf(stderr, "%s: bad color input data\n",
					progname);
			exit(1);
		}
		if (n == 0) {				/* calibration white */
			xyY2RGB(whiteclr, xyYin);
			gotwhite++;
		} else {				/* standard color */
			n--;
			xyY2RGB(inpRGB[n], xyYin);
			inpflags |= 1L<<n;
		}
	}
					/* normalize colors */
	if (!gotwhite) {
		if (!(inpflags & 1L<<White)) {
			fprintf(stderr, "%s: missing input for White\n",
					progname);
			exit(1);
		}
		setcolor(whiteclr,
			colval(inpRGB[White],RED)/colval(mbRGB[White],RED),
			colval(inpRGB[White],GRN)/colval(mbRGB[White],GRN),
			colval(inpRGB[White],BLU)/colval(mbRGB[White],BLU));
	}
	for (n = 0; n < 24; n++)
		if (inpflags & 1L<<n)
			setcolor(inpRGB[n],
				colval(inpRGB[n],RED)/colval(whiteclr,RED),
				colval(inpRGB[n],GRN)/colval(whiteclr,GRN),
				colval(inpRGB[n],BLU)/colval(whiteclr,BLU));
}


static void
bresp(		/* piecewise linear interpolation of primaries */
	COLOR	y,
	COLOR	x
)
{
	register int	i, n;

	for (i = 0; i < 3; i++) {
		for (n = 0; n < NMBNEU-2; n++)
			if (colval(x,i) < colval(bramp[n+1][0],i))
				break;
		colval(y,i) = ((colval(bramp[n+1][0],i) - colval(x,i)) *
						colval(bramp[n][1],i) +
				(colval(x,i) - colval(bramp[n][0],i)) *
						colval(bramp[n+1][1],i)) /
			(colval(bramp[n+1][0],i) - colval(bramp[n][0],i));
	}
}


static void
compute(void)			/* compute color mapping */
{
	COLOR	clrin[24], clrout[24];
	long	cflags;
	COLOR	ctmp;
	register int	i, n;
					/* did we get what we need? */
	if ((inpflags & REQFLGS) != REQFLGS) {
		fprintf(stderr, "%s: missing required input colors\n",
				progname);
		exit(1);
	}
					/* compute piecewise luminance curve */
	for (i = 0; i < NMBNEU; i++) {
		copycolor(bramp[i][0], inpRGB[mbneu[i]]);
		for (n = i ? 3 : 0; n--; )
			if (colval(bramp[i][0],n) <=
					colval(bramp[i-1][0],n)+1e-7) {
				fprintf(stderr,
		"%s: non-increasing neutral patch\n", progname);
				exit(1);
			}
		copycolor(bramp[i][1], mbRGB[mbneu[i]]);
	}
					/* compute color space gamut */
	if (scanning) {
		copycolor(colmin, cblack);
		copycolor(colmax, cwhite);
		scalecolor(colmax, irrad);
	} else
		for (i = 0; i < 3; i++) {
			colval(colmin,i) = colval(bramp[0][0],i) -
				colval(bramp[0][1],i) *
				(colval(bramp[1][0],i)-colval(bramp[0][0],i)) /
				(colval(bramp[1][1],i)-colval(bramp[1][0],i));
			colval(colmax,i) = colval(bramp[NMBNEU-2][0],i) +
				(1.-colval(bramp[NMBNEU-2][1],i)) *
				(colval(bramp[NMBNEU-1][0],i) -
					colval(bramp[NMBNEU-2][0],i)) /
				(colval(bramp[NMBNEU-1][1],i) -
					colval(bramp[NMBNEU-2][1],i));
		}
					/* compute color mapping */
	do {
		cflags = inpflags & ~gmtflags;
		n = 0;				/* compute transform matrix */
		for (i = 0; i < 24; i++)
			if (cflags & 1L<<i) {
				bresp(clrin[n], inpRGB[i]);
				copycolor(clrout[n], mbRGB[i]);
				n++;
			}
		compsoln(clrin, clrout, n);
		if (irrad > 0.99 && irrad < 1.01)	/* check gamut */
			for (i = 0; i < 24; i++)
				if (cflags & 1L<<i && cvtcolor(ctmp, mbRGB[i]))
					gmtflags |= 1L<<i;
	} while (cflags & gmtflags);
	if (gmtflags & MODFLGS)
		fprintf(stderr,
		"%s: warning - some moderate colors are out of gamut\n",
				progname);
}


static void
putmapping(void)			/* put out color mapping */
{
	static char	cchar[3] = {'r', 'g', 'b'};
	register int	i, j;
					/* print brightness mapping */
	for (j = 0; j < 3; j++) {
		printf("%cxa(i) : select(i", cchar[j]);
		for (i = 0; i < NMBNEU; i++)
			printf(",%g", colval(bramp[i][0],j));
		printf(");\n");
		printf("%cya(i) : select(i", cchar[j]);
		for (i = 0; i < NMBNEU; i++)
			printf(",%g", colval(bramp[i][1],j));
		printf(");\n");
		printf("%cfi(n) = if(n-%g, %d, if(%cxa(n+1)-%c, n, %cfi(n+1)));\n",
				cchar[j], NMBNEU-1.5, NMBNEU-1, cchar[j],
				cchar[j], cchar[j]);
		printf("%cndx = %cfi(1);\n", cchar[j], cchar[j]);
		printf("%c%c = ((%cxa(%cndx+1)-%c)*%cya(%cndx) + ",
				cchar[j], scanning?'n':'o', cchar[j],
				cchar[j], cchar[j], cchar[j], cchar[j]);
		printf("(%c-%cxa(%cndx))*%cya(%cndx+1)) /\n",
				cchar[j], cchar[j], cchar[j],
				cchar[j], cchar[j]);
		printf("\t\t(%cxa(%cndx+1) - %cxa(%cndx)) ;\n",
				cchar[j], cchar[j], cchar[j], cchar[j]);
	}
					/* print color mapping */
	if (scanning) {
		printf("r = ri(1); g = gi(1); b = bi(1);\n");
		printf("ro = %g*rn + %g*gn + %g*bn ;\n",
				solmat[0][0], solmat[0][1], solmat[0][2]);
		printf("go = %g*rn + %g*gn + %g*bn ;\n",
				solmat[1][0], solmat[1][1], solmat[1][2]);
		printf("bo = %g*rn + %g*gn + %g*bn ;\n",
				solmat[2][0], solmat[2][1], solmat[2][2]);
	} else {
		printf("r1 = ri(1); g1 = gi(1); b1 = bi(1);\n");
		printf("r = %g*r1 + %g*g1 + %g*b1 ;\n",
				solmat[0][0], solmat[0][1], solmat[0][2]);
		printf("g = %g*r1 + %g*g1 + %g*b1 ;\n",
				solmat[1][0], solmat[1][1], solmat[1][2]);
		printf("b = %g*r1 + %g*g1 + %g*b1 ;\n",
				solmat[2][0], solmat[2][1], solmat[2][2]);
	}
}


static void
compsoln(		/* solve 3xN system using least-squares */
	COLOR	cin[],
	COLOR	cout[],
	int	n
)
{
	extern double	mx3d_adjoint(), fabs();
	double	mat[3][3], invmat[3][3];
	double	det;
	double	colv[3], rowv[3];
	register int	i, j, k;

	if (n < 3) {
		fprintf(stderr, "%s: too few colors to match!\n", progname);
		exit(1);
	}
	if (n == 3)
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				mat[i][j] = colval(cin[j],i);
	else {				/* compute A^t A */
		for (i = 0; i < 3; i++)
			for (j = i; j < 3; j++) {
				mat[i][j] = 0.;
				for (k = 0; k < n; k++)
					mat[i][j] += colval(cin[k],i) *
							colval(cin[k],j);
			}
		for (i = 1; i < 3; i++)		/* using symmetry */
			for (j = 0; j < i; j++)
				mat[i][j] = mat[j][i];
	}
	det = mx3d_adjoint(mat, invmat);
	if (fabs(det) < 1e-4) {
		fprintf(stderr, "%s: cannot compute color mapping\n",
				progname);
		solmat[0][0] = solmat[1][1] = solmat[2][2] = 1.;
		solmat[0][1] = solmat[0][2] = solmat[1][0] =
		solmat[1][2] = solmat[2][0] = solmat[2][1] = 0.;
		return;
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			invmat[i][j] /= det;
	for (i = 0; i < 3; i++) {
		if (n == 3)
			for (j = 0; j < 3; j++)
				colv[j] = colval(cout[j],i);
		else
			for (j = 0; j < 3; j++) {
				colv[j] = 0.;
				for (k = 0; k < n; k++)
					colv[j] += colval(cout[k],i) *
							colval(cin[k],j);
			}
		mx3d_transform(colv, invmat, rowv);
		for (j = 0; j < 3; j++)
			solmat[i][j] = rowv[j];
	}
}


static void
cwarp(void)				/* compute color warp map */
{
	register int	i;

	if ((wcor = new3dw(W3EXACT)) == NULL)
		goto memerr;
	for (i = 0; i < 24; i++)
		if (!add3dpt(wcor, inpRGB[i], mbRGB[i]))
			goto memerr;
	return;
memerr:
	perror(progname);
	exit(1);
}


static int
cvtcolor(		/* convert color according to our mapping */
	COLOR	cout,
	COLOR	cin
)
{
	COLOR	ctmp;
	int	clipped;

	if (wcor != NULL) {
		clipped = warp3d(cout, cin, wcor);
		clipped |= clipgamut(cout,bright(cout),CGAMUT,colmin,colmax);
	} else if (scanning) {
		bresp(ctmp, cin);
		clipped = cresp(cout, ctmp);
	} else {
		clipped = cresp(ctmp, cin);
		bresp(cout, ctmp);
	}
	return(clipped);
}


static int
cresp(		/* transform color according to matrix */
	COLOR	cout,
	COLOR	cin
)
{
	colortrans(cout, solmat, cin);
	return(clipgamut(cout, bright(cout), CGAMUT, colmin, colmax));
}


static void
xyY2RGB(		/* convert xyY to RGB */
	COLOR	rgbout,
	register float	xyYin[3]
)
{
	COLOR	ctmp;
	double	d;

	d = xyYin[2] / xyYin[1];
	ctmp[0] = xyYin[0] * d;
	ctmp[1] = xyYin[2];
	ctmp[2] = (1. - xyYin[0] - xyYin[1]) * d;
				/* allow negative values */
	colortrans(rgbout, xyz2rgbmat, ctmp);
}


static void
picdebug(void)			/* put out debugging picture */
{
	static COLOR	blkcol = BLKCOLOR;
	COLOR	*scan;
	int	y, i;
	register int	x, rg;

	if (fseek(stdin, 0L, 0) == EOF) {
		fprintf(stderr, "%s: cannot seek on input picture\n", progname);
		exit(1);
	}
	getheader(stdin, NULL, NULL);		/* skip input header */
	fgetresolu(&xmax, &ymax, stdin);
						/* allocate scanline */
	scan = (COLOR *)malloc(xmax*sizeof(COLOR));
	if (scan == NULL) {
		perror(progname);
		exit(1);
	}
						/* finish debug header */
	fputformat(COLRFMT, debugfp);
	putc('\n', debugfp);
	fprtresolu(xmax, ymax, debugfp);
						/* write debug picture */
	for (y = ymax-1; y >= 0; y--) {
		if (freadscan(scan, xmax, stdin) < 0) {
			fprintf(stderr, "%s: error rereading input picture\n",
					progname);
			exit(1);
		}
		for (x = 0; x < xmax; x++) {
			rg = chartndx(x, y, &i);
			if (rg == RG_CENT) {
				if (!(1L<<i & gmtflags) || (x+y)&07) {
					copycolor(scan[x], mbRGB[i]);
					clipgamut(scan[x], bright(scan[x]),
						CGAMUT, colmin, colmax);
				} else
					copycolor(scan[x], blkcol);
			} else if (rg == RG_CORR)
				cvtcolor(scan[x], scan[x]);
			else if (rg != RG_ORIG)
				copycolor(scan[x], blkcol);
		}
		if (fwritescan(scan, xmax, debugfp) < 0) {
			fprintf(stderr, "%s: error writing debugging picture\n",
					progname);
			exit(1);
		}
	}
						/* clean up */
	fclose(debugfp);
	free((void *)scan);
}


static void
clrdebug(void)			/* put out debug picture from color input */
{
	static COLR	blkclr = BLKCOLR;
	COLR	mbclr[24], cvclr[24], orclr[24];
	COLR	*scan;
	COLOR	ctmp, ct2;
	int	y, i;
	register int	x, rg;
						/* convert colors */
	for (i = 0; i < 24; i++) {
		copycolor(ctmp, mbRGB[i]);
		clipgamut(ctmp, bright(ctmp), CGAMUT, cblack, cwhite);
		setcolr(mbclr[i], colval(ctmp,RED),
				colval(ctmp,GRN), colval(ctmp,BLU));
		if (inpflags & 1L<<i) {
			copycolor(ctmp, inpRGB[i]);
			clipgamut(ctmp, bright(ctmp), CGAMUT, cblack, cwhite);
			setcolr(orclr[i], colval(ctmp,RED),
					colval(ctmp,GRN), colval(ctmp,BLU));
			if (rawmap)
				copycolr(cvclr[i], mbclr[i]);
			else {
				bresp(ctmp, inpRGB[i]);
				colortrans(ct2, solmat, ctmp);
				clipgamut(ct2, bright(ct2), CGAMUT,
						cblack, cwhite);
				setcolr(cvclr[i], colval(ct2,RED),
						colval(ct2,GRN),
						colval(ct2,BLU));
			}
		}
	}
						/* allocate scanline */
	scan = (COLR *)malloc(xmax*sizeof(COLR));
	if (scan == NULL) {
		perror(progname);
		exit(1);
	}
						/* finish debug header */
	fputformat(COLRFMT, debugfp);
	putc('\n', debugfp);
	fprtresolu(xmax, ymax, debugfp);
						/* write debug picture */
	for (y = ymax-1; y >= 0; y--) {
		for (x = 0; x < xmax; x++) {
			rg = chartndx(x, y, &i);
			if (rg == RG_CENT) {
				if (!(1L<<i & gmtflags) || (x+y)&07)
					copycolr(scan[x], mbclr[i]);
				else
					copycolr(scan[x], blkclr);
			} else if (rg == RG_BORD || !(1L<<i & inpflags))
				copycolr(scan[x], blkclr);
			else if (rg == RG_ORIG)
				copycolr(scan[x], orclr[i]);
			else /* rg == RG_CORR */
				copycolr(scan[x], cvclr[i]);
		}
		if (fwritecolrs(scan, xmax, debugfp) < 0) {
			fprintf(stderr, "%s: error writing debugging picture\n",
					progname);
			exit(1);
		}
	}
						/* clean up */
	fclose(debugfp);
	free((void *)scan);
}


static void
getpos(		/* get boundary position */
	char	*name,
	int	bnds[2],
	FILE	*fp
)
{
	char	buf[64];

	fprintf(stderr, "\tSelect corner: %s\n", name);
	if (fgets(buf, sizeof(buf), fp) == NULL ||
			sscanf(buf, "%d %d", &bnds[0], &bnds[1]) != 2) {
		fprintf(stderr, "%s: read error from display process\n",
				progname);
		exit(1);
	}
}


static void
pickchartpos(		/* display picture and pick chart location */
	char	*pfn
)
{
	char	combuf[PATH_MAX];
	FILE	*pfp;

	sprintf(combuf, DISPCOM, pfn);
	if ((pfp = popen(combuf, "r")) == NULL) {
		perror(combuf);
		exit(1);
	}
	fputs("Use middle mouse button to select chart corners:\n", stderr);
	getpos("upper left (dark skin)", bounds[0], pfp);
	getpos("upper right (bluish green)", bounds[1], pfp);
	getpos("lower left (white)", bounds[2], pfp);
	getpos("lower right (black)", bounds[3], pfp);
	fputs("Got it -- quit display program.\n", stderr);
	pclose(pfp);
}
