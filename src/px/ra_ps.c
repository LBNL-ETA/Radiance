#ifndef lint
static const char	RCSid[] = "$Id: ra_ps.c,v 2.28 2004/03/28 20:33:14 schorsch Exp $";
#endif
/*
 *  Radiance picture to PostScript file translator -- one way!
 */

#include  <stdio.h>
#include  <string.h>
#include  <math.h>
#include  <ctype.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

#define UPPER(c)	((c)&~0x20)		/* ASCII trick */

#define CODE6GAM	1.47			/* gamma for 6-bit codes */
#define DEFGGAM		1.0			/* greyscale device gamma */
#define DEFCGAM		1.8			/* color device gamma */

#define GRY		-1			/* artificial index for grey */

#define DEFMARG		0.5			/* default margin (inches) */
#define DEFWIDTH	8.5			/* default page width */
#define DEFHEIGHT	11			/* default page height */
#define HMARGIN		(hmarg*72)		/* horizontal margin */
#define VMARGIN		(vmarg*72)		/* vertical margin */
#define PWIDTH		(width*72-2*HMARGIN)	/* width of device */
#define PHEIGHT		(height*72-2*VMARGIN)	/* height of device */

#define RUNCHR		'*'			/* character to start rle */
#define MINRUN		4			/* minimum run-length */
#define RSTRT		'!'			/* character for MINRUN */
#define MAXRUN		(MINRUN+'~'-RSTRT)	/* maximum run-length */

char  code[] =			/* 6-bit code lookup table */
	"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ@+";

int  wrongformat = 0;			/* input in wrong format? */
double	pixaspect = 1.0;		/* pixel aspect ratio */

double  devgam = 0.;			/* device gamma response */
double  hmarg = DEFMARG,
	vmarg = DEFMARG;		/* horizontal and vertical margins */
double	width = DEFWIDTH,
	height = DEFHEIGHT;		/* default paper width and height */
double  dpi = 0;			/* print density (0 if unknown) */
int  docolor = 1;			/* produce color image? */
int  bradj = 0;				/* brightness adjustment */
int  ncopies = 1;			/* number of copies */

char  *progname;
int  xmax, ymax;			/* input image dimensions */

typedef void putprimf_t(COLR *scn, int pri);

static gethfunc headline;
static putprimf_t Aputprim, Bputprim, Cputprim;

static double unit2inch(register char *s);
static int matchid(char *name, char *id);
static void parsepaper(char *ps);
static void quiterr(char *err);
static void PSheader(int ac, char **av);
static void PStrailer(void);
static void PSprocdef(char *nam);
static void ra2ps(void);
static void putrle(int cnt, int cod);


putprimf_t *putprim = Aputprim;		/* function for writing scanline */


static int
headline(		/* check header line */
	char	*s,
	void	*p
)
{
	char  fmt[32];

	if (isformat(s)) {
		formatval(fmt, s);
		wrongformat = strcmp(fmt, COLRFMT);
	} else if (isaspect(s))
		pixaspect *= aspectval(s);
	return(0);
}

int
main(int  argc, char  *argv[])
{
	int  i;
	double	d;
	
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'b':		/* produce b&w PostScript */
				docolor = 0;
				break;
			case 'c':		/* produce color PostScript */
				docolor = 1;
				break;
			case 'A':		/* standard ASCII encoding */
				putprim = Aputprim;
				break;
			case 'B':		/* standard binary encoding */
				putprim = Bputprim;
				break;
			case 'C':		/* compressed ASCII encoding */
				putprim = Cputprim;
				break;
			case 'd':		/* print density */
				dpi = atof(argv[++i]);
				break;
			case 'g':		/* device gamma adjustment */
				devgam = atof(argv[++i]);
				break;
			case 'p':		/* paper size */
				parsepaper(argv[++i]);
				break;
			case 'm':		/* margin */
				d = atof(argv[i+1]);
				d *= unit2inch(argv[i+1]);
				switch (argv[i][2]) {
				case '\0':
					hmarg = vmarg = d;
					break;
				case 'h':
					hmarg = d;
					break;
				case 'v':
					vmarg = d;
					break;
				default:
					goto userr;
				}
				i++;
				break;
			case 'e':		/* exposure adjustment */
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'n':		/* number of copies */
				ncopies = atoi(argv[++i]);
				break;
			default:
				goto userr;
			}
		else
			break;

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
	SET_FILE_BINARY(stdin);
				/* get our header */
	getheader(stdin, headline, NULL);
	if (wrongformat || fgetresolu(&xmax, &ymax, stdin) < 0)
		quiterr("bad picture format");
				/* gamma compression */
	if (devgam <= 0.05)
		devgam = docolor ? DEFCGAM : DEFGGAM;
	if (putprim == Cputprim)
		setcolrgam(CODE6GAM);
	else if (devgam != 1.)
		setcolrgam(devgam);
				/* write header */
	PSheader(argc, argv);
				/* convert file */
	ra2ps();
				/* write trailer */
	PStrailer();
	exit(0);
userr:
	fprintf(stderr,
"Usage: %s [-b|c][-A|B|C][-e +/-stops][-p paper][-m[h|v] margin][-d dpi][-g gamma] [input [output]]\n",
			progname);
	exit(1);
}


static double
unit2inch(		/* determine unit */
	register char	*s
)
{
	static struct unit {char n; float f;} u[] = {
		{'i', 1.},
		{'m', 1./25.4},
		{'c', 1./2.54},
		{'\0',0} };
	register struct unit	*up;

	while (*s && !isalpha(*s))
		s++;
	for (up = u; up->n; up++)
		if (up->n == *s)
			return(up->f);
	return(1.);
}


static int
matchid(	/* see if name matches id (case insensitive) */
	char	*name,
	register char	*id
)
{
	register char	*s = name;

	while (*s) {
		if (isalpha(*s)) {
			if (!isalpha(*id) || UPPER(*s) != UPPER(*id))
				return(0);
		} else if (*s != *id)
			return(0);
		s++; id++;
	}
	return(!*id || s-name >= 3);	/* substrings >= 3 chars OK */
}


static void
parsepaper(		/* determine paper size from name */
	char	*ps
)
{
	static struct psize {char n[12]; float w,h;} p[] = {
		{"envelope", 4.12, 9.5},
		{"executive", 7.25, 10.5},
		{"letter", 8.5, 11.},
		{"lettersmall", 7.68, 10.16},
		{"legal", 8.5, 14.},
		{"monarch", 3.87, 7.5},
		{"statement", 5.5, 8.5},
		{"tabloid", 11., 17.},
		{"A3", 11.69, 16.54},
		{"A4", 8.27, 11.69},
		{"A4small", 7.47, 10.85},
		{"A5", 6.00, 8.27},
		{"A6", 4.13, 6.00},
		{"B4", 10.12, 14.33},
		{"B5", 7.17, 10.12},
		{"C5", 6.38, 9.01},
		{"C6", 4.49, 6.38},
		{"DL", 4.33, 8.66},
		{"hagaki", 3.94, 5.83},
		{"",0.0,0.0} };
	register struct psize	*pp;
	register char	*s = ps;
	double	d;

	if (isdigit(*s)) {		/* check for WWxHH specification */
		width = atof(s);
		while (*s && !isalpha(*s))
			s++;
		d = unit2inch(s);
		height = atof(++s);
		width *= d;
		height *= d;
		if ((width >= 1.) & (height >= 1.))
			return;
	} else				/* check for match to standard size */
		for (pp = p; pp->n[0]; pp++)
			if (matchid(s, pp->n)) {
				width = pp->w;
				height = pp->h;
				return;
			}
	fprintf(stderr, "%s: unknown paper size \"%s\" -- known sizes:\n",
			progname, ps);
	fprintf(stderr, "_Name________Width_Height_(inches)\n");
	for (pp = p; pp->n[0]; pp++)
		fprintf(stderr, "%-11s  %5.2f  %5.2f\n", pp->n, pp->w, pp->h);
	fprintf(stderr, "Or use WWxHH size specification\n");
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
PSheader(		/* print PostScript header */
	int  ac,
	char  **av
)
{
	char  *rstr;
	int  landscape, rotate, n;
	double	pwidth, pheight;
	double	iwidth, iheight;
					/* EPS comments */
	puts("%!PS-Adobe-2.0 EPSF-2.0");
	printf("%%%%Title: "); printargs(ac, av, stdout);
	printf("%%%%Creator: %s\n", progname);
	printf("%%%%Pages: %d\n", ncopies);
	if ( (landscape = xmax > pixaspect*ymax) )
		puts("%%Orientation: Landscape");
	else
		puts("%%Orientation: Portrait");
	if ( (rotate = (PWIDTH > PHEIGHT) ^ landscape) ) {
		pwidth = PHEIGHT;
		pheight = PWIDTH;
	} else {
		pwidth = PWIDTH;
		pheight = PHEIGHT;
	}
	if (dpi > 100 && (pixaspect >= 0.99) & (pixaspect <= 1.01))
		if (pheight/pwidth > ymax/xmax) {
			n = pwidth*dpi/xmax;	/* floor */
			iwidth = n > 0 ? (double)(n*xmax)/dpi : pwidth;
			iheight = iwidth*ymax/xmax;
		} else {
			n = pheight*dpi/ymax;	/* floor */
			iheight = n > 0 ? (double)(n*ymax)/dpi : pheight;
			iwidth = iheight*xmax/ymax;
		}
	else
		if (pheight/pwidth > pixaspect*ymax/xmax) {
			iwidth = pwidth;
			iheight = iwidth*pixaspect*ymax/xmax;
		} else {
			iheight = pheight;
			iwidth = iheight*xmax/(pixaspect*ymax);
		}
	if (rotate)
		printf("%%%%BoundingBox: %.0f %.0f %.0f %.0f\n",
				HMARGIN+(pheight-iheight)*.5,
				VMARGIN+(pwidth-iwidth)*.5,
				HMARGIN+(pheight-iheight)*.5+iheight,
				VMARGIN+(pwidth-iwidth)*.5+iwidth);
	else
		printf("%%%%BoundingBox: %.0f %.0f %.0f %.0f\n",
				HMARGIN+(pwidth-iwidth)*.5,
				VMARGIN+(pheight-iheight)*.5,
				HMARGIN+(pwidth-iwidth)*.5+iwidth,
				VMARGIN+(pheight-iheight)*.5+iheight);
	puts("%%EndComments");
	puts("gsave save");
	puts("17 dict begin");
					/* define image reader */
	if (docolor) {
		printf("/redline %d string def\n", xmax);
		printf("/grnline %d string def\n", xmax);
		printf("/bluline %d string def\n", xmax);
	} else
		printf("/gryline %d string def\n", xmax);
					/* use compressed encoding? */
	if (putprim == Cputprim)
		PSprocdef("read6bitRLE");
					/* set up transformation matrix */
	printf("%f %f translate\n", HMARGIN, VMARGIN);
	if (rotate) {
		printf("%f 0 translate\n", PWIDTH);
		puts("90 rotate");
	}
	printf("%f %f translate\n", (pwidth-iwidth)*.5, (pheight-iheight)*.5);
	printf("%f %f scale\n", iwidth, iheight);
	puts("%%EndProlog");
					/* start image procedure */
	printf("%d %d 8 [%d 0 0 %d 0 %d]\n", xmax, ymax, xmax, -ymax, ymax);
	if (putprim == Cputprim) {
		if (docolor) {
			puts("{redline read6bitRLE}");
			puts("{grnline read6bitRLE}");
			puts("{bluline read6bitRLE}");
			puts("true 3 colorimage");
		} else
			puts("{gryline read6bitRLE} image");
	} else {
		rstr = putprim==Aputprim ? "readhexstring" : "readstring";
		if (docolor) {
			printf("{currentfile redline %s pop}\n", rstr);
			printf("{currentfile grnline %s pop}\n", rstr);
			printf("{currentfile bluline %s pop}\n", rstr);
			puts("true 3 colorimage");
		} else
			printf("{currentfile gryline %s pop} image\n", rstr);
	}
}


static void
PStrailer(void)			/* print PostScript trailer */
{
	puts("%%Trailer");
	if (ncopies > 1)
		printf("/#copies %d def\n", ncopies);
	puts("showpage");
	puts("end");
	puts("restore grestore");
	puts("%%EOF");
}


static void
PSprocdef(			/* define PS procedure to read image */
	char  *nam
)
{
	short  itab[128];
	register int  i;
				/* assign code values */
	for (i = 0; i < 128; i++)	/* clear */
		itab[i] = -1;
	for (i = 1; i < 63; i++)	/* assign greys */
		itab[(int)code[i]] = 256.0*pow((i+.5)/64.0, CODE6GAM/devgam);
	itab[(int)code[0]] = 0;		/* black is black */
	itab[(int)code[63]] = 255;		/* and white is white */
	printf("/codetab [");
	for (i = 0; i < 128; i++) {
		if (!(i & 0xf))
			putchar('\n');
		printf(" %3d", itab[i]);
	}
	printf("\n] def\n");
	printf("/nrept 0 def\n");
	printf("/readbyte { currentfile read not {stop} if } bind def\n");
	printf("/decode { codetab exch get } bind def\n");
	printf("/%s {\t%% scanbuffer\n", nam);
	printf("\t/scanline exch def\n");
	printf("\t{ 0 1 %d { scanline exch\n", xmax-1);
	printf("\t\tnrept 0 le\n");
	printf("\t\t\t{ { readbyte dup %d eq\n", RUNCHR);
	printf("\t\t\t\t\t{ pop /nrept readbyte %d sub def\n", RSTRT-MINRUN+1);
	printf("\t\t\t\t\t\t/reptv readbyte decode def\n");
	printf("\t\t\t\t\t\treptv exit }\n");
	printf("\t\t\t\t\t{ decode dup 0 lt {pop} {exit} ifelse }\n");
	printf("\t\t\t\tifelse } loop }\n");
	printf("\t\t\t{ /nrept nrept 1 sub def reptv }\n");
	printf("\t\tifelse put\n");
	printf("\t\t} for\n");
	printf("\t} stopped {pop pop 0 string} {scanline} ifelse\n");
	printf("} bind def\n");
}


static void
ra2ps(void)				/* convert Radiance scanlines to 6-bit */
{
	register COLR	*scanin;
	int	y;
						/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in ra2ps");
						/* convert image */
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(scanin, xmax, stdin) < 0)
			quiterr("error reading Radiance picture");
		if (putprim == Cputprim || devgam != 1.) {
			if (bradj)			/* adjust exposure */
				shiftcolrs(scanin, xmax, bradj);
			colrs_gambs(scanin, xmax);	/* gamma compression */
		} else
			normcolrs(scanin, xmax, bradj);
		if (docolor) {
			(*putprim)(scanin, RED);
			(*putprim)(scanin, GRN);
			(*putprim)(scanin, BLU);
		} else
			(*putprim)(scanin, GRY);
		if (ferror(stdout))
			quiterr("error writing PostScript file");
	}
	putchar('\n');
						/* free scanline */
	free((void *)scanin);
}


static void
Aputprim(		/* put out hex ASCII primary from scanline */
	COLR	*scn,
	int	pri
)
{
	static char	hexdigit[] = "0123456789ABCDEF";
	static int	col = 0;
	register int	x, c;

	for (x = 0; x < xmax; x++) {
		if (pri == GRY)
			c = normbright(scn[x]);
		else
			c = scn[x][pri];
		if (c > 255) c = 255;
		putchar(hexdigit[c>>4]);
		putchar(hexdigit[c&0xf]);
		if ((col += 2) >= 72) {
			putchar('\n');
			col = 0;
		}
	}
}


static void
Bputprim(		/* put out binary primary from scanline */
	COLR	*scn,
	int	pri
)
{
	register int	x, c;

	for (x = 0; x < xmax; x++) {
		if (pri == GRY)
			c = normbright(scn[x]);
		else
			c = scn[x][pri];
		if (c > 255) c = 255;
		putchar(c);
	}
}


static void
Cputprim(		/* put out compressed primary from scanline */
	COLR	*scn,
	int	pri
)
{
	register int	c;
	register int	x;
	int	lastc, cnt;

	lastc = -1; cnt = 0;
	for (x = 0; x < xmax; x++) {
		if (pri == GRY)
			c = normbright(scn[x]) + 2;
		else
			c = scn[x][pri] + 2;
		if (c > 255) c = 255;
		c = code[c>>2];
		if (c == lastc && cnt < MAXRUN)
			cnt++;
		else {
			putrle(cnt, lastc);
			lastc = c;
			cnt = 1;
		}
	}
	putrle(cnt, lastc);
}


static void
putrle(		/* put out cnt of cod */
	register int	cnt,
	register int	cod
)
{
	static int	col = 0;

	if (cnt >= MINRUN) {
		col += 3;
		putchar(RUNCHR);
		putchar(RSTRT-MINRUN+cnt);
		putchar(cod);
	} else {
		col += cnt;
		while (cnt-- > 0)
			putchar(cod);
	}
	if (col >= 72) {
		putchar('\n');
		col = 0;
	}
}
