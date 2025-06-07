#ifndef lint
static const char	RCSid[] = "$Id: ra_tiff.c,v 2.40 2025/06/07 05:09:46 greg Exp $";
#endif
/*
 *  Program to convert between RADIANCE and TIFF files.
 *  Added LogLuv encodings 7/97 (GWL).
 *  Added white-balance adjustment 10/01 (GW).
 */

#include  <math.h>
#include  <ctype.h>
#include  "rtio.h"
#include  "platform.h"
#include  "tiffio.h"
#include  "color.h"
#include  "resolu.h"

#define  GAMCOR		2.2		/* default gamma */

				/* conversion flags */
#define C_CXFM		0x1		/* needs color transformation */
#define C_GAMUT		0x2		/* needs gamut mapping */
#define C_GAMMA		0x4		/* needs gamma correction */
#define C_GRY		0x8		/* TIFF is greyscale */
#define C_XYZE		0x10		/* Radiance is XYZE */
#define C_SPEC		0x20		/* Radiance is spectral */
#define C_RFLT		0x40		/* Radiance data is float */
#define C_TFLT		0x80		/* TIFF data is float */
#define C_TWRD		0x100		/* TIFF data is 16-bit */
#define C_PRIM		0x200		/* has assigned primaries */

typedef void colcvf_t(uint32 y);

struct CONVST {
	uint16	flags;		/* conversion flags (defined above) */
	char	capdate[20];	/* capture date/time */
	char	owner[256];	/* content owner */
	uint16	comp;		/* TIFF compression type */
	uint16	phot;		/* TIFF photometric type */
	uint16	pconf;		/* TIFF planar configuration */
	float	gamcor;		/* gamma correction value */
	short	bradj;		/* Radiance exposure adjustment (stops) */
	uint16	orient;		/* visual orientation (TIFF spec.) */
	double	stonits;	/* input conversion to nits */
	float	pixrat;		/* pixel aspect ratio */
	int	ncc;		/* # color components for spectral */
	float	wpt[4];		/* wavelength partitions */
	FILE	*rfp;		/* Radiance stream pointer */
	TIFF	*tif;		/* TIFF pointer */
	uint32	xmax, ymax;	/* image dimensions */
	COLORMAT	cmat;	/* color transformation matrix */
	RGBPRIMS	prims;	/* RGB primaries */
	union {
		COLR	*colrs;		/* 4-byte ???E pointer */
		COLOR	*colors;	/* float array pointer */
		void	*p;		/* generic pointer */
	} r;			/* Radiance scanline */
	union {
		uint8	*bp;		/* byte pointer */
		uint16	*wp;		/* word pointer */
		float	*fp;		/* float pointer */
		void	*p;		/* generic pointer */
	} t;			/* TIFF scanline */
	colcvf_t *tf;	/* translation procedure */
}	cvts = {	/* conversion structure */
	0, "", "", COMPRESSION_NONE, PHOTOMETRIC_RGB,
	PLANARCONFIG_CONTIG, GAMCOR, 0, 1, 1., 1., 3,
};

#define CHK(f)		(cvts.flags & (f))
#define SET(f)		(cvts.flags |= (f))
#define CLR(f)		(cvts.flags &= ~(f))
#define TGL(f)		(cvts.flags ^= (f))

static colcvf_t Luv2Color, L2Color, RGB2Colr, Gry2Colr;
static colcvf_t Color2Luv, Color2L, Colr2RGB, Colr2Gry;
static colcvf_t RRGGBB2Color, GGry2Color, Color2RRGGBB, Color2GGry;

static gethfunc headline;
static void quiterr(char *err);
static void allocbufs(void);
static void initfromtif(void);
static void tiff2ra(int  ac, char  *av[]);
static void initfromrad(void);
static void ra2tiff(int  ac, char  *av[]);
static void ReadRadScan(struct CONVST *cvp);


#define RfGfBf2Color	Luv2Color
#define Gryf2Color	L2Color
#define	Color2Gryf	Color2L
#define Color2RfGfBf	Color2Luv

short	ortab[8] = {		/* orientation conversion table */
	YMAJOR|YDECR,
	YMAJOR|YDECR|XDECR,
	YMAJOR|XDECR,
	YMAJOR,
	YDECR,
	XDECR|YDECR,
	XDECR,
	0
};

#define pixorder()	ortab[cvts.orient-1]

extern char	TMSTR[];	/* "CAPDATE=" from header.c */
char		OWNSTR[] = "OWNER=";


int
main(int  argc, char  *argv[])
{
	int  reverse = 0;
	int  i;
					/* set global progname */
	fixargv0(argv[0]);

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':		/* gamma correction */
				cvts.gamcor = atof(argv[++i]);
				break;
			case 'x':		/* XYZE Radiance output */
				TGL(C_XYZE);
				break;
			case 'z':		/* LZW compressed output */
				cvts.comp = COMPRESSION_LZW;
				break;
			case 'L':		/* LogLuv 32-bit output */
				cvts.comp = COMPRESSION_SGILOG;
				cvts.phot = PHOTOMETRIC_LOGLUV;
				break;
			case 'l':		/* LogLuv 24-bit output */
				cvts.comp = COMPRESSION_SGILOG24;
				cvts.phot = PHOTOMETRIC_LOGLUV;
				break;
			case 'w':		/* 16-bit/primary output? */
				TGL(C_TWRD);
				break;
			case 'f':		/* IEEE float output? */
				TGL(C_TFLT);
				break;
			case 'b':		/* greyscale output? */
				TGL(C_GRY);
				break;
			case 'e':		/* exposure adjustment */
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				cvts.bradj = atoi(argv[++i]);
				break;
			case 'r':		/* reverse conversion? */
				reverse = !reverse;
				break;
			case '\0':
				goto doneopts;
			default:
				goto userr;
			}
		else
			break;
doneopts:
	if (reverse) {

		if (i != argc-2 && i != argc-1)
			goto userr;

		tiff2ra(i, argv);

	} else {

		if (i != argc-2)
			goto userr;
						/* consistency checks */
		if (CHK(C_GRY)) {
			if (cvts.phot == PHOTOMETRIC_RGB)
				cvts.phot = PHOTOMETRIC_MINISBLACK;
			else {
				cvts.phot = PHOTOMETRIC_LOGL;
				cvts.comp = COMPRESSION_SGILOG;
			}
		}
		if (CHK(C_TWRD|C_TFLT) == (C_TWRD|C_TFLT))
			goto userr;

		ra2tiff(i, argv);
	}

	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s [-z|-L|-l|-f|-w][-b][-e +/-stops][-g gamma] {in.hdr|-} out.tif\n",
			progname);
	fprintf(stderr,
	"   Or: %s -r [-x][-e +/-stops][-g gamma] in.tif [out.hdr|-]\n",
			progname);
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
allocbufs(void)			/* allocate scanline buffers */
{
	int	rsiz, tsiz;

	rsiz = CHK(C_RFLT) ? sizeof(COLOR) : sizeof(COLR);
	tsiz = (CHK(C_TFLT) ? sizeof(float) : 
			CHK(C_TWRD) ? sizeof(uint16) : sizeof(uint8)) *
			(CHK(C_GRY) ? 1 : 3);
	cvts.r.p = malloc(rsiz*cvts.xmax);
	cvts.t.p = malloc(tsiz*cvts.xmax);
	if ((cvts.r.p == NULL) | (cvts.t.p == NULL))
		quiterr("no memory to allocate scanline buffers");
}


static void
initfromtif(void)		/* initialize conversion from TIFF input */
{
	uint16	hi;
	char	*cp;
	float	*fa, f1, f2;

	TIFFGetFieldDefaulted(cvts.tif, TIFFTAG_PLANARCONFIG, &cvts.pconf);

	if (TIFFGetField(cvts.tif, TIFFTAG_PRIMARYCHROMATICITIES, &fa)) {
		cvts.prims[RED][CIEX] = fa[0];
		cvts.prims[RED][CIEY] = fa[1];
		cvts.prims[GRN][CIEX] = fa[2];
		cvts.prims[GRN][CIEY] = fa[3];
		cvts.prims[BLU][CIEX] = fa[4];
		cvts.prims[BLU][CIEY] = fa[5];
		cvts.prims[WHT][CIEX] = 1./3.;
		cvts.prims[WHT][CIEY] = 1./3.;
		if (TIFFGetField(cvts.tif, TIFFTAG_WHITEPOINT, &fa)) {
			cvts.prims[WHT][CIEX] = fa[0];
			cvts.prims[WHT][CIEY] = fa[1];
		}
		SET(C_PRIM);
	}

	if (!TIFFGetField(cvts.tif, TIFFTAG_COMPRESSION, &cvts.comp))
		cvts.comp = COMPRESSION_NONE;

	if (TIFFGetField(cvts.tif, TIFFTAG_XRESOLUTION, &f1) &&
			TIFFGetField(cvts.tif, TIFFTAG_YRESOLUTION, &f2))
		cvts.pixrat = f1/f2;

	TIFFGetFieldDefaulted(cvts.tif, TIFFTAG_ORIENTATION, &cvts.orient);

	if (!TIFFGetFieldDefaulted(cvts.tif, TIFFTAG_PHOTOMETRIC, &cvts.phot))
		quiterr("TIFF has unspecified photometric type");

	switch (cvts.phot) {
	case PHOTOMETRIC_LOGLUV:
		SET(C_RFLT|C_TFLT);
		if (!CHK(C_XYZE)) {
			cpcolormat(cvts.cmat, xyz2rgbmat);
			SET(C_CXFM|C_GAMUT);
		} else if (cvts.comp == COMPRESSION_SGILOG)
			SET(C_GAMUT);		/* may be outside XYZ gamut */
		if (cvts.pconf != PLANARCONFIG_CONTIG)
			quiterr("cannot handle separate Luv planes");
		TIFFSetField(cvts.tif, TIFFTAG_SGILOGDATAFMT,
				SGILOGDATAFMT_FLOAT);
		cvts.tf = Luv2Color;
		break;
	case PHOTOMETRIC_LOGL:
		SET(C_GRY|C_RFLT|C_TFLT|C_GAMUT);
		cvts.pconf = PLANARCONFIG_CONTIG;
		TIFFSetField(cvts.tif, TIFFTAG_SGILOGDATAFMT,
				SGILOGDATAFMT_FLOAT);
		cvts.tf = L2Color;
		break;
	case PHOTOMETRIC_YCBCR:
		if (cvts.comp == COMPRESSION_JPEG &&
				cvts.pconf == PLANARCONFIG_CONTIG) {
			TIFFSetField(cvts.tif, TIFFTAG_JPEGCOLORMODE,
					JPEGCOLORMODE_RGB);
			cvts.phot = PHOTOMETRIC_RGB;
		} else
			quiterr("unsupported photometric type");
		/* fall through */
	case PHOTOMETRIC_RGB:
		SET(C_GAMMA);
		setcolrgam(cvts.gamcor);
		if (CHK(C_XYZE)) {
			comprgb2xyzWBmat(cvts.cmat,
					CHK(C_PRIM) ? cvts.prims : stdprims);
			SET(C_CXFM);
		}
		if (!TIFFGetField(cvts.tif, TIFFTAG_SAMPLESPERPIXEL, &hi) ||
				hi != 3)
			quiterr("unsupported samples per pixel for RGB");
		if (!TIFFGetField(cvts.tif, TIFFTAG_BITSPERSAMPLE, &hi))
			hi = -1;
		switch (hi) {
		case 8:
			cvts.tf = RGB2Colr;
			break;
		case 16:
			cvts.tf = RRGGBB2Color;
			SET(C_RFLT|C_TWRD);
			break;
		case 32:
			cvts.tf = RfGfBf2Color;
			SET(C_RFLT|C_TFLT);
			break;
		default:
			quiterr("unsupported bits per sample for RGB");
		}
		break;
	case PHOTOMETRIC_MINISBLACK:
		SET(C_GRY|C_GAMMA);
		setcolrgam(cvts.gamcor);
		cvts.pconf = PLANARCONFIG_CONTIG;
		if (!TIFFGetField(cvts.tif, TIFFTAG_SAMPLESPERPIXEL, &hi) ||
				hi != 1)
			quiterr("unsupported samples per pixel for greyscale");
		if (!TIFFGetField(cvts.tif, TIFFTAG_BITSPERSAMPLE, &hi))
			hi = -1;
		switch (hi) {
		case 8:
			cvts.tf = Gry2Colr;
			break;
		case 16:
			cvts.tf = GGry2Color;
			SET(C_RFLT|C_TWRD);
			break;
		case 32:
			cvts.tf = Gryf2Color;
			SET(C_RFLT|C_TFLT);
			break;
		default:
			quiterr("unsupported bits per sample for Gray");
		}
		break;
	default:
		quiterr("unsupported photometric type");
		break;
	}

	if (!TIFFGetField(cvts.tif, TIFFTAG_IMAGEWIDTH, &cvts.xmax) ||
		!TIFFGetField(cvts.tif, TIFFTAG_IMAGELENGTH, &cvts.ymax))
		quiterr("unknown input image resolution");

	if (!TIFFGetField(cvts.tif, TIFFTAG_STONITS, &cvts.stonits))
		cvts.stonits = -1.;

	if (!TIFFGetField(cvts.tif, TIFFTAG_DATETIME, &cp))
		cvts.capdate[0] = '\0';
	else {
		strncpy(cvts.capdate, cp, 19);
		cvts.capdate[19] = '\0';
	}
	if (!TIFFGetField(cvts.tif, TIFFTAG_ARTIST, &cp))
		cvts.owner[0] = '\0';
	else {
		strncpy(cvts.owner, cp, sizeof(cvts.owner));
		cvts.owner[sizeof(cvts.owner)-1] = '\0';
	}
					/* add to Radiance header */
	if (cvts.pixrat < .99 || cvts.pixrat > 1.01)
		fputaspect(cvts.pixrat, cvts.rfp);
	if (CHK(C_XYZE)) {
		if (cvts.stonits > .0)
			fputexpos(pow(2.,(double)cvts.bradj)/cvts.stonits, cvts.rfp);
		fputformat(CIEFMT, cvts.rfp);
	} else {
		if (CHK(C_PRIM))
			fputprims(cvts.prims, cvts.rfp);
		if (cvts.stonits > .0)
			fputexpos(WHTEFFICACY*pow(2.,(double)cvts.bradj)/cvts.stonits,
					cvts.rfp);
		fputformat(COLRFMT, cvts.rfp);
	}
	if (cvts.capdate[0])
		fprintf(cvts.rfp, "%s %s\n", TMSTR, cvts.capdate);
	if (cvts.owner[0])
		fprintf(cvts.rfp, "%s %s\n", OWNSTR, cvts.owner);

	allocbufs();			/* allocate scanline buffers */
}


static void
tiff2ra(		/* convert TIFF image to Radiance picture */
	int  ac,
	char  *av[]
)
{
	int32	y;
					/* open TIFF input */
	if ((cvts.tif = TIFFOpen(av[ac], "r")) == NULL)
		quiterr("cannot open TIFF input");
					/* open Radiance output */
	if (av[ac+1] == NULL || !strcmp(av[ac+1], "-"))
		cvts.rfp = stdout;
	else if ((cvts.rfp = fopen(av[ac+1], "w")) == NULL)
		quiterr("cannot open Radiance output picture");
	SET_FILE_BINARY(cvts.rfp);
					/* start output header */
	newheader("RADIANCE", cvts.rfp);
	printargs(ac, av, cvts.rfp);

	initfromtif();			/* initialize conversion */

	fputc('\n', cvts.rfp);		/* finish Radiance header */
	fputresolu(pixorder(), (int)cvts.xmax, (int)cvts.ymax, cvts.rfp);

	for (y = 0; y < cvts.ymax; y++)		/* convert image */
		(*cvts.tf)(y);
						/* clean up */
	fclose(cvts.rfp);
	TIFFClose(cvts.tif);
}


static int
headline(			/* process Radiance input header line */
	char	*s,
	void	*p
)
{
	static int	tmstrlen = 0;
	static int	ownstrlen = 0;
	struct CONVST	*cvp = (struct CONVST *)p;
	char		fmt[MAXFMTLEN];

	if (!tmstrlen) {
		tmstrlen = strlen(TMSTR);
		ownstrlen = strlen(OWNSTR);
	}
	if (formatval(fmt, s)) {
		CLR(C_XYZE|C_SPEC);
		if (!strcmp(fmt, CIEFMT))
			SET(C_XYZE);
		else if (!strcmp(fmt, SPECFMT))
			SET(C_SPEC);
		else if (strcmp(fmt, COLRFMT))
			return(-1);
		return(1);
	}
	if (isexpos(s)) {
		cvp->stonits /= exposval(s);
		return(1);
	}
	if (isaspect(s)) {
		cvp->pixrat *= aspectval(s);
		return(1);
	}
	if (isprims(s)) {
		primsval(cvp->prims, s);
		SET(C_PRIM);
		return(1);
	}
	if (isdate(s)) {
		if (s[tmstrlen] == ' ')
			strncpy(cvp->capdate, s+tmstrlen+1, 19);
		else
			strncpy(cvp->capdate, s+tmstrlen, 19);
		cvp->capdate[19] = '\0';
		return(1);
	}
	if (!strncmp(s, OWNSTR, ownstrlen)) {
		char	*cp = s + ownstrlen;

		while (isspace(*cp))
			++cp;
		strncpy(cvp->owner, cp, sizeof(cvp->owner));
		cvp->owner[sizeof(cvp->owner)-1] = '\0';
		for (cp = cvp->owner; *cp; cp++)
			;
		while (cp > cvp->owner && isspace(cp[-1]))
			*--cp = '\0';
		return(1);
	}
	if (isncomp(s)) {
		cvp->ncc = ncompval(s);
		return((3 <= cvp->ncc) & (cvp->ncc < MAXCSAMP) ? 1 : -1);
	}
	if (iswlsplit(s))
		return(wlsplitval(cvp->wpt, s) ? 1 : -1);
	return(0);
}


static void
initfromrad(void)			/* initialize input from a Radiance picture */
{
	int	i1, i2, po;
						/* read Radiance header */
	memcpy(cvts.wpt, WLPART, sizeof(cvts.wpt));
	if (getheader(cvts.rfp, headline, &cvts) < 0 ||
			(po = fgetresolu(&i1, &i2, cvts.rfp)) < 0)
		quiterr("bad Radiance picture");
	cvts.xmax = i1; cvts.ymax = i2;
	for (i1 = 0; i1 < 8; i1++)		/* interpret orientation */
		if (ortab[i1] == po) {
			cvts.orient = i1 + 1;
			break;
		}
	if (i1 >= 8)
		quiterr("internal error 1 in initfromrad");
	if (!(po & YMAJOR))
		cvts.pixrat = 1./cvts.pixrat;
	if (!CHK(C_XYZE))
		cvts.stonits *= WHTEFFICACY;
						/* set up conversion */
	TIFFSetField(cvts.tif, TIFFTAG_COMPRESSION, cvts.comp);
	TIFFSetField(cvts.tif, TIFFTAG_PHOTOMETRIC, cvts.phot);

	switch (cvts.phot) {
	case PHOTOMETRIC_LOGLUV:
		SET(C_RFLT|C_TFLT);
		CLR(C_GRY|C_TWRD);
		if (!CHK(C_XYZE|C_SPEC)) {
			comprgb2xyzWBmat(cvts.cmat,
					CHK(C_PRIM) ? cvts.prims : stdprims);
			SET(C_CXFM);
		}
		if (cvts.comp != COMPRESSION_SGILOG &&
				cvts.comp != COMPRESSION_SGILOG24)
			quiterr("internal error 2 in initfromrad");
		TIFFSetField(cvts.tif, TIFFTAG_SGILOGDATAFMT,
				SGILOGDATAFMT_FLOAT);
		cvts.tf = Color2Luv;
		break;
	case PHOTOMETRIC_LOGL:
		SET(C_GRY|C_RFLT|C_TFLT);
		CLR(C_TWRD);
		if (cvts.comp != COMPRESSION_SGILOG)	
			quiterr("internal error 3 in initfromrad");
		TIFFSetField(cvts.tif, TIFFTAG_SGILOGDATAFMT,
				SGILOGDATAFMT_FLOAT);
		cvts.tf = Color2L;
		break;
	case PHOTOMETRIC_RGB:
		SET(C_GAMMA|C_GAMUT);
		CLR(C_GRY);
		setcolrgam(cvts.gamcor);
		if (CHK(C_TWRD|C_TFLT) ? CHK(C_XYZE|C_SPEC) : CHK(C_XYZE)) {
			compxyz2rgbWBmat(cvts.cmat,
					CHK(C_PRIM) ? cvts.prims : stdprims);
			SET(C_CXFM);
		}
		if (CHK(C_PRIM)) {
			TIFFSetField(cvts.tif, TIFFTAG_PRIMARYCHROMATICITIES,
					(float *)cvts.prims);
			TIFFSetField(cvts.tif, TIFFTAG_WHITEPOINT,
					(float *)cvts.prims[WHT]);
		}
		if (CHK(C_TWRD)) {
			cvts.tf = Color2RRGGBB;
			SET(C_RFLT);
		} else if (CHK(C_TFLT)) {
			TIFFSetField(cvts.tif, TIFFTAG_SAMPLEFORMAT,
					SAMPLEFORMAT_IEEEFP);
			cvts.tf = Color2RfGfBf;
			SET(C_RFLT);
			CLR(C_GAMUT);
		} else
			cvts.tf = Colr2RGB;
		break;
	case PHOTOMETRIC_MINISBLACK:
		SET(C_GRY|C_GAMMA|C_GAMUT);
		setcolrgam(cvts.gamcor);
		if (CHK(C_TWRD)) {
			cvts.tf = Color2GGry;
			SET(C_RFLT);
		} else if (CHK(C_TFLT)) {
			TIFFSetField(cvts.tif, TIFFTAG_SAMPLEFORMAT,
					SAMPLEFORMAT_IEEEFP);
			cvts.tf = Color2Gryf;
			SET(C_RFLT);
		} else
			cvts.tf = Colr2Gry;
		break;
	default:
		quiterr("internal error 4 in initfromrad");
		break;
	}
						/* set other TIFF fields */
	TIFFSetField(cvts.tif, TIFFTAG_IMAGEWIDTH, cvts.xmax);
	TIFFSetField(cvts.tif, TIFFTAG_IMAGELENGTH, cvts.ymax);
	TIFFSetField(cvts.tif, TIFFTAG_SAMPLESPERPIXEL, CHK(C_GRY) ? 1 : 3);
	TIFFSetField(cvts.tif, TIFFTAG_BITSPERSAMPLE,
			CHK(C_TFLT) ? 32 : CHK(C_TWRD) ? 16 : 8);
	TIFFSetField(cvts.tif, TIFFTAG_XRESOLUTION, 72.);
	TIFFSetField(cvts.tif, TIFFTAG_YRESOLUTION, 72./cvts.pixrat);
	TIFFSetField(cvts.tif, TIFFTAG_ORIENTATION, cvts.orient);
	TIFFSetField(cvts.tif, TIFFTAG_RESOLUTIONUNIT, 2);
	TIFFSetField(cvts.tif, TIFFTAG_PLANARCONFIG, cvts.pconf);
	TIFFSetField(cvts.tif, TIFFTAG_STONITS,
			cvts.stonits/pow(2.,(double)cvts.bradj));
	if (cvts.capdate[0])
		TIFFSetField(cvts.tif, TIFFTAG_DATETIME, cvts.capdate);
	if (cvts.owner[0])
		TIFFSetField(cvts.tif, TIFFTAG_ARTIST, cvts.owner);
	if (cvts.comp == COMPRESSION_NONE)
		i1 = TIFFScanlineSize(cvts.tif);
	else
		i1 = 3*cvts.xmax;	/* conservative guess */
	i2 = 8192/i1;				/* compute good strip size */
	if (i2 < 1) i2 = 1;
	TIFFSetField(cvts.tif, TIFFTAG_ROWSPERSTRIP, (uint32)i2);

	allocbufs();				/* allocate scanline buffers */
}


static void
ra2tiff(		/* convert Radiance picture to TIFF image */
	int  ac,
	char  *av[]
)
{
	uint32	y;
						/* open Radiance file */
	if (!strcmp(av[ac], "-"))
		cvts.rfp = stdin;
	else if ((cvts.rfp = fopen(av[ac], "r")) == NULL)
		quiterr("cannot open Radiance input picture");
	SET_FILE_BINARY(cvts.rfp);
						/* open TIFF file */
	if ((cvts.tif = TIFFOpen(av[ac+1], "w")) == NULL)
		quiterr("cannot open TIFF output");

	initfromrad();				/* initialize conversion */

	for (y = 0; y < cvts.ymax; y++)		/* convert image */
		(*cvts.tf)(y);
						/* clean up */
	TIFFClose(cvts.tif);
	fclose(cvts.rfp);
}


static void
ReadRadScan(struct CONVST *cvp)
{
	static struct CONVST	*lastcvp = NULL;
	static COLOR		scomp[MAXCSAMP];

	if (cvp->ncc > 3) {			/* convert spectral pixels */
		SCOLR		sclr;
		SCOLOR		scol;
		COLOR		xyz;
		int		x, n;
		if (lastcvp != cvp) {
		    for (n = cvp->ncc; n--; )
			spec_cie(scomp[n],
				cvp->wpt[0] + (cvp->wpt[3] - cvp->wpt[0])*(n+1)/cvp->ncc,
				cvp->wpt[0] + (cvp->wpt[3] - cvp->wpt[0])*n/cvp->ncc);
		    lastcvp = cvp;
		}
		for (x = 0; x < cvp->xmax; x++) {
			if (getbinary(sclr, cvp->ncc+1, 1, cvp->rfp) != 1)
				goto readerr;
			scolr2scolor(scol, sclr, cvp->ncc);
			setcolor(cvp->r.colors[x], 0, 0, 0);
			for (n = cvp->ncc; n--; ) {
				copycolor(xyz, scomp[n]);
				scalecolor(xyz, scol[n]);
				addcolor(cvp->r.colors[x], xyz);
			}
		}
		return;
	}
	if (freadscan(cvp->r.colors, cvp->xmax, cvp->rfp) >= 0)
		return;
readerr:
	quiterr("error reading Radiance picture");
}


static void
Luv2Color(			/* read/convert/write Luv->COLOR scanline */
	uint32	y
)
{
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != (C_RFLT|C_TFLT))
		quiterr("internal error 1 in Luv2Color");

	if (cvts.pconf != PLANARCONFIG_CONTIG)
		quiterr("cannot handle separate 32-bit color planes");

	if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error reading TIFF input");
					/* also works for float RGB */
	for (x = cvts.xmax; x--; ) {
		setcolor(cvts.r.colors[x], 
				cvts.t.fp[3*x],
				cvts.t.fp[3*x + 1],
				cvts.t.fp[3*x + 2]);
		if (CHK(C_CXFM))
			colortrans(cvts.r.colors[x], cvts.cmat,
					cvts.r.colors[x]);
		if (CHK(C_GAMUT))
			clipgamut(cvts.r.colors[x], cvts.t.fp[3*x + 1],
					CGAMUT_LOWER, cblack, cwhite);
	}
	if (cvts.bradj) {
		double	m = pow(2.,(double)cvts.bradj);
		for (x = cvts.xmax; x--; )
			scalecolor(cvts.r.colors[x], m);
	}

	if (fwritescan(cvts.r.colors, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
}


static void
RRGGBB2Color(			/* read/convert/write RGB16->COLOR scanline */
	uint32	y
)
{
	int	dogamma = (cvts.gamcor < 0.99) | (cvts.gamcor > 1.01);
	double	d;
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != (C_TWRD|C_RFLT))
		quiterr("internal error 1 in RRGGBB2Color");

	if (cvts.pconf != PLANARCONFIG_CONTIG)
		quiterr("cannot handle separate 16-bit color planes");

	if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error reading TIFF input");
	
	for (x = cvts.xmax; x--; ) {
		d = (cvts.t.wp[3*x] + 0.5)*(1./(1L<<16));
		if (dogamma) d = pow(d, cvts.gamcor);
		colval(cvts.r.colors[x],RED) = d;
		d = (cvts.t.wp[3*x + 1] + 0.5)*(1./(1L<<16));
		if (dogamma) d = pow(d, cvts.gamcor);
		colval(cvts.r.colors[x],GRN) = d;
		d = (cvts.t.wp[3*x + 2] + 0.5)*(1./(1L<<16));
		if (dogamma) d = pow(d, cvts.gamcor);
		colval(cvts.r.colors[x],BLU) = d;
		if (CHK(C_CXFM))
			colortrans(cvts.r.colors[x], cvts.cmat,
					cvts.r.colors[x]);
	}
	if (cvts.bradj) {
		d = pow(2.,(double)cvts.bradj);
		for (x = cvts.xmax; x--; )
			scalecolor(cvts.r.colors[x], d);
	}

	if (fwritescan(cvts.r.colors, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
}


static void
L2Color(			/* read/convert/write Lfloat->COLOR scanline */
	uint32	y
)
{
	float	m = pow(2., (double)cvts.bradj);
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != (C_RFLT|C_TFLT|C_GRY))
		quiterr("internal error 1 in L2Color");

	if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error reading TIFF input");
					/* also works for float greyscale */
	for (x = cvts.xmax; x--; ) {
		float	f = cvts.t.fp[x];
		if (cvts.bradj) f *= m;
		setcolor(cvts.r.colors[x], f, f, f);
	}
	if (fwritescan(cvts.r.colors, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
}


static void
RGB2Colr(			/* read/convert/write RGB->COLR scanline */
	uint32	y
)
{
	COLOR	ctmp;
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY))
		quiterr("internal error 1 in RGB2Colr");

	if (cvts.pconf == PLANARCONFIG_CONTIG) {
		if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
			goto readerr;
		for (x = cvts.xmax; x--; ) {
			cvts.r.colrs[x][RED] = cvts.t.bp[3*x];
			cvts.r.colrs[x][GRN] = cvts.t.bp[3*x + 1];
			cvts.r.colrs[x][BLU] = cvts.t.bp[3*x + 2];
		}
	} else {
		if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
			goto readerr;
		if (TIFFReadScanline(cvts.tif,
				(tdata_t)(cvts.t.bp + cvts.xmax), y, 1) < 0)
			goto readerr;
		if (TIFFReadScanline(cvts.tif,
				(tdata_t)(cvts.t.bp + 2*cvts.xmax), y, 2) < 0)
			goto readerr;
		for (x = cvts.xmax; x--; ) {
			cvts.r.colrs[x][RED] = cvts.t.bp[x];
			cvts.r.colrs[x][GRN] = cvts.t.bp[cvts.xmax + x];
			cvts.r.colrs[x][BLU] = cvts.t.bp[2*cvts.xmax + x];
		}
	}

	gambs_colrs(cvts.r.colrs, cvts.xmax);
	if (CHK(C_CXFM))
		for (x = cvts.xmax; x--; ) {
			colr_color(ctmp, cvts.r.colrs[x]);
			colortrans(ctmp, cvts.cmat, ctmp);
			if (CHK(C_GAMUT))	/* !CHK(C_XYZE) */
				clipgamut(ctmp, bright(ctmp), CGAMUT_LOWER,
						cblack, cwhite);
			setcolr(cvts.r.colrs[x], colval(ctmp,RED),
					colval(ctmp,GRN), colval(ctmp,BLU));
		}
	if (cvts.bradj)
		shiftcolrs(cvts.r.colrs, cvts.xmax, cvts.bradj);

	if (fwritecolrs(cvts.r.colrs, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
	return;
readerr:
	quiterr("error reading TIFF input");
}


static void
Gry2Colr(			/* read/convert/write G8->COLR scanline */
	uint32	y
)
{
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != C_GRY)
		quiterr("internal error 1 in Gry2Colr");

	if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error reading TIFF input");

	for (x = cvts.xmax; x--; )
		cvts.r.colrs[x][RED] =
		cvts.r.colrs[x][GRN] =
		cvts.r.colrs[x][BLU] = cvts.t.bp[x];

	gambs_colrs(cvts.r.colrs, cvts.xmax);
	if (cvts.bradj)
		shiftcolrs(cvts.r.colrs, cvts.xmax, cvts.bradj);

	if (fwritecolrs(cvts.r.colrs, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
}


static void
GGry2Color(			/* read/convert/write G16->COLOR scanline */
	uint32	y
)
{
	int	dogamma = (cvts.gamcor < 0.99) | (cvts.gamcor > 1.01);
	double	m;
	double	d;
	int	x;

	if (CHK(C_TFLT|C_TWRD|C_GRY|C_RFLT) != (C_GRY|C_RFLT|C_TWRD))
		quiterr("internal error 1 in GGry2Color");

	if (TIFFReadScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error reading TIFF input");

	if (cvts.bradj)
		m = pow(2., (double)cvts.bradj);
	for (x = cvts.xmax; x--; ) {
		d = (cvts.t.wp[x] + 0.5)*(1./(1L<<16));
		if (dogamma) d = pow(d, cvts.gamcor);
		if (cvts.bradj) d *= m;
		colval(cvts.r.colors[x],RED) =
		colval(cvts.r.colors[x],GRN) =
		colval(cvts.r.colors[x],BLU) = d;
	}
	if (fwritescan(cvts.r.colors, cvts.xmax, cvts.rfp) < 0)
		quiterr("error writing Radiance picture");
}


static void
Color2GGry(			/* read/convert/write COLOR->G16 scanline */
	uint32	y
)
{
	int	dogamma = (cvts.gamcor < 0.99) | (cvts.gamcor > 1.01);
	float	m = pow(2.,(double)cvts.bradj);
	int	x;

	if (CHK(C_RFLT|C_TFLT|C_TWRD|C_GRY) != (C_RFLT|C_TWRD|C_GRY))
		quiterr("internal error 1 in Color2GGry");

	ReadRadScan(&cvts);

	for (x = cvts.xmax; x--; ) {
		float	f = m*( CHK(C_XYZE) ?
						colval(cvts.r.colors[x],CIEY)
						: bright(cvts.r.colors[x]) );
		if (f <= 0)
			cvts.t.wp[x] = 0;
		else if (f >= 1)
			cvts.t.wp[x] = 0xffff;
		else if (dogamma)
			cvts.t.wp[x] = (int)((float)(1L<<16) *
						pow(f, 1./cvts.gamcor));
		else
			cvts.t.wp[x] = (int)((float)(1L<<16) * f);
	}

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}


static void
Color2L(			/* read/convert/write COLOR->Lfloat scanline */
	uint32	y
)
{
	float	m = pow(2.,(double)cvts.bradj);
	int	x;

	if (CHK(C_RFLT|C_TFLT|C_TWRD|C_GRY) != (C_RFLT|C_TFLT|C_GRY))
		quiterr("internal error 1 in Color2L");

	ReadRadScan(&cvts);

	for (x = cvts.xmax; x--; )
		cvts.t.fp[x] = m*( CHK(C_XYZE) ? colval(cvts.r.colors[x],CIEY)
						: bright(cvts.r.colors[x]) );

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}


static void
Color2Luv(			/* read/convert/write COLOR->Luv scanline */
	uint32	y
)
{
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != (C_RFLT|C_TFLT))
		quiterr("internal error 1 in Color2Luv");

	ReadRadScan(&cvts);

	if (CHK(C_CXFM))
		for (x = cvts.xmax; x--; )
			colortrans(cvts.r.colors[x], cvts.cmat,
					cvts.r.colors[x]);
	if (cvts.bradj) {
		double	m = pow(2.,(double)cvts.bradj);
		for (x = cvts.xmax; x--; )
			scalecolor(cvts.r.colors[x], m);
	}
					/* also works for float RGB */
	for (x = cvts.xmax; x--; ) {
		cvts.t.fp[3*x] = colval(cvts.r.colors[x],CIEX);
		cvts.t.fp[3*x+1] = colval(cvts.r.colors[x],CIEY);
		cvts.t.fp[3*x+2] = colval(cvts.r.colors[x],CIEZ);
	}

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}


static void
Color2RRGGBB(			/* read/convert/write COLOR->RGB16 scanline */
	uint32	y
)
{
	int	dogamma = (cvts.gamcor < 0.99) | (cvts.gamcor > 1.01);
	float	m = pow(2.,(double)cvts.bradj);
	int	x, i;

	if (CHK(C_RFLT|C_TFLT|C_TWRD|C_GRY) != (C_RFLT|C_TWRD))
		quiterr("internal error 1 in Color2RRGGBB");

	ReadRadScan(&cvts);

	for (x = cvts.xmax; x--; ) {
	    if (CHK(C_CXFM)) {
			colortrans(cvts.r.colors[x], cvts.cmat,
					cvts.r.colors[x]);
		if (CHK(C_GAMUT))
			clipgamut(cvts.r.colors[x], bright(cvts.r.colors[x]),
					CGAMUT_LOWER, cblack, cwhite);
	    }
	    for (i = 3; i--; ) {
		float	f = m*colval(cvts.r.colors[x],i);
		if (f <= 0)
			cvts.t.wp[3*x + i] = 0;
		else if (f >= 1)
			cvts.t.wp[3*x + i] = 0xffff;
		else if (dogamma)
			cvts.t.wp[3*x + i] = (int)((float)(1L<<16) *
						pow(f, 1./cvts.gamcor));
		else
			cvts.t.wp[3*x + i] = (int)((float)(1L<<16)*f);
	    }
	}

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}


static void
Colr2Gry(			/* read/convert/write COLR->RGB scanline */
	uint32	y
)
{
	int	x;

	if (CHK(C_RFLT|C_TWRD|C_TFLT|C_GRY) != C_GRY)
		quiterr("internal error 1 in Colr2Gry");

	if (fread2colrs(cvts.r.colrs, cvts.xmax, cvts.rfp,
				cvts.ncc, cvts.wpt) < 0)
		quiterr("error reading Radiance picture");

	if (cvts.bradj)
		shiftcolrs(cvts.r.colrs, cvts.xmax, cvts.bradj);
	for (x = cvts.xmax; x--; )
		colval(cvts.r.colrs[x],CIEY) = normbright(cvts.r.colrs[x]);
	colrs_gambs(cvts.r.colrs, cvts.xmax);

	for (x = cvts.xmax; x--; )
		cvts.t.bp[x] = colval(cvts.r.colrs[x],CIEY);

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}


static void
Colr2RGB(			/* read/convert/write COLR->RGB scanline */
	uint32	y
)
{
	COLOR	ctmp;
	int	x;

	if (CHK(C_RFLT|C_TFLT|C_TWRD|C_GRY))
		quiterr("internal error 1 in Colr2RGB");

	if (fread2colrs(cvts.r.colrs, cvts.xmax, cvts.rfp,
				cvts.ncc, cvts.wpt) < 0)
		quiterr("error reading Radiance picture");

	if (cvts.bradj)
		shiftcolrs(cvts.r.colrs, cvts.xmax, cvts.bradj);
	if (CHK(C_CXFM))
		for (x = cvts.xmax; x--; ) {
			colr_color(ctmp, cvts.r.colrs[x]);
			colortrans(ctmp, cvts.cmat, ctmp);
			if (CHK(C_GAMUT))
				clipgamut(ctmp, bright(ctmp), CGAMUT,
						cblack, cwhite);
			setcolr(cvts.r.colrs[x], colval(ctmp,RED),
					colval(ctmp,GRN), colval(ctmp,BLU));
		}
	colrs_gambs(cvts.r.colrs, cvts.xmax);

	for (x = cvts.xmax; x--; ) {
		cvts.t.bp[3*x] = cvts.r.colrs[x][RED];
		cvts.t.bp[3*x+1] = cvts.r.colrs[x][GRN];
		cvts.t.bp[3*x+2] = cvts.r.colrs[x][BLU];
	}

	if (TIFFWriteScanline(cvts.tif, cvts.t.p, y, 0) < 0)
		quiterr("error writing TIFF output");
}
