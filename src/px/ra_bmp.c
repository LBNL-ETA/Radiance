#ifndef lint
static const char RCSid[] = "$Id: ra_bmp.c,v 2.6 2004/04/30 17:56:06 greg Exp $";
#endif
/*
 *  program to convert between RADIANCE and Windows BMP file
 */

#include  <stdio.h>
#include  <string.h>

#include  "platform.h"
#include  "color.h"
#include  "tonemap.h"
#include  "resolu.h"
#include  "bmpfile.h"

int  bradj = 0;				/* brightness adjustment */

double	gamcor = 2.2;			/* gamma correction value */

char  *progname;

static void quiterr(const char *err);
static void tmap2bmp(char *fnin, char *fnout, char *expec,
				RGBPRIMP monpri, double gamval);
static void rad2bmp(FILE *rfp, BMPWriter *bwr, int inv, int gry);
static void bmp2rad(BMPReader *brd, FILE *rfp, int inv);


int
main(int argc, char *argv[])
{
	char		*inpfile=NULL, *outfile=NULL;
	char		*expec = NULL;
	int		reverse = 0;
	RGBPRIMP	rgbp = stdprims;
	RGBPRIMS	myprims;
	RESOLU		rs;
	int		i;
	
	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-' && argv[i][1])
			switch (argv[i][1]) {
			case 'b':
				rgbp = NULL;
				break;
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					expec = argv[++i];
				else
					bradj = atoi(argv[++i]);
				break;
			case 'p':
				if (argc-i < 9)
					goto userr;
				myprims[RED][CIEX] = atof(argv[++i]);
				myprims[RED][CIEY] = atof(argv[++i]);
				myprims[GRN][CIEX] = atof(argv[++i]);
				myprims[GRN][CIEY] = atof(argv[++i]);
				myprims[BLU][CIEX] = atof(argv[++i]);
				myprims[BLU][CIEY] = atof(argv[++i]);
				myprims[WHT][CIEX] = atof(argv[++i]);
				myprims[WHT][CIEY] = atof(argv[++i]);
				if (rgbp == stdprims)
					rgbp = myprims;
				break;
			case 'r':
				reverse = !reverse;
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;

	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	SET_DEFAULT_BINARY();

	if (i <= argc-1 && strcmp(argv[i], "-"))
		inpfile = argv[i];

	if (i == argc-2 && strcmp(argv[i+1], "-"))
		outfile = argv[i+1];
					/* check for tone-mapping */
	if (expec != NULL) {
		if (reverse)
			goto userr;
		tmap2bmp(inpfile, outfile, expec, rgbp, gamcor);
		return(0);
	}

	setcolrgam(gamcor);		/* set up conversion */

	if (reverse) {
		BMPReader       *rdr;
					/* open BMP file or stream */
		if (inpfile != NULL)
			rdr = BMPopenInputFile(inpfile);
		else
			rdr = BMPopenInputStream(stdin);
			
		if (rdr == NULL) {
			fprintf(stderr, "%s: cannot open or recognize BMP\n",
				inpfile != NULL ? inpfile : "<stdin>");
			exit(1);
		}
					/* open Radiance output */
		if (outfile != NULL && freopen(outfile, "w", stdout) == NULL) {
			fprintf(stderr, "%s: cannot open for output\n",
					outfile);
			exit(1);
		}
					/* put Radiance header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		rs.xr = rdr->hdr->width;
		rs.yr = rdr->hdr->height;
		rs.rt = YMAJOR;
					/* write scans downward if we can */
		if (rdr->hdr->yIsDown || inpfile != NULL)
			rs.rt |= YDECR;
		fputsresolu(&rs, stdout);
					/* convert file */
		bmp2rad(rdr, stdout, !rdr->hdr->yIsDown && inpfile!=NULL);
					/* flush output */
		BMPcloseInput(rdr);
		if (fflush(stdout) < 0)
			quiterr("error writing Radiance output");
	} else {
		BMPHeader       *hdr;
		BMPWriter       *wtr;
					/* open Radiance input */
		if (inpfile != NULL && freopen(inpfile, "r", stdin) == NULL) {
			fprintf(stderr, "%s: cannot open input file\n",
					inpfile);
			exit(1);
		}
					/* get header info. */
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				!fgetsresolu(&rs, stdin))
			quiterr("bad Radiance picture format");
					/* initialize BMP header */
		if (rgbp == NULL) {
			hdr = BMPmappedHeader(scanlen(&rs),
						numscans(&rs), 0, 256);
			if (outfile != NULL)
				hdr->compr = BI_RLE8;
		} else
			hdr = BMPtruecolorHeader(scanlen(&rs),
						numscans(&rs), 0);
		if (hdr == NULL)
			quiterr("cannot initialize BMP header");
					/* set up output direction */
		hdr->yIsDown = (rs.rt & YDECR) &&
				((outfile == NULL) | (hdr->compr == BI_RLE8));
					/* open BMP output */
		if (outfile != NULL)
			wtr = BMPopenOutputFile(outfile, hdr);
		else
			wtr = BMPopenOutputStream(stdout, hdr);
		if (wtr == NULL)
			quiterr("cannot allocate writer structure");
					/* convert file */
		rad2bmp(stdin, wtr, !hdr->yIsDown && rs.rt&YDECR, rgbp==NULL);
					/* flush output */
		if (fflush((FILE *)wtr->c_data) < 0)
			quiterr("error writing BMP output");
		BMPcloseOutput(wtr);
	}
	return(0);			/* success */
userr:
	fprintf(stderr,
"Usage: %s [-b][-g gamma][-e spec][-p xr yr xg yg xb yb xw yw] [input|- [output]]\n",
			progname);
	fprintf(stderr,
		"   or: %s -r [-g gamma][-e +/-stops] [input|- [output]]\n",
			progname);
	return(1);
}

/* print message and exit */
static void
quiterr(const char *err)
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}

/* convert Radiance picture to BMP */
static void
rad2bmp(FILE *rfp, BMPWriter *bwr, int inv, int gry)
{
	COLR	*scanin;
	int	y, yend, ystp;
	int     x;
						/* allocate scanline */
	scanin = (COLR *)malloc(bwr->hdr->width*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in rad2bmp");
						/* convert image */
	if (inv) {
		y = bwr->hdr->height - 1;
		ystp = -1; yend = -1;
	} else {
		y = 0;
		ystp = 1; yend = bwr->hdr->height;
	}
						/* convert each scanline */
	for ( ; y != yend; y += ystp) {
		if (freadcolrs(scanin, bwr->hdr->width, rfp) < 0)
			quiterr("error reading Radiance picture");
		if (bradj)
			shiftcolrs(scanin, bwr->hdr->width, bradj);
		for (x = gry ? bwr->hdr->width : 0; x--; )
			scanin[x][GRN] = normbright(scanin[x]);
		colrs_gambs(scanin, bwr->hdr->width);
		if (gry)
			for (x = bwr->hdr->width; x--; )
				bwr->scanline[x] = scanin[x][GRN];
		else
			for (x = bwr->hdr->width; x--; ) {
				bwr->scanline[3*x] = scanin[x][BLU];
				bwr->scanline[3*x+1] = scanin[x][GRN];
				bwr->scanline[3*x+2] = scanin[x][RED];
			}
		bwr->yscan = y;
		x = BMPwriteScanline(bwr);
		if (x != BIR_OK)
			quiterr(BMPerrorMessage(x));
	}
						/* free scanline */
	free((void *)scanin);
}

/* convert BMP file to Radiance */
static void
bmp2rad(BMPReader *brd, FILE *rfp, int inv)
{
	COLR	*scanout;
	int	y, yend, ystp;
	int     x;
						/* allocate scanline */
	scanout = (COLR *)malloc(brd->hdr->width*sizeof(COLR));
	if (scanout == NULL)
		quiterr("out of memory in bmp2rad");
						/* convert image */
	if (inv) {
		y = brd->hdr->height - 1;
		ystp = -1; yend = -1;
	} else {
		y = 0;
		ystp = 1; yend = brd->hdr->height;
	}
						/* convert each scanline */
	for ( ; y != yend; y += ystp) {
		x = BMPseekScanline(y, brd);
		if (x != BIR_OK)
			quiterr(BMPerrorMessage(x));
		for (x = brd->hdr->width; x--; ) {
			RGBquad		rgbq = BMPdecodePixel(x, brd);
			scanout[x][RED] = rgbq.r;
			scanout[x][GRN] = rgbq.g;
			scanout[x][BLU] = rgbq.b;
		}
		gambs_colrs(scanout, brd->hdr->width);
		if (bradj)
			shiftcolrs(scanout, brd->hdr->width, bradj);
		if (fwritecolrs(scanout, brd->hdr->width, rfp) < 0)
			quiterr("error writing Radiance picture");
	}
						/* clean up */
	free((void *)scanout);
}

/* Tone-map and convert Radiance picture */
static void
tmap2bmp(char *fnin, char *fnout, char *expec, RGBPRIMP monpri, double gamval)
{
	int		tmflags;
	BMPHeader       *hdr;
	BMPWriter       *wtr;
	RESOLU		rs;
	FILE		*fp;
	int		xr, yr;
	BYTE		*pa;
	int		i;
					/* check tone-mapping spec */
	i = strlen(expec);
	if (i && !strncmp(expec, "auto", i))
		tmflags = TM_F_CAMERA;
	else if (i && !strncmp(expec, "human", i))
		tmflags = TM_F_HUMAN & ~TM_F_UNIMPL;
	else if (i && !strncmp(expec, "linear", i))
		tmflags = TM_F_LINEAR;
	else
		quiterr("illegal exposure specification (auto|human|linear)");
	if (monpri == NULL) {
		tmflags |= TM_F_BW;
		monpri = stdprims;
	}
					/* open Radiance input */
	if (fnin == NULL)
		fp = stdin;
	else if ((fp = fopen(fnin, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fnin);
		exit(1);
	}
					/* get picture orientation */
	if (fnin != NULL) {
		if (getheader(fp, NULL, NULL) < 0 || !fgetsresolu(&rs, fp))
			quiterr("bad Radiance picture format");
		rewind(fp);
	} else				/* assume stdin has normal orient */
		rs.rt = PIXSTANDARD;
					/* tone-map picture */
	if (tmMapPicture(&pa, &xr, &yr, tmflags, monpri, gamval,
			0., 0., fnin, fp) != TM_E_OK)
		exit(1);
					/* initialize BMP header */
	if (tmflags & TM_F_BW) {
		hdr = BMPmappedHeader(xr, yr, 0, 256);
		if (fnout != NULL)
			hdr->compr = BI_RLE8;
	} else
		hdr = BMPtruecolorHeader(xr, yr, 0);
	if (hdr == NULL)
		quiterr("cannot initialize BMP header");
					/* open BMP output */
	if (fnout != NULL)
		wtr = BMPopenOutputFile(fnout, hdr);
	else
		wtr = BMPopenOutputStream(stdout, hdr);
	if (wtr == NULL)
		quiterr("cannot allocate writer structure");
					/* write to BMP file */
	while (wtr->yscan < yr) {
		BYTE    *scn = pa + xr*((tmflags & TM_F_BW) ? 1 : 3)*
					((rs.rt & YDECR) ?
						(yr-1 - wtr->yscan) :
						wtr->yscan);
		if (tmflags & TM_F_BW)
			memcpy((void *)wtr->scanline, (void *)scn, xr);
		else
			for (i = xr; i--; ) {
				wtr->scanline[3*i] = scn[3*i+BLU];
				wtr->scanline[3*i+1] = scn[3*i+GRN];
				wtr->scanline[3*i+2] = scn[3*i+RED];
			}
		if ((i = BMPwriteScanline(wtr)) != BIR_OK)
			quiterr(BMPerrorMessage(i));
	}
					/* flush output */
	if (fflush((FILE *)wtr->c_data) < 0)
		quiterr("error writing BMP output");
					/* clean up */
	free((void *)pa);
	BMPcloseOutput(wtr);
}
