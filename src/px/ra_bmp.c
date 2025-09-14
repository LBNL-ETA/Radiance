#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 *  program to convert between RADIANCE and Windows BMP file
 */

#include  <math.h>

#include  "rtio.h"
#include  "platform.h"
#include  "color.h"
#include  "tonemap.h"
#include  "resolu.h"
#include  "bmpfile.h"

int		bradj = 0;		/* brightness adjustment */

double		gamcor = 2.2;		/* gamma correction value */

char		*info = "";		/* information header string */
int		infolen = 0;		/* information header length */

extern void quiterr(const char *err);
extern void addBMPcspace(RGBPRIMP pp, double gamma);
extern void tmap2bmp(char *fnin, char *fnout, char *expec,
				RGBPRIMP monpri, double gamval);
extern void rad2bmp(FILE *rfp, BMPWriter *bwr, int inv, RGBPRIMP monpri);
extern void bmp2rad(BMPReader *brd, FILE *rfp, int inv);
extern void info2rad(char *infs, int len, FILE *fout);
extern char *growInfo(int n);
extern gethfunc headline;

#define	add2info(s)	{int _n=strlen(s); memcpy(growInfo(_n), s, _n+1);}
#define clearInfo()	growInfo(-infolen)

RGBPRIMP	rgbinp = stdprims;	/* RGB input primitives */
RGBPRIMS	myinprims;		/* custom primitives holder */


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
	
	fixargv0(argv[0]);		/* assigns progname */

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

	if (expec != NULL) {		/* check for tone-mapping */
		if (reverse)
			goto userr;
		tmap2bmp(inpfile, outfile, expec, rgbp, gamcor);
		return(0);
	}
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
		info2rad(BMPinfo(rdr->hdr), rdr->hdr->infoSiz, stdout);
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
					/* set up conversion */
		setcolrgam(gamcor);
					/* convert file */
		bmp2rad(rdr, stdout, !rdr->hdr->yIsDown & (inpfile!=NULL));
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
					/* get/save header info. */
		if (getheader(stdin, headline, NULL) < 0 ||
					!fgetsresolu(&rs, stdin))
			quiterr("bad Radiance picture format");
					/* record color space */
		addBMPcspace(rgbp, gamcor);
					/* open output/write BMP header */
		if (rgbp == NULL) {
			hdr = BMPmappedHeader(scanlen(&rs),
						numscans(&rs), infolen+1, 256);
			/*
			if (outfile != NULL)
				hdr->compr = BI_RLE8;
			*/
		} else
			hdr = BMPtruecolorHeader(scanlen(&rs),
						numscans(&rs), infolen+1);
		if (hdr == NULL)
			quiterr("cannot create BMP output");
					/* copy info to BMP header */
		strcpy(BMPinfo(hdr), info);
		clearInfo();
					/* set up output direction */
		hdr->yIsDown = ((outfile == NULL) | (hdr->compr == BI_RLE8));
					/* open BMP output */
		if (outfile != NULL)
			wtr = BMPopenOutputFile(outfile, hdr);
		else
			wtr = BMPopenOutputStream(stdout, hdr);
		if (wtr == NULL)
			quiterr("cannot allocate writer structure");
					/* set up conversion */
		setcolrgam(gamcor);
					/* convert file */
		rad2bmp(stdin, wtr, !hdr->yIsDown, rgbp);
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
void
quiterr(const char *err)
{
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}

/* grow (or shrink) saved info header string */
char *
growInfo(int n)
{
	char	*ns = NULL;

	if (infolen + n <= 0) {
		if (info) free(info);
		info = "";
		infolen = 0;
		return(NULL);
	}
	if (infolen)
		info = (char *)realloc(info, infolen+n+1);
	else
		info = (char *)malloc(n+1);

	if (info == NULL)
		quiterr("out of memory in growInfo()");

	if (n > 0) memset(ns = info+infolen, 0, n+1);

	infolen += n;
	return(ns);
}

/* process header line (don't echo) */
int
headline(char *s, void *p)
{
	char	fmt[MAXFMTLEN];

	if (isheadid(s))		/* skip header magic ID */
		return(0);
	if (formatval(fmt, s)) {	/* check if format string */
		if (!strcmp(fmt,COLRFMT))
			return(0);
		if (!strcmp(fmt,CIEFMT)) {
			rgbinp = TM_XYZPRIM;
			return(0);
		}
		if (!strcmp(fmt,SPECFMT))
			return(0);
		return(-1);
	}
	if (isprims(s)) {		/* get input primaries */
		primsval(myinprims, s);
		rgbinp = myinprims;
		return(0);
	}
	if (isexpos(s))
		return(0);		/* ignore this on input */
	if (!strncmp(s, "GAMMA=", 6))
		return(0);		/* should not be here! */
	if (isncomp(s)) {
		NCSAMP = ncompval(s);
		return(0);
	}
	if (iswlsplit(s)) {
		wlsplitval(WLPART, s);
		return(0);
	}
	add2info(s);			/* else save info string */
	return(1);
}

/* add BMP output color space to info string */
void
addBMPcspace(RGBPRIMP pp, double gamma)
{
	char	ibuf[196];
	char	*cp = ibuf;

	sprintf(cp, "GAMMA=%.2f\n", gamma);
	cp += strlen(cp);
	if (pp != NULL) {
		sprintf(cp,
			"%s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
				PRIMARYSTR,
				pp[RED][CIEX],pp[RED][CIEY],
				pp[GRN][CIEX],pp[GRN][CIEY],
				pp[BLU][CIEX],pp[BLU][CIEY],
				pp[WHT][CIEX],pp[WHT][CIEY]);
		/* cp += strlen(cp); */
	}
	add2info(ibuf);
}

/* write out Radiance header from BMP info string */
void
info2rad(char *infs, int len, FILE *fout)
{
	char	*cp;
					/* must fit metadata profile */
	if (len < 3 || infs[0] == '\n' ||
			infs[--len] != '\0' || infs[len-1] != '\n')
		return;			/* not what we expected */
	if (strlen(infs) < len || strstr(infs, "\n\n") != NULL)
		return;			/* also not cool */
					/* check for gamma */
	if ((cp = strstr(infs, "GAMMA=")) != NULL) {
					/* copy what came before */
		fwrite(infs, 1, cp-infs, fout);
		cp += 6;
		gamcor = atof(cp);	/* record setting */
		while (*cp++ != '\n')
			;
		len -= cp - infs;
		infs = cp;		/* & elide from output */
	}
	fputs(infs, fout);		/* copy the remainder */
}

/* convert Radiance picture to BMP */
void
rad2bmp(FILE *rfp, BMPWriter *bwr, int inv, RGBPRIMP monpri)
{
	int	usexfm = 0;
	COLORMAT	xfm;
	COLR	*scanin;
	COLOR	cval;
	int	y, yend, ystp;
	int     x;
						/* allocate scanline */
	scanin = (COLR *)malloc(bwr->hdr->width*sizeof(COLR));
	if (scanin == NULL)
		quiterr("out of memory in rad2bmp");
						/* set up color conversion */
	usexfm = (monpri != NULL) ? (rgbinp != monpri) :
			((rgbinp != TM_XYZPRIM) & (rgbinp != stdprims));
	if (usexfm) {
		RGBPRIMP	destpri = monpri != NULL ? monpri : stdprims;
		double		expcomp = pow(2.0, (double)bradj);
		if (rgbinp == TM_XYZPRIM)
			compxyz2rgbWBmat(xfm, destpri);
		else
			comprgb2rgbWBmat(xfm, rgbinp, destpri);
		for (y = 0; y < 3; y++)
			for (x = 0; x < 3; x++)
				xfm[y][x] *= expcomp;
	}
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
		if (fread2colrs(scanin, bwr->hdr->width, rfp, NCSAMP, WLPART) < 0)
			quiterr("error reading Radiance picture");
		if (usexfm)
			for (x = bwr->hdr->width; x--; ) {
				colr_color(cval, scanin[x]);
				colortrans(cval, xfm, cval);
				setcolr(scanin[x], colval(cval,RED),
						colval(cval,GRN),
						colval(cval,BLU));
			}
		else if (bradj)
			shiftcolrs(scanin, bwr->hdr->width, bradj);
		if ((monpri == NULL) & (rgbinp != TM_XYZPRIM))
			for (x = bwr->hdr->width; x--; )
				scanin[x][GRN] = normbright(scanin[x]);
		colrs_gambs(scanin, bwr->hdr->width);
		if (monpri == NULL)
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
	free(scanin);				/* free scanline */
}

/* convert BMP file to Radiance */
void
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
	free(scanout);				/* clean up */
}

/* Tone-map and convert Radiance picture */
void
tmap2bmp(char *fnin, char *fnout, char *expec, RGBPRIMP monpri, double gamval)
{
	int		tmflags;
	BMPHeader       *hdr;
	BMPWriter       *wtr;
	FILE		*fp;
	int		xr, yr;
	uby8		*pa;
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

	tmflags |= (monpri == NULL)*TM_F_BW;
					/* open Radiance input */
	if (fnin == NULL)
		fp = stdin;
	else if ((fp = fopen(fnin, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fnin);
		exit(1);
	}
					/* tone-map picture */
	if (tmMapPicture(&pa, &xr, &yr, tmflags,
			tmflags&TM_F_BW ? stdprims : monpri, gamval,
			0., 0., fnin, fp) != TM_E_OK)
		exit(1);
					/* try to retrieve info */
	if (fseek(fp, 0L, SEEK_SET) == 0)
		getheader(fp, headline, NULL);
					/* add output color space */
	addBMPcspace(monpri, gamval);
					/* initialize BMP header */
	if (tmflags & TM_F_BW) {
		hdr = BMPmappedHeader(xr, yr, infolen+1, 256);
		if (fnout != NULL)
			hdr->compr = BI_RLE8;
	} else
		hdr = BMPtruecolorHeader(xr, yr, infolen+1);
	if (hdr == NULL)
		quiterr("cannot initialize BMP header");

	strcpy(BMPinfo(hdr), info);	/* copy info if any */
	clearInfo();
					/* open BMP output */
	if (fnout != NULL)
		wtr = BMPopenOutputFile(fnout, hdr);
	else
		wtr = BMPopenOutputStream(stdout, hdr);
	if (wtr == NULL)
		quiterr("cannot allocate writer structure");
					/* write to BMP file */
	while (wtr->yscan < yr) {
		uby8    *scn = pa + xr*((tmflags & TM_F_BW) ? 1 : 3)*
						(yr-1 - wtr->yscan);
		if (tmflags & TM_F_BW)
			memcpy(wtr->scanline, scn, xr);
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
	if (fnin != NULL)
		fclose(fp);
	free(pa);
	BMPcloseOutput(wtr);
}
