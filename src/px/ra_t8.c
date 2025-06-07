#ifndef lint
static const char	RCSid[] = "$Id: ra_t8.c,v 2.19 2025/06/07 05:09:46 greg Exp $";
#endif
/*
 *  ra_t8.c - program to convert between RADIANCE and
 *		Targa 8-bit color-mapped images.
 *
 *	8/22/88		Adapted from ra_pr.c
 */

#include  <math.h>
#include  "platform.h"
#include  "rtio.h"
#include  "rtmisc.h"
#include  "color.h"
#include  "resolu.h"
#include  "clrtab.h"
#include  "targa.h"


#define	 goodpic(h)	(my_imType(h) && my_mapType(h))
#define	 my_imType(h)	(((h)->dataType==IM_CMAP || (h)->dataType==IM_CCMAP) \
				&& (h)->dataBits==8 && (h)->imType==0)
#define	 my_mapType(h)	((h)->mapType==CM_HASMAP && \
				((h)->CMapBits==24 || (h)->CMapBits==32))

#define	 taralloc(h)	(uby8 *)emalloc((h)->x*(h)->y)

uby8  clrtab[256][3];
extern int	samplefac;
double	gamv = 2.2;			/* gamv correction */
int  bradj = 0;				/* brightness adjustment */
char  errmsg[128];
COLR	*inl;
uby8	*tarData;
int  xmax, ymax;

static int getint2(FILE	*fp);
static void putint2(int  i, FILE	*fp);
static void quiterr(char  *err);
static int getthead(struct hdStruct	 *hp, char  *ip, FILE  *fp);
static int putthead(struct hdStruct	 *hp, char  *ip, FILE  *fp);
static int getrhead(struct hdStruct  *h, FILE  *fp);
static void tg2ra(struct hdStruct	 *hp);
static void getmapped(int  nc, int  dith);
static void getgrey(int  nc);
static void writetarga(struct hdStruct	 *h, uby8  *d, FILE  *fp);
static void readtarga(struct hdStruct	 *h, uby8  *data, FILE  *fp);


int
main(
	int  argc,
	char  *argv[]
)
{
	struct hdStruct	 head;
	int  dither = 1;
	int  reverse = 0;
	int  ncolors = 256;
	int  greyscale = 0;
	int  i;
	
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	fixargv0(argv[0]);
	samplefac = 0;

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'd':
				dither = !dither;
				break;
			case 'g':
				gamv = atof(argv[++i]);
				break;
			case 'r':
				reverse = !reverse;
				break;
			case 'b':
				greyscale = 1;
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'c':
				ncolors = atoi(argv[++i]);
				break;
			case 'n':
				samplefac = atoi(argv[++i]);
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		sprintf(errmsg, "cannot open input \"%s\"", argv[i]);
		quiterr(errmsg);
	}
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		sprintf(errmsg, "cannot open output \"%s\"", argv[i+1]);
		quiterr(errmsg);
	}
	if (reverse) {
					/* get header */
		if (getthead(&head, NULL, stdin) < 0)
			quiterr("bad targa file");
		if (!goodpic(&head))
			quiterr("incompatible format");
		xmax = head.x;
		ymax = head.y;
					/* put header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fprtresolu(xmax, ymax, stdout);
					/* convert file */
		tg2ra(&head);
	} else {
		if (getrhead(&head, stdin) < 0)
			quiterr("bad Radiance input");
					/* write header */
		putthead(&head, NULL, stdout);
					/* convert file */
		if (greyscale)
			getgrey(ncolors);
		else
			getmapped(ncolors, dither);
					/* write data */
		writetarga(&head, tarData, stdout);
	}
	quiterr(NULL);
userr:
	fprintf(stderr,
	"Usage: %s [-d][-n samp][-c ncolors][-b][-g gamv][-e +/-stops] input [output]\n",
			progname);
	fprintf(stderr, "   Or: %s -r [-g gamv][-e +/-stops] [input [output]]\n",
			progname);
	exit(1);
}


static int
getint2(			/* get a 2-byte positive integer */
	FILE	*fp
)
{
	int  b1, b2;

	if ((b1 = getc(fp)) == EOF || (b2 = getc(fp)) == EOF)
		quiterr("read error");

	return(b1 | b2<<8);
}


static void
putint2(			/* put a 2-byte positive integer */
	int  i,
	FILE	*fp
)
{
	putc(i&0xff, fp);
	putc(i>>8&0xff, fp);
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


void
eputs(const char *s)
{
	fputs(s, stderr);
}


void
quit(int code)
{
	exit(code);
}


static int
getthead(		/* read header from input */
	struct hdStruct	 *hp,
	char  *ip,
	FILE  *fp
)
{
	int	nidbytes;

	if ((nidbytes = getc(fp)) == EOF)
		return(-1);
	hp->mapType = getc(fp);
	hp->dataType = getc(fp);
	hp->mapOrig = getint2(fp);
	hp->mapLength = getint2(fp);
	hp->CMapBits = getc(fp);
	hp->XOffset = getint2(fp);
	hp->YOffset = getint2(fp);
	hp->x = getint2(fp);
	hp->y = getint2(fp);
	hp->dataBits = getc(fp);
	hp->imType = getc(fp);

	if (ip != NULL)
		if (nidbytes)
			fread((char *)ip, nidbytes, 1, fp);
		else
			*ip = '\0';
	else if (nidbytes)
		fseek(fp, (long)nidbytes, 1);

	return(feof(fp) || ferror(fp) ? -1 : 0);
}


static int
putthead(		/* write header to output */
	struct hdStruct	 *hp,
	char  *ip,
	FILE  *fp
)
{
	if (ip != NULL)
		putc(strlen(ip), fp);
	else
		putc(0, fp);
	putc(hp->mapType, fp);
	putc(hp->dataType, fp);
	putint2(hp->mapOrig, fp);
	putint2(hp->mapLength, fp);
	putc(hp->CMapBits, fp);
	putint2(hp->XOffset, fp);
	putint2(hp->YOffset, fp);
	putint2(hp->x, fp);
	putint2(hp->y, fp);
	putc(hp->dataBits, fp);
	putc(hp->imType, fp);

	if (ip != NULL)
		fputs(ip, fp);

	return(ferror(fp) ? -1 : 0);
}


static int
getrhead(			/* load RADIANCE input file header */
	struct hdStruct  *h,
	FILE  *fp
)
{
					/* get header info. */
	if (checkheader(fp, COLRFMT, NULL) < 0 ||
			fgetresolu(&xmax, &ymax, fp) < 0)
		return(-1);
					/* assign header */
	h->textSize = 0;
	h->mapType = CM_HASMAP;
	h->dataType = IM_CMAP;
	h->mapOrig = 0;
	h->mapLength = 256;
	h->CMapBits = 24;
	h->XOffset = 0;
	h->YOffset = 0;
	h->x = xmax;
	h->y = ymax;
	h->dataBits = 8;
	h->imType = 0;
					/* allocate scanline */
	inl = (COLR *)emalloc(xmax*sizeof(COLR));
					/* allocate targa data */
	tarData = taralloc(h);

	return(0);
}


static void
tg2ra(			/* targa file to RADIANCE file */
	struct hdStruct	 *hp
)
{
	union {
		uby8  c3[256][3];
		uby8  c4[256][4];
	} map;
	COLR  ctab[256];
	COLR  *scanline;
	int  i, j;

					/* get color table */
	if ((hp->CMapBits==24 ? fread((char *)(map.c3+hp->mapOrig),
				3*hp->mapLength,1,stdin) :
			fread((char *)(map.c4+hp->mapOrig),
				4*hp->mapLength,1,stdin)) != 1)
		quiterr("error reading color table");
					/* convert table */
	for (i = hp->mapOrig; i < hp->mapOrig+hp->mapLength; i++)
		if (hp->CMapBits == 24)
			setcolr(ctab[i],
					pow((map.c3[i][2]+.5)/256.,gamv),
					pow((map.c3[i][1]+.5)/256.,gamv),
					pow((map.c3[i][0]+.5)/256.,gamv));
		else
			setcolr(ctab[i],
					pow((map.c4[i][3]+.5)/256.,gamv),
					pow((map.c4[i][2]+.5)/256.,gamv),
					pow((map.c4[i][1]+.5)/256.,gamv));
	if (bradj)
		shiftcolrs(ctab, 256, bradj);
					/* allocate targa data */
	tarData = taralloc(hp);
					/* get data */
	readtarga(hp, tarData, stdin);
					/* allocate input scanline */
	scanline = (COLR *)emalloc(xmax*sizeof(COLR));
					/* convert file */
	for (i = ymax-1; i >= 0; i--) {
		for (j = 0; j < xmax; j++)
			copycolr(scanline[j], ctab[tarData[i*xmax+j]]);
		if (fwritecolrs(scanline, xmax, stdout) < 0)
			quiterr("error writing RADIANCE file");
	}
	free((void *)scanline);
	free((void *)tarData);
}


static void
getmapped(		/* read in and quantize image */
	int  nc,		/* number of colors to use */
	int  dith		/* use dithering? */
)
{
	long  fpos;
	int  y;

	setcolrgam(gamv);
	fpos = ftell(stdin);
	if ((samplefac ? neu_init(xmax*ymax) : new_histo(xmax*ymax)) == -1)
		quiterr("cannot initialized histogram");
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(inl, xmax, stdin) < 0)
			quiterr("error reading Radiance input");
		if (bradj)
			shiftcolrs(inl, xmax, bradj);
		colrs_gambs(inl, xmax);
		if (samplefac)
			neu_colrs(inl, xmax);
		else
			cnt_colrs(inl, xmax);
	}
	if (fseek(stdin, fpos, 0) == EOF)
		quiterr("Radiance input must be from a file");
	if (samplefac)			/* map colors */
		neu_clrtab(nc);
	else
		new_clrtab(nc);
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(inl, xmax, stdin) < 0)
			quiterr("error reading Radiance input");
		if (bradj)
			shiftcolrs(inl, xmax, bradj);
		colrs_gambs(inl, xmax);
		if (samplefac)
			if (dith)
				neu_dith_colrs(tarData+y*xmax, inl, xmax);
			else
				neu_map_colrs(tarData+y*xmax, inl, xmax);
		else
			if (dith)
				dith_colrs(tarData+y*xmax, inl, xmax);
			else
				map_colrs(tarData+y*xmax, inl, xmax);
	}
}


static void
getgrey(			/* read in and convert to greyscale image */
	int  nc		/* number of colors to use */
)
{
	int  y;
	uby8  *dp;
	int  x;

	setcolrgam(gamv);
	dp = tarData+xmax*ymax;;
	for (y = ymax-1; y >= 0; y--) {
		if (freadcolrs(inl, xmax, stdin) < 0)
			quiterr("error reading Radiance input");
		if (bradj)
			shiftcolrs(inl, xmax, bradj);
		x = xmax;
		while (x--)
			inl[x][GRN] = normbright(inl[x]);
		colrs_gambs(inl, xmax);
		x = xmax;
		if (nc < 256)
			while (x--)
				*--dp = ((long)inl[x][GRN]*nc+nc/2)>>8;
		else
			while (x--)
				*--dp = inl[x][GRN];
	}
	for (x = 0; x < nc; x++)
		clrtab[x][RED] = clrtab[x][GRN] =
			clrtab[x][BLU] = ((long)x*256+128)/nc;
}


static void
writetarga(		/* write out targa data */
	struct hdStruct	 *h,
	uby8  *d,
	FILE  *fp
)
{
	int  i, j;

	for (i = 0; i < h->mapLength; i++)	/* write color map */
		for (j = 2; j >= 0; j--)
			putc(clrtab[i][j], fp);
	if (h->dataType == IM_CMAP) {		/* uncompressed */
		if (fwrite((char *)d,h->x*sizeof(uby8),h->y,fp) != h->y)
			quiterr("error writing targa file");
		return;
	}
	quiterr("unsupported output type");
}


static void
readtarga(		/* read in targa data */
	struct hdStruct	 *h,
	uby8  *data,
	FILE  *fp
)
{
	int  cnt, c;
	uby8	*dp;

	if (h->dataType == IM_CMAP) {		/* uncompressed */
		if (fread((char *)data,h->x*sizeof(uby8),h->y,fp) != h->y)
			goto readerr;
		return;
	}
	for (dp = data; dp < data+h->x*h->y; ) {
		if ((c = getc(fp)) == EOF)
			goto readerr;
		cnt = (c & 0x7f) + 1;
		if (c & 0x80) {			/* repeated pixel */
			if ((c = getc(fp)) == EOF)
				goto readerr;
			while (cnt--)
				*dp++ = c;
		} else				/* non-repeating pixels */
			while (cnt--) {
				if ((c = getc(fp)) == EOF)
					goto readerr;
				*dp++ = c;
			}
	}
	return;
readerr:
	quiterr("error reading targa file");
}
