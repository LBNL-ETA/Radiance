#ifndef lint
static const char	RCSid[] = "$Id: ra_t16.c,v 2.12 2025/06/07 05:09:46 greg Exp $";
#endif
/*
 *  ra_t16.c - program to convert between RADIANCE and
 *		Targa 16, 24 and 32-bit images.
 *
 *	11/2/88		Adapted from ra_t8.c
 */

#include  <math.h>
#include  "platform.h"
#include  "rtio.h"
#include  "rtmisc.h"
#include  "color.h"
#include  "resolu.h"
#include  "random.h"
#include  "targa.h"

#define  goodpic(h)	(((h)->dataType==IM_RGB || (h)->dataType==IM_CRGB) \
				&& ((h)->dataBits==16 || (h)->dataBits==24))

#define  taralloc(h)	(unsigned char *)emalloc((h)->x*(h)->y*(h)->dataBits/8)

#define  readtarga(h,d,f)	((h)->dataBits==16 ? \
			readt16(h,(unsigned short*)d,f) : readt24(h,d,f))

#define  writetarga(h,d,f)	((h)->dataBits==16 ? \
			writet16(h,(unsigned short*)d,f) : writet24(h,d,f))

double	gamcor = 2.2;			/* gamma correction */
int  bradj = 0;				/* brightness adjustment */
char  msg[128];

static int getint2(FILE *fp);
static void putint2(int i, FILE *fp);
static void quiterr(char *err);
static int getthead(struct hdStruct *hp, char *ip, FILE *fp);
static int putthead(struct hdStruct *hp, char *ip, FILE *fp);
static void tg2ra(struct hdStruct *hp);
static void ra2tg(struct hdStruct *hp);
static void writet24(struct hdStruct *h, unsigned char *d, FILE *fp);
static void writet16(struct hdStruct *h, unsigned short *d, FILE *fp);
static void readt24(struct hdStruct *h, unsigned char *data, FILE *fp);
static void readt16(struct hdStruct *h, unsigned short *data, FILE *fp);


int
main(int  argc, char  *argv[])
{
	struct hdStruct  head;
	int  reverse = 0;
	int  i;
	
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	fixargv0(argv[0]);

	head.dataBits = 16;
	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':
				gamcor = atof(argv[++i]);
				break;
			case 'r':
				reverse = !reverse;
				break;
			case '2':
				head.dataBits = 16;
				break;
			case '3':
				head.dataBits = 24;
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;
					/* open input file */
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		sprintf(msg, "can't open input \"%s\"", argv[i]);
		quiterr(msg);
	}
					/* open output file */
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		sprintf(msg, "can't open output \"%s\"", argv[i+1]);
		quiterr(msg);
	}
					/* set gamma */
	setcolrgam(gamcor);
					/* convert */
	if (reverse) {
					/* get header */
		if (getthead(&head, NULL, stdin) < 0)
			quiterr("bad targa file");
		if (!goodpic(&head))
			quiterr("incompatible format");
					/* put header */
		newheader("RADIANCE", stdout);
		printargs(i, argv, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fprtresolu(head.x, head.y, stdout);
					/* convert file */
		tg2ra(&head);
	} else {
		if (checkheader(stdin, COLRFMT, NULL) < 0 ||
				fgetresolu(&head.x, &head.y, stdin) < 0)
			quiterr("bad picture file");
					/* assign header */
		head.textSize = 0;
		head.mapType = CM_NOMAP;
		head.dataType = IM_RGB;
		head.XOffset = 0;
		head.YOffset = 0;
		head.imType = 0;
					/* write header */
		putthead(&head, NULL, stdout);
					/* convert file */
		ra2tg(&head);
	}
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-2|-3|-r][-g gamma][-e +/-stops] [input [output]]\n",
			progname);
	exit(1);
}


static int
getint2(			/* get a 2-byte positive integer */
	register FILE	*fp
)
{
	register int  b1, b2;

	if ((b1 = getc(fp)) == EOF || (b2 = getc(fp)) == EOF)
		quiterr("read error");

	return(b1 | b2<<8);
}


static void
putint2(			/* put a 2-byte positive integer */
	register int  i,
	register FILE	*fp
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
	fprintf(stderr, "%s: %s\n", progname, err);
	exit(1);
}


void
eputs(
	char *s
)
{
	fputs(s, stderr);
}


void
quit(
	int code
)
{
	exit(code);
}


static int
getthead(		/* read header from input */
	struct hdStruct  *hp,
	char  *ip,
	register FILE  *fp
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
	struct hdStruct  *hp,
	char  *ip,
	register FILE  *fp
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


static void
tg2ra(			/* targa file to RADIANCE file */
	struct hdStruct  *hp
)
{
	COLR  *scanline;
	unsigned char  *tarData;
	register int  i, j;
					/* skip color table */
	if (hp->mapType == CM_HASMAP)
		fseek(stdin, (long)hp->mapLength*hp->CMapBits/8, 1);
					/* allocate targa data */
	tarData = taralloc(hp);
					/* get data */
	readtarga(hp, tarData, stdin);
					/* allocate input scanline */
	scanline = (COLR *)emalloc(hp->x*sizeof(COLR));
					/* convert file */
	for (i = hp->y-1; i >= 0; i--) {
		if (hp->dataBits == 16) {
			register unsigned short  *dp;
			dp = (unsigned short *)tarData + i*hp->x;
			for (j = 0; j < hp->x; j++) {
				scanline[j][RED] = *dp>>7 & 0xf8;
				scanline[j][GRN] = *dp>>2 & 0xf8;
				scanline[j][BLU] = *dp<<3 & 0xf8;
				dp++;
			}
		} else {	/* hp->dataBits == 24 */
			register unsigned char  *dp;
			dp = (unsigned char *)tarData + i*3*hp->x;
			for (j = 0; j < hp->x; j++) {
				scanline[j][RED] = dp[2];
				scanline[j][GRN] = dp[1];
				scanline[j][BLU] = dp[0];
				dp += 3;
			}
		}
		gambs_colrs(scanline, hp->x);
		if (bradj)
			shiftcolrs(scanline, hp->x, bradj);
		if (fwritecolrs(scanline, hp->x, stdout) < 0)
			quiterr("error writing RADIANCE file");
	}
	free((void *)scanline);
	free((void *)tarData);
}


static void
ra2tg(			/* convert radiance to targa file */
	struct hdStruct  *hp
)
{
	register int	i, j;
	unsigned char  *tarData;
	COLR	*inl;
					/* allocate space for data */
	inl = (COLR *)emalloc(hp->x*sizeof(COLR));
	tarData = taralloc(hp);
					/* convert file */
	for (j = hp->y-1; j >= 0; j--) {
		if (freadcolrs(inl, hp->x, stdin) < 0)
			quiterr("error reading RADIANCE file");
		if (bradj)
			shiftcolrs(inl, hp->x, bradj);
		colrs_gambs(inl, hp->x);
		if (hp->dataBits == 16) {
			register unsigned short  *dp;
			register int  v;
			dp = (unsigned short *)tarData + j*hp->x;
			for (i = 0; i < hp->x; i++) {
				v = inl[i][RED] + (random()&7);
				if (v > 255) v = 255;
				*dp = (v&0xf8)<<7;
				v = inl[i][GRN] + (random()&7);
				if (v > 255) v = 255;
				*dp |= (v&0xf8)<<2;
				v = inl[i][BLU] + (random()&7);
				if (v > 255) v = 255;
				*dp++ |= v>>3;
			}
		} else {	/* hp->dataBits == 24 */
			register unsigned char  *dp;
			dp = (unsigned char *)tarData + j*3*hp->x;
			for (i = 0; i < hp->x; i++) {
				*dp++ = inl[i][BLU];
				*dp++ = inl[i][GRN];
				*dp++ = inl[i][RED];
			}
		}
	}
						/* write out targa data */
	writetarga(hp, tarData, stdout);

	free((void *)inl);
	free((void *)tarData);
}


static void
writet24(		/* write out 24-bit targa data */
	struct hdStruct  *h,
	unsigned char  *d,
	FILE  *fp
)
{
	if (h->dataType == IM_RGB) {		/* uncompressed */
		if (fwrite((char *)d, 3*h->x, h->y, fp) != h->y)
			quiterr("error writing targa file");
		return;
	}
	quiterr("unsupported output type");
}


static void
writet16(		/* write out 16-bit targa data */
	struct hdStruct  *h,
	register unsigned short  *d,
	FILE  *fp
)
{
	register int  cnt;

	if (h->dataType == IM_RGB) {		/* uncompressed */
		for (cnt = h->x*h->y; cnt-- > 0; )
			putint2(*d++, fp);
		if (ferror(fp))
			quiterr("error writing targa file");
		return;
	}
	quiterr("unsupported output type");
}


static void
readt24(		/* read in 24-bit targa data */
	register struct hdStruct  *h,
	unsigned char  *data,
	FILE  *fp
)
{
	register int  cnt, c;
	register unsigned char  *dp;
	int  r, g, b;

	if (h->dataType == IM_RGB) {		/* uncompressed */
		if (fread((char *)data, 3*h->x, h->y, fp) != h->y)
			goto readerr;
		return;
	}
	for (dp = data; dp < data+3*h->x*h->y; ) {
		if ((c = getc(fp)) == EOF)
			goto readerr;
		cnt = (c & 0x7f) + 1;
		if (c & 0x80) {			/* repeated pixel */
			b = getc(fp); g = getc(fp);
			if ((r = getc(fp)) == EOF)
				goto readerr;
			while (cnt--) {
				*dp++ = b;
				*dp++ = g;
				*dp++ = r;
			}
		} else				/* non-repeating pixels */
			while (cnt--) {
				*dp++ = getc(fp); *dp++ = getc(fp);
				if ((r = getc(fp)) == EOF)
					goto readerr;
				*dp++ = r;
			}
	}
	return;
readerr:
	quiterr("error reading targa file");
}


static void
readt16(		/* read in 16-bit targa data */
	register struct hdStruct  *h,
	unsigned short  *data,
	FILE  *fp
)
{
	register int  cnt, c;
	register unsigned short  *dp;

	if (h->dataType == IM_RGB) {		/* uncompressed */
		dp = data;
		for (cnt = h->x*h->y; cnt-- > 0; )
			*dp++ = getint2(fp);
		return;
	}
	for (dp = data; dp < data+h->x*h->y; ) {
		if ((c = getc(fp)) == EOF)
			goto readerr;
		cnt = (c & 0x7f) + 1;
		if (c & 0x80) {			/* repeated pixel */
			c = getint2(fp);
			while (cnt--)
				*dp++ = c;
		} else				/* non-repeating pixels */
			while (cnt--)
				*dp++ = getint2(fp);
	}
	return;
readerr:
	quiterr("error reading targa file");
}
