#ifndef lint
static const char RCSid[] = "$Id: bmpfile.c,v 2.10 2004/04/30 17:00:29 greg Exp $";
#endif
/*
 *  Windows and OS/2 BMP file support
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmpfile.h"

#ifdef getc_unlocked		/* avoid horrendous overhead of flockfile */
#define getc    getc_unlocked
#define putc    putc_unlocked
#endif

/* get corresponding error message */
const char *
BMPerrorMessage(int ec)
{
	switch (ec) {
	case BIR_OK:
		return "No error";
	case BIR_EOF:
		return "End of BMP image";
	case BIR_TRUNCATED:
		return "Truncated BMP image";
	case BIR_UNSUPPORTED:
		return "Unsupported BMP feature";
	case BIR_RLERROR:
		return "BMP runlength encoding error";
	case BIR_SEEKERR:
		return "BMP seek error";
	}
	return "Unknown BMP error";
}

/* check than header is sensible */
static int
BMPheaderOK(const BMPHeader *hdr)
{
	if (!hdr)
		return 0;
	if ((hdr->width <= 0) | (hdr->height <= 0))
		return 0;
	switch (hdr->bpp) {		/* check compression */
	case 1:
	case 24:
		if (hdr->compr != BI_UNCOMPR)
			return 0;
		break;
	case 16:
	case 32:
		if ((hdr->compr != BI_UNCOMPR) & (hdr->compr != BI_BITFIELDS))
			return 0;
		break;
	case 4:
		if ((hdr->compr != BI_UNCOMPR) & (hdr->compr != BI_RLE4))
			return 0;
		break;
	case 8:
		if ((hdr->compr != BI_UNCOMPR) & (hdr->compr != BI_RLE8))
			return 0;
		break;
	default:
		return 0;
	}
	if (hdr->compr == BI_BITFIELDS && (BMPbitField(hdr)[0] &
				BMPbitField(hdr)[1] & BMPbitField(hdr)[2]))
		return 0;
	if ((hdr->nColors < 0) | (hdr->impColors < 0))
		return 0;
	if (hdr->impColors > hdr->nColors)
		return 0;
	if (hdr->nColors > BMPpalLen(hdr))
		return 0;
	return 1;
}

					/* compute uncompressed scan size */
#define getScanSiz(h)   ( ((((h)->bpp*(h)->width+7) >>3) + 3) & ~03 )

					/* get next byte from reader */
#define rdbyte(c,br)    ((br)->fpos += (c=(*(br)->cget)((br)->c_data))!=EOF, c)

/* read n bytes */
static int
rdbytes(char *bp, uint32 n, BMPReader *br)
{
	int     c;

	while (n--) {
		if (rdbyte(c, br) == EOF)
			return BIR_TRUNCATED;
		*bp++ = c;
	}
	return BIR_OK;
}

/* read 32-bit integer in littlendian order */
static int32
rdint32(BMPReader *br)
{
	int32   i;
	int     c;

	i = rdbyte(c, br);
	i |= rdbyte(c, br) << 8;
	i |= rdbyte(c, br) << 16;
	i |= rdbyte(c, br) << 24;
	return i;			/* -1 on EOF */
}

/* read 16-bit unsigned integer in littlendian order */
static int
rduint16(BMPReader *br)
{
	int     i;
	int     c;
	
	i = rdbyte(c, br);
	i |= rdbyte(c, br) << 8;
	return i;			/* -1 on EOF */
}

/* seek on reader or return 0 (BIR_OK) on success */
static int
rdseek(uint32 pos, BMPReader *br)
{
	if (pos == br->fpos)
		return BIR_OK;
	if (br->seek == NULL || (*br->seek)(pos, br->c_data) != 0)
		return BIR_SEEKERR;
	br->fpos = pos;
	return BIR_OK;
}

/* open BMP stream for reading and get first scanline */
BMPReader *
BMPopenReader(int (*cget)(void *), int (*seek)(uint32, void *), void *c_data)
{
	BMPReader       *br;
	uint32		bmPos, hdrSiz;
	int		palLen;

	if (cget == NULL)
		return NULL;
	int     magic[2];		/* check magic number */
	magic[0] = (*cget)(c_data);
	if (magic[0] != 'B')
		return NULL;
	magic[1] = (*cget)(c_data);
	if (magic[1] != 'M' && magic[1] != 'A')
		return NULL;
	br = (BMPReader *)calloc(1, sizeof(BMPReader));
	if (br == NULL)
		return NULL;
	br->cget = cget;
	br->seek = seek;
	br->c_data = c_data;
	br->hdr = (BMPHeader *)malloc(sizeof(BMPHeader));
	if (br->hdr == NULL)
		goto err;
	br->fpos = 2;
					/* read & verify file header */
	(void)rdint32(br);			/* file size */
	(void)rdint32(br);			/* reserved word */
	bmPos = rdint32(br);			/* offset to bitmap */
	hdrSiz = 2 + 3*4 + rdint32(br);		/* header size */
	if (hdrSiz < 2 + 6*4 + 2*2 + 6*4)
		goto err;
	br->hdr->width = rdint32(br);		/* bitmap width */
	br->hdr->height = rdint32(br);		/* bitmap height */
	if (((br->hdr->width <= 0) | (br->hdr->height == 0)))
		goto err;
	if ((br->hdr->yIsDown = br->hdr->height < 0))
		br->hdr->height = -br->hdr->height;
	if (rduint16(br) != 1)			/* number of planes */
		goto err;
	br->hdr->bpp = rduint16(br);		/* bits per pixel */
	br->hdr->compr = rdint32(br);		/* compression mode */
	(void)rdint32(br);			/* bitmap size */
	br->hdr->hRes = rdint32(br);		/* horizontal resolution */
	br->hdr->vRes = rdint32(br);		/* vertical resolution */
	br->hdr->nColors = rdint32(br);		/* # colors used */
	br->hdr->impColors = rdint32(br);       /* # important colors */
	if (br->hdr->impColors < 0)
		goto err;			/* catch premature EOF */
	if (!BMPheaderOK(br->hdr))
		goto err;
	palLen = BMPpalLen(br->hdr);
	if (br->hdr->bpp <= 8) {		/* normalize color counts */
		if (br->hdr->nColors <= 0)
			br->hdr->nColors = palLen;
		if (br->hdr->impColors <= 0)
			br->hdr->impColors = br->hdr->nColors;
	}
						/* extend header */
	if (bmPos < hdrSiz + sizeof(RGBquad)*palLen)
		goto err;
	br->hdr->infoSiz = bmPos - (hdrSiz + sizeof(RGBquad)*palLen);
	if (palLen > 0 || br->hdr->infoSiz > 0) {
		br->hdr = (BMPHeader *)realloc((void *)br->hdr,
					sizeof(BMPHeader) +
					sizeof(RGBquad)*palLen +
					br->hdr->infoSiz);
		if (br->hdr == NULL)
			goto err;
	}
						/* read colors or fields */
	if (br->hdr->compr == BI_BITFIELDS) {
		BMPbitField(br->hdr)[0] = (uint32)rdint32(br);
		BMPbitField(br->hdr)[1] = (uint32)rdint32(br);
		BMPbitField(br->hdr)[2] = (uint32)rdint32(br);
	} else if (rdbytes((char *)br->hdr->palette,
			sizeof(RGBquad)*palLen, br) != BIR_OK)
		goto err;
						/* read add'l information */
	if (rdbytes(BMPinfo(br->hdr), br->hdr->infoSiz, br) != BIR_OK)
		goto err;
						/* read first scanline */
	br->scanline = (uint8 *)calloc(getScanSiz(br->hdr), sizeof(uint8));
	if (br->scanline == NULL)
		goto err;
	br->yscan = -1;
	if (seek != NULL && ((br->hdr->compr == BI_RLE8) |
					(br->hdr->compr == BI_RLE4))) {
		BMPReader       *newbr = (BMPReader *)realloc((void *)br,
						sizeof(BMPReader) +
						sizeof(br->scanpos[0]) *
							br->hdr->height);
		if (newbr == NULL)
			goto err;
		br = newbr;
		memset((void *)(br->scanpos + 1), 0,
				sizeof(br->scanpos[0])*br->hdr->height);
	}
	br->scanpos[0] = br->fpos;
	if (BMPreadScanline(br) == BIR_OK)
		return br;
err:
	if (br->hdr != NULL)
		free((void *)br->hdr);
	if (br->scanline != NULL)
		free((void *)br->scanline);
	free((void *)br);
	return NULL;
}

/* determine if image is grayscale */
int
BMPisGrayscale(const BMPHeader *hdr)
{
	const RGBquad   *rgbp;
	int		n;

	if (hdr == NULL)
		return -1;
	if (hdr->bpp > 8)		/* assume they had a reason for it */
		return 0;
	for (rgbp = hdr->palette, n = hdr->nColors; n-- > 0; rgbp++)
		if (((rgbp->r != rgbp->g) | (rgbp->g != rgbp->b)))
			return 0;
	return 1;			/* all colors neutral in map */
}

/* read and decode next BMP scanline */
int
BMPreadScanline(BMPReader *br)
{
	int     n;
	int8    *sp;

	if (br->yscan + 1 >= br->hdr->height)
		return BIR_EOF;
	br->yscan++;			/* prepare to read */
	n = getScanSiz(br->hdr);	/* reading uncompressed data? */
	if (br->hdr->compr == BI_UNCOMPR || br->hdr->compr == BI_BITFIELDS)
		return rdbytes((char *)br->scanline, n, br);
	/*
	 * RLE4/RLE8 Decoding
	 *
	 * Certain aspects of this scheme are completely insane, so
	 * we don't support them.  Fortunately, they rarely appear.
	 * One is the mid-file EOD (0x0001) and another is the ill-conceived
	 * "delta" (0x0002), which is like a "goto" statement for bitmaps.
	 * Whoever thought this up should be wrestled to the ground and told
	 * why it's impossible to support such a scheme in any reasonable way.
	 * Also, RLE4 mode allows runs to stop halfway through a byte,
	 * which is likewise uncodeable, so we don't even try.
	 * Finally, the scanline break is ambiguous -- we assume here that
	 * it is required at the end of each scanline, though I haven't
	 * found anywhere this is written.  Otherwise, we would read to
	 * the end of the scanline, assuming the next bit of data belongs
	 * the following scan.  If a break follows the last pixel, as it
	 * seems to in the files I've tested out of Photoshop, you end up
	 * painting every other line black.  Also, I assume any skipped
	 * pixels are painted with color 0, which is often black.  Nowhere
	 * is it specified what we should assume for missing pixels.  This
	 * is undoubtedly the most brain-dead format I've ever encountered.
	 */
	sp = br->scanline;
	n = br->hdr->width;
	if (br->hdr->compr == BI_RLE4)
		n = (n + 1) >> 1;
	while (n > 0) {
		int     skipOdd, len, val;

		if (rdbyte(len, br) == EOF)
			return BIR_TRUNCATED;
		if (len > 0) {		/* got a run */
			if (br->hdr->compr == BI_RLE4) {
				if (len & 1)
					return BIR_UNSUPPORTED;
				len >>= 1;
			}
			if (len > n)
				return BIR_RLERROR;
			if (rdbyte(val, br) == EOF)
				return BIR_TRUNCATED;
			n -= len;
			while (len--)
				*sp++ = val;
			continue;
		}
					/* check for escape */
		switch (rdbyte(len, br)) {
		case EOF:
			return BIR_TRUNCATED;
		case 0:			/* end of line */
			while (n--)
				*sp++ = 0;
			/* leaves n == -1 as flag for test after loop */
			continue;
		case 1:			/* end of bitmap */
		case 2:			/* delta */
			return BIR_UNSUPPORTED;
		}
					/* absolute mode */
		if (br->hdr->compr == BI_RLE4) {
			if (len & 1)
				return BIR_UNSUPPORTED;
			len >>= 1;
		}
		skipOdd = len & 1;
		if (len > n)
			return BIR_RLERROR;
		n -= len;
		while (len--) {
			if (rdbyte(val, br) == EOF)
				return BIR_TRUNCATED;
			*sp++ = val;
		}
		if (skipOdd && rdbyte(val, br) == EOF)
			return BIR_TRUNCATED;
	}
					/* verify break at end of line */
	if (!n && (rdbyte(n, br) != 0 || (rdbyte(n, br) != 0 &&
				(n != 1 || br->yscan != br->hdr->height-1))))
		return BIR_RLERROR;
	if (br->seek != NULL)		/* record next scanline position */
		br->scanpos[br->yscan + 1] = br->fpos;
	return BIR_OK;
}

/* read a specific scanline */
int
BMPseekScanline(int y, BMPReader *br)
{
	int     rv;
					/* check arguments */
	if (br == NULL)
		return BIR_EOF;
	if (y < 0)
		return BIR_SEEKERR;
	if (y >= br->hdr->height)
		return BIR_EOF;
					/* already read? */
	if (y == br->yscan)
		return BIR_OK;
					/* shall we seek? */
	if (y != br->yscan + 1 && br->seek != NULL) {
		int     yseek;
		uint32  seekp;
		if (br->hdr->compr == BI_UNCOMPR ||
					br->hdr->compr == BI_BITFIELDS) {
			yseek = y;
			seekp = br->scanpos[0] + y*getScanSiz(br->hdr);
		} else {
			yseek = br->yscan + 1;
			while (yseek < y && br->scanpos[yseek+1] != 0)
				++yseek;
			if (y < yseek && br->scanpos[yseek=y] == 0)
				return BIR_SEEKERR;
			seekp = br->scanpos[yseek];
		}
		if ((rv = rdseek(seekp, br)) != BIR_OK)
			return rv;
		br->yscan = yseek - 1;
	} else if (y < br->yscan)       /* else we can't back up */
		return BIR_SEEKERR;
					/* read until we get there */
	while (br->yscan < y)
		if ((rv = BMPreadScanline(br)) != BIR_OK)
			return rv;
	return BIR_OK;
}

/* get ith pixel from last scanline */
RGBquad
BMPdecodePixel(int i, const BMPReader *br)
{
	static const uint32     std16mask[3] = {0x7c00, 0x3e0, 0x1f};
	static const RGBquad    black = {0, 0, 0, 0};
	const uint32		*mask;
	const uint8		*pp;
	uint32			pval, v;
	RGBquad			cval;
	
	if (((br == NULL) | (i < 0)) || i >= br->hdr->width)
		return black;

	cval.padding = 0;

	switch (br->hdr->bpp) {
	case 24:
		pp = br->scanline + 3*i;
		cval.b = *pp++;
		cval.g = *pp++;
		cval.r = *pp;
		return cval;
	case 32:
		if (br->hdr->compr == BI_UNCOMPR)
			return ((RGBquad *)br->scanline)[i];
						/* convert bit fields */
		pp = br->scanline + 4*i;
		pval = *pp++;
		pval |= *pp++ << 8;
		pval |= *pp++ << 16;
		pval |= *pp << 24;
		mask = BMPbitField(br->hdr);
		v = pval & mask[0];
		while (v & ~0xff) v >>= 8;
		cval.r = v;
		v = pval & mask[1];
		while (v & ~0xff) v >>= 8;
		cval.g = v;
		v = pval & mask[2];
		while (v & ~0xff) v >>= 8;
		cval.b = v;
		return cval;
	case 8:
		return br->hdr->palette[br->scanline[i]];
	case 1:
		return br->hdr->palette[br->scanline[i>>3]>>((7-i)&7) & 1];
	case 4:
		return br->hdr->palette[br->scanline[i>>1]>>(i&1?4:0) & 0xf];
	case 16:
		pp = br->scanline + 2*i;
		pval = *pp++;
		pval |= *pp++ << 8;
		mask = std16mask;
		if (br->hdr->compr == BI_BITFIELDS)
			mask = BMPbitField(br->hdr);
		cval.r = ((pval & mask[0]) << 8) / (mask[0] + 1);
		cval.g = ((pval & mask[1]) << 8) / (mask[1] + 1);
		cval.b = ((pval & mask[2]) << 8) / (mask[2] + 1);
		return cval;
	}
	return black;				/* should never happen */
}

/* free BMP reader resources */
void
BMPfreeReader(BMPReader *br)
{
	if (br == NULL)
		return;
	free((void *)br->hdr);
	free((void *)br->scanline);
	free((void *)br);
}

/* stdio getc() callback */
int
stdio_getc(void *p)
{
	if (!p)
		return EOF;
	return getc((FILE *)p);
}

/* stdio putc() callback */
void
stdio_putc(int c, void *p)
{
	if (p)
		putc(c, (FILE *)p);
}

/* stdio fseek() callback */
int
stdio_fseek(uint32 pos, void *p)
{
	if (!p)
		return -1;
	return fseek((FILE *)p, (long)pos, 0);
}

/* allocate uncompressed (24-bit) RGB header */
BMPHeader *
BMPtruecolorHeader(int xr, int yr, int infolen)
{
	BMPHeader       *hdr;
	
	if (xr <= 0 || yr <= 0 || infolen < 0)
		return NULL;
	hdr = (BMPHeader *)malloc(sizeof(BMPHeader) - sizeof(hdr->palette) +
					infolen);
	if (hdr == NULL)
		return NULL;
	hdr->width = xr;
	hdr->height = yr;
	hdr->yIsDown = 0;			/* default to upwards order */
	hdr->bpp = 24;
	hdr->compr = BI_UNCOMPR;
	hdr->hRes = hdr->vRes = 2835;		/* default to 72 ppi */
	hdr->nColors = hdr->impColors = 0;
	hdr->infoSiz = infolen;
	return hdr;
}

/* allocate color-mapped header (defaults to minimal grayscale) */
BMPHeader *
BMPmappedHeader(int xr, int yr, int infolen, int ncolors)
{
	int		n;
	BMPHeader       *hdr;
	
	if (xr <= 0 || yr <= 0 || infolen < 0 || ncolors < 2)
		return NULL;
	if (ncolors <= 2)
		n = 1;
	else if (ncolors <= 16)
		n = 4;
	else if (ncolors <= 256)
		n = 8;
	else
		return NULL;
	hdr = (BMPHeader *)malloc(sizeof(BMPHeader) +
					sizeof(RGBquad)*(1<<n) -
					sizeof(hdr->palette) +
					infolen);
	if (hdr == NULL)
		return NULL;
	hdr->width = xr;
	hdr->height = yr;
	hdr->yIsDown = 0;			/* default to upwards order */
	hdr->bpp = n;
	hdr->compr = BI_UNCOMPR;		/* compression needs seek */
	hdr->hRes = hdr->vRes = 2835;		/* default to 72 ppi */
	hdr->nColors = ncolors;
	hdr->impColors = 0;			/* says all colors important */
	hdr->infoSiz = infolen;
	memset((void *)hdr->palette, 0, sizeof(RGBquad)*(1<<n) + infolen);
	for (n = ncolors; n--; )
		hdr->palette[n].r = hdr->palette[n].g =
			hdr->palette[n].b = n*255/(ncolors-1);
	return hdr;
}

					/* put byte to writer */
#define wrbyte(c,bw)    ( (*(bw)->cput)(c,(bw)->c_data), \
				++(bw)->fpos > (bw)->flen ? \
					((bw)->flen = (bw)->fpos) : \
					(bw)->fpos )

/* write out a string of bytes */
static void
wrbytes(char *bp, uint32 n, BMPWriter *bw)
{
	while (n--)
		wrbyte(*bp++, bw);
}

/* write 32-bit integer in littlendian order */
static void
wrint32(int32 i, BMPWriter *bw)
{
	wrbyte(i& 0xff, bw);
	wrbyte(i>>8 & 0xff, bw);
	wrbyte(i>>16 & 0xff, bw);
	wrbyte(i>>24 & 0xff, bw);
}

/* write 16-bit unsigned integer in littlendian order */
static void
wruint16(uint16 ui, BMPWriter *bw)
{
	wrbyte(ui & 0xff, bw);
	wrbyte(ui>>8 & 0xff, bw);
}

/* seek to the specified file position, returning 0 (BIR_OK) on success */
static int
wrseek(uint32 pos, BMPWriter *bw)
{
	if (pos == bw->fpos)
		return BIR_OK;
	if (bw->seek == NULL)
		return BIR_SEEKERR;
	if ((*bw->seek)(pos, bw->c_data) != 0)
		return BIR_SEEKERR;
	bw->fpos = pos;
	if (pos > bw->flen)
		bw->flen = pos;
	return BIR_OK;
}

/* open BMP stream for writing */
BMPWriter *
BMPopenWriter(void (*cput)(int, void *), int (*seek)(uint32, void *),
			void *c_data, BMPHeader *hdr)
{
	BMPWriter       *bw;
	uint32		hdrSiz, palSiz, scanSiz, bmSiz;
						/* check arguments */
	if (cput == NULL)
		return NULL;
	if (!BMPheaderOK(hdr))
		return NULL;
	if ((hdr->bpp == 16) | (hdr->compr == BI_RLE4))
		return NULL;			/* unsupported */
/* no seek means we may have the wrong file length, but most app's don't care
	if (seek == NULL && ((hdr->compr == BI_RLE8) | (hdr->compr == BI_RLE4)))
		return NULL;
*/
						/* compute sizes */
	hdrSiz = 2 + 6*4 + 2*2 + 6*4;
	if (hdr->compr == BI_BITFIELDS)
		hdrSiz += sizeof(uint32)*3;
	palSiz = sizeof(RGBquad)*BMPpalLen(hdr);
	scanSiz = getScanSiz(hdr);
	bmSiz = hdr->height*scanSiz;		/* wrong if compressed */
						/* initialize writer */
	bw = (BMPWriter *)malloc(sizeof(BMPWriter));
	if (bw == NULL)
		return NULL;
	bw->hdr = hdr;
	bw->yscan = 0;
	bw->scanline = (uint8 *)calloc(scanSiz, sizeof(uint8));
	if (bw->scanline == NULL) {
		free((void *)bw);
		return NULL;
	}
	bw->fbmp = hdrSiz + palSiz + hdr->infoSiz;
	bw->fpos = bw->flen = 0;
	bw->cput = cput;
	bw->seek = seek;
	bw->c_data = c_data;
						/* write out header */
	wrbyte('B', bw); wrbyte('M', bw);       /* magic number */
	wrint32(bw->fbmp + bmSiz, bw);		/* file size */
	wrint32(0, bw);				/* reserved word */
	wrint32(bw->fbmp, bw);			/* offset to bitmap */
	wrint32(hdrSiz - bw->fpos, bw);		/* info header size */
	wrint32(hdr->width, bw);		/* bitmap width */
	if (hdr->yIsDown)			/* bitmap height */
		wrint32(-hdr->height, bw);
	else
		wrint32(hdr->height, bw);
	wruint16(1, bw);			/* number of planes */
	wruint16(hdr->bpp, bw);			/* bits per pixel */
	wrint32(hdr->compr, bw);		/* compression mode */
	wrint32(bmSiz, bw);			/* bitmap size */
	wrint32(hdr->hRes, bw);			/* horizontal resolution */
	wrint32(hdr->vRes, bw);			/* vertical resolution */
	wrint32(hdr->nColors, bw);		/* # colors used */
	wrint32(hdr->impColors, bw);		/* # important colors */
						/* write out color palette */
	wrbytes((char *)hdr->palette, palSiz, bw);
						/* write add'l information */
	wrbytes(BMPinfo(hdr), hdr->infoSiz, bw);
#ifndef NDEBUG
	if (bw->fpos != bw->fbmp) {
		fputs("Coding error 1 in BMPopenWriter\n", stderr);
		exit(1);
	}
#endif
	return bw;
}

/* find position of next run of 5 or more identical bytes, or 255 if none */
static int
findNextRun(const int8 *bp, int len)
{
	int     pos, cnt;
						/* look for run */
	for (pos = 0; (len > 0) & (pos < 255); pos++, bp++, len--) {
		if (len < 5)			/* no hope left? */
			continue;
		cnt = 1;			/* else let's try it */
		while (bp[cnt] == bp[0])
			if (++cnt >= 5)
				return pos;     /* long enough */
	}
	return pos;				/* didn't find any */
}
				
/* write the current scanline */
int
BMPwriteScanline(BMPWriter *bw)
{
	const int8      *sp;
	int		n;

	if (bw->yscan >= bw->hdr->height)
		return BIR_EOF;
						/* writing uncompressed? */
	if (bw->hdr->compr == BI_UNCOMPR || bw->hdr->compr == BI_BITFIELDS) {
		uint32  scanSiz = getScanSiz(bw->hdr);
		uint32  slpos = bw->fbmp + bw->yscan*scanSiz;
		if (wrseek(slpos, bw) != BIR_OK)
			return BIR_SEEKERR;
		wrbytes((char *)bw->scanline, scanSiz, bw);
		bw->yscan++;
		return BIR_OK;
	}
	/*
	 * RLE8 Encoding
	 *
	 * See the notes in BMPreadScanline() on this encoding.  Needless
	 * to say, we avoid the nuttier aspects of this specification.
	 * We also assume that every scanline ends in a line break
	 * (0x0000) except for the last, which ends in a bitmap break
	 * (0x0001).  We don't support RLE4 at all; it's too awkward.
	 */
	sp = bw->scanline;
	n = bw->hdr->width;
	while (n > 0) {
		int     cnt, val;
		cnt = findNextRun(sp, n);       /* 0-255 < n */
		if (cnt >= 3) {			/* output absolute */
			int     skipOdd = cnt & 1;
			wrbyte(0, bw);
			wrbyte(cnt, bw);
			n -= cnt;
			while (cnt--)
				wrbyte(*sp++, bw);
			if (skipOdd)
				wrbyte(0, bw);
		}
		if (n <= 0)			/* was that it? */
			break;
		val = *sp;			/* output run */
		for (cnt = 1; --n && cnt < 255; cnt++)
			if (*++sp != val)
				break;
		wrbyte(cnt, bw);
		wrbyte(val, bw);
	}
	bw->yscan++;				/* write line break or EOD */
	if (bw->yscan == bw->hdr->height) {
		wrbyte(0, bw); wrbyte(1, bw);   /* end of bitmap marker */
		if (wrseek(2, bw) != BIR_OK)
			return BIR_OK;		/* no one may care */
		wrint32(bw->flen, bw);		/* correct file length */
		if (wrseek(34, bw) != BIR_OK)
			return BIR_OK;
		wrint32(bw->flen-bw->fbmp, bw); /* correct bitmap length */
	} else {
		wrbyte(0, bw); wrbyte(0, bw);   /* end of line marker */
	}
	return BIR_OK;
}

/* free BMP writer resources */
void
BMPfreeWriter(BMPWriter *bw)
{
	if (bw == NULL)
		return;
	free((void *)bw->hdr);
	free((void *)bw->scanline);
	free((void *)bw);
}
