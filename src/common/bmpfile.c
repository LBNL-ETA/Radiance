#ifndef lint
static const char RCSid[] = "$Id: bmpfile.c,v 2.2 2004/03/26 21:29:19 schorsch Exp $";
#endif
/*
 *  Windows and OS/2 BMP file support
 */

/* XXX reading and writing of compressed files currently unsupported */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmpfile.h"

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
	case BIR_SEEKERR:
		return "BMP seek error";
	}
	return "Unknown BMP error";
}

/* compute uncompressed scanline size */
static uint32
getScanSiz(BMPHeader *hdr)
{
	uint32  scanSiz;
						/* bytes per scanline */
	scanSiz = (hdr->bpp*hdr->width + 7) >> 3;
	if (scanSiz & 3)			/* 32-bit word boundary */
		scanSiz += 4 - (scanSiz & 3);

	return scanSiz;
}

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
	hdrSiz = rdint32(br) + 14;		/* header size */
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
	if (br->hdr->bpp != 1 && br->hdr->bpp != 4 &&
			br->hdr->bpp != 8 && br->hdr->bpp != 24 &&
			br->hdr->bpp != 32)
		goto err;
	br->hdr->compr = rdint32(br);		/* compression mode */
	if (br->hdr->compr != BI_UNCOMPR && br->hdr->compr != BI_RLE8)
		goto err;			/* don't support the rest */
	(void)rdint32(br);			/* bitmap size */
	br->hdr->hRes = rdint32(br);		/* horizontal resolution */
	br->hdr->vRes = rdint32(br);		/* vertical resolution */
	br->hdr->nColors = rdint32(br);		/* # colors used */
	br->hdr->impColors = rdint32(br);       /* # important colors */
	if (br->hdr->impColors < 0)
		goto err;			/* catch premature EOF */
	palLen = 0;
	if (br->hdr->bpp <= 8) {
		palLen = 1 << br->hdr->bpp;
		if (br->hdr->nColors <= 0)
			br->hdr->nColors = palLen;
		else if (br->hdr->nColors > palLen)
			goto err;
		if (br->hdr->impColors <= 0)
			br->hdr->impColors = br->hdr->nColors;
		else if (br->hdr->impColors > br->hdr->nColors)
			goto err;
	} else if (br->hdr->nColors > 0 || br->hdr->compr != BI_UNCOMPR)
		goto err;
						/* extend header */
	if (bmPos < hdrSiz + sizeof(RGBquad)*palLen)
		goto err;
	br->hdr->infoSiz = bmPos - (hdrSiz + sizeof(RGBquad)*palLen);
	if (palLen > 0 || br->hdr->infoSiz > 0) {
		br->hdr = (BMPHeader *)realloc((void *)br->hdr,
					sizeof(BMPHeader) +
					sizeof(RGBquad)*palLen -
					sizeof(br->hdr->palette) +
					br->hdr->infoSiz);
		if (br->hdr == NULL)
			goto err;
	}
						/* read color palette */
	if (rdbytes((char *)br->hdr->palette,
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
	if (seek != NULL && br->hdr->compr != BI_UNCOMPR) {
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
	if (BMPreadScanline(br) != BIR_OK)
		goto err;
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
	for (rgbp = hdr->palette, n = hdr->impColors; n-- > 0; rgbp++)
		if (((rgbp->r != rgbp->g) | (rgbp->g != rgbp->b)))
			return 0;
	return 1;			/* all colors neutral in map */
}

/* read and decode next BMP scanline */
int
BMPreadScanline(BMPReader *br)
{
	if (br->yscan + 1 >= br->hdr->height)
		return BIR_EOF;
	br->yscan++;			/* prepare to read */
	if (br->hdr->compr == BI_UNCOMPR)
		return rdbytes((char *)br->scanline, getScanSiz(br->hdr), br);
/* XXX need to perform actual decompression */
	return BIR_UNSUPPORTED;
	if (br->seek != NULL)
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
		if (br->hdr->compr == BI_UNCOMPR) {
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
BMPdecodePixel(int i, BMPReader *br)
{
	static const RGBquad    black = {0, 0, 0, 0};
	
	if (((br == NULL) | (i < 0)) || i >= br->hdr->width)
		return black;

	switch (br->hdr->bpp) {
	case 24: {
		uint8   *bgr = br->scanline + 3*i;
		RGBquad cval;
		cval.b = bgr[0]; cval.g = bgr[1]; cval.r = bgr[2];
		cval.padding = 0;
		return cval;
		}
	case 32:
		return ((RGBquad *)br->scanline)[i];
	case 8:
		return br->hdr->palette[br->scanline[i]];
	case 1:
		return br->hdr->palette[br->scanline[i>>3]>>(i&7) & 1];
	case 4:
		return br->hdr->palette[br->scanline[i>>1]>>((i&1)<<2) & 0xf];
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
	hdr->compr = BI_UNCOMPR;
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
				++((bw)->fpos) > (bw)->flen ? \
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
	wrbyte(i>>24  & 0xff, bw);
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
	if (bw->seek == NULL) {
		if (pos < bw->fpos)
			return BIR_SEEKERR;
		while (bw->fpos < pos)
			wrbyte(0, bw);
		return BIR_OK;
	}
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
	if (cput == NULL || hdr == NULL)
		return NULL;
	if (hdr->width <= 0 || hdr->height <= 0)
		return NULL;
	if (hdr->compr != BI_UNCOMPR && hdr->compr != BI_RLE8)
		return NULL;			/* unsupported */
	if (hdr->bpp > 8 && hdr->compr != BI_UNCOMPR)
		return NULL;
	palSiz = BMPpalLen(hdr);
	if (hdr->nColors > palSiz)
		return NULL;
	if (hdr->compr != BI_UNCOMPR && seek == NULL)
		return NULL;
						/* compute sizes */
	hdrSiz = 2 + 6*4 + 2*2 + 6*4;
	palSiz *= sizeof(RGBquad);
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
				
/* write the current scanline */
int
BMPwriteScanline(BMPWriter *bw)
{
	if (bw->yscan >= bw->hdr->height)
		return BIR_EOF;
						/* writing uncompressed? */
	if (bw->hdr->compr == BI_UNCOMPR) {
		uint32  scanSiz = getScanSiz(bw->hdr);
		uint32  slpos = bw->fbmp + bw->yscan*scanSiz;
		if (wrseek(slpos, bw) != BIR_OK)
			return BIR_SEEKERR;
		wrbytes((char *)bw->scanline, scanSiz, bw);
		bw->yscan++;
		return BIR_OK;
	}
						/* else write compressed */
/* XXX need to do actual compression */
	return BIR_UNSUPPORTED;
						/* write file length at end */
	if (bw->yscan == bw->hdr->height) {
		if (bw->seek == NULL || (*bw->seek)(2, bw->c_data) != 0)
			return BIR_SEEKERR;
		bw->fpos = 2;
		wrint32(bw->flen, bw);
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
