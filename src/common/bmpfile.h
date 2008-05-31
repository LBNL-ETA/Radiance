/* RCSid $Id: bmpfile.h,v 2.6 2008/05/31 19:38:36 greg Exp $ */
/*
 *  Windows and OS/2 BMP file support
 */

#ifndef _RAD_BMPFILE_H_
#define _RAD_BMPFILE_H_

#include "tiff.h"			/* needed for int32, etc. */

#ifdef __cplusplus
extern "C" {
#endif

					/* compression modes */
#define BI_UNCOMPR		0		/* no compression */
#define BI_RLE8			1		/* runlength on bytes */
#define BI_RLE4			2		/* runlength on half-bytes */
#define BI_BITFIELDS		3		/* uncompressed w/ bitfields */

/* A Color table entry */
typedef struct {
	uint8		b, g, r;	/* blue, green, & red bytes */
	uint8		padding;	/* padding to 32-bit boundary */
} RGBquad;

/* Allocated BMP header data */
typedef struct {
	/* the following fields may be altered before the open call */
	int		yIsDown;	/* scanlines proceed downward? */
	int32		hRes;		/* horizontal resolution pixels/meter */
	int32		vRes;		/* vertical resolution pixels/meter */
	int		nColors;	/* total color palette size */
	int		impColors;      /* number of colors actually used */
	int		compr;		/* compression */
	int32		width;		/* bitmap width (pixels) */
	int32		height;		/* bitmap height (pixels) */
	/* the following fields must not be altered after allocation */
	int		bpp;		/* bits per sample (1,4,8,16,24,32) */
	uint32		infoSiz;	/* info buffer size (bytes) */
	/* but the color table should be filled by writer before open call */
	RGBquad		palette[3];     /* color palette (extends struct) */
} BMPHeader;
				
					/* access to bit field triple */
#define BMPbitField(h)  ((uint32 *)(h)->palette)

					/* info buffer access */
#define BMPinfo(h)      ((char *)((h)->palette + (h)->nColors))

					/* function return values */
#define BIR_OK			0		/* all is well */
#define BIR_EOF			(-1)		/* reached end of file */
#define BIR_TRUNCATED		1		/* unexpected EOF */
#define BIR_UNSUPPORTED		2		/* unsupported encoding */
#define BIR_RLERROR		3		/* RLE error */
#define BIR_SEEKERR		4		/* could not seek */

/* A BMP reader structure */
typedef struct BMPReader {
	/* members in this structure should be considered read-only */
	uint8		*scanline;      /* unpacked scanline data */
	int		yscan;		/* last scanline read */
	BMPHeader       *hdr;		/* pointer to allocated header */
	uint32		fpos;		/* current file position */
					/* get character callback */
	int		(*cget)(void *);
					/* absolute seek callback (or NULL) */
	int		(*seek)(uint32, void *);
	void		*c_data;	/* client's opaque data */
	uint32		scanpos[1];     /* recorded scanline position(s) */
} BMPReader;

/* A BMP writer structure */
typedef struct {
	/* the scanline data is filled in by caller before each write */
	uint8		*scanline;      /* caller-prepared scanline data */
	/* modify yscan only if seek is defined & data is uncompressed */
	int		yscan;		/* scanline for next write */
	/* the following fields should not be altered directly */
	BMPHeader       *hdr;		/* allocated header */
	uint32		fbmp;		/* beginning of bitmap data */
	uint32		fpos;		/* current file position */
	uint32		flen;		/* last character written */
					/* put character callback */
	void		(*cput)(int, void *);
					/* absolute seek callback (or NULL) */
	int		(*seek)(uint32, void *);
	void		*c_data;	/* client's opaque data */
} BMPWriter;

					/* open BMP stream for reading */
BMPReader       *BMPopenReader(int (*cget)(void *),
				int (*seek)(uint32, void *), void *c_data);

					/* determine if image is grayscale */
int		BMPisGrayscale(const BMPHeader *hdr);

					/* read next BMP scanline */
int		BMPreadScanline(BMPReader *br);

					/* read a specific scanline */
int		BMPseekScanline(int y, BMPReader *br);

					/* get ith pixel from last scanline */
RGBquad		BMPdecodePixel(int i, const BMPReader *br);

					/* free BMP reader resources */
void		BMPfreeReader(BMPReader *br);

					/* allocate uncompressed RGB header */
BMPHeader       *BMPtruecolorHeader(int xr, int yr, int infolen);

					/* allocate color-mapped header */
BMPHeader       *BMPmappedHeader(int xr, int yr, int infolen, int ncolors);

					/* open BMP stream for writing */
BMPWriter	*BMPopenWriter(void (*cput)(int, void *),
				int (*seek)(uint32, void *), void *c_data,
					BMPHeader *hdr);
				
					/* write the prepared scanline */
int		BMPwriteScanline(BMPWriter *bw);

					/* free BMP writer resources */
void		BMPfreeWriter(BMPWriter *bw);

					/* get corresponding error message */
const char      *BMPerrorMessage(int ec);

					/* stdio callback functions */
int		stdio_getc(void *p);
void		stdio_putc(int c, void *p);
int		stdio_fseek(uint32 pos, void *p);

					/* open stdio input stream */
#define BMPopenInputStream(fp)  BMPopenReader(&stdio_getc, NULL, (void *)fp)

					/* open input file */
#define BMPopenInputFile(fn)    BMPopenReader(&stdio_getc, &stdio_fseek, \
					(void *)fopen(fn, "rb"))

					/* close stdio input file or stream */
#define BMPcloseInput(br)       ( fclose((FILE *)(br)->c_data), \
					BMPfreeReader(br) )

					/* open stdio output stream */
#define BMPopenOutputStream(fp,hdr) \
			BMPopenWriter(&stdio_putc, NULL, (void *)fp, hdr)

					/* open stdio output file */
#define BMPopenOutputFile(fn,hdr) \
			BMPopenWriter(&stdio_putc, &stdio_fseek, \
					(void *)fopen(fn, "wb"), hdr)

					/* close stdio output file or stream */
#define BMPcloseOutput(bw)      ( fclose((FILE *)(bw)->c_data), \
					BMPfreeWriter(bw) )

#ifdef __cplusplus
}
#endif
#endif  /* ! _RAD_BMPFILE_H_ */
