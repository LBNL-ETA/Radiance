#ifndef lint
static const char	RCSid[] = "$Id: ra_xim.c,v 3.2 2003/06/30 14:59:12 schorsch Exp $";
#endif
/***********************************************************************
*   Author: Philip Thompson (phils@athena.mit.edu)
*   - convert from RADIANCE image file format to Xim and vice-versa.
*  Copyright (c) 1990, Massachusetts Institute of Technology
*   			 Philip R. Thompson (phils@athena.mit.edu)
*                Computer Resource Laboratory (CRL)
*                Dept. of Architecture and Planning
*                M.I.T., Rm 9-526
*                Cambridge, MA  02139
*   This  software and its documentation may be used, copied, modified,
*   and distributed for any purpose without fee, provided:
*       --  The above copyright notice appears in all copies.
*       --  This disclaimer appears in all source code copies.
*       --  The names of M.I.T. and the CRL are not used in advertising
*           or publicity pertaining to distribution of the software
*           without prior specific written permission from me or CRL.
*   I provide this software freely as a public service.  It is NOT
*   public domain nor a commercial product, and therefore is not subject
*   to an an implied warranty of merchantability or fitness for a
*   particular purpose.  I provide it as is, without warranty.
*   This software is furnished  only on the basis that any party who
*   receives it indemnifies and holds harmless the parties who furnish
*   it against any claims, demands, or liabilities connected with using
*   it, furnishing it to others, or providing it to a third party.
***********************************************************************/
#if (!defined(lint) && !defined(SABER))
static char  radiance_rcsid[] =
    "$Header: /home/cvsd/radiance/ray/src/px/Attic/ra_xim.c,v 3.2 2003/06/30 14:59:12 schorsch Exp $";
#endif

#include <stdio.h>
#include <math.h>
#include <pwd.h>
#include <sys/types.h>
#include <string.h>
#include <sys/time.h>

/* XimHeader.h */

/*
*  Author: Philip R. Thompson
*  Address:  phils@athena.mit.edu, 9-526 
*  $Header: /home/cvsd/radiance/ray/src/px/Attic/ra_xim.c,v 3.2 2003/06/30 14:59:12 schorsch Exp $
*  $Date: 2003/06/30 14:59:12 $
*  $Source: /home/cvsd/radiance/ray/src/px/Attic/ra_xim.c,v $
*  Copyright (c) 1989, 1990 Massachusetts Institute of Technology
*                1988, Philip R. Thompson (phils@athena.mit.edu)
*                Computer Resource Laboratory (CRL)
*                Dept. of Architecture and Planning
*                M.I.T., Rm 9-526
*                Cambridge, MA  02139
*   This  software and its documentation may be used, copied, modified,
*   and distributed for any purpose without fee, provided:
*       --  The above copyright notice appears in all copies.
*       --  This disclaimer appears in all source code copies.
*       --  The names of M.I.T. and the CRL are not used in advertising
*           or publicity pertaining to distribution of the software
*           without prior specific written permission from me or CRL.
*   I provide this software freely as a public service.  It is NOT
*   public domain nor a commercial product, and therefore is not subject
*   to an an implied warranty of merchantability or fitness for a
*   particular purpose.  I provide it as is, without warranty.
*
*   This software is furnished  only on the basis that any party who
*   receives it indemnifies and holds harmless the parties who furnish
*   it against any claims, demands, or liabilities connected with using
*   it, furnishing it to others, or providing it to a third party.
*/

#define IMAGE_VERSION    3
#ifndef _BYTE
typedef unsigned char  byte;
#define _BYTE  1
#endif
#include <sys/types.h>


/* External ascii file format
*  Note:  size of header should be 1024 (1K) bytes.
*/
typedef struct ImageHeader {
    char file_version[8];   /* header version */
    char header_size[8];    /* Size of file header in bytes (1024) */
    char image_width[8];    /* Width of the raster image */
    char image_height[8];   /* Height of the raster imgage */
    char num_colors[4];     /* Actual number of entries in c_map */
    char greyscale[2];      /* Flag if no c_map needed */
	char flag[2];     		/* Not used */
    char num_channels[3];   /* 0 or 1 = pixmap, 3 = RG&B buffers */
    char bytes_per_line[5]; /* bytes per scanline */
    char num_pictures[4];   /* Number of pictures in file */
    char bits_per_channel[4]; /* usually 1 or 8 */
    char alpha_channel[4];  /* Alpha channel flag */
    char runlength[4];      /* Runlength encoded flag */
    char author[48];        /* Name of who made it */
    char date[32];          /* Date and time image was made */
    char program[16];       /* Program that created this file */
    char comment[96];       /* other viewing info. for this image */
    unsigned char c_map[256][3]; /* RGB values of the pixmap indices */
} XimAsciiHeader;


/* Internal binary format
*/
typedef struct Color {
    byte red, grn, blu, pixel;
} Color;

typedef struct XimImage {
    int  width;             /* width of the image in pixels */
    int  height;            /* height of the image in pixels */
    int  bytes_line;        /* bytes to hold one scanline */
    unsigned  datasize;     /* size of one channel of data */
    short  nchannels;       /* data channels per pixel, usually 1 or3 */
    short  bits_channel;    /* depth of a channel, usually 1 or 8 */
    byte  *data;            /* pixmap or red channel data */
    byte  *grn_data;        /* green channel data */
    byte  *blu_data;        /* blue  channel data */
    byte  *other;           /* other (alpha) data */
    short  ncolors;         /* number of colors in the color table */
    Color  *colors;         /* colortable for pixmaps */
    short  curr_pic, npics; /* current & total file number of images */
    unsigned alpha_flag: 1; /* alpha channel flag - "other" is used */
    unsigned grey_flag : 1; /* is pixmap grey scaled */
    unsigned packed_flag:1; /* is data packed in one chunk of memory */
    unsigned runlen_flag:1; /* is data compressed flag (not used) */
    unsigned : 0;           /* future flags space and word alignment */
    char   *author;         /* author credit, copyright, etc */
    char   *date;           /* date image was made, grabbed, etc. */
    char   *program;        /* program used to make this */
    char   *name;           /* name of the image - usually filename */
    short  ncomments;       /* number of comments strings */
    char   **comments;      /* pointers to null terminated comments */
    char   *offset;         /* original offset in memory when packed */
    float  gamcor;           /* image storage gamma */
    float  chroma_red[2];   /* x, y image chromacity coords */
    float  chroma_grn[2];
    float  chroma_blu[2];
    float  chroma_wht[2];
} XimImage;


/* Future external ascii variable length header - under review
*/
#if (VERSION == 4)
typedef struct XimAsciiHeader {
    char file_version[4];   /* header version */
    char header_size[8];    /* Size of file header (fixed part only) */
    char image_height[8];   /* Height of the raster imgage in pixels */
    char image_width[8];    /* Width of the raster image in pixels */
    char bytes_line[8];     /* Actual # of bytes separating scanlines */
    char bits_channel[4];   /* Bits per channel (usually 1 or 8) */
    char num_channels[4];   /* 1 = pixmap, 3 = RG&B buffers */
    char alpha_channel[2];  /* Alpha channel flag */
    char num_colors[4];     /* Number of entries in c_map (if any) */
    char num_pictures[4];   /* Number of images in file */
    char compress_flag[2];  /* Runlength encoded flag */
    char greyscale_flag[2]; /* Greyscale flag */
    char future_flags[4];
    char author[48];        /* Name of who made it, from passwd entry */
    char date[32];          /* Unix format date */
    char program[32];       /* Program that created this */
    char gamcor[12];         /* image storage gamma */
    char chroma_red[24];    /* image red primary chromaticity coords. */
    char chroma_grn[24];    /*   "   green "          "         "     */
    char chroma_blu[24];    /*   "   blue  "          "         "     */
    char chroma_wht[24];    /*   "   white point      "         "     */
    char comment_length[8]  /* Total length of all comments */
    /* char *comment;           Null separated comments  */
    /* byte c_map[];   RGB Colortable, (ncolors * 3 bytes) */
} XimAsciiHeader;
#endif VERSION 4

#ifndef rnd
#define rnd(x)  ((int)((float)(x) + 0.5)) /* round a float to an int */
#endif

/* Note:
* - All external data is in char's in order to maintain easily
*   portability across machines and some human readibility.
* - Images may be stored as bitmaps (8 pixels/byte), pixmaps
*   (1 byte /pixel) or as seperate red, green, blue channel data
*   (3 bytes/pixel).
* - An alpha channel is optional and is found after every num_channels
*   of data.
* - Pixmaps or RG&B (and alpha) channel data are stored respectively
*   after the header.
* - If num_channels == 1 a pixmap is assumed and the colormap in the
*   header is used unless its greyscale.
* - Datasize is for one channel.
*/

/*** end ImageHeader.h ***/

#define XDECR   1
#define YDECR   2
#define YMAJOR  4
#define RED     0
#define GRN     1
#define BLU     2
#define EXP     3
#define MAXLINE 512

#define setcolor(col,r,g,b) (col[RED]=(r), col[GRN]=(g), col[BLU]=(b))
#define  copycolr(c1,c2)    (c1[0]=c2[0],c1[1]=c2[1], \
    	c1[2]=c2[2],c1[3]=c2[3])
#define  copycolor(c1,c2)   ((c1)[0]=(c2)[0],(c1)[1]=(c2)[1],\
		(c1)[2]=(c2)[2])
#define COLXS      128      /* excess used for exponent */
typedef float COLOR[3];     /* red, green, blue */
typedef byte  COLR[4];      /* red, green, blue, exponent */

#define	FREE(a) if ((a) != NULL) free((char*)(a)), (a) = NULL

char *AllocAndCopy();
extern debug_flag;
static byte RadiancePixel();


XimReadRadiance(ifp, xim, gamcor)
FILE *ifp;
XimImage *xim;
double gamcor;
{
    register int x, y, i;
    int width, height;
    COLOR *scanline;		/* RADIANCE variable */
    struct timeval  time_val;
    struct timezone  time_zone;
	struct passwd  *getpwnam(), *getpwuid(), *pwd = NULL;
	char  *login_ptr = NULL, *ctime(), *getlogin();
	byte gmap[256];

	/* get width, height from radiance file */
    memset(xim, '\0', sizeof(XimImage));
    if (!RadianceGetSize(ifp, &width, &height))
		return 0;
	xim->width = xim->bytes_line = width;
	xim->height = height;
	xim->nchannels = 3;
	xim->bits_channel = 8;
	xim->datasize = xim->height * xim->bytes_line;
	xim->npics = 1;
    xim->data = (byte *)malloc(xim->datasize);
    xim->grn_data = (byte *)malloc(xim->datasize);
    xim->blu_data = (byte *)malloc(xim->datasize);
	if (!xim->data || !xim->grn_data || !xim->blu_data) {
		FREE(xim->data);
		FREE(xim->grn_data);
		FREE(xim->blu_data);
        XimWarn("XimReadRadiance: can't malloc xim buffers");
		return 0;
	}
    if (!(scanline = (COLOR*)malloc(width * sizeof(COLOR)))) {
        XimWarn("XimReadRadiance: can't malloc scanline");
        return(0);
    }
	for (i=0; i < 256; i++) {
		gmap[i] = (byte)(0.5 + 255.0 * pow(i / 255.0, 1.0/gamcor));
		if (debug_flag)
			fprintf(stderr, "%2d: %u\n", i, gmap[i]);
	}
    for (y=0, i=0; y < height; y++) {
    	if (!freadscan(scanline, width, ifp)) {
        	XimWarn("freadscan: read error");
			free((char*)scanline);
        	return 0;
    	}
        for (x=0; x < width; x++, i++) {
        	xim->data[i] = gmap[RadiancePixel(scanline[x], RED)];
            xim->grn_data[i] = gmap[RadiancePixel(scanline[x], GRN)];
            xim->blu_data[i] = gmap[RadiancePixel(scanline[x], BLU)];
        }
    }
	free((char*)scanline);
	xim->comments = (char**)calloc(2, sizeof(char*));
    xim->comments[0] =AllocAndCopy("converted from RADIANCE format",40);
	xim->ncomments = 1;
	(void)gettimeofday(&time_val, &time_zone);
	xim->date = ctime(&time_val.tv_sec);
	xim->date[strlen(xim->date)-1] = '\0';
	/* date points to static memory, so copy it to make freeable */
	xim->date = AllocAndCopy(xim->date, strlen(xim->date)+1);
    if ((login_ptr = getlogin()) != NULL)
        pwd = getpwnam(login_ptr);
    else
        pwd = getpwuid(getuid());
    if (pwd && !(xim->author = AllocAndCopy(pwd->pw_gecos,
            strlen(pwd->pw_gecos)+1))) {
        XimWarn("XimReadXwdFile: Can't AllocAndCopy author's name");
        return 0;
    } else
        xim->author = NULL;
	return 1;
}


XimWriteRadiance(ofp, xim, gamcor)
FILE *ofp;
XimImage *xim;
double gamcor;
{
    COLOR *scanline;
    int width, height;
    register int x, y;
	register byte *red, *grn, *blu;
	float gmap[256];

	/* get width, height from xim */
    width = xim->width;
    height = xim->height;

	/* allocate scan line space */
    if (!(scanline = (COLOR *)malloc(width * sizeof(COLOR)))) {
        XimWarn("XimWriteRadiance: cant' alloc scanline");
        return(0);
    }
	/* write out the RADIANCE header */
    radiance_writeheader(ofp, width, height);
	red = xim->data;
	grn = xim->grn_data;
	blu = xim->blu_data;
	for (x=0; x < 256; x++) {
		gmap[x] = (float)pow((x / 255.0), gamcor);
		if (debug_flag)
			fprintf(stderr, "%2d: %f\n", x, gmap[x]);
	}
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++)
            setcolor(scanline[x], gmap[red[x]], gmap[grn[x]], gmap[blu[x]]);
        radiance_writescanline(ofp, scanline, width);
		red += xim->bytes_line;
		grn += xim->bytes_line;
		blu += xim->bytes_line;
    }
	return 1;
}

static RadianceGetSize(fp, w, h)
FILE *fp;
int *w, *h;
{
    if (!GetHeader(fp, NULL))
		return 0;
    if (fgetresolu(w, h, fp) != (YMAJOR | YDECR)) {
        fprintf(stderr, "bad RADIANCE format\n");
        return 0;
    }
    return 1;
}

static byte RadiancePixel(pixel, chan)
COLOR pixel;
int chan;
{
    float value;

    value = pixel[chan];
    if (value <= 0.0)
        return (byte)0;
    else if (value >= 1.0)
        return (byte)255;
    return (byte)((int)((value * 256.0)));
}

static radiance_writeheader(fp, x, y)
FILE *fp;
int x;
int y;
{
    fputc('\n', fp);
    fputresolu(YMAJOR | YDECR, x, y, fp);
    fflush(fp);
}


radiance_writescanline(fp, buf, x)
FILE *fp;
COLOR *buf;
int x;
{
    if (fwritescan(buf, x, fp) < 0) {
        fprintf(stderr, "write error?\n");
        perror("fwrite");
        exit(1);
    }
}


freadscan(scanline, len, fp)           /* read in a scanline */
register COLOR *scanline;
int len;
register FILE *fp;
{
    COLR thiscolr;
    int rshift;
    register int i;

    rshift = 0;
    while (len > 0) {
        thiscolr[RED] = getc(fp);
        thiscolr[GRN] = getc(fp);
        thiscolr[BLU] = getc(fp);
        thiscolr[EXP] = getc(fp);
        if (feof(fp) || ferror(fp))
            return 0;
        if (thiscolr[RED] == 1 &&
            thiscolr[GRN] == 1 &&
            thiscolr[BLU] == 1) {
            for (i = thiscolr[EXP] << rshift; i > 0; i--) {
                copycolor(scanline[0], scanline[-1]);
                scanline++;
                len--;
            }
            rshift += 8;
        } else {
            colr_color(scanline[0], thiscolr);
            scanline++;
            len--;
            rshift = 0;
        }
    }
    return 1;
}


colr_color(col, clr)               /* convert short to float color */
register COLOR col;
register COLR clr;
{
    double f;

    if (clr[EXP] == 0)
        col[RED] = col[GRN] = col[BLU] = 0.0;
    else {
        f = ldexp(1.0, (int) clr[EXP] - (COLXS + 8));
        col[RED] = (clr[RED] + 0.5) * f;
        col[GRN] = (clr[GRN] + 0.5) * f;
        col[BLU] = (clr[BLU] + 0.5) * f;
    }
}

static GetHeader(fp, f)               /* get header from file */
FILE *fp;
int (*f)();
{
    char buf[MAXLINE];

    for (;;) {
        buf[MAXLINE - 2] = '\n';
        if (fgets(buf, sizeof(buf), fp) == NULL)
            return (-1);
        if (buf[0] == '\n')
            return 1;
        if (buf[MAXLINE - 2] != '\n') {
            ungetc(buf[MAXLINE - 2], fp);   /* prevent false end */
            buf[MAXLINE - 2] = '\0';
        }
        if (f != NULL)
            (*f)(buf);
    }
}


static fgetresolu(xrp, yrp, fp)           /* get x and y resolution */
int *xrp, *yrp;
FILE *fp;
{
    char buf[64], *xndx, *yndx;
    register char *cp;
    register int ord;

    if (fgets(buf, sizeof(buf), fp) == NULL)
        return (-1);
    xndx = yndx = NULL;
    for (cp = buf + 1; *cp; cp++)
        if (*cp == 'X')
            xndx = cp;
        else if (*cp == 'Y')
            yndx = cp;
    if (xndx == NULL || yndx == NULL)
        return (-1);
    ord = 0;
    if (xndx > yndx)
        ord |= YMAJOR;
    if (xndx[-1] == '-')
        ord |= XDECR;
    if (yndx[-1] == '-')
        ord |= YDECR;
    if ((*xrp = atoi(xndx + 1)) <= 0)
        return (-1);
    if ((*yrp = atoi(yndx + 1)) <= 0)
        return (-1);
    return(ord);
}


fputresolu(ord, xres, yres, fp)        /* put x and y resolution */
register int ord;
int xres, yres;
FILE *fp;
{
    if (ord & YMAJOR)
        fprintf(fp, "%cY %d %cX %d\n", ord & YDECR ? '-' : '+', yres,
            ord & XDECR ? '-' : '+', xres);
    else
        fprintf(fp, "%cX %d %cY %d\n", ord & XDECR ? '-' : '+', xres,
            ord & YDECR ? '-' : '+', yres);
}


fwritescan(scanline, len, fp)          /* write out a scanline */
register COLOR *scanline;
int len;
register FILE *fp;
{
    int rept;
    COLR lastcolr, thiscolr;

    lastcolr[RED] = lastcolr[GRN] = lastcolr[BLU] = 1;
    lastcolr[EXP] = 0;
    rept = 0;
    while (len > 0) {
        setcolr(thiscolr, scanline[0][RED],
            scanline[0][GRN],
            scanline[0][BLU]);
        if (thiscolr[EXP] == lastcolr[EXP] &&
            thiscolr[RED] == lastcolr[RED] &&
            thiscolr[GRN] == lastcolr[GRN] &&
            thiscolr[BLU] == lastcolr[BLU])
            rept++;
        else {
            while (rept) {  /* write out count */
                putc(1, fp);
                putc(1, fp);
                putc(1, fp);
                putc(rept & 255, fp);
                rept >>= 8;
            }
            putc(thiscolr[RED], fp);    /* new color */
            putc(thiscolr[GRN], fp);
            putc(thiscolr[BLU], fp);
            putc(thiscolr[EXP], fp);
            copycolr(lastcolr, thiscolr);
            rept = 0;
        }
        scanline++;
        len--;
    }
    while (rept) {  /* write out count */
        putc(1, fp);
        putc(1, fp);
        putc(1, fp);
        putc(rept & 255, fp);
        rept >>= 8;
    }
    return (ferror(fp) ? -1 : 0);
}

static setcolr(clr, r, g, b)          /* assign a short color value */
register COLR clr;
double r, g, b;
{
    double d, frexp();
    int e;

    d = r > g ? r : g;
    if (b > d)
        d = b;
    if (d <= 1e-32) {
        clr[RED] = clr[GRN] = clr[BLU] = 0;
        clr[EXP] = 0;
        return;
    }
    d = frexp(d, &e) * 256.0 / d;

    clr[RED] = (byte)(r * d);
    clr[GRN] = (byte)(g * d);
    clr[BLU] = (byte)(b * d);
    clr[EXP] = (byte)(e + COLXS);
}

/*** radiance_io.c ***/
