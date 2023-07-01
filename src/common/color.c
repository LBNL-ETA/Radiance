#ifndef lint
static const char	RCSid[] = "$Id: color.c,v 2.26 2023/07/01 01:31:17 greg Exp $";
#endif
/*
 *  color.c - routines for color calculations.
 *
 *  Externals declared in color.h
 */

#include "copyright.h"

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "color.h"

#ifdef getc_unlocked		/* avoid horrendous overhead of flockfile */
#undef getc
#undef putc
#undef ferror
#define getc    getc_unlocked
#define putc    putc_unlocked
#define ferror	ferror_unlocked
#endif

#define  MINELEN	8	/* minimum scanline length for encoding */
#define  MAXELEN	0x7fff	/* maximum scanline length for encoding */
#define  MINRUN		4	/* minimum run length */


void *
tempbuffer(			/* get a temporary buffer */
	unsigned int  len
)
{
	static void		*tempbuf = NULL;
	static unsigned int	tempbuflen = 0;

	if (!len) {		/* call to free */
		if (tempbuflen) {
			free(tempbuf);
			tempbuf = NULL;
			tempbuflen = 0;
		}
		return(NULL);
	}
	if (len <= tempbuflen)	/* big enough already? */
		return(tempbuf);
				/* else free & reallocate */
	if (tempbuflen)
		free(tempbuf);
	tempbuf = malloc(len);
	tempbuflen = len*(tempbuf != NULL);
	return(tempbuf);
}


int
fwritecolrs(			/* write out a colr scanline */
	COLR  *scanline,
	int  len,
	FILE  *fp
)
{
	int  i, j, beg, cnt = 1;
	int  c2;
	
	if ((len < MINELEN) | (len > MAXELEN))	/* OOBs, write out flat */
		return(fwrite((char *)scanline,sizeof(COLR),len,fp) - len);
					/* put magic header */
	putc(2, fp);
	putc(2, fp);
	putc(len>>8, fp);
	putc(len&0xff, fp);
					/* put components seperately */
	for (i = 0; i < 4; i++) {
	    for (j = 0; j < len; j += cnt) {	/* find next run */
		for (beg = j; beg < len; beg += cnt) {
		    for (cnt = 1; (cnt < 127) & (beg+cnt < len) &&
			    scanline[beg+cnt][i] == scanline[beg][i]; cnt++)
			;
		    if (cnt >= MINRUN)
			break;			/* long enough */
		}
		if ((beg-j > 1) & (beg-j < MINRUN)) {
		    c2 = j+1;
		    while (scanline[c2++][i] == scanline[j][i])
			if (c2 == beg) {	/* short run */
			    putc(128+beg-j, fp);
			    putc(scanline[j][i], fp);
			    j = beg;
			    break;
			}
		}
		while (j < beg) {		/* write out non-run */
		    if ((c2 = beg-j) > 128) c2 = 128;
		    putc(c2, fp);
		    while (c2--)
			putc(scanline[j++][i], fp);
		}
		if (cnt >= MINRUN) {		/* write out run */
		    putc(128+cnt, fp);
		    putc(scanline[beg][i], fp);
		} else
		    cnt = 0;
	    }
	}
	return(ferror(fp) ? -1 : 0);
}

/*
 * An old-format scanline is either a stream of valid RGBE or XYZE real
 * pixels or at least one real pixel followed by some number of
 * invalid real pixels of the form (1,1,1,n), where n is a count.
 * These can themselves be repeated to create a multibyte repeat
 * count, with the least significant byte first (little-endian order.)
 * Repeat counts are limited by the size of an int; if a repetition
 * leads to an overrun, the rest of the the repetition will be
 * silently ignored.
 */
static int
oldreadcolrs(			/* read in an old-style colr scanline */
	COLR  *scanline,
	int  len,
	FILE  *fp
)
{
	int  rshift = 0;
	int  i;
	
	while (len > 0) {
		scanline[0][RED] = getc(fp);
		scanline[0][GRN] = getc(fp);
		scanline[0][BLU] = getc(fp);
		scanline[0][EXP] = i = getc(fp);
		if (i == EOF)
			return(-1);
		if (scanline[0][GRN] == 1 &&
				(scanline[0][RED] == 1) &
				(scanline[0][BLU] == 1)) {
			i = scanline[0][EXP] << rshift;
			while (i--) {
				copycolr(scanline[0], scanline[-1]);
				if (--len <= 0)
					return(0);
				scanline++;
			}
			rshift += 8;
		} else {
			scanline++;
			len--;
			rshift = 0;
		}
	}
	return(0);
}

/* 
 * There are two scanline formats: old and new.  The old format
 * compresses runs of RGBE or XYZE four-byte real pixels; the new
 * format breaks the pixels into R, G, B, and E lines (or XYZE lines)
 * which are individually run-length encoded.
 *
 * An old-format scanline always begins with a valid real pixel; at
 * least one of the RGB (or XYZ) values will have its high-order bit
 * set.  A new-format scanline begins with four bytes which are not a
 * valid real pixel: (2, 2, lenhigh, lenlow) where lenhigh is always
 * less than 128 and hence never has a high-order bit set.
 *
 * A new-format scanline is broken into its RGBE or XYZE components.
 * Each is output and run-length encoded separately so that a scanline
 * is broken into four records.  In turn, each record is organized
 * into chunks of up to 128 characters, which begin with a count byte.
 * If the count byte is greater than 128, the following data byte is
 * repeated (count-128) times.  If not, the count byte is followed by
 * that many data bytes.
 */
int
freadcolrs(			/* read in an encoded colr scanline */
	COLR  *scanline,
	int  len,
	FILE  *fp
)
{
	int  i, j;
	int  code, val;
					/* determine scanline type */
	if (len <= 0)
		return(0);
	if ((i = getc(fp)) == EOF)
		return(-1);
	scanline[0][RED] = i;
	scanline[0][GRN] = getc(fp);
	scanline[0][BLU] = getc(fp);
	if ((i = getc(fp)) == EOF)
		return(-1);
	if ((scanline[0][RED] != 2) | (scanline[0][GRN] != 2) |
			(scanline[0][BLU] & 0x80)) {
		scanline[0][EXP] = i;
		return(oldreadcolrs(scanline+1, len-1, fp));
	}
	if ((scanline[0][BLU]<<8 | i) != len)
		return(-1);		/* length mismatch! */
					/* read each component */
	for (i = 0; i < 4; i++)
	    for (j = 0; j < len; ) {
		if ((code = getc(fp)) == EOF)
		    return(-1);
		if (code > 128) {	/* run */
		    code &= 127;
		    if ((val = getc(fp)) == EOF)
			return -1;
		    if (j + code > len)
		    	return -1;	/* overrun */
		    while (code--)
			scanline[j++][i] = val;
		} else {		/* non-run */
		    if (j + code > len)
		    	return -1;	/* overrun */
		    while (code--) {
			if ((val = getc(fp)) == EOF)
			    return -1;
			scanline[j++][i] = val;
		    }
		}
	    }
	return(0);
}


int
fwritescan(			/* write out a scanline */
	COLOR  *scanline,
	int  len,
	FILE  *fp
)
{
	COLR  *clrscan;
	int  n;
	COLR  *sp;
					/* get scanline buffer */
	if ((sp = (COLR *)tempbuffer(len*sizeof(COLR))) == NULL)
		return(-1);
	clrscan = sp;
					/* convert scanline */
	n = len;
	while (n-- > 0) {
		setcolr(sp[0], scanline[0][RED],
				  scanline[0][GRN],
				  scanline[0][BLU]);
		scanline++;
		sp++;
	}
	return(fwritecolrs(clrscan, len, fp));
}


int
freadscan(			/* read in a scanline */
	COLOR  *scanline,
	int  len,
	FILE  *fp
)
{
	COLR  *clrscan;

	if ((clrscan = (COLR *)tempbuffer(len*sizeof(COLR))) == NULL)
		return(-1);
	if (freadcolrs(clrscan, len, fp) < 0)
		return(-1);
					/* convert scanline */
	colr_color(scanline[0], clrscan[0]);
	while (--len > 0) {
		scanline++; clrscan++;
		if (clrscan[0][GRN] == clrscan[-1][GRN] &&
			    (clrscan[0][RED] == clrscan[-1][RED]) &
			    (clrscan[0][BLU] == clrscan[-1][BLU]) &
			    (clrscan[0][EXP] == clrscan[-1][EXP]))
			copycolor(scanline[0], scanline[-1]);
		else
			colr_color(scanline[0], clrscan[0]);
	}
	return(0);
}


void
setcolr(			/* assign a short color value */
	COLR  clr,
	double  r,
	double  g,
	double  b
)
{
	double  d;
	int  e;
	
	d = r > g ? r : g;
	if (b > d) d = b;

	if (d <= 1e-32) {
		clr[RED] = clr[GRN] = clr[BLU] = 0;
		clr[EXP] = 0;
		return;
	}

	d = frexp(d, &e) * 256.0 / d;

	if (r > 0.0)
		clr[RED] = r * d;
	else
		clr[RED] = 0;
	if (g > 0.0)
		clr[GRN] = g * d;
	else
		clr[GRN] = 0;
	if (b > 0.0)
		clr[BLU] = b * d;
	else
		clr[BLU] = 0;

	clr[EXP] = e + COLXS;
}


void
colr_color(			/* convert short to float color */
	COLOR  col,
	COLR  clr
)
{
	double  f;
	
	if (clr[EXP] == 0)
		col[RED] = col[GRN] = col[BLU] = 0.0;
	else {
		f = ldexp(1.0, (int)clr[EXP]-(COLXS+8));
		col[RED] = (clr[RED] + 0.5)*f;
		col[GRN] = (clr[GRN] + 0.5)*f;
		col[BLU] = (clr[BLU] + 0.5)*f;
	}
}


int
bigdiff(				/* c1 delta c2 > md? */
	COLOR  c1,
	COLOR  c2,
	double  md
)
{
	int  i;

	for (i = 0; i < 3; i++)
		if (colval(c1,i)-colval(c2,i) > md*colval(c2,i) ||
				colval(c2,i)-colval(c1,i) > md*colval(c1,i))
			return(1);
	return(0);
}
