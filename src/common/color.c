#ifndef lint
static const char	RCSid[] = "$Id: color.c,v 2.21 2019/05/27 15:44:26 greg Exp $";
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


char *
tempbuffer(			/* get a temporary buffer */
	unsigned int  len
)
{
	static char  *tempbuf = NULL;
	static unsigned  tempbuflen = 0;

	if (!len | (len > tempbuflen)) {
		if (tempbuflen)
			free(tempbuf);
		tempbuf = len ? (char *)malloc(len) : (char *)NULL;
		tempbuflen = len*(tempbuf != NULL);
	}
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
	if ((len < MINELEN) | (len > MAXELEN))
		return(oldreadcolrs(scanline, len, fp));
	if ((i = getc(fp)) == EOF)
		return(-1);
	if (i != 2) {
		ungetc(i, fp);
		return(oldreadcolrs(scanline, len, fp));
	}
	scanline[0][GRN] = getc(fp);
	scanline[0][BLU] = getc(fp);
	if ((i = getc(fp)) == EOF)
		return(-1);
	if ((scanline[0][GRN] != 2) | ((scanline[0][BLU] & 0x80) != 0)) {
		scanline[0][RED] = 2;
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
