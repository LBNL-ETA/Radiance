#ifndef lint
static const char	RCSid[] = "$Id: color.c,v 2.31 2023/11/27 21:00:14 greg Exp $";
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


int  CNDX[4] = {0,1,2,3};	/* RGBE indices for SCOLOR, SCOLR */
float  WLPART[4] = {780,588,480,380};	/* RGB wavelength limits+partitions (nm) */


int
setspectrsamp(			/* assign spectral sampling, 1 if good, -1 if bad */
	int cn[4],		/* input cn[3]=nsamps */
	float wlpt[4]		/* input wlpt[0],wlpt[3]=extrema */
)
{
	static const float	PKWL[3] = {607, 553, 469};
	int			i, j;

	if (cn[3] < 3)
		return(-1);		/* reject this */

	if (wlpt[0] < wlpt[3]) {
		float	tf = wlpt[0];
		wlpt[0] = wlpt[3]; wlpt[3] = tf;
	}
	if (wlpt[0] - wlpt[3] < 50.f)
		return(-1);		/* also reject */

	if (cn[3] > MAXCSAMP)
		cn[3] = MAXCSAMP;

	if ((wlpt[3] >= PKWL[2]) | (wlpt[0] <= PKWL[0])) {
		wlpt[1] = wlpt[0] + 0.333333f*(wlpt[3]-wlpt[0]);
		wlpt[2] = wlpt[0] + 0.666667f*(wlpt[3]-wlpt[0]);
		cn[0] = 0; cn[1] = cn[3]/3; cn[2] = cn[3]*2/3;
		return(0);		/* unhappy but non-fatal return value */
	}
	wlpt[1] = 588.f;		/* tuned for standard green channel */
	wlpt[2] = 480.f;
	if (cn[3] == 3) {		/* nothing to tune? */
		cn[0] = 0; cn[1] = 1; cn[2] = 2;
	} else {			/* else find nearest color indices */
		double	curwl[3];
		memset(curwl, 0, sizeof(curwl));
		for (i = cn[3]; i--; ) {
			const float	cwl = (i+.5f)/cn[3]*(wlpt[3]-wlpt[0]) + wlpt[0];
			for (j = 3; j--; )
				if (fabs(cwl - PKWL[j]) < fabs(curwl[j] - PKWL[j])) {
					curwl[j] = cwl;
					cn[j] = i;
				}
		}
	}
	return(1);			/* happy return value */
}


void
setscolor(			/* assign spectral color from RGB */
	SCOLOR scol,
	double r,
	double g,
	double b
)
{
	const double	step = (WLPART[3] - WLPART[0])/(double)NCSAMP;
	double		cwl = WLPART[0] + .5*step;
	int		i;

	for (i = 0; i < NCSAMP; i++) {
		if (cwl >= WLPART[1])
			scol[i] = r;
		else if (cwl >= WLPART[2])
			scol[i] = g;
		else
			scol[i] = b;
		cwl += step;
	}
}


void
scolor2color(			/* assign RGB color from spectrum */
	COLOR col,
	SCOLOR scol,		/* uses average over bands */
	int ncs,
	const float wlpt[4]
)
{
	const double	step = (wlpt[3] - wlpt[0])/(double)ncs;
	double		cwl = wlpt[0] + .5*step;
	int		i, j=0, n=0;

	setcolor(col, 0, 0, 0);
	for (i = 0; i < ncs; i++) {
		if (cwl < wlpt[j+1]) {
			if (n > 1) col[j] /= (COLORV)n;
			j++;
			n = 0;
		}
		col[j] += scol[i];
		n++;
		cwl += step;
	}
	if (n > 1) col[j] /= (COLORV)n;
}


void
scolor2colr(			/* assign RGBE from spectral color */
	COLR clr,
	SCOLOR scol,		/* uses average over bands */
	int ncs,
	const float wlpt[4]
)
{
	COLOR	col;

	scolor2color(col, scol, ncs, wlpt);
	setcolr(clr, col[RED], col[GRN], col[BLU]);
}


void
scolor2scolr(			/* float spectrum to common exponent */
	SCOLR sclr,
	SCOLOR scol,
	int ncs
)
{
	int	i = ncs;
	COLORV	p = scol[--i];

	while (i)
		if (scol[--i] > p)
			p = scol[i];
	if (p <= 1e-32) {
		memset(sclr, 0, ncs+1);
		return;
	}
	p = frexp(p, &i) * 256.0 / p;
	sclr[ncs] = i + COLXS;
	for (i = ncs; i--; )
		sclr[i] = (scol[i] > 0) * (int)(scol[i]*p);
}


void
scolr2scolor(			/* common exponent to float spectrum */
	SCOLOR scol,
	SCOLR sclr,
	int ncs
)
{
	double	f;
	int	i;

	if (sclr[ncs] == 0) {
		memset(scol, 0, sizeof(COLORV)*ncs);
		return;
	}
	f = ldexp(1.0, (int)sclr[ncs]-(COLXS+8));

	for (i = ncs; i--; )
		scol[i] = (sclr[i] + 0.5)*f;
}


double
scolor_mean(			/* compute average for spectral color */
	SCOLOR  scol
)
{
	int	i = NCSAMP;
	double	sum = 0;

	while (i--)
		sum += scol[i];

	return sum/(double)NCSAMP;
}


double
sintens(			/* find maximum value from spectrum */
	SCOLOR  scol
)
{
	int	i = NCSAMP;
	COLORV	peak = scol[--i];

	while (i)
		if (scol[--i] > peak)
			peak = scol[i];

	return peak;
}


void
convertscolor(			/* spectrum conversion, zero-fill ends */
	SCOLOR dst,		/* destination spectrum */
	int dnc,		/* destination # of spectral samples/intervals */
	double dwl0,		/* starting destination wavelength (longer) */
	double dwl1,		/* ending destination wavelength (shorter) */
	const COLORV src[],	/* source spectrum array */
	int snc,
	double swl0,		/* long/short wavelengths may be reversed */
	double swl1
)
{
	const int	sdir = 1 - 2*(swl0 < swl1);
	const double	sstp = (swl1 - swl0)/(double)snc;
	const double	dstp = (dwl1 - dwl0)/(double)dnc;
	const double	rdstp = 1./dstp;
	int		si, ssi, di;
	double		wl;

	if ((dnc < 3) | (dwl0 <= dwl1) | (dst == src))
		return;		/* invalid destination */

	if (dnc == snc && (dwl0-swl0)*(dwl0-swl0) + (dwl1-swl1)*(dwl1-swl1) <= .5) {
		memcpy(dst, src, sizeof(COLORV)*dnc);
		return;		/* same spectral sampling */
	}
	memset(dst, 0, sizeof(COLORV)*dnc);
				/* set starting positions */
	if ((sdir>0 ? swl0 : swl1) <= dwl0) {
		if (sdir > 0) {
			wl = swl0;
			ssi = 0;
		} else {
			wl = swl1;
			ssi = snc-1;
		}
		si = 0;
		di = (wl - dwl0)*rdstp;
	} else {
		wl = dwl0;
		if (sdir > 0) {
			ssi = si = (wl - swl0)/sstp;
		} else {
			si = (wl - swl1)/sstp;
			ssi = snc-1 - si;
		}
		di = 0;
	}
	swl0 += (sdir < 0)*sstp;
				/* step through intervals */
	while ((si < snc) & (di < dnc)) {
		double	intvl;
		if (swl0 + (ssi+sdir)*sstp < dwl0 + (di+1)*dstp) {
			intvl = dwl0 + (di+1)*dstp - wl;
			dst[di++] += src[ssi]*intvl*rdstp;
		} else {
			intvl = swl0 + (ssi+sdir)*sstp - wl;
			dst[di] += src[ssi]*intvl*rdstp;
			ssi += sdir;
			si++;
		}
		wl += intvl;
	}
}


void *
tempbuffer(			/* get a temporary buffer */
	size_t  len
)
{
	static void	*tempbuf = NULL;
	static size_t	tempbuflen = 0;

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


/* read an nc-component common-exponent color scanline */
int
freadscolrs(uby8 *scanline, int nc, int len, FILE *fp)
{
	if (fread(scanline, nc+1, len, fp) != len)
		return(-1);
	return(0);
}


/* write an common-exponent spectral color scanline */
int
fwritescolrs(uby8 *sscanline, int nc, int len, FILE *fp)
{
	if (fwrite(sscanline, nc+1, len, fp) != len)
		return(-1);
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


/* read an nc-component color scanline */
int
freadsscan(COLORV *sscanline, int nc, int len, FILE *fp)
{
	uby8	*tscn = (uby8 *)tempbuffer((nc+1)*len);
	int	i;

	if (tscn == NULL || freadscolrs(tscn, nc, len, fp) < 0)
		return(-1);
	for (i = len; i-- > 0; ) {
		scolr2scolor(sscanline, tscn, nc);
		sscanline += nc;
		tscn += nc+1;
	}
	return(0);
}


/* write an spectral color scanline (NCSAMP) */
int
fwritesscan(COLORV *sscanline, int nc, int len, FILE *fp)
{
	uby8	*tscn = (uby8 *)tempbuffer((nc+1)*len);
	int	i;

	if (tscn == NULL)
		return(-1);
	for (i = 0; i < len; i++) {
		scolor2scolr(tscn+i*(nc+1), sscanline, nc);
		sscanline += nc;
	}
	return(fwritescolrs(tscn, nc, len, fp));
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

	clr[RED] = (r > 0) * (int)(r*d);
	clr[GRN] = (g > 0) * (int)(g*d);
	clr[BLU] = (b > 0) * (int)(b*d);
	clr[EXP] = e + COLXS;
}


void
colr_color(			/* convert short to float color */
	COLOR  col,
	COLR  clr
)
{
	double  f;
	
	if (clr[EXP] == 0) {
		col[RED] = col[GRN] = col[BLU] = 0.0;
		return;
	}
	f = ldexp(1.0, (int)clr[EXP]-(COLXS+8));
	col[RED] = (clr[RED] + 0.5)*f;
	col[GRN] = (clr[GRN] + 0.5)*f;
	col[BLU] = (clr[BLU] + 0.5)*f;
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
		if ((colval(c1,i)-colval(c2,i) > md*colval(c2,i)) |
				(colval(c2,i)-colval(c1,i) > md*colval(c1,i)))
			return(1);
	return(0);
}


int
sbigsdiff(				/* sc1 delta sc2 > md? */
	SCOLOR  c1,
	SCOLOR  c2,
	double  md
)
{
	int  i = NCSAMP;

	while (i--)
		if ((c1[i]-c2[i] > md*c2[i]) | (c2[i]-c1[i] > md*c1[i]))
			return(1);
	return(0);
}
