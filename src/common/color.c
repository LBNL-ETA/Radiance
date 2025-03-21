#ifndef lint
static const char	RCSid[] = "$Id: color.c,v 2.41 2025/03/21 20:04:15 greg Exp $";
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

#define  MINELEN	17	/* minimum scanline length for encoding */
#define  MAXELEN	0x7fff	/* maximum scanline length for encoding */
#define  MINRUN		4	/* minimum run length */

const float	cxponent[256] = {	/* exponent look-up */
.0f, 2.2958874e-41f, 4.59177481e-41f, 9.18354962e-41f,
1.83670992e-40f, 3.67341985e-40f, 7.34683969e-40f, 1.46936794e-39f,
2.93873588e-39f, 5.87747175e-39f, 1.17549435e-38f, 2.3509887e-38f,
4.7019774e-38f, 9.40395481e-38f, 1.88079096e-37f, 3.76158192e-37f,
7.52316385e-37f, 1.50463277e-36f, 3.00926554e-36f, 6.01853108e-36f,
1.20370622e-35f, 2.40741243e-35f, 4.81482486e-35f, 9.62964972e-35f,
1.92592994e-34f, 3.85185989e-34f, 7.70371978e-34f, 1.54074396e-33f,
3.08148791e-33f, 6.16297582e-33f, 1.23259516e-32f, 2.46519033e-32f,
4.93038066e-32f, 9.86076132e-32f, 1.97215226e-31f, 3.94430453e-31f,
7.88860905e-31f, 1.57772181e-30f, 3.15544362e-30f, 6.31088724e-30f,
1.26217745e-29f, 2.5243549e-29f, 5.04870979e-29f, 1.00974196e-28f,
2.01948392e-28f, 4.03896783e-28f, 8.07793567e-28f, 1.61558713e-27f,
3.23117427e-27f, 6.46234854e-27f, 1.29246971e-26f, 2.58493941e-26f,
5.16987883e-26f, 1.03397577e-25f, 2.06795153e-25f, 4.13590306e-25f,
8.27180613e-25f, 1.65436123e-24f, 3.30872245e-24f, 6.6174449e-24f,
1.32348898e-23f, 2.64697796e-23f, 5.29395592e-23f, 1.05879118e-22f,
2.11758237e-22f, 4.23516474e-22f, 8.47032947e-22f, 1.69406589e-21f,
3.38813179e-21f, 6.77626358e-21f, 1.35525272e-20f, 2.71050543e-20f,
5.42101086e-20f, 1.08420217e-19f, 2.16840434e-19f, 4.33680869e-19f,
8.67361738e-19f, 1.73472348e-18f, 3.46944695e-18f, 6.9388939e-18f,
1.38777878e-17f, 2.77555756e-17f, 5.55111512e-17f, 1.11022302e-16f,
2.22044605e-16f, 4.4408921e-16f, 8.8817842e-16f, 1.77635684e-15f,
3.55271368e-15f, 7.10542736e-15f, 1.42108547e-14f, 2.84217094e-14f,
5.68434189e-14f, 1.13686838e-13f, 2.27373675e-13f, 4.54747351e-13f,
9.09494702e-13f, 1.8189894e-12f, 3.63797881e-12f, 7.27595761e-12f,
1.45519152e-11f, 2.91038305e-11f, 5.82076609e-11f, 1.16415322e-10f,
2.32830644e-10f, 4.65661287e-10f, 9.31322575e-10f, 1.86264515e-09f,
3.7252903e-09f, 7.4505806e-09f, 1.49011612e-08f, 2.98023224e-08f,
5.96046448e-08f, 1.1920929e-07f, 2.38418579e-07f, 4.76837158e-07f,
9.53674316e-07f, 1.90734863e-06f, 3.81469727e-06f, 7.62939453e-06f,
1.52587891e-05f, 3.05175781e-05f, 6.10351562e-05f, 0.000122070312f,
0.000244140625f, 0.00048828125f, 0.0009765625f, 0.001953125f,
0.00390625f, 0.0078125f, 0.015625f, 0.03125f, 0.0625f, 0.125f,
0.25f, 0.5f, 1.f, 2.f, 4.f, 8.f, 16.f, 32.f, 64.f, 128.f, 256.f,
512.f, 1024.f, 2048.f, 4096.f, 8192.f, 16384.f, 32768.f, 65536.f,
131072.f, 262144.f, 524288.f, 1048576.f, 2097152.f, 4194304.f,
8388608.f, 16777216.f, 33554432.f, 67108864.f, 134217728.f,
268435456.f, 536870912.f, 1.07374182e+09f, 2.14748365e+09f,
4.2949673e+09f, 8.58993459e+09f, 1.71798692e+10f, 3.43597384e+10f,
6.87194767e+10f, 1.37438953e+11f, 2.74877907e+11f, 5.49755814e+11f,
1.09951163e+12f, 2.19902326e+12f, 4.39804651e+12f, 8.79609302e+12f,
1.7592186e+13f, 3.51843721e+13f, 7.03687442e+13f, 1.40737488e+14f,
2.81474977e+14f, 5.62949953e+14f, 1.12589991e+15f, 2.25179981e+15f,
4.50359963e+15f, 9.00719925e+15f, 1.80143985e+16f, 3.6028797e+16f,
7.2057594e+16f, 1.44115188e+17f, 2.88230376e+17f, 5.76460752e+17f,
1.1529215e+18f, 2.30584301e+18f, 4.61168602e+18f, 9.22337204e+18f,
1.84467441e+19f, 3.68934881e+19f, 7.37869763e+19f, 1.47573953e+20f,
2.95147905e+20f, 5.9029581e+20f, 1.18059162e+21f, 2.36118324e+21f,
4.72236648e+21f, 9.44473297e+21f, 1.88894659e+22f, 3.77789319e+22f,
7.55578637e+22f, 1.51115727e+23f, 3.02231455e+23f, 6.0446291e+23f,
1.20892582e+24f, 2.41785164e+24f, 4.83570328e+24f, 9.67140656e+24f,
1.93428131e+25f, 3.86856262e+25f, 7.73712525e+25f, 1.54742505e+26f,
3.0948501e+26f, 6.1897002e+26f, 1.23794004e+27f, 2.47588008e+27f,
4.95176016e+27f, 9.90352031e+27f, 1.98070406e+28f, 3.96140813e+28f,
7.92281625e+28f, 1.58456325e+29f, 3.1691265e+29f, 6.338253e+29f,
1.2676506e+30f, 2.5353012e+30f, 5.0706024e+30f, 1.01412048e+31f,
2.02824096e+31f, 4.05648192e+31f, 8.11296384e+31f, 1.62259277e+32f,
3.24518554e+32f, 6.49037107e+32f, 1.29807421e+33f, 2.59614843e+33f,
5.19229686e+33f, 1.03845937e+34f, 2.07691874e+34f, 4.15383749e+34f,
8.30767497e+34f, 1.66153499e+35f, 3.32306999e+35f, 6.64613998e+35f
};

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
	double		step, cwl;
	int		i;

	if (NCSAMP == 3) {
		setcolor(scol, r, g, b);
		return;
	}
	step = (WLPART[3] - WLPART[0])/(double)NCSAMP;
	cwl = WLPART[0] + .5*step;
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
setscolr(			/* assign common-exponent spectral from RGB */
	SCOLR sclr,
	double r,
	double g,
	double b
)
{
	SCOLOR	scol;

	setscolor(scol, r, g, b);
	scolor2scolr(sclr, scol, NCSAMP);
}


void
scolor2color(			/* assign RGB color from spectrum */
	COLOR col,
	const SCOLOR scol,		/* uses average over bands */
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
scolr2colr(			/* assign RGBE from common-exponent spectrum */
	COLR clr,
	const SCOLR sclr,
	int ncs,
	const float wlpt[4]
)
#if 0			/* fancier method seems to be slower(!) */
{
	const double	step = (wlpt[3] - wlpt[0])/(double)ncs;
	double		cwl;
	int		csum[3], cnt[3], eshft;
	int		i, j;

	csum[0] = csum[1] = csum[2] = 0;
	cnt[0] = cnt[1] = cnt[2] = 0;
	cwl = wlpt[j=0] + .5*step;
	for (i = 0; i < ncs; i++) {
		csum[j] += sclr[i];
		++cnt[j];
		j += ((cwl += step) < wlpt[j+1]);
	}
	eshft = 7;		/* compute exponent shift */
	for (j = 3; (eshft > 0) & (j-- > 0); ) {
		i = 0;
		while (csum[j] < 128*cnt[j] >> i)
			if (++i >= eshft)
				break;
		if (eshft > i)
			eshft = i;
	}
	if (sclr[ncs] <= eshft) {
		clr[RED] = clr[GRN] = clr[BLU] = 0;
		clr[EXP] = 0;
		return;
	}
	for (j = 3; j--; )
		clr[j] = (csum[j]<<eshft)/cnt[j];

	clr[EXP] = sclr[ncs] - eshft;
}
#else
{
	SCOLOR	scol;
	COLOR	col;

	scolr2scolor(scol, sclr, ncs);
	scolor2color(col, scol, ncs, wlpt);
	setcolr(clr, col[RED], col[GRN], col[BLU]);
}
#endif


void
scolor2colr(			/* assign RGBE from spectral color */
	COLR clr,
	const SCOLOR scol,		/* uses average over bands */
	int ncs,
	const float wlpt[4]
)
{
	COLOR	col;

	scolor2color(col, scol, ncs, wlpt);
	setcolr(clr, col[RED], col[GRN], col[BLU]);
}


void
scolr2color(			/* assign RGB from common exponent */
	COLOR col,
	const SCOLR sclr,
	int ncs,
	const float wlpt[4]
)
{
	SCOLOR	scol;

	scolr2scolor(scol, sclr, ncs);
	scolor2color(col, scol, ncs, wlpt);
}


void
scolor2scolr(			/* float spectrum to common exponent */
	SCOLR sclr,
	const SCOLOR scol,
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
	const SCOLR sclr,
	int ncs
)
{
	const COLORV	f = cxponent[sclr[ncs]];
	int		i = ncs;

	while (i--)
		scol[i] = (sclr[i] + (COLORV)0.5)*f;
}


double
scolor_mean(			/* compute average for spectral color */
	const SCOLOR  scol
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
	const SCOLOR  scol
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
freadscolrs(COLRV *scanline, int nc, int len, FILE *fp)
{
	if (nc < 3)
		return(-1);
	if (nc == 3)
		return(freadcolrs((COLR *)scanline, len, fp));

	if (fread(scanline, nc+1, len, fp) != len)
		return(-1);
	return(0);
}


/* read nc-component common-exponent color scan and convert to COLR's */
int
fread2colrs(COLR *scanline, int len, FILE *fp, int nc, const float wlpt[4])
{
	COLRV  *sclrscan;
	int  n;

	if (nc < 3)
		return(-1);
	if (nc == 3)
		return(freadcolrs(scanline, len, fp));

	sclrscan = (COLRV *)tempbuffer(sizeof(COLRV)*(nc+1)*len);
	if (sclrscan == NULL || freadscolrs(sclrscan, nc, len, fp) < 0)
		return(-1);
	for (n = len; n--; ) {
		scolr2colr(*scanline++, sclrscan, nc, wlpt);
		sclrscan += nc+1;
	}
	return(0);
}


/* write an common-exponent spectral color scanline */
int
fwritescolrs(const COLRV *sscanline, int nc, int len, FILE *fp)
{
	if (nc < 3)
		return(-1);
	if (nc == 3)
		return(fwritecolrs((COLR *)sscanline, len, fp));

	if (fwrite(sscanline, nc+1, len, fp) != len)
		return(-1);
	return(0);
}


int
fwritescan(		/* write out an RGB or XYZ scanline */
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
freadscan(		/* read in an RGB or XYZ scanline */
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
	while (len-- > 0) {
		colr_color(scanline[0], clrscan[0]);
		scanline++; clrscan++;
	}
	return(0);
}


/* read an nc-component color scanline */
int
freadsscan(COLORV *sscanline, int nc, int len, FILE *fp)
{
	COLRV	*tscn = (COLRV *)tempbuffer((nc+1)*len);
	int	i, j;

	if (tscn == NULL || freadscolrs(tscn, nc, len, fp) < 0)
		return(-1);
	for (i = len; i-- > 0; ) {	/* inline scolr2scolor() */
		const COLORV	f = cxponent[tscn[nc]];
		for (j = nc; j--; )
			sscanline[j] = (tscn[j] + (COLORV)0.5)*f;
		sscanline += nc;
		tscn += nc+1;
	}
	return(0);
}


/* read an nc-component color scanline and return as RGB */
int
fread2scan(COLOR *scanline, int len, FILE *fp, int nc, const float wlpt[4])
{
	COLRV	*tscn;
	int	i;

	if (nc < 3)
		return(-1);
	if (nc == 3)
		return(freadscan(scanline, len, fp));

	tscn = (COLRV *)tempbuffer((nc+1)*len);
	if (tscn == NULL || freadscolrs(tscn, nc, len, fp) < 0)
		return(-1);
	for (i = len; i-- > 0; ) {
		scolr2color(*scanline++, tscn, nc, wlpt);
		tscn += nc+1;
	}
	return(0);
}


/* write an nc-component spectral color scanline */
int
fwritesscan(const COLORV *sscanline, int nc, int len, FILE *fp)
{
	COLRV	*tscn = (COLRV *)tempbuffer((nc+1)*len);
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


int
bigdiff(				/* c1 delta c2 > md? */
	const COLOR  c1,
	const COLOR  c2,
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
	const SCOLOR  c1,
	const SCOLOR  c2,
	double  md
)
{
	int  i = NCSAMP;

	while (i--)
		if ((c1[i]-c2[i] > md*c2[i]) | (c2[i]-c1[i] > md*c1[i]))
			return(1);
	return(0);
}
