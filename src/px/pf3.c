#ifndef lint
static const char RCSid[] = "$Id: pf3.c,v 2.17 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
 *  pf3.c - routines for gaussian and box filtering
 *
 *     8/13/86
 */

#include  "standard.h"

#include  <string.h>

#include  "color.h"

#define  RSCA		1.13	/* square-radius multiplier: sqrt(4/PI) */
#define  TEPS		0.2	/* threshold proximity goal */
#define  REPS		0.1	/* radius proximity goal */

extern double  CHECKRAD;	/* radius over which gaussian is summed */

extern double  rad;		/* output pixel radius for filtering */

extern double  thresh;		/* maximum contribution for subpixel */

extern int  nrows;		/* number of rows for output */
extern int  ncols;		/* number of columns for output */

extern int  xres, yres;		/* resolution of input */

extern double  x_c, y_r;	/* conversion factors */

extern int  xrad;		/* x search radius */
extern int  yrad;		/* y search radius */
extern int  xbrad;		/* x box size */
extern int  ybrad;		/* y box size */

extern int  barsize;		/* size of input scan bar */
extern COLOR  **scanin;		/* input scan bar */
extern COLOR  *scanout;		/* output scan line */
extern COLOR  **scoutbar;	/* output scan bar (if thresh > 0) */
extern float  **greybar;	/* grey-averaged input values */
extern int  obarsize;		/* size of output scan bar */
extern int  orad;		/* output window radius */

extern int  wrapfilt;		/* wrap filter horizontally? */

extern char  *progname;

float  *gausstable;		/* gauss lookup table */

float  *ringsum;		/* sum of ring values */
short  *ringwt;			/* weight (count) of ring values */
short  *ringndx;		/* ring index table */
float  *warr;			/* array of pixel weights */

extern double  (*ourbright)();	/* brightness computation function */

double  pickfilt();

#define	 lookgauss(x)		gausstable[(int)(10.*(x)+.5)]


initmask()			/* initialize gaussian lookup table */
{
	int  gtabsiz;
	double  gaussN;
	double  d;
	register int  x;

	gtabsiz = 111*CHECKRAD*CHECKRAD;
	gausstable = (float *)malloc(gtabsiz*sizeof(float));
	if (gausstable == NULL)
		goto memerr;
	d = x_c*y_r*0.25/(rad*rad);
	gausstable[0] = exp(-d);
	for (x = 1; x < gtabsiz; x++)
		if (x*0.1 <= d)
			gausstable[x] = gausstable[0];
		else
			gausstable[x] = exp(-x*0.1);
	if (obarsize == 0)
		return;
					/* compute integral of filter */
	gaussN = PI*d*exp(-d);			/* plateau */
	for (d = sqrt(d)+0.05; d <= RSCA*CHECKRAD; d += 0.1)
		gaussN += 0.1*2.0*PI*d*exp(-d*d);
					/* normalize filter */
	gaussN = x_c*y_r/(rad*rad*gaussN);
	for (x = 0; x < gtabsiz; x++)
		gausstable[x] *= gaussN;
					/* create ring averages table */
	ringndx = (short *)malloc((2*orad*orad+1)*sizeof(short));
	if (ringndx == NULL)
		goto memerr;
	for (x = 2*orad*orad+1; --x > orad*orad; )
		ringndx[x] = -1;
	do
		ringndx[x] = sqrt((double)x);
	while (x--);
	ringsum = (float *)malloc((orad+1)*sizeof(float));
	ringwt = (short *)malloc((orad+1)*sizeof(short));
	warr = (float *)malloc(obarsize*obarsize*sizeof(float));
	if ((ringsum == NULL) | (ringwt == 0) | (warr == NULL))
		goto memerr;
	return;
memerr:
	fprintf(stderr, "%s: out of memory in initmask\n", progname);
	quit(1);
}


dobox(csum, xcent, ycent, c, r)			/* simple box filter */
COLOR  csum;
int  xcent, ycent;
int  c, r;
{
	int  wsum;
	double  d;
	int  y;
	register int  x, offs;
	register COLOR	*scan;
	
	wsum = 0;
	setcolor(csum, 0.0, 0.0, 0.0);
	for (y = ycent+1-ybrad; y <= ycent+ybrad; y++) {
		if (y < 0) continue;
		if (y >= yres) break;
		d = y_r < 1.0 ? y_r*y - (r+.5) : (double)(y - ycent);
		if (d < -0.5) continue;
		if (d >= 0.5) break;
		scan = scanin[y%barsize];
		for (x = xcent+1-xbrad; x <= xcent+xbrad; x++) {
			offs = x < 0 ? xres : x >= xres ? -xres : 0;
			if (offs && !wrapfilt)
				continue;
			d = x_c < 1.0 ? x_c*x - (c+.5) : (double)(x - xcent);
			if (d < -0.5) continue;
			if (d >= 0.5) break;
			wsum++;
			addcolor(csum, scan[x+offs]);
		}
	}
	if (wsum > 1) {
		d = 1.0/wsum;
		scalecolor(csum, d);
	}
}


dogauss(csum, xcent, ycent, c, r)		/* gaussian filter */
COLOR  csum;
int  xcent, ycent;
int  c, r;
{
	double  dy, dx, weight, wsum;
	COLOR  ctmp;
	int  y;
	register int  x, offs;
	register COLOR	*scan;

	wsum = FTINY;
	setcolor(csum, 0.0, 0.0, 0.0);
	for (y = ycent-yrad; y <= ycent+yrad; y++) {
		if (y < 0) continue;
		if (y >= yres) break;
		dy = (y_r*(y+.5) - (r+.5))/rad;
		scan = scanin[y%barsize];
		for (x = xcent-xrad; x <= xcent+xrad; x++) {
			offs = x < 0 ? xres : x >= xres ? -xres : 0;
			if (offs && !wrapfilt)
				continue;
			dx = (x_c*(x+.5) - (c+.5))/rad;
			weight = lookgauss(dx*dx + dy*dy);
			wsum += weight;
			copycolor(ctmp, scan[x+offs]);
			scalecolor(ctmp, weight);
			addcolor(csum, ctmp);
		}
	}
	weight = 1.0/wsum;
	scalecolor(csum, weight);
}


dothresh(xcent, ycent, ccent, rcent)	/* gaussian threshold filter */
int  xcent, ycent;
int  ccent, rcent;
{
	double  d;
	int  r, y, offs;
	register int  c, x;
	register float  *gscan;
					/* compute ring sums */
	memset((char *)ringsum, '\0', (orad+1)*sizeof(float));
	memset((char *)ringwt, '\0', (orad+1)*sizeof(short));
	for (r = -orad; r <= orad; r++) {
		if (rcent+r < 0) continue;
		if (rcent+r >= nrows) break;
		gscan = greybar[(rcent+r)%obarsize];
		for (c = -orad; c <= orad; c++) {
			offs = ccent+c < 0 ? ncols :
					ccent+c >= ncols ? -ncols : 0;
			if (offs && !wrapfilt)
				continue;
			x = ringndx[c*c + r*r];
			if (x < 0) continue;
			ringsum[x] += gscan[ccent+c+offs];
			ringwt[x]++;
		}
	}
					/* filter each subpixel */
	for (y = ycent+1-ybrad; y <= ycent+ybrad; y++) {
		if (y < 0) continue;
		if (y >= yres) break;
		d = y_r < 1.0 ? y_r*y - (rcent+.5) : (double)(y - ycent);
		if (d < -0.5) continue;
		if (d >= 0.5) break;
		for (x = xcent+1-xbrad; x <= xcent+xbrad; x++) {
			offs = x < 0 ? xres : x >= xres ? -xres : 0;
			if (offs && !wrapfilt)
				continue;
			d = x_c < 1.0 ? x_c*x - (ccent+.5) : (double)(x - xcent);
			if (d < -0.5) continue;
			if (d >= 0.5) break;
			sumans(x, y, rcent, ccent,
			pickfilt((*ourbright)(scanin[y%barsize][x+offs])));
		}
	}
}


double
pickfilt(p0)			/* find filter multiplier for p0 */
double  p0;
{
	double  m = 1.0;
	double  t, num, denom, avg, wsum;
	double  mlimit[2];
	int  ilimit = 4.0/TEPS;
	register int  i;
				/* iterative search for m */
	mlimit[0] = 1.0; mlimit[1] = orad/rad/CHECKRAD;
	do {
					/* compute grey weighted average */
		i = RSCA*CHECKRAD*rad*m + .5;
		if (i > orad) i = orad;
		avg = wsum = 0.0;
		while (i--) {
			t = (i+.5)/(m*rad);
			t = lookgauss(t*t);
			avg += t*ringsum[i];
			wsum += t*ringwt[i];
		}
		if (avg < 1e-20)		/* zero inclusive average */
			return(1.0);
		avg /= wsum;
					/* check threshold */
		denom = m*m/gausstable[0] - p0/avg;
		if (denom <= FTINY) {		/* zero exclusive average */
			if (m >= mlimit[1]-REPS)
				break;
			m = mlimit[1];
			continue;
		}
		num = p0/avg - 1.0;
		if (num < 0.0) num = -num;
		t = num/denom;
		if (t <= thresh) {
			if (m <= mlimit[0]+REPS || (thresh-t)/thresh <= TEPS)
				break;
		} else {
			if (m >= mlimit[1]-REPS || (t-thresh)/thresh <= TEPS)
				break;
		}
		t = m;			/* remember current m */
					/* next guesstimate */
		m = sqrt(gausstable[0]*(num/thresh + p0/avg));
		if (m < t) {		/* bound it */
			if (m <= mlimit[0]+FTINY)
				m = 0.5*(mlimit[0] + t);
			mlimit[1] = t;
		} else {
			if (m >= mlimit[1]-FTINY)
				m = 0.5*(mlimit[1] + t);
			mlimit[0] = t;
		}
	} while (--ilimit > 0);
	return(m);
}


sumans(px, py, rcent, ccent, m)		/* sum input pixel to output */
int  px, py;
int  rcent, ccent;
double  m;
{
	double  dy2, dx;
	COLOR  pval, ctmp;
	int  ksiz, r, offs;
	double  pc, pr, norm;
	register int  i, c;
	register COLOR	*scan;
	/*
	 * This normalization method fails at the picture borders because
	 * a different number of input pixels contribute there.
	 */
	scan = scanin[py%barsize] + (px < 0 ? xres : px >= xres ? -xres : 0);
	copycolor(pval, scan[px]);
	pc = x_c*(px+.5);
	pr = y_r*(py+.5);
	ksiz = CHECKRAD*m*rad + 1;
	if (ksiz > orad) ksiz = orad;
						/* compute normalization */
	norm = 0.0;
	i = 0;
	for (r = rcent-ksiz; r <= rcent+ksiz; r++) {
		if (r < 0) continue;
		if (r >= nrows) break;
		dy2 = (pr - (r+.5))/(m*rad);
		dy2 *= dy2;
		for (c = ccent-ksiz; c <= ccent+ksiz; c++) {
			if (!wrapfilt) {
				if (c < 0) continue;
				if (c >= ncols) break;
			}
			dx = (pc - (c+.5))/(m*rad);
			norm += warr[i++] = lookgauss(dx*dx + dy2);
		}
	}
	norm = 1.0/norm;
	if (x_c < 1.0) norm *= x_c;
	if (y_r < 1.0) norm *= y_r;
						/* sum pixels */
	i = 0;
	for (r = rcent-ksiz; r <= rcent+ksiz; r++) {
		if (r < 0) continue;
		if (r >= nrows) break;
		scan = scoutbar[r%obarsize];
		for (c = ccent-ksiz; c <= ccent+ksiz; c++) {
			offs = c < 0 ? ncols : c >= ncols ? -ncols : 0;
			if (offs && !wrapfilt)
				continue;
			copycolor(ctmp, pval);
			dx = norm*warr[i++];
			scalecolor(ctmp, dx);
			addcolor(scan[c+offs], ctmp);
		}
	}
}
