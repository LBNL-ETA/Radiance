#ifndef lint
static const char	RCSid[] = "$Id: mgraph.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  mgraph.c - routines for plotting graphs from variables.
 *
 *     6/23/86
 *
 *     Greg Ward Larson
 */

#include  <stdio.h>
#include  <string.h>

#include  "meta.h"
#include  "mgvars.h"
#include  "mgraph.h"

extern char  *progname;			/* argv[0] */

extern double  goodstep(), floor(), ceil(), sin(), cos();

static BOUNDS  xbounds, ybounds;	/* the boundaries for the graph */

static double  period = DEFPERIOD;	/* period for polar plot */

static double  axbegin, axsize;		/* the mapped x axis boundaries */
static double  aybegin, aysize;		/* the mapped y axis boundaries */

static int  npltbl[MAXCUR];		/* plottable points per curve */

static double  lastx, lasty;		/* last curve postion */
static int  nplottable;			/* number of plottable points */
static int  nplotted;			/* number of plotted points */

static void getbounds(void);
static void polaraxis(void);
static void makeaxis(void);
static void plotcurves(void);
static void cartaxis(void);
static void stretchbounds(int  c, double  x, double  y);
static void boxstring(int  xmin, int ymin, int xmax, int ymax,
	char  *s, int  d, int width, int color);
static void drawcircle(int  x, int y, int r, int  typ, int wid, int col);
static int rconv(double  r);
static int xconv(double  x);
static int yconv(double  y);
static void csymbol(int  c, double  u, double v);
static int cmline(int  c, int  x, int y);
static void cmsymbol(int  c, int  x, int y);
static int inbounds(double  x, double y);
static void climline(int  c, double  x, double y, double xout, double yout);
static void nextpoint(register int  c, double  x, double  y);
static void cline(int  c, double  u1, double v1, double u2, double v2);


void
mgraph(void)			/* plot the current graph */
{
					/* load the symbol file */
	if (gparam[SYMFILE].flags & DEFINED)
		minclude(gparam[SYMFILE].v.s);
					/* check for polar plot */
	if (gparam[PERIOD].flags & DEFINED)
		period = varvalue(gparam[PERIOD].name);

	getbounds();			/* get the boundaries */
	makeaxis();			/* draw the coordinate axis */
	plotcurves();			/* plot the curves */

	mendpage();			/* advance page */
}


void
getbounds(void)			/* compute the boundaries */
{
	int  i;

	xbounds.min = gparam[XMIN].flags & DEFINED ?
			 varvalue(gparam[XMIN].name) - FTINY :
			 FHUGE ;
	xbounds.max = gparam[XMAX].flags & DEFINED ?
			 varvalue(gparam[XMAX].name) + FTINY :
			 -FHUGE ;
	ybounds.min = gparam[YMIN].flags & DEFINED ?
			 varvalue(gparam[YMIN].name) - FTINY :
			 FHUGE ;
	ybounds.max = gparam[YMAX].flags & DEFINED ?
			 varvalue(gparam[YMAX].name) + FTINY :
			 -FHUGE ;

	nplottable = 0;
	for (i = 0; i < MAXCUR; i++) {
		npltbl[i] = 0;
		mgcurve(i, stretchbounds);
		nplottable += npltbl[i];
	}
	if (nplottable == 0) {
		fprintf(stderr, "%s: no plottable data\n", progname);
		quit(1);
	}

	xbounds.step = gparam[XSTEP].flags & DEFINED ?
			  varvalue(gparam[XSTEP].name) :
			  period > FTINY ?
				DEFPLSTEP*period :
				goodstep(xbounds.max - xbounds.min) ;
	if (!(gparam[XMIN].flags & DEFINED))
		xbounds.min = floor(xbounds.min/xbounds.step) *
					xbounds.step;
	if (!(gparam[XMAX].flags & DEFINED))
		xbounds.max = ceil(xbounds.max/xbounds.step) *
					xbounds.step;
	ybounds.step = gparam[YSTEP].flags & DEFINED ?
			  varvalue(gparam[YSTEP].name) :
			  period > FTINY ?
				goodstep((ybounds.max - ybounds.min)*1.75) :
				goodstep(ybounds.max - ybounds.min) ;
	if (!(gparam[YMIN].flags & DEFINED))
		ybounds.min = floor(ybounds.min/ybounds.step) *
					ybounds.step;
	if (!(gparam[YMAX].flags & DEFINED))
		ybounds.max = ceil(ybounds.max/ybounds.step) *
					ybounds.step;
	if (gparam[XMAP].flags & DEFINED) {
		axbegin = funvalue(gparam[XMAP].name, 1, &xbounds.min);
		axsize = funvalue(gparam[XMAP].name, 1, &xbounds.max);
		axsize -= axbegin;
	} else {
		axbegin = xbounds.min;
		axsize = xbounds.max;
		axsize -= axbegin;
	}
	if (gparam[YMAP].flags & DEFINED) {
		aybegin = funvalue(gparam[YMAP].name, 1, &ybounds.min);
		aysize = funvalue(gparam[YMAP].name, 1, &ybounds.max);
		aysize -= aybegin;
	} else {
		aybegin = ybounds.min;
		aysize = ybounds.max;
		aysize -= aybegin;
	}
}


void
makeaxis(void)				/* draw the coordinate axis */
{
	char  stmp[64];

	if (period > FTINY)
		polaraxis();
	else
		cartaxis();
						/* x axis label */
	if (gparam[XLABEL].flags & DEFINED)
		boxstring(XL_L,XL_D,XL_R,XL_U,gparam[XLABEL].v.s,'r',0,0);
						/* x mapping */
	if (gparam[XMAP].flags & DEFINED) {
		mgtoa(stmp, &gparam[XMAP]);
		boxstring(XM_L,XM_D,XM_R,XM_U,stmp,'r',0,0);
	}
						/* y axis label */
	if (gparam[YLABEL].flags & DEFINED)
		boxstring(YL_L,YL_D,YL_R,YL_U,gparam[YLABEL].v.s,'u',0,0);
						/* y mapping */
	if (gparam[YMAP].flags & DEFINED) {
		mgtoa(stmp, &gparam[YMAP]);
		boxstring(YM_L,YM_D,YM_R,YM_U,stmp,'u',0,0);
	}
						/* title */
	if (gparam[TITLE].flags & DEFINED)
		boxstring(TI_L,TI_D,TI_R,TI_U,gparam[TITLE].v.s,'r',2,0);
						/* subtitle */
	if (gparam[SUBTITLE].flags & DEFINED)
		boxstring(ST_L,ST_D,ST_R,ST_U,gparam[SUBTITLE].v.s,'r',1,0);
						/* legend */
	if (gparam[LEGEND].flags & DEFINED)
		mtext(LT_X, LT_Y, gparam[LEGEND].v.s, CPI, 0);
}


void
polaraxis(void)				/* print polar coordinate axis */
{
	int  lw, tstyle, t0, t1;
	double  d, d1, ybeg;
	char  stmp[64], *fmt, *goodformat();
						/* get tick style */
	if (gparam[TSTYLE].flags & DEFINED)
		tstyle = varvalue(gparam[TSTYLE].name) + 0.5;
	else
		tstyle = DEFTSTYLE;
						/* start of numbering */
	ybeg = ceil(ybounds.min/ybounds.step)*ybounds.step;
						/* theta (x) numbering */
	fmt = goodformat(xbounds.step);
	for (d = 0.0; d < period-FTINY; d += xbounds.step) {
		sprintf(stmp, fmt, d);
		d1 = d*(2*PI)/period;
		t0 = TN_X + TN_R*cos(d1) + .5;
		if (t0 < TN_X)
			t0 -= strlen(stmp)*CWID;
		mtext(t0,(int)(TN_Y+TN_R*sin(d1)+.5),stmp,CPI,0);
	}
						/* radius (y) numbering */
	fmt = goodformat(ybounds.step);
	lw = PL_R+RN_S;
	for (d = ybeg; d <= ybounds.max+FTINY; d += ybounds.step) {
		t0 = rconv(d);
		if (t0 >= lw+RN_S || t0 <= lw-RN_S) {
			sprintf(stmp, fmt, d);
			mtext(RN_X+t0-strlen(stmp)*(CWID/2),RN_Y,stmp,CPI,0);
			lw = t0;
		}
	}
						/* frame */
	if (gparam[FTHICK].flags & DEFINED)
		lw = varvalue(gparam[FTHICK].name) + 0.5;
	else
		lw = DEFFTHICK;
	if (lw-- > 0) {
		drawcircle(PL_X,PL_Y,PL_R,0,lw,0);
		switch (tstyle) {
		case 1:		/* outside */
			t0 = 0; t1 = TLEN; break;
		case 2:		/* inside */
			t0 = TLEN; t1 = 0; break;
		case 3:		/* accross */
			t0 = TLEN/2; t1 = TLEN/2; break;
		default:	/* none */
			t0 = t1 = 0; break;
		}
		if (t0 + t1) {
			for (d = 0.0; d < 2*PI-FTINY;
					d += xbounds.step*(2*PI)/period) {
				mline((int)(PL_X+(PL_R-t0)*cos(d)+.5),
					(int)(PL_Y+(PL_R-t0)*sin(d)+.5),
					0, lw, 0);
				mdraw((int)(PL_X+(PL_R+t1)*cos(d)+.5),
					(int)(PL_Y+(PL_R+t1)*sin(d)+.5));
			}
		}
	}
						/* origin */
	if (gparam[OTHICK].flags & DEFINED)
		lw = varvalue(gparam[OTHICK].name) + 0.5;
	else
		lw = DEFOTHICK;
	if (lw-- > 0) {
		mline(PL_X-PL_R,PL_Y,0,lw,0);
		mdraw(PL_X+PL_R,PL_Y);
		mline(PL_X,PL_Y-PL_R,0,lw,0);
		mdraw(PL_X,PL_Y+PL_R);
		if (tstyle > 0)
			for (d = ybeg; d <= ybounds.max+FTINY;
					d += ybounds.step) {
				t0 = rconv(d);
				mline(PL_X+t0,PL_Y-TLEN/2,0,lw,0);
				mdraw(PL_X+t0,PL_Y+TLEN/2);
				mline(PL_X-TLEN/2,PL_Y+t0,0,lw,0);
				mdraw(PL_X+TLEN/2,PL_Y+t0);
				mline(PL_X-t0,PL_Y-TLEN/2,0,lw,0);
				mdraw(PL_X-t0,PL_Y+TLEN/2);
				mline(PL_X-TLEN/2,PL_Y-t0,0,lw,0);
				mdraw(PL_X+TLEN/2,PL_Y-t0);
			}
	}
						/* grid */
	if (gparam[GRID].flags & DEFINED)
		lw = varvalue(gparam[GRID].name);
	else
		lw = DEFGRID;
	if (lw-- > 0) {
		for (d = 0.0; d < PI-FTINY; d += xbounds.step*(2*PI)/period) {
			mline((int)(PL_X+PL_R*cos(d)+.5),
				(int)(PL_Y+PL_R*sin(d)+.5),2,0,0);
			mdraw((int)(PL_X-PL_R*cos(d)+.5),
				(int)(PL_Y-PL_R*sin(d)+.5));
		}
		for (d = ybeg; d <= ybounds.max + FTINY; d += ybounds.step)
			drawcircle(PL_X,PL_Y,rconv(d),2,0,0);
	}
}

void
cartaxis(void)				/* print Cartesian coordinate axis */
{
	int  lw, t0, t1, tstyle;
	double  d, xbeg, ybeg;
	char  stmp[64], *fmt, *goodformat();
	register int  i;
						/* get tick style */
	if (gparam[TSTYLE].flags & DEFINED)
		tstyle = varvalue(gparam[TSTYLE].name) + 0.5;
	else
		tstyle = DEFTSTYLE;
						/* start of numbering */
	xbeg = ceil(xbounds.min/xbounds.step)*xbounds.step;
	ybeg = ceil(ybounds.min/ybounds.step)*ybounds.step;
	
						/* x numbering */
	fmt = goodformat(xbounds.step);
	lw = 2*AX_L-AX_R;
	for (d = xbeg;
			d <= xbounds.max + FTINY;
			d += xbounds.step)
		if ((i = xconv(d)) >= lw+XN_S || i <= lw-XN_S) {
			sprintf(stmp, fmt, d);
			mtext(i-strlen(stmp)*(CWID/2)+XN_X,XN_Y,stmp,CPI,0);
			lw = i;
		}
						/* y numbering */
	fmt = goodformat(ybounds.step);
	lw = 2*AX_D-AX_U;
	for (d = ybeg;
			d <= ybounds.max + FTINY;
			d += ybounds.step)
		if ((i = yconv(d)) >= lw+YN_S || i <= lw-YN_S) {
			sprintf(stmp, fmt, d);
			mtext(YN_X-strlen(stmp)*CWID,i+YN_Y,stmp,CPI,0);
			lw = i;
		}
						/* frame */
	if (gparam[FTHICK].flags & DEFINED)
		lw = varvalue(gparam[FTHICK].name) + 0.5;
	else
		lw = DEFFTHICK;
	if (lw-- > 0) {
		mline(AX_L,AX_D,0,lw,0);
		mdraw(AX_R,AX_D);
		mdraw(AX_R,AX_U);
		mdraw(AX_L,AX_U);
		mdraw(AX_L,AX_D);
		switch (tstyle) {
		case 1:		/* outside */
			t0 = 0; t1 = TLEN; break;
		case 2:		/* inside */
			t0 = TLEN; t1 = 0; break;
		case 3:		/* accross */
			t0 = TLEN/2; t1 = TLEN/2; break;
		default:	/* none */
			t0 = t1 = 0; break;
		}
		if (t0 + t1) {
			for (d = xbeg;
					d <= xbounds.max + FTINY;
					d += xbounds.step) {
				i = xconv(d);
				mline(i,AX_D+t0,0,lw,0);
				mdraw(i,AX_D-t1);
				mline(i,AX_U-t0,0,lw,0);
				mdraw(i,AX_U+t1);
			}
			for (d = ybeg;
					d <= ybounds.max + FTINY;
					d += ybounds.step) {
				i = yconv(d);
				mline(AX_L+t0,i,0,lw,0);
				mdraw(AX_L-t1,i);
				mline(AX_R-t0,i,0,lw,0);
				mdraw(AX_R+t1,i);
			}
		}
	}
						/* origin */
	if (gparam[OTHICK].flags & DEFINED)
		lw = varvalue(gparam[OTHICK].name) + 0.5;
	else
		lw = DEFOTHICK;
	if (lw-- > 0) {
		i = yconv(0.0);
		if (i >= AX_D && i <= AX_U) {
			mline(AX_L,i,0,lw,0);
			mdraw(AX_R,i);
			if (tstyle > 0)
				for (d = xbeg; d <= xbounds.max+FTINY;
						d += xbounds.step) {
					mline(xconv(d),i+TLEN/2,0,lw,0);
					mdraw(xconv(d),i-TLEN/2);
				}
		}
		i = xconv(0.0);
		if (i >= AX_L && i <= AX_R) {
			mline(i,AX_D,0,lw,0);
			mdraw(i,AX_U);
			if (tstyle > 0)
				for (d = ybeg; d <= ybounds.max+FTINY;
						d += ybounds.step) {
					mline(i+TLEN/2,yconv(d),0,lw,0);
					mdraw(i-TLEN/2,yconv(d));
				}
		}
	}
						/* grid */
	if (gparam[GRID].flags & DEFINED)
		lw = varvalue(gparam[GRID].name);
	else
		lw = DEFGRID;
	if (lw-- > 0) {
		for (d = xbeg;
				d <= xbounds.max + FTINY;
				d += xbounds.step) {
			i = xconv(d);
			mline(i,AX_D,2,0,0);
			mdraw(i,AX_U);
		}
		for (d = ybeg;
				d <= ybounds.max + FTINY;
				d += ybounds.step) {
			i = yconv(d);
			mline(AX_L,i,2,0,0);
			mdraw(AX_R,i);
		}
	}
}

void
plotcurves(void)				/* plot the curves */
{
	int  i, j, k;

	for (i = 0; i < MAXCUR; i++) {
		nplottable = nplotted = 0;
		lastx = FHUGE;
		if (mgcurve(i, nextpoint) > 0 &&
				cparam[i][CLABEL].flags & DEFINED) {
			j = (LE_U-LE_D)/MAXCUR;
			k = LE_U - i*j;
			mtext(LE_L+(LE_R-LE_L)/8,k+j/3,
					cparam[i][CLABEL].v.s,CPI,0);
			cmsymbol(i,LE_L,k);
			if (cmline(i,LE_L,k) == 0)
			    mdraw(LE_R-(LE_R-LE_L)/4,k);
		}
	}
}

void
nextpoint(			/* plot the next point for c */
	register int  c,
	double  x,
	double  y
)
{
	if (inbounds(x, y)) {

		if (!(cparam[c][CNPOINTS].flags & DEFINED) ||
				nplotted * npltbl[c] <= nplottable *
				(int)varvalue(cparam[c][CNPOINTS].name) ) {
			csymbol(c, x, y);
			nplotted++;
		}
		nplottable++;
		if (lastx != FHUGE)
			climline(c, x, y, lastx, lasty);

	} else if (inbounds(lastx, lasty)) {

		climline(c, lastx, lasty, x, y);

	}
	lastx = x;
	lasty = y;
}


void
stretchbounds(			/* stretch our boundaries */
	int  c,
	double  x,
	double  y
)
{
	if (gparam[XMIN].flags & DEFINED &&
			x < xbounds.min)
		return;
	if (gparam[XMAX].flags & DEFINED &&
			x > xbounds.max)
		return;
	if (gparam[YMIN].flags & DEFINED &&
			y < ybounds.min)
		return;
	if (gparam[YMAX].flags & DEFINED &&
			y > ybounds.max)
		return;

	if (x < xbounds.min)
		xbounds.min = x;
	if (x > xbounds.max)
		xbounds.max = x;
	if (y < ybounds.min)
		ybounds.min = y;
	if (y > ybounds.max)
		ybounds.max = y;
		
	npltbl[c]++;
}


#define  exp10(x)	exp((x)*2.3025850929940456)

double
goodstep(		/* determine a good step for the interval */
	double  interval
)
{
	static int  steps[] = {50, 20, 10, 5, 2, 1};
	double  fact, exp(), log10(), floor();
	int  i;

	if (interval <= FTINY)
		return(1.0);
	fact = exp10(floor(log10(interval)))/10;
	interval /= fact * MINDIVS;
	for (i = 0; interval < steps[i]; i++)
		;
	return(steps[i] * fact);
}

#undef  exp10


int
xconv(			/* convert x to meta coords */
	double  x
)
{
	if (gparam[XMAP].flags & DEFINED)
		x = funvalue(gparam[XMAP].name, 1, &x);
	x = (x - axbegin)/axsize;
	return( AX_L + (int)(x*(AX_R-AX_L)) );
}


int
yconv(			/* convert y to meta coords */
	double  y
)
{
	if (gparam[YMAP].flags & DEFINED)
		y = funvalue(gparam[YMAP].name, 1, &y);
	y = (y - aybegin)/aysize;
	return( AX_D + (int)(y*(AX_U-AX_D)) );
}


void
pconv(		/* convert theta and radius to meta coords */
	int  *xp,
	int  *yp,
	double  t,
	double  r
)
{
	t *= (2.*PI)/period;
	r = rconv(r);
	*xp = r*cos(t) + (PL_X+.5);
	*yp = r*sin(t) + (PL_Y+.5);
}


int
rconv(			/* convert radius to meta coords */
	double  r
)
{
	if (gparam[YMAP].flags & DEFINED)
		r = funvalue(gparam[YMAP].name, 1, &r);

	return((r - aybegin)*PL_R/aysize + .5);
}
	

void
boxstring(	/* put string in box */
	int  xmin,
	int ymin,
	int xmax,
	int ymax,
	char  *s,
	int  d,
	int width,
	int color
)
{
	register long  size;

	if (d == 'u' || d == 'd') {		/* up or down */
		size = strlen(s)*(xmax-xmin)/ASPECT;
		size -= ymax-ymin;
		size /= 2;
		if (size < 0) {				/* center */
			ymin -= size;
			ymax += size;
		}
	} else {				/* left or right */
		size = strlen(s)*(ymax-ymin)/ASPECT;
		size -= xmax-xmin;
		size /= 2;
		if (size < 0) {				/* center */
			xmin -= size;
			xmax += size;
		}
	}
	mvstr(xmin, ymin, xmax, ymax, s, d, width, color);	/* print */
}


char *
goodformat(			/* return a suitable format string for d */
	double  d
)
{
	static char  *f[5] = {"%.0f", "%.1f", "%.2f", "%.3f", "%.4f"};
	double  floor();
	register int  i;

	if (d < 0.0)
		d = -d;
	if (d > 1e-4 && d < 1e6)
		for (i = 0; i < 5; i++) {
			if (d - floor(d+FTINY) <= FTINY)
				return(f[i]);
			d *= 10.0;
		}
	return("%.1e");
}


void
drawcircle(	/* draw a circle */
	int  x,
	int y,
	int r,
	int  typ,
	int wid,
	int col
)
{
	double  d;

	if (r <= 0)
		return;
	mline(x+r, y, typ, wid, col);
	for (d = 2*PI*PL_F; d <= 2*PI+FTINY; d += 2*PI*PL_F)
		mdraw((int)(x+r*cos(d)+.5), (int)(y+r*sin(d)+.5));
}


void
climline(	/* print line from/to out of bounds */
	int  c,
	double  x,
	double y,
	double xout,
	double yout
)
{
	for ( ; ; )
		if (xout < xbounds.min) {
			yout = y + (yout - y)*(xbounds.min - x)/(xout - x);
			xout = xbounds.min;
		} else if (yout < ybounds.min) {
			xout = x + (xout - x)*(ybounds.min - y)/(yout - y);
			yout = ybounds.min;
		} else if (xout > xbounds.max) {
			yout = y + (yout - y)*(xbounds.max - x)/(xout - x);
			xout = xbounds.max;
		} else if (yout > ybounds.max) {
			xout = x + (xout - x)*(ybounds.max - y)/(yout - y);
			yout = ybounds.max;
		} else {
			cline(c, x, y, xout, yout);
			break;
		}        
}


void
cline(		/* print a curve line */
	int  c,
	double  u1,
	double v1,
	double u2,
	double v2
)
{
	int  x, y;
	double  ustep, vstep;
	
	if (period > FTINY) {		/* polar */
		if (u1 > u2) {
			ustep = u1; u1 = u2; u2 = ustep;
			vstep = v1; v1 = v2; v2 = vstep;
		}
		pconv(&x, &y, u1, v1);
		if (cmline(c, x, y) < 0)
			return;
		ustep = period*PL_F;
		if (u2-u1 > ustep) {
			vstep = ustep*(v2-v1)/(u2-u1);
			while ((u1 += ustep) < u2) {
				v1 += vstep;
				pconv(&x, &y, u1, v1);
				mdraw(x, y);
			}
		}
		pconv(&x, &y, u2, v2);
		mdraw(x, y);
	} else if (cmline(c, xconv(u1), yconv(v1)) == 0)
		mdraw(xconv(u2), yconv(v2));
}


int
cmline(			/* start curve line in meta coords */
	int  c,
	int  x,
	int y
)
{
	int  lw, lt, col;
	register VARIABLE  *cv;

	cv = cparam[c];
	if (cv[CLINTYPE].flags & DEFINED)
		lt = varvalue(cv[CLINTYPE].name);
	else
		lt = DEFLINTYPE;
	if (lt-- <= 0)
		return(-1);
	if (cv[CTHICK].flags & DEFINED)
		lw = varvalue(cv[CTHICK].name);
	else
		lw = DEFTHICK;
	if (lw-- <= 0)
		return(-1);
	if (cv[CCOLOR].flags & DEFINED)
		col = varvalue(cv[CCOLOR].name);
	else
		col = DEFCOLOR;
	if (col-- <= 0)
		return(-1);
	mline(x, y, lt, lw, col);
	return(0);
}


void
csymbol(		/* plot curve symbol */
	int  c,
	double  u,
	double v
)
{
	int  x, y;
	
	if (period > FTINY) {
		pconv(&x, &y, u, v);
		cmsymbol(c, x, y);
	} else
		cmsymbol(c, xconv(u), yconv(v));
}


void
cmsymbol(		/* print curve symbol in meta coords */
	int  c,
	int  x,
	int y
)
{
	int  col, ss;
	register VARIABLE  *cv;

	cv = cparam[c];
	if (!(cv[CSYMTYPE].flags & DEFINED))
		return;
	if (cv[CSYMSIZE].flags & DEFINED)
		ss = varvalue(cv[CSYMSIZE].name);
	else
		ss = DEFSYMSIZE;
	if (ss <= 0)
		return;
	if (cv[CCOLOR].flags & DEFINED)
		col = varvalue(cv[CCOLOR].name);
	else
		col = DEFCOLOR;
	if (col-- <= 0)
		return;
	msegment(x-ss,y-ss,x+ss,y+ss,
			cv[CSYMTYPE].v.s,'r',0,col);
}


int
inbounds(			/* determine if x and y are within gbounds */
	double  x,
	double y
)
{
	if (x < xbounds.min || x > xbounds.max)
		return(0);
	if (y < ybounds.min || y > ybounds.max)
		return(0);
	return(1);
}
