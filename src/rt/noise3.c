#ifndef lint
static const char	RCSid[] = "$Id: noise3.c,v 2.11 2010/09/03 21:16:50 greg Exp $";
#endif
/*
 *  noise3.c - noise functions for random textures.
 *
 *     Credit for the smooth algorithm goes to Ken Perlin.
 *     (ref. SIGGRAPH Vol 19, No 3, pp 287-96)
 */

#include "copyright.h"

#include  <math.h>

#include  "calcomp.h"
#include  "func.h"

#define  A		0
#define  B		1
#define  C		2
#define  D		3

#define  rand3a(x,y,z)	frand(67*(x)+59*(y)+71*(z))
#define  rand3b(x,y,z)	frand(73*(x)+79*(y)+83*(z))
#define  rand3c(x,y,z)	frand(89*(x)+97*(y)+101*(z))
#define  rand3d(x,y,z)	frand(103*(x)+107*(y)+109*(z))

#define  hpoly1(t)	((2.0*t-3.0)*t*t+1.0)
#define  hpoly2(t)	(-2.0*t+3.0)*t*t
#define  hpoly3(t)	((t-2.0)*t+1.0)*t
#define  hpoly4(t)	(t-1.0)*t*t

#define  hermite(p0,p1,r0,r1,t)  (	p0*hpoly1(t) + \
					p1*hpoly2(t) + \
					r0*hpoly3(t) + \
					r1*hpoly4(t) )

static char  noise_name[4][8] = {"noise3x", "noise3y", "noise3z", "noise3"};
static char  fnoise_name[] = "fnoise3";
static char  hermite_name[] = "hermite";

static long  xlim[3][2];
static double  xarg[3];

#define  EPSILON	.001		/* error allowed in fractal */

#define  frand3(x,y,z)	frand(17*(x)+23*(y)+29*(z))

static double l_noise3(char  *nam);
static double l_hermite(char *nm);
static double * noise3(double  xnew[3]);
static void interpolate(double  f[4], int  i, int  n);
static double frand(long  s);
static double fnoise3(double  p[3]);


static double
l_noise3(			/* compute a noise function */
	register char  *nam
)
{
	register int  i;
	double  x[3];
					/* get point */
	x[0] = argument(1);
	x[1] = argument(2);
	x[2] = argument(3);
					/* make appropriate call */
	if (nam == fnoise_name)
		return(fnoise3(x));
	i = 4;
	while (i--)
		if (nam == noise_name[i])
			return(noise3(x)[i]);
	eputs(nam);
	eputs(": called l_noise3!\n");
	quit(1);
	return 1; /* pro forma return */
}


static double
l_hermite(char *nm)		/* library call for hermite interpolation */
{
	double  t;
	
	t = argument(5);
	return( hermite(argument(1), argument(2),
			argument(3), argument(4), t) );
}


extern void
setnoisefuncs(void)			/* add noise functions to library */
{
	register int  i;

	funset(hermite_name, 5, ':', l_hermite);
	funset(fnoise_name, 3, ':', l_noise3);
	i = 4;
	while (i--)
		funset(noise_name[i], 3, ':', l_noise3);
}


static double *
noise3(			/* compute the noise function */
	register double  xnew[3]
)
{
	static double  x[3] = {-100000.0, -100000.0, -100000.0};
	static double  f[4];

	if (x[0]==xnew[0] && x[1]==xnew[1] && x[2]==xnew[2])
		return(f);
	x[0] = xnew[0]; x[1] = xnew[1]; x[2] = xnew[2];
	xlim[0][0] = floor(x[0]); xlim[0][1] = xlim[0][0] + 1;
	xlim[1][0] = floor(x[1]); xlim[1][1] = xlim[1][0] + 1;
	xlim[2][0] = floor(x[2]); xlim[2][1] = xlim[2][0] + 1;
	xarg[0] = x[0] - xlim[0][0];
	xarg[1] = x[1] - xlim[1][0];
	xarg[2] = x[2] - xlim[2][0];
	interpolate(f, 0, 3);
	return(f);
}


static void
interpolate(
	double  f[4],
	register int  i,
	register int  n
)
{
	double  f0[4], f1[4], hp1, hp2;

	if (n == 0) {
		f[A] = rand3a(xlim[0][i&1],xlim[1][i>>1&1],xlim[2][i>>2]);
		f[B] = rand3b(xlim[0][i&1],xlim[1][i>>1&1],xlim[2][i>>2]);
		f[C] = rand3c(xlim[0][i&1],xlim[1][i>>1&1],xlim[2][i>>2]);
		f[D] = rand3d(xlim[0][i&1],xlim[1][i>>1&1],xlim[2][i>>2]);
	} else {
		n--;
		interpolate(f0, i, n);
		interpolate(f1, i | 1<<n, n);
		hp1 = hpoly1(xarg[n]); hp2 = hpoly2(xarg[n]);
		f[A] = f0[A]*hp1 + f1[A]*hp2;
		f[B] = f0[B]*hp1 + f1[B]*hp2;
		f[C] = f0[C]*hp1 + f1[C]*hp2;
		f[D] = f0[D]*hp1 + f1[D]*hp2 +
				f0[n]*hpoly3(xarg[n]) + f1[n]*hpoly4(xarg[n]);
	}
}


static double
frand(			/* get random number from seed */
	register long  s
)
{
	s = s<<13 ^ s;
	return(1.0-((s*(s*s*15731+789221)+1376312589)&0x7fffffff)/1073741824.0);
}


static double
fnoise3(			/* compute fractal noise function */
	double  p[3]
)
{
	long  t[3], v[3], beg[3];
	double  fval[8], fc;
	int  branch;
	register long  s;
	register int  i, j;
						/* get starting cube */
	s = (long)(1.0/EPSILON);
	for (i = 0; i < 3; i++) {
		t[i] = s*p[i];
		beg[i] = s*floor(p[i]);
	}
	for (j = 0; j < 8; j++) {
		for (i = 0; i < 3; i++) {
			v[i] = beg[i];
			if (j & 1<<i)
				v[i] += s;
		}
		fval[j] = frand3(v[0],v[1],v[2]);
	}
						/* compute fractal */
	for ( ; ; ) {
		fc = 0.0;
		for (j = 0; j < 8; j++)
			fc += fval[j];
		fc *= 0.125;
		if ((s >>= 1) == 0)
			return(fc);		/* close enough */
		branch = 0;
		for (i = 0; i < 3; i++) {	/* do center */
			v[i] = beg[i] + s;
			if (t[i] > v[i]) {
				branch |= 1<<i;
			}
		}
		fc += s*EPSILON*frand3(v[0],v[1],v[2]);
		fval[~branch & 7] = fc;
		for (i = 0; i < 3; i++) {	/* do faces */
			if (branch & 1<<i)
				v[i] += s;
			else
				v[i] -= s;
			fc = 0.0;
			for (j = 0; j < 8; j++)
				if (~(j^branch) & 1<<i)
					fc += fval[j];
			fc = 0.25*fc + s*EPSILON*frand3(v[0],v[1],v[2]);
			fval[~(branch^1<<i) & 7] = fc;
			v[i] = beg[i] + s;
		}
		for (i = 0; i < 3; i++) {	/* do edges */
			if ((j = i+1) == 3) j = 0;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			if (++j == 3) j = 0;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			fc = fval[branch & ~(1<<i)];
			fc += fval[branch | 1<<i];
			fc = 0.5*fc + s*EPSILON*frand3(v[0],v[1],v[2]);
			fval[branch^1<<i] = fc;
			if ((j = i+1) == 3) j = 0;
			v[j] = beg[j] + s;
			if (++j == 3) j = 0;
			v[j] = beg[j] + s;
		}
		for (i = 0; i < 3; i++)		/* new cube */
			if (branch & 1<<i)
				beg[i] += s;
	}
}
