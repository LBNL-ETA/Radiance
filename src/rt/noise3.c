#ifndef lint
static const char	RCSid[] = "$Id: noise3.c,v 2.6 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 *  noise3.c - noise functions for random textures.
 *
 *     Credit for the smooth algorithm goes to Ken Perlin.
 *     (ref. SIGGRAPH Vol 19, No 3, pp 287-96)
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include  <math.h>

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

double  *noise3(), fnoise3(), argument(), frand();
static  interpolate();

static long  xlim[3][2];
static double  xarg[3];

#define  EPSILON	.001		/* error allowed in fractal */

#define  frand3(x,y,z)	frand(17*(x)+23*(y)+29*(z))


static double
l_noise3(nam)			/* compute a noise function */
register char  *nam;
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
}


double
l_hermite()			/* library call for hermite interpolation */
{
	double  t;
	
	t = argument(5);
	return( hermite(argument(1), argument(2),
			argument(3), argument(4), t) );
}


setnoisefuncs()			/* add noise functions to library */
{
	register int  i;

	funset(hermite_name, 5, ':', l_hermite);
	funset(fnoise_name, 3, ':', l_noise3);
	i = 4;
	while (i--)
		funset(noise_name[i], 3, ':', l_noise3);
}


double *
noise3(xnew)			/* compute the noise function */
register double  xnew[3];
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


static
interpolate(f, i, n)
double  f[4];
register int  i, n;
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


double
frand(s)			/* get random number from seed */
register long  s;
{
	s = s<<13 ^ s;
	return(1.0-((s*(s*s*15731+789221)+1376312589)&0x7fffffff)/1073741824.0);
}


double
fnoise3(p)			/* compute fractal noise function */
double  p[3];
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
			j = (i+1)%3;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			j = (i+2)%3;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			fc = fval[branch & ~(1<<i)];
			fc += fval[branch | 1<<i];
			fc = 0.5*fc + s*EPSILON*frand3(v[0],v[1],v[2]);
			fval[branch^1<<i] = fc;
			j = (i+1)%3;
			v[j] = beg[j] + s;
			j = (i+2)%3;
			v[j] = beg[j] + s;
		}
		for (i = 0; i < 3; i++)		/* new cube */
			if (branch & 1<<i)
				beg[i] += s;
	}
}
