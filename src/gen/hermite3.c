#ifndef lint
static const char	RCSid[] = "$Id: hermite3.c,v 2.2 2003/02/22 02:07:23 greg Exp $";
#endif
/*
 *  hermite.c - routines for 3D hermite curves.
 *
 *     10/29/85
 */

#include  <stdio.h>


hermite3(hp, p0, p1, r0, r1, t)		/* compute point on hermite curve */
double  hp[3];		/* returned hermite point */
double  p0[3];		/* first endpoint */
double  p1[3];		/* second endpoint */
double  r0[3];		/* tangent at p0 */
double  r1[3];		/* tangent at p1 */
double  t;		/* position parameter */
{
	register int  i;
	double  tmh[4];
	
	tmh[0] = (2.0*t - 3.0)*t*t + 1.0;
	tmh[1] = (-2.0*t + 3.0)*t*t;
	tmh[2] = ((t - 2.0)*t + 1.0)*t;
	tmh[3] = (t - 1.0)*t*t;
	
	for (i = 0; i < 3; i++)
		hp[i] = p0[i]*tmh[0] + p1[i]*tmh[1] +
			r0[i]*tmh[2] + r1[i]*tmh[3];
}


htan3(ht, p0, p1, r0, r1, t)		/* compute tangent on hermite curve */
double  ht[3];		/* returned hermite tangent */
double  p0[3];		/* first endpoint */
double  p1[3];		/* second endpoint */
double  r0[3];		/* tangent at p0 */
double  r1[3];		/* tangent at p1 */
double  t;		/* position parameter */
{
	register int  i;
	double  tpmh[4];
	
	tpmh[0] = (6.0*t - 6.0)*t;
	tpmh[1] = (-6.0*t + 6.0)*t;
	tpmh[2] = (3.0*t - 4.0)*t + 1.0;
	tpmh[3] = (3.0*t - 2.0)*t;
	
	for (i = 0; i < 3; i++)
		ht[i] = p0[i]*tpmh[0] + p1[i]*tpmh[1] +
			r0[i]*tpmh[2] + r1[i]*tpmh[3];
}
