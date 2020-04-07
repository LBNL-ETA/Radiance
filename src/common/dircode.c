#ifndef lint
static const char RCSid[] = "$Id: dircode.c,v 2.13 2020/04/07 18:20:52 greg Exp $";
#endif
/*
 * Compute a 4-byte direction code (externals defined in rtmath.h).
 *
 * Mean accuracy is 0.0022 degrees, with a maximum error of 0.0058 degrees.
 */

#include "rtmath.h"

#define	DCSCALE		11584.5		/* ((1<<13)-.5)*sqrt(2) */
#define FXNEG		01
#define FYNEG		02
#define FZNEG		04
#define F1X		010
#define F2Z		020
#define F1SFT		5
#define F2SFT		18
#define FMASK		0x1fff

int32
encodedir(FVECT dv)		/* encode a normalized direction vector */
{
	int32	dc = 0;
	int	cd[3], cm;
	int	i;

	for (i = 0; i < 3; i++)
		if (dv[i] < 0.) {
			cd[i] = (int)(dv[i] * -DCSCALE + .5);
			dc |= FXNEG<<i;
		} else
			cd[i] = (int)(dv[i] * DCSCALE + .5);
	if (!(cd[0] | cd[1] | cd[2]))
		return(0);		/* zero normal */
	if (cd[0] <= cd[1]) {
		dc |= F1X | cd[0] << F1SFT;
		cm = cd[1];
	} else {
		dc |= cd[1] << F1SFT;
		cm = cd[0];
	}
	if (cd[2] <= cm)
		dc |= F2Z | cd[2] << F2SFT;
	else
		dc |= cm << F2SFT;
	if (!dc)	/* don't generate 0 code normally */
		dc = F1X;
	return(dc);
}

#if 0		/* original version for reference */

void
decodedir(FVECT dv, int32 dc)	/* decode a normalized direction vector */
{
	double	d1, d2, der;

	if (!dc) {		/* special code for zero normal */
		dv[0] = dv[1] = dv[2] = 0.;
		return;
	}
	d1 = (dc>>F1SFT & FMASK)*(1./DCSCALE);
	d2 = (dc>>F2SFT & FMASK)*(1./DCSCALE);
	der = sqrt(1. - d1*d1 - d2*d2);
	if (dc & F1X) {
		dv[0] = d1;
		if (dc & F2Z) { dv[1] = der; dv[2] = d2; }
		else { dv[1] = d2; dv[2] = der; }
	} else {
		dv[1] = d1;
		if (dc & F2Z) { dv[0] = der; dv[2] = d2; }
		else { dv[0] = d2; dv[2] = der; }
	}
	if (dc & FXNEG) dv[0] = -dv[0];
	if (dc & FYNEG) dv[1] = -dv[1];
	if (dc & FZNEG) dv[2] = -dv[2];
}

#else	

void
decodedir(FVECT dv, int32 dc)	/* decode a normalized direction vector */
{
	static const short	itab[4][3] = {
					{1,0,2},{0,1,2},{1,2,0},{0,2,1}
				};
	static const RREAL	neg[2] = {1., -1.};
	const int		ndx = ((dc & F2Z) != 0)<<1 | ((dc & F1X) != 0);
	double			d1, d2, der;

	if (!dc) {		/* special code for zero normal */
		dv[0] = dv[1] = dv[2] = 0.;
		return;
	}
	d1 = (dc>>F1SFT & FMASK)*(1./DCSCALE);
	d2 = (dc>>F2SFT & FMASK)*(1./DCSCALE);
	der = sqrt(1. - d1*d1 - d2*d2);
	dv[itab[ndx][0]] = d1;
	dv[itab[ndx][1]] = d2;
	dv[itab[ndx][2]] = der;
	dv[0] *= neg[(dc&FXNEG)!=0];
	dv[1] *= neg[(dc&FYNEG)!=0];
	dv[2] *= neg[(dc&FZNEG)!=0];
}

#endif

double
dir2diff(int32 dc1, int32 dc2)	/* approx. radians^2 between directions */
{
	FVECT	v1, v2;

	if (dc1 == dc2)
		return 0.;

	decodedir(v1, dc1);
	decodedir(v2, dc2);

	return(2. - 2.*DOT(v1,v2));
}

double
fdir2diff(int32 dc1, FVECT v2)	/* approx. radians^2 between directions */
{
	FVECT	v1;

	decodedir(v1, dc1);

	return(2. - 2.*DOT(v1,v2));
}
