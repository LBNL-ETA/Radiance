#ifndef lint
static const char	RCSid[] = "$Id: dircode.c,v 2.1 2003/03/04 05:49:21 greg Exp $";
#endif
/*
 * Compute 4-byte direction code (assume this fits into int)
 */

#include "standard.h"

#define	DCSCALE		11585.2		/* (1<<13)*sqrt(2) */
#define FXNEG		01
#define FYNEG		02
#define FZNEG		04
#define F1X		010
#define F2Z		020
#define F1SFT		5
#define F2SFT		18
#define FMASK		0x1fff

int4
encodedir(dv)		/* encode a normalized direction vector */
FVECT	dv;
{
	register int4	dc = 0;
	int	cd[3], cm;
	register int	i;

	for (i = 0; i < 3; i++)
		if (dv[i] < 0.) {
			cd[i] = dv[i] * -DCSCALE;
			dc |= FXNEG<<i;
		} else
			cd[i] = dv[i] * DCSCALE;
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
	return(dc);
}


void
decodedir(dv, dc)	/* decode a normalized direction vector */
register FVECT	dv;	/* returned */
register int4	dc;
{
	double	d1, d2, der;

	d1 = ((dc>>F1SFT & FMASK)+.5)*(1./DCSCALE);
	d2 = ((dc>>F2SFT & FMASK)+.5)*(1./DCSCALE);
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


double
dir2diff(dc1, dc2)		/* approx. radians^2 between directions */
int4	dc1, dc2;
{
	FVECT	v1, v2;

	decodedir(v1, dc1);
	decodedir(v2, dc2);

	return(2. - 2.*DOT(v1,v2));
}


double
fdir2diff(dc1, v2)		/* approx. radians^2 between directions */
int4	dc1;
register FVECT	v2;
{
	FVECT	v1;

	decodedir(v1, dc1);

	return(2. - 2.*DOT(v1,v2));
}
