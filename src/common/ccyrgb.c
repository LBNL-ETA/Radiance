#ifndef lint
static const char	RCSid[] = "$Id: ccyrgb.c,v 3.2 2011/06/28 21:11:04 greg Exp $";
#endif
/*
 * Convert MGF color to RGB representation
 */

#include <stdio.h>
#include "color.h"
#include "ccolor.h"

void
ccy2rgb(			/* convert MGF color to RGB */
	C_COLOR	*cin,	/* input MGF chrominance */
	double	cieY,	/* input luminance or reflectance */
	COLOR	cout	/* output RGB color */
)
{
	double	d;
	COLOR	xyz;
					/* get CIE XYZ representation */
	c_ccvt(cin, C_CSXY);
	d = cin->cx/cin->cy;
	xyz[CIEX] = d * cieY;
	xyz[CIEY] = cieY;
	xyz[CIEZ] = (1./cin->cy - d - 1.) * cieY;
	cie_rgb(cout, xyz);
}

/* Convert RGB to MGF color and value */
double
rgb2ccy(COLOR cin, C_COLOR *cout)
{
	COLOR	xyz;
	double	df;

	rgb_cie(xyz, cin);
	*cout = c_dfcolor;
	df = xyz[CIEX] + xyz[CIEY] + xyz[CIEZ];
	if (df <= .0)
		return(.0);
	df = 1./df;
	cout->cx = xyz[CIEX]*df;
	cout->cy = xyz[CIEZ]*df;
	cout->flags = (C_CSXY|C_CDXY);

	return(xyz[CIEY]);
}
