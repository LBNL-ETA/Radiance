#ifndef lint
static const char	RCSid[] = "$Id: ccyrgb.c,v 3.5 2023/11/15 18:02:52 greg Exp $";
#endif
/*
 * Convert MGF color to Radiance RGB and spectral representations
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

void
ccy2scolor(			/* convert MGF color, Y to spectral color */
	C_COLOR	*cin,	/* input MGF color or spectrum */
	double	cieY,	/* input luminance or reflectance */
	SCOLOR	sco	/* output spectral color */
)
{
	COLORV	nf, sorig[C_CNSS];
	int	i;

	if (cin->flags & C_CDXY || NCSAMP <= 3) {
		COLOR	col;		/* tristimulus in or out, so... */
		ccy2rgb(cin, cieY, col);
		setscolor(sco, colval(col,RED), colval(col,GRN), colval(col,BLU));
		return;
	}
	c_ccvt(cin, C_CSXY);
	nf = cieY*C_CNSS/(cin->ssum*3*cin->cy);
	for (i = C_CNSS; i--; )		/* may as well reverse order */
		sorig[C_CNSS-1-i] = nf*(COLORV)cin->ssamp[i];

	convertscolor(sco, NCSAMP, WLPART[0], WLPART[3],
			sorig, C_CNSS, C_CMAXWL+.5*C_CWLI, C_CMINWL-.5*C_CWLI);
}

/* Convert spectral color to MGF color and Y value */
double
scolor2ccy(SCOLOR sci, C_COLOR *cout)
{
	float		hinc, cspec[MAXCSAMP];
	int		i = NCSAMP;

	if (i <= 3)		/* defined as RGB, so... */
		return(rgb2ccy(sci, cout));

	while (i-- > 0)		/* need to reverse array */
		cspec[i] = sci[NCSAMP-1-i];

	hinc = (WLPART[0]-WLPART[3])/(float)(2*NCSAMP);

	return(c_sset(cout, WLPART[3]+hinc, WLPART[0]-hinc, cspec, NCSAMP));
}
