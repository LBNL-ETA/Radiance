/* RCSid $Id: ccolor.h,v 3.3 2012/05/17 17:10:23 greg Exp $ */
/*
 *  Header file for spectral colors.
 *
 */

#ifndef _MGF_COLOR_H_
#define _MGF_COLOR_H_
#ifdef __cplusplus
extern "C" {
#endif

#define C_CMINWL	380		/* minimum wavelength */
#define C_CMAXWL	780		/* maximum wavelength */
#define C_CNSS		41		/* number of spectral samples */
#define C_CWLI		((C_CMAXWL-C_CMINWL)/(C_CNSS-1))
#define C_CMAXV		10000		/* nominal maximum sample value */
#define C_CLPWM		(683./C_CMAXV)	/* peak lumens/watt multiplier */

#define C_CSSPEC	01		/* flag if spectrum is set */
#define C_CDSPEC	02		/* flag if defined w/ spectrum */
#define C_CSXY		04		/* flag if xy is set */
#define C_CDXY		010		/* flag if defined w/ xy */
#define C_CSEFF		020		/* flag if efficacy set */

typedef struct {
	int	clock;		/* incremented each change */
	void	*client_data;	/* pointer to private client-owned data */
	short	flags;		/* what's been set and how */
	short	ssamp[C_CNSS];	/* spectral samples, min wl to max */
	long	ssum;		/* straight sum of spectral values */
	float	cx, cy;		/* xy chromaticity value */
	float	eff;		/* efficacy (lumens/watt) */
} C_COLOR;

#define C_DEFCOLOR	{ 1, NULL, C_CDXY|C_CSXY|C_CSSPEC|C_CSEFF,\
			{C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,\
			C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV,C_CMAXV},\
			(long)C_CNSS*C_CMAXV, 1./3., 1./3., 178.006 }

#define c_cval(c,l)	((double)(c)->ssamp[((l)-C_MINWL)/C_CWLI] / (c)->ssum)

extern C_COLOR		c_dfcolor;		/* default color */

						/* assign arbitrary spectrum */
extern double	c_sset(C_COLOR *clr, double wlmin, double wlmax,
				const float spec[], int nwl);
extern void	c_ccvt(C_COLOR *clr, int fl);	/* fix color representation */
extern int	c_isgrey(C_COLOR *clr);		/* check if color is grey */
						/* mix two colors */
extern void	c_cmix(C_COLOR *cres, double w1, C_COLOR *c1,
				double w2, C_COLOR *c2);
						/* multiply two colors */
extern double	c_cmult(C_COLOR *cres, C_COLOR *c1, double y1,
				C_COLOR *c2, double y2);
						/* set black body spectrum */
extern int	c_bbtemp(C_COLOR *clr, double tk);
						/* convert to RGB color */
extern void	ccy2rgb(C_COLOR *cin, double cieY, float cout[3]);
						/* convert from RGB color */
extern double	rgb2ccy(float cin[3], C_COLOR *cout);

#ifdef __cplusplus
}
#endif
#endif	/* _MGF_COLOR_H_ */
