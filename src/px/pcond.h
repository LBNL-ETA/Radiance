/* RCSid: $Id: pcond.h,v 3.15 2017/11/30 18:43:05 greg Exp $ */
/*
 * Header file for picture file conditioning.
 */
#ifndef _RAD_PCOND_H_
#define _RAD_PCOND_H_

#include "standard.h"
#include "color.h"
#include "view.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef ADJ_VEIL
#define ADJ_VEIL	0		/* adjust veil to preserve contrast? */
#endif

#define SWNORM		2.26		/* scotopic/photopic ratio for white */
#define WHTSEFFICACY	(SWNORM*WHTEFFICACY)

#define BotMesopic	5.62e-3		/* top of scotopic range */
#define TopMesopic	5.62		/* bottom of photopic range */

#define FOVDIA		(1.0*PI/180.0)	/* foveal diameter (radians) */

#define	HISTRES		100		/* histogram resolution */
#define MAXPREHIST	1024		/* maximum precomputed histogram */

#define LMIN		1e-7		/* minimum visible world luminance */
#define LMAX		1e7		/* maximum visible world luminance */

#define Bl(Lw)		log(Lw)		/* brightness function */
#define Bl1(Lw)		(1.0/(Lw))	/* first derivative of Bl(Lw) */
#define Lb(Bw)		exp(Bw)		/* inverse of brightness function */
#define Lb1(Bw)		exp(Bw)		/* first derivative of Lb(Bw) */

				/* Flags of what to do */
#define DO_ACUITY	01
#define DO_VEIL		02
#define DO_HSENS	04
#define DO_COLOR	010
#define DO_CWEIGHT	020
#define DO_FIXHIST	040
#define DO_PREHIST	0100
#define DO_LINEAR	0200

#define DO_HUMAN	(DO_ACUITY|DO_VEIL|DO_HSENS|DO_COLOR)

extern int	what2do;		/* desired adjustments */

extern double	ldmax;			/* maximum output luminance */
extern double	lddyn;			/* display dynamic range */
extern double	Bldmin, Bldmax;		/* Bl(ldmax/lddyn) and Bl(ldmax) */

extern char	*progname;		/* global argv[0] */

extern char	*infn;			/* input file name */
extern FILE	*infp;			/* input stream */
extern double	(*lumf)();		/* input luminance function */
extern double	inpexp;			/* input exposure value */

#define	plum(clr)	((*lumf)(clr,0)/inpexp)
#define slum(clr)	((*lumf)(clr,1)/inpexp)

#define ldmin		(ldmax/lddyn)

extern COLOR	*fovimg;		/* foveal (1 degree) averaged image */
extern int	fvxr, fvyr;		/* foveal image resolution */
extern float	*crfimg;		/* contrast reduction factors */

#define fovscan(y)	(fovimg+(y)*fvxr)
#define crfscan(y)	(crfimg+(y)*fvxr)

extern double	fixfrac;		/* histogram share due to fixations */
extern short	(*fixlst)[2];		/* fixation history list */
extern int	nfixations;		/* number of fixation points */

extern double	bwhist[HISTRES];	/* luminance histogram */
extern double	histot;			/* total count of histogram */
extern double	bwmin, bwmax;		/* histogram limits */
extern double	bwavg;			/* mean brightness */

#define bwhi(B)		(int)(HISTRES*((B)-bwmin)/(bwmax-bwmin))

extern RGBPRIMP	inprims;		/* input primaries */
extern COLORMAT	inrgb2xyz;		/* convert input RGB to XYZ */

extern RGBPRIMP	outprims;		/* output primaries */

extern double	scalef;			/* linear scaling factor */

extern VIEW	ourview;		/* picture view */
extern double	pixaspect;		/* pixel aspect ratio */
extern RESOLU	inpres;			/* input picture resolution */

extern char	*mbcalfile;		/* macbethcal mapping file */
extern char	*cwarpfile;		/* color warp mapping file */



	/* defined in pcond.c */
extern void syserror(char *s);

	/* defined in pcond2.c */
extern double rgblum(COLOR clr, int scotopic);	/* compute (scotopic) luminance of RGB color */
extern double cielum(COLOR xyz, int scotopic);	/* compute (scotopic) luminance of CIE color */
extern COLOR	*nextscan(void);		/* next processed scanline */
extern COLOR	*firstscan(void);		/* first processed scanline */

	/* defined in pcond3.c */
extern void getfixations(FILE *fp);	/* load fixation history list */
extern void gethisto(FILE *fp);		/* load precomputed luminance histogram */
extern void comphist(void);		/* create foveal sampling histogram */
extern double htcontrs(double La);	/* human contrast sens. func. */
extern double clampf(double La);	/* histogram clamping function */
extern double crfactor(double La);	/* contrast reduction factor */
extern int mkbrmap(void);		/* make dynamic range map */
extern void putmapping(FILE	*fp);	/* put out mapping function */
extern void scotscan(COLOR *scan, int xres);	/* apply scotopic color sensitivity loss */
extern void mapscan(COLOR *scan, int xres);	/* apply tone mapping operator to scanline */

	/* defined in pcond4.c */
extern void compveil(void);		/* compute veiling image */
#if ADJ_VEIL
extern void adjveil(void);		/* adjust veil image */
#endif
extern void acuscan(COLOR *scln, int y);	/* get acuity-sampled scanline */
extern void addveil(COLOR *sl, int y);	/* add veil to scanline */
extern void initacuity(void);		/* initialize variable acuity sampling */
extern double hacuity(double La);	/* human acuity func. (cycles/deg.) */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PCOND_H_ */
