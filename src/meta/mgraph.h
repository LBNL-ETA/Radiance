/* RCSid: $Id: mgraph.h,v 1.3 2003/11/15 02:13:37 schorsch Exp $ */
/*
 *  mgraph.h - header file for graphing routines.
 *
 *     6/23/86
 *
 *     Greg Ward Larson
 */
#ifndef _RAD_MGRAPH_H_
#define _RAD_MGRAPH_H_

#ifdef __cplusplus
extern "C" {
#endif

#define  PI		3.14159265358979323846
#define  FHUGE		1e20
#define  FTINY		1e-7
					/* text density */
#define  CPI		10
#define  CWID		(2048/CPI)
					/* character aspect ratio */
#define  ASPECT		1.67
					/* fraction of period for polar arc */
#define  PL_F		0.02
					/* minimum # of axis divisions */
#define  MINDIVS	4
					/* coordinate axis box */
#define  AX_L		2048
#define  AX_R		12288
#define  AX_D		2048
#define  AX_U		12288
					/* polar plot center and radius */
#define  PL_X		6912
#define  PL_Y		7168
#define  PL_R		5120
					/* axis tick length */
#define  TLEN		200
					/* x numbering offsets */
#define  XN_X		0
#define  XN_Y		(AX_D-TLEN-300)
#define  XN_S		1000
					/* y numbering offsets */
#define  YN_X		(AX_L-TLEN-200)
#define  YN_Y		64
#define  YN_S		500
					/* polar numbering offsets */
#define  TN_X		PL_X
#define  TN_Y		(PL_Y+150)
#define  TN_R		(PL_R+TLEN+250)
#define  RN_X		PL_X
#define  RN_Y		(PL_Y-TLEN/2-200)
#define  RN_S		XN_S
					/* title box */
#define  TI_L		0
#define  TI_R		16383
#define  TI_D		14500
#define  TI_U		15300
					/* subtitle box */
#define  ST_L		0
#define  ST_R		16383
#define  ST_D		13675
#define  ST_U		14170
					/* x label box */
#define  XL_L		2048
#define  XL_R		9850
#define  XL_D		500
#define  XL_U		1000
					/* x mapping box */
#define  XM_L		10000
#define  XM_R		12288
#define  XM_D		670
#define  XM_U		1000
					/* y label box */
#define  YL_L		0
#define  YL_R		500
#define  YL_D		2048
#define  YL_U		9850
					/* y mapping box */
#define  YM_L		170
#define  YM_R		500
#define  YM_D		10000
#define  YM_U		12288
					/* legend box */
#define  LE_L		13100
#define  LE_R		16383
#define  LE_D		4000
#define  LE_U		10000
					/* legend title */
#define  LT_X		LE_L
#define  LT_Y		(LE_U+800)
					/* parameter defaults */
#define  DEFFTHICK	3			/* frame thickness */
#define  DEFTSTYLE	1			/* tick mark style */
#define  DEFOTHICK	0			/* origin thickness */
#define  DEFGRID	0			/* default grid */
#define  DEFSYMSIZE	100			/* symbol size */
#define  DEFLINTYPE	1			/* line type */
#define  DEFTHICK	2			/* line thickness */
#define  DEFCOLOR	1			/* color */
#define  DEFPERIOD	0.0			/* period for polar plot */
#define  DEFPLSTEP	(1./12.)		/* default angular step */

/*
 *  Bounds are used to hold the axis specifications.
 */

typedef struct {
	double  min, max, step;
}  BOUNDS;

extern void mgraph(void);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_MGRAPH_H_ */

