/* RCSid $Id: falsecolor.h,v 3.2 2005/11/15 06:52:38 greg Exp $ */
/*
 * Header file for false color tone-mapping.
 *
 * Include after "color.h" and "tonemap.h"
 */
#ifndef _RAD_FALSECOLOR_H_
#define _RAD_FALSECOLOR_H_

#ifdef __cplusplus
extern "C" {
#endif
				/* false color mapping data structure */
typedef struct {
	TMbright	mbrmin, mbrmax;	/* mapped min. & max. brightnesses */
	BYTE		*lumap;		/* false color luminance map */
	BYTE		(*scale)[3];	/* false color ordinal scale */
} FCstruct;

extern BYTE	fcDefaultScale[256][3];		/* default false color scale */

extern FCstruct *
fcInit(BYTE fcscale[256][3]);
/*
	Allocate and initialize false color mapping data structure.
	
	fcscale	-	false color ordinal scale.
	
	returns	-	new false color mapping pointer, or NULL if no memory.
 */

extern int
fcFixedLinear(FCstruct *fcs, double Lwmax);
/*
	Assign fixed linear false color map.
	
	fcs	-	false color structure pointer.
	Lwmax	-	maximum luminance for scaling
	
	returns	-	0 on success, TM_E_* on failure.
 */

extern int
fcFixedLog(FCstruct *fcs, double Lwmin, double Lwmax);
/*
	Assign fixed logarithmic false color map.
	
	fcs	-	false color structure pointer.
	Lwmin	-	minimum luminance for scaling
	Lwmax	-	maximum luminance for scaling
	
	returns	-	0 on success, TM_E_* on failure.
 */

extern int
fcLinearMapping(FCstruct *fcs, TMstruct *tms, int pctile);
/*
	Compute linear false color map.
	
	fcs	-	false color structure pointer.
	tms	-	tone mapping structure pointer, with histogram.
	pctile	-	percentile to ignore on top
	
	returns	-	0 on success, TM_E_* on failure.
 */

extern int
fcLogMapping(FCstruct *fcs, TMstruct *tms, int pctile);
/*
	Compute logarithmic false color map.
	
	fcs	-	false color structure pointer.
	tms	-	tone mapping structure pointer, with histogram.
	pctile	-	percentile to ignore on top and bottom
	
	returns	-	0 on success, TM_E_* on failure.
 */
 
extern int
fcMapPixels(FCstruct *fcs, BYTE *ps, TMbright *ls, int len);
/*
	Apply false color mapping to pixel values.
	
	fcs	-	false color structure pointer.
	ps	-	returned RGB pixel values.
	ls	-	encoded luminance values.
	len	-	number of pixels.
	
	returns	-	0 on success, TM_E_* on failure.
 */

extern int
fcIsLogMap(FCstruct *fcs);
/*
	Determine if false color mapping is logarithmic.

	fcs     -       false color structure pointer.

	returns	-	1 if map follows logarithmic mapping, -1 on error.
*/

extern FCstruct *
fcDup(FCstruct *fcs);
/*
	Duplicate a false color structure.

	fcs     -       false color structure pointer.

	returns	-	duplicate structure, or NULL if no memory.
*/

extern void
fcDone(FCstruct *fcs);
/*
	Free data associated with the given false color mapping structure.

	fcs	-	false color mapping structure to free.
*/

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_FALSECOLOR_H_ */
