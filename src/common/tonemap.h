/* RCSid $Id: tonemap.h,v 3.17 2005/01/07 20:33:02 greg Exp $ */
/*
 * Header file for tone mapping functions.
 *
 * Include after "color.h"
 */
#ifndef _RAD_TONEMAP_H_
#define _RAD_TONEMAP_H_

#include	"tifftypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/****    Argument Macros    ****/
				/* Flags of what to do */
#define	TM_F_HCONTR	01		/* human contrast sensitivity */
#define	TM_F_MESOPIC	02		/* mesopic color sensitivity */
#define	TM_F_LINEAR	04		/* linear brightness mapping */
#define	TM_F_ACUITY	010		/* acuity adjustment */
#define	TM_F_VEIL	020		/* veiling glare */
#define	TM_F_CWEIGHT	040		/* center weighting */
#define	TM_F_FOVEAL	0100		/* use foveal sample size */
#define	TM_F_BW		0200		/* luminance only */
#define	TM_F_NOSTDERR	0400		/* don't report errors to stderr */
				/* combined modes */
#define	TM_F_CAMERA	0
#define	TM_F_HUMAN	(TM_F_HCONTR|TM_F_MESOPIC|TM_F_VEIL|\
				TM_F_ACUITY|TM_F_FOVEAL)
#define	TM_F_UNIMPL	(TM_F_CWEIGHT|TM_F_VEIL|TM_F_ACUITY|TM_F_FOVEAL)

				/* special pointer values */
#define	TM_XYZPRIM	(RGBPRIMP)NULL	/* indicate XYZ primaries (Note 1) */
#define	TM_NOCHROM	(BYTE *)NULL	/* indicate no chrominance */
#define	TM_NOCHROMP	(BYTE **)NULL	/* indicate no chrominances */
#define	TM_GETFILE	(FILE *)NULL	/* indicate file must be opened */


/****    Error Return Values    ****/

#define TM_E_OK		0		/* normal return status */
#define	TM_E_NOMEM	1		/* out of memory */
#define	TM_E_ILLEGAL	2		/* illegal argument value */
#define	TM_E_TMINVAL	3		/* no valid tone mapping */
#define	TM_E_TMFAIL	4		/* cannot compute tone mapping */
#define	TM_E_BADFILE	5		/* cannot open or understand file */
#define TM_E_CODERR1	6		/* code consistency error 1 */
#define TM_E_CODERR2	7		/* code consistency error 2 */


/****    Conversion Constants and Table Sizes    ****/

#define	TM_BRTSCALE	128		/* brightness scale factor (integer) */

#define TM_NOBRT	(-1<<15)	/* bogus brightness value */
#define TM_NOLUM	(1e-17)		/* ridiculously small luminance */

#define	TM_MAXPKG	8		/* maximum number of color formats */


/****    Global Data Types and Structures    ****/

#ifndef	MEM_PTR
#define	MEM_PTR		void *
#endif

extern char	*tmErrorMessage[];	/* error messages */
extern int	tmLastError;		/* last error incurred by library */
extern char	*tmLastFunction;	/* error-generating function name */

typedef short	TMbright;		/* encoded luminance type */

				/* basic tone mapping data structure */
typedef struct tmStruct {
	int		flags;		/* flags of what to do */
	RGBPRIMP	monpri;		/* monitor RGB primaries */
	double		mongam;		/* monitor gamma value (approx.) */
	COLOR		clf;		/* computed luminance coefficients */
	int		cdiv[3];	/* computed color divisors */
	RGBPRIMP	inppri;		/* current input primaries */
	double		inpsf;		/* current input scalefactor */
	COLORMAT	cmat;		/* color conversion matrix */
	TMbright	hbrmin, hbrmax;	/* histogram brightness limits */	
	int		*histo;		/* input histogram */
	TMbright	mbrmin, mbrmax;	/* mapped brightness limits */
	unsigned short	*lumap;		/* computed luminance map */
	MEM_PTR		pd[TM_MAXPKG];	/* pointers to private data */
} TMstruct;

				/* conversion package functions */
struct tmPackage {
	MEM_PTR		(*Init)(TMstruct *tms);
	void		(*NewSpace)(TMstruct *tms);
	void		(*Free)(MEM_PTR pp);
};
				/* our list of conversion packages */
extern struct tmPackage	*tmPkg[TM_MAXPKG];
extern int	tmNumPkgs;	/* number of registered packages */


/****    Useful Macros    ****/

				/* compute luminance from encoded value */
#define	tmLuminance(li)	exp((li)*(1./TM_BRTSCALE))

				/* does tone mapping need color matrix? */
#define tmNeedMatrix(t)	((t)->monpri != (t)->inppri)

				/* register a conversion package */
#define tmRegPkg(pf)	( tmNumPkgs >= TM_MAXPKG ? -1 : \
				(tmPkg[tmNumPkgs] = (pf), tmNumPkgs++) )

				/* get the specific package's data */
#define tmPkgData(t,i)	((t)->pd[i]!=NULL ? (t)->pd[i] : (*tmPkg[i]->Init)(t))


/****    Library Function Calls    ****/


extern TMstruct *
tmInit(int flags, RGBPRIMP monpri, double gamval);
/*
	Allocate and initialize new tone mapping.

	flags	-	TM_F_* flags indicating what is to be done.
	monpri	-	display monitor primaries (Note 1).
	gamval	-	display gamma response (can be approximate).

	returns	-	new tone-mapping pointer, or NULL if no memory.
*/

extern int
tmSetSpace(TMstruct *tms, RGBPRIMP pri, double sf);
/*
	Set color primaries and scale factor for incoming scanlines.

	tms	-	tone mapping structure pointer.
	pri	-	RGB color input primaries (Note 1).
	sf	-	scale factor to get to luminance in cd/m^2.

	returns	-	0 on success, TM_E_* code on failure.
*/

extern void
tmClearHisto(TMstruct *tms);
/*
	Clear histogram for current tone mapping.

	tms	-	tone mapping structure pointer.
*/

extern int
tmAddHisto(TMstruct *tms, TMbright *ls, int len, int wt);
/*
	Add brightness values to current histogram.

	tms	-	tone mapping structure pointer.
	ls	-	encoded luminance values.
	len	-	number of luminance values.
	wt	-	integer weight to use for each value (usually 1 or -1).

	returns	-	0 on success, TM_E_* on error.
*/

extern int
tmFixedMapping(TMstruct *tms, double expmult, double gamval);
/*
	Assign a fixed, linear tone-mapping using the given multiplier,
	which is the ratio of maximum output to uncalibrated input.
	This mapping will be used in subsequent calls to tmMapPixels()
	until a new tone mapping is computed.
	Only the min. and max. values are used from the histogram.
	
	tms	-	tone mapping structure pointer.
	expmult	-	the fixed exposure multiplier to use.
	gamval	-	display gamma response (0. for default).
	returns -	0 on success, TM_E_* on error.
*/

extern int
tmComputeMapping(TMstruct *tms, double gamval, double Lddyn, double Ldmax);
/*
	Compute tone mapping function from the current histogram.
	This mapping will be used in subsequent calls to tmMapPixels()
	until a new tone mapping is computed.
	I.e., calls to tmAddHisto() have no immediate effect.

	tms	-	tone mapping structure pointer.
	gamval	-	display gamma response (0. for default).
	Lddyn	-	the display's dynamic range (0. for default).
	Ldmax	-	maximum display luminance in cd/m^2 (0. for default).

	returns	-	0 on success, TM_E_* on failure.
*/

extern int
tmMapPixels(TMstruct *tms, BYTE *ps, TMbright *ls, BYTE *cs, int len);
/*
	Apply tone mapping function to pixel values.

	tms	-	tone mapping structure pointer.
	ps	-	returned pixel values (Note 2).
	ls	-	encoded luminance values.
	cs	-	encoded chrominance values (Note 2).
	len	-	number of pixels.

	returns	-	0 on success, TM_E_* on failure.
*/

extern TMstruct *
tmDup(TMstruct *orig);
/*
	Duplicate the given tone mapping into a new struct.

	orig	-	tone mapping structure to duplicate.
	returns	-	pointer to new struct, or NULL on error.
*/

extern void
tmDone(TMstruct *tms);
/*
	Free data associated with the given tone mapping structure.

	tms	-	tone mapping structure to free.
*/

extern int
tmCvColors(TMstruct *tms, TMbright *ls, BYTE *cs, COLOR *scan, int len);
/*
	Convert RGB/XYZ float scanline to encoded luminance and chrominance.

	tms	-	tone mapping structure pointer.
	ls	-	returned encoded luminance values.
	cs	-	returned encoded chrominance values (Note 2).
	scan	-	input scanline.
	len	-	scanline length.

	returns	-	0 on success, TM_E_* on error.
*/

extern int
tmCvGrays(TMstruct *tms, TMbright *ls, float *scan, int len);
/*
	Convert gray float scanline to encoded luminance.

	tms	-	tone mapping structure pointer.
	ls	-	returned encoded luminance values.
	scan	-	input scanline.
	len	-	scanline length.

	returns	-	0 on success, TM_E_* on error.
*/

extern int
tmCvColrs(TMstruct *tms, TMbright *ls, BYTE *cs, COLR *scan, int len);
/*
	Convert RGBE/XYZE scanline to encoded luminance and chrominance.

	tms	-	tone mapping structure pointer.
	ls	-	returned encoded luminance values.
	cs	-	returned encoded chrominance values (Note 2).
	scan	-	input scanline.
	len	-	scanline length.

	returns	-	0 on success, TM_E_* on error.
*/

extern int
tmLoadPicture(TMstruct *tms, TMbright **lpp, BYTE **cpp, int *xp, int *yp,
		char *fname, FILE *fp);
/*
	Load Radiance picture and convert to tone mapping representation.
	Memory for the luminance and chroma arrays is allocated using
	malloc(3), and should be freed with free(3) when no longer needed.
	Calls tmSetSpace() to calibrate input color space.

	tms	-	tone mapping structure pointer.
	lpp	-	returned array of encoded luminances, picture ordering.
	cpp	-	returned array of encoded chrominances (Note 2).
	xp, yp	-	returned picture dimensions.
	fname	-	picture file name.
	fp	-	pointer to open file (Note 3).

	returns	-	0 on success, TM_E_* on failure.
*/

extern int
tmMapPicture(BYTE **psp, int *xp, int *yp, int flags,
		RGBPRIMP monpri, double gamval, double Lddyn, double Ldmax,
		char *fname, FILE *fp);
/*
	Load and apply tone mapping to Radiance picture.
	If fp is TM_GETFILE and (flags&TM_F_UNIMPL)!=0, tmMapPicture()
	calls pcond to perform the actual conversion, which takes
	longer but gives access to all the TM_F_* features.
	Memory for the final pixel array is allocated using malloc(3),
	and should be freed with free(3) when it is no longer needed.

	psp	-	returned array of tone mapped pixels, picture ordering.
	xp, yp	-	returned picture dimensions.
	flags	-	TM_F_* flags indicating what is to be done.
	monpri	-	display monitor primaries (Note 1).
	gamval	-	display gamma response.
	Lddyn	-	the display's dynamic range (0. for default).
	Ldmax	-	maximum display luminance in cd/m^2 (0. for default).
	fname	-	picture file name.
	fp	-	pointer to open file (Note 3).

	returns	-	0 on success, TM_E_* on failure.
*/

extern int
tmCvRGB48(TMstruct *tms, TMbright *ls, BYTE *cs, uint16 (*scan)[3],
		int len, double gv);
/*
	Convert 48-bit RGB scanline to encoded luminance and chrominance.

	tms	-	tone mapping structure pointer.
	ls	-	returned encoded luminance values.
	cs	-	returned encoded chrominance values (Note 2).
	scan	-	input scanline.
	len	-	scanline length.
	gv	-	input gamma value.

	returns	-	0 on success, TM_E_* on error.
*/

extern int
tmCvGray16(TMstruct *tms, TMbright *ls, uint16 *scan, int len, double gv);
/*
	Convert 16-bit gray scanline to encoded luminance.

	tms	-	tone mapping structure pointer.
	ls	-	returned encoded luminance values.
	scan	-	input scanline.
	len	-	scanline length.
	gv	-	input gamma value.

	returns	-	0 on success, TM_E_* on error.
*/

/****    Notes    ****/
/*
	General:

	The usual sequence after calling tmInit() is to convert some
	pixel values to chroma and luminance encodings, which can
	be passed to tmAddHisto() to put into the tone mapping histogram.
	This histogram is then used along with the display parameters
	by tmComputeMapping() to compute the luminance mapping function.
	(Colors are tone-mapped as they are converted if TM_F_MESOPIC
	is set.)  The encoded chroma and luminance values may then be
	passed to tmMapPixels() to apply the computed tone mapping in
	a second pass.

	Especially for RGBE colors in the same color space as the
	target display, the conversion to chroma and luminance values
	is fast enough that it may be recomputed on the second pass
	if memory is at a premium.  (Since the chroma values are only
	used during final mapping, setting the cs parameter to TM_NOCHROM
	on the first pass will save time.)  Another memory saving option
	if third and subsequent passes are not needed is to use the
	same array to store the mapped pixels as used to store
	the encoded chroma values.  This way, only two extra bytes
	for storing encoded luminances are required per pixel.  This
	is the method employed by tmMapPicture(), for example.
*/
/*
	Note 1:

	All RGBPRIMP values should be pointers to static structures.
	I.e., the same pointer should be used for the same primaries,
	and those primaries should not change from one call to the next.
	This is because the routines use the pointer value to identify
	computed conversion matrices, so the matrices do not have to be
	rebuilt each time.  If you are using the standard primaries,
	use the standard primary pointer, stdprims.

	To indicate a set of input primaries are actually CIE XYZ,
	the TM_XYZPRIM macro may be passed.  If the input primaries
	are unknown, it is wisest to pass tmTop->monpri to tmSetSpace().
	By default, the input primaries equal the monitor primaries
	and the scalefactor is set to WHTEFFICACY in tmInit().
*/
/*
	Note 2:

	Chrominance is optional in most of these routines, whether
	TM_F_BW is set or not.  Passing a value of TM_NOCHROM (or
	TM_NOCHROMP for picture reading) indicates that returned
	color values are not desired.
	
	If TM_NOCHROM is passed to a mapping routine expecting chroma
	input, it will do without it and return luminance channel
	rather than RGB pixel values.  Otherwise, the array for
	returned pixel values may be the same array used to pass
	encoded chrominances if no other mappings are desired.
*/
/*
	Note 3:

	Picture files may either be opened by the calling routine or by
	the library function.  If opened by the calling routine, the
	pointer to the stream is passed, which may be a pipe or other
	non-seekable input as it is only read once.  It must be
	positioned at the beginning when the function is called, and it
	will be positioned at the end on return.  It will not be closed.
	If the passed stream pointer is TM_GETFILE, then the library
	function will open the file, read its contents and close it
	before returning, whether or not an error was encountered.
*/

#ifdef __cplusplus
}
#endif
#endif /* _RAD_TONEMAP_H_ */

