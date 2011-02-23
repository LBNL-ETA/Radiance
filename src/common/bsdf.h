/* RCSid $Id: bsdf.h,v 2.9 2011/02/23 21:58:31 greg Exp $ */
/*
 *  bsdf.h
 *  
 *  Declarations for bidirectional scattering distribution functions.
 *  Assumes <stdio.h> already included.
 *
 *  A material is oriented in right-hand coordinate system with X-axis
 *	in the surface plane pointed to the right as seen from the front.
 *	This means the Y-axis is "up" and the Z-axis is the surface normal.
 *
 *  Note that we reverse the identification of "front" and "back" from
 *	the conventions used in WINDOW 6.  "Front" in our library points
 *	in the +Z direction, towards the interior of the space rather
 *	than the exterior.
 *
 *  BSDF vectors always oriented away from surface, even when "incident."
 *
 *  Created by Greg Ward on 1/10/11.
 *
 */

#ifndef _BSDF_H_
#define	_BSDF_H_

#include "fvect.h"
#include "ccolor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define	SDnameLn	128		/* maximum BSDF name length */
#define	SDmaxCh		3		/* maximum # spectral channels */

/* Component flags for SDsampBSDF() and SDdirectHemi() */
#define	SDsampR		0x1		/* include reflection */
#define	SDsampT		0x2		/* include transmission */
#define	SDsampS		0x3		/* include scattering (R+T) */
#define	SDsampSp	0x4		/* include non-diffuse portion */
#define	SDsampDf	0x8		/* include diffuse portion */
#define	SDsampSpR	0x5		/* include non-diffuse reflection */
#define	SDsampSpT	0x6		/* include non-diffuse transmission */
#define	SDsampSpS	0x7		/* include non-diffuse scattering */
#define	SDsampAll	0xF		/* include everything */

/* Projected solid angle query flags for SDsizeBSDF() */
#define SDqueryVal	0x0		/* query single value */
#define	SDqueryMin	0x1		/* query minimum proj. solid angle */
#define	SDqueryMax	0x2		/* query maximum proj. solid angle */

/* Error codes: normal return, out of memory, file i/o, file format, bad argument,
		 bad data, unsupported feature, internal error, unknown error */
typedef enum {SDEnone=0, SDEmemory, SDEfile, SDEformat, SDEargument,
		SDEdata, SDEsupport, SDEinternal, SDEunknown} SDError;

/* English ASCII strings corresponding to ennumerated errors */
extern const char	*SDerrorEnglish[];

/* Additional information on last error (ASCII English) */
extern char		SDerrorDetail[];

/* Holder for BSDF value and spectral color */
typedef struct {
	double		cieY;		/* photopic BSDF (Y) value */
	C_COLOR		spec;		/* spectral and (x,y) color */
} SDValue;

/* Cached, encoded, cumulative distribution for one incident (solid) angle */
#define SD_CDIST_BASE	double		cTotal;	\
			struct SDCDst_s	*next
typedef struct SDCDst_s {
	SD_CDIST_BASE;			/* base fields first */
	/* ...encoded distribution extends struct */
} SDCDst;

/* Forward declaration of BSDF component */
typedef struct SDComp_s	SDComponent;

/* Methods needed to handle BSDF components (nothing is optional) */
typedef const struct {
					/* return non-diffuse BSDF */
	int		(*getBSDFs)(float coef[SDmaxCh], const FVECT outVec,
				    const FVECT inVec, const void *dist);
					/* query non-diffuse PSA for vector */
	SDError		(*queryProjSA)(double *psa, const FVECT vec,
						int qflags, const void *dist);
					/* get cumulative distribution */
	const SDCDst	*(*getCDist)(const FVECT inVec, SDComponent *sdc);
					/* sample cumulative distribution */
	SDError		(*sampCDist)(FVECT outVec, double randX,
						const SDCDst *cdp);
					/* free a spectral BSDF component */
	void		(*freeSC)(void *dist);
} SDFunc;

/* Structure to hold a spectral BSDF component (typedef SDComponent above) */
struct SDComp_s {
	C_COLOR		cspec[SDmaxCh];	/* component spectral bases */
	SDFunc		*func;		/* methods for this component */
	void		*dist;		/* loaded distribution data */
	SDCDst		*cdList;	/* cumulative distribution cache */
};

/* Container for non-diffuse BSDF components */
typedef struct {
	double		minProjSA;	/* minimum projected solid angle */
	double		maxHemi;	/* maximum directional hemispherical */
	int		ncomp;		/* number of separate components */
	SDComponent	comp[1];	/* BSDF components (extends struct) */
} SDSpectralDF;

/* Loaded BSDF data */
typedef struct {
	char		name[SDnameLn];	/* BSDF name (derived from file) */
	char		*mgf;		/* geometric description (if any) */
	float		dim[3];		/* width, height, thickness (meters) */
	SDValue		rLambFront;	/* diffuse front reflectance */
	SDValue		rLambBack;	/* diffuse rear reflectance */
	SDValue		tLamb;		/* diffuse transmission */
	SDSpectralDF	*rf, *rb, *tf;	/* non-diffuse BSDF components */
} SDData;

/* List of loaded BSDFs */
extern struct SDCache_s {
	SDData		bsdf;		/* BSDF data */
	unsigned	refcnt;		/* how many callers are using us? */
	struct SDCache_s		/* next in cache list */
			*next;
} *SDcacheList;		/* Global BSDF cache */

/* BSDF cache retention preference */
#define SDretainNone	0		/* free unreferenced data (default) */
#define	SDretainBSDFs	1		/* keep loaded BSDFs in cache */
#define SDretainAll	2		/* also keep cumulative cache data */

extern int		SDretainSet;	/* set to SDretainNone by default */

/*****************************************************************
 * The following routines are less commonly used by applications.
 */

#define	SDisLoaded(sd)	((sd)->rLambFront.spec.flags != 0)

/* Report an error to the indicated stream (in English) */
extern SDError		SDreportEnglish(SDError ec, FILE *fp);

/* Shorten file path to useable BSDF name, removing suffix */
extern void		SDclipName(char res[SDnameLn], const char *fname);

/* Allocate new spectral distribution function */
extern SDSpectralDF	*SDnewSpectralDF(int nc);

/* Free a spectral distribution function */
extern void		SDfreeSpectralDF(SDSpectralDF *df);

/* Initialize an unused BSDF struct (clears to zeroes) */
extern void		SDclearBSDF(SDData *sd);

/* Load a BSDF struct from the given file (keeps name unchanged) */
extern SDError		SDloadFile(SDData *sd, const char *fname);

/* Free data associated with BSDF struct */
extern void		SDfreeBSDF(SDData *sd);

/* Find writeable BSDF by name, or allocate new cache entry if absent */
extern SDData		*SDgetCache(const char *bname);

/* Free cached cumulative distributions for BSDF component */
extern void		SDfreeCumulativeCache(SDSpectralDF *df);

/* Sample an individual BSDF component */
extern SDError		SDsampComponent(SDValue *sv, FVECT outVec,
					const FVECT inVec, double randX,
					SDComponent *sdc);

/* Convert 1-dimensional random variable to N-dimensional */
extern void		SDmultiSamp(double t[], int n, double randX);

/* Map a [0,1]^2 square to a unit radius disk */
void			SDsquare2disk(double ds[2], double seedx, double seedy);

/* Map point on unit disk to a unit square in [0,1]^2 range */
void			SDdisk2square(double sq[2], double diskx, double disky);

/*****************************************************************
 * The calls below are the ones most applications require.
 * All directions are assumed to be unit vectors.
 */

/* Get BSDF from cache (or load and cache it on first call) */
/* Report any problems to stderr and return NULL on failure */
extern const SDData	*SDcacheFile(const char *fname);

/* Free a BSDF from our cache (clear all if NULL) */
extern void		SDfreeCache(const SDData *sd);

/* Query projected solid angle resolution for non-diffuse BSDF direction */
extern SDError		SDsizeBSDF(double *projSA, const FVECT vec,
					int qflags, const SDData *sd);

/* Return BSDF for the given incident and scattered ray vectors */
extern SDError		SDevalBSDF(SDValue *sv, const FVECT outVec,
					const FVECT inVec, const SDData *sd);

/* Compute directional hemispherical scattering at given incident angle */
extern double		SDdirectHemi(const FVECT inVec,
					int sflags, const SDData *sd);

/* Sample BSDF direction based on the given random variable */
extern SDError		SDsampBSDF(SDValue *sv, FVECT outVec,
					const FVECT inVec, double randX,
					int sflags, const SDData *sd);

/*****************************************************************
 * Vector math for getting between world and local BSDF coordinates.
 * Directions may be passed unnormalized to these routines.
 */

/* Compute World->BSDF transform from surface normal and up (Y) vector */
extern SDError		SDcompXform(RREAL vMtx[3][3], const FVECT sNrm,
					const FVECT uVec);

/* Compute inverse transform */
extern SDError		SDinvXform(RREAL iMtx[3][3], RREAL vMtx[3][3]);

/* Transform and normalize direction (column) vector */
extern SDError		SDmapDir(FVECT resVec, RREAL vMtx[3][3],
					const FVECT inpVec);

/* System-specific BSDF loading routine (not part of our library) */
extern SDData		*loadBSDF(char *name);

/* System-specific BSDF error translator (not part of our library) */
extern char		*transSDError(SDError ec);

/*################################################################*/
/*######### DEPRECATED DEFINITIONS AWAITING PERMANENT REMOVAL #######*/
/*
 * Header for BSDF i/o and access routines
 */

#include "mat4.h"
				/* up directions */
typedef enum {
	UDzneg=-3,
	UDyneg=-2,
	UDxneg=-1,
	UDunknown=0,
	UDxpos=1,
	UDypos=2,
	UDzpos=3
} UpDir;
				/* BSDF coordinate calculation routines */
				/* vectors always point away from surface */

typedef int	b_vecf2(FVECT v, int n, void *cd);
typedef int	b_ndxf2(FVECT v, void *cd);
typedef double	b_radf2(int n, void *cd);

/* Bidirectional Scattering Distribution Function */
struct BSDF_data {
	int	ninc;			/* number of incoming directions */
	int	nout;			/* number of outgoing directions */
	float	dim[3];			/* width, height, thickness (meters) */
	char	*mgf;			/* geometric description (if any) */
	void	*ib_priv;		/* input basis private data */
	b_vecf2	*ib_vec;		/* get input vector from index */
	b_ndxf2	*ib_ndx;		/* get input index from vector */
	b_radf2	*ib_ohm;		/* get input radius for index */
	void	*ob_priv;		/* output basis private data */
	b_vecf2	*ob_vec;		/* get output vector from index */
	b_ndxf2	*ob_ndx;		/* get output index from vector */
	b_radf2	*ob_ohm;		/* get output radius for index */
	float	*bsdf;			/* scattering distribution data */
};				/* bidirectional scattering distrib. func. */

#define getBSDF_incvec(v,b,i)	(*(b)->ib_vec)(v,i,(b)->ib_priv)
#define getBSDF_incndx(b,v)	(*(b)->ib_ndx)(v,(b)->ib_priv)
#define getBSDF_incohm(b,i)	(*(b)->ib_ohm)(i,(b)->ib_priv)
#define getBSDF_outvec(v,b,o)	(*(b)->ob_vec)(v,o,(b)->ob_priv)
#define getBSDF_outndx(b,v)	(*(b)->ob_ndx)(v,(b)->ob_priv)
#define getBSDF_outohm(b,o)	(*(b)->ob_ohm)(o,(b)->ob_priv)
#define BSDF_value(b,i,o)	(b)->bsdf[(o)*(b)->ninc + (i)]

extern struct BSDF_data *load_BSDF(char *fname);
extern void free_BSDF(struct BSDF_data *b);
extern int r_BSDF_incvec(FVECT v, struct BSDF_data *b, int i,
				double rv, MAT4 xm);
extern int r_BSDF_outvec(FVECT v, struct BSDF_data *b, int o,
				double rv, MAT4 xm);
extern int getBSDF_xfm(MAT4 xm, FVECT nrm, UpDir ud, char *xfbuf);

/*######### END DEPRECATED DEFINITIONS #######*/
/*################################################################*/

#ifdef __cplusplus
}
#endif
#endif	/* ! _BSDF_H_ */
