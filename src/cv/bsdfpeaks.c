#ifndef lint
static const char RCSid[] = "$Id: bsdfpeaks.c,v 2.6 2025/06/07 05:09:45 greg Exp $";
#endif
/*
 *  Compute minimum FWHM peak for each incident direction in SIR input.
 *  Report FWHM of corresponding peaks in XML representations if provided.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bsdfrep.h"

typedef struct {
	float		peakv;		/* peak BSDF value */
	float		width;		/* smallest FWHM (deg) */
	const RBFNODE	*rbs;		/* incident system */
	int		ndx;		/* peak index for RBFVAL */
} FWHM;			/* struct to hold peak value */

typedef double	eval_f(const FVECT vin, const FVECT vout, const void *p);

/* Comparison function to put larger peaks first */
int
cmpFWHM(const void *p0, const void *p1)
{
	float	diff = (*(const FWHM *)p0).peakv -
			(*(const FWHM *)p1).peakv;

	if (diff > 0) return(-1);
	if (diff < 0) return(1);
	return(0);
}

/* BSDF evaluation function for RBF system */
double
rbf_eval(const FVECT vin, const FVECT vout, const void *p)
{
	/* XXX verify vin == p->invec ? */
	return(eval_rbfrep((const RBFNODE *)p, vout));
}

/* BSDF evaluation for XML input */
double
bsdf_eval(const FVECT vin, const FVECT vout, const void *p)
{
	SDValue	sv;

	if (SDreportError(
			SDevalBSDF(&sv, vin, vout, (const SDData *)p),
			stderr))
		exit(1);

	return(sv.cieY);
}

/* Find full-width, half-maximum in radians around BSDF direction */
double
getFWHM(const FVECT vin, const FVECT vc, double rad0, eval_f *ev, const void *p)
{
	const double	peakv = (*ev)(vin, vc, p);
	double		rad1 = rad0;			/* current radii */

	while (rad0 < M_PI/2.) {			/* look for FWHM */
	    FVECT	v0, vt;
	    double	phi;
	    v0[0] = 1; v0[1] = v0[2] = 0;
	    geodesic(v0, vc, v0, rad0, GEOD_RAD);	/* use vc as pivot */
	    for (phi = 0; phi < 2.*M_PI; phi += M_PI/18.) {
	    	spinvector(vt, v0, vc, phi);
	    	if ((*ev)(vin, vt, p) <= .5*peakv) {	/* found one side? */
		    FVECT	vt1;
		    while (rad1 < M_PI/2.) {		/* bracket peak */
			geodesic(vt1, vt, vc, rad0+rad1, GEOD_RAD);
			if ((*ev)(vin, vt1, p) <= .5*peakv)
			    return(rad0+rad1);		/* got both! */
			rad1 *= 1.05;			/* else bump rad1 */
		    }
	    	}
	    }
	    rad1 = rad0 *= 1.05;			/* or expand search */
	}
	return(M_PI);			/* failure return */
}

/* Get outgoing direction for the given FWHM record */
void
getOutDir(FVECT vo, FWHM *dp)
{
	const RBFVAL	*vp = dp->rbs->rbfa + dp->ndx;

	ovec_from_pos(vo, vp->gx, vp->gy);
}

/* Assign FWHM record for specified RBF system */
void
assignFWHM(FWHM *dp, const RBFNODE *rbf)
{
	FVECT	vo;
	int	j;
	double	rad;

	dp->rbs = rbf;
	dp->ndx = 0;		/* find peak outgoing */
	for (j = rbf->nrbf; --j; )
		if (rbf->rbfa[j].peak > rbf->rbfa[dp->ndx].peak)
			dp->ndx = j;
				/* record peak */
	getOutDir(vo, dp);
	dp->peakv = eval_rbfrep(rbf, vo);
				/* get FWHM angle in degrees */
	dp->width = 180./M_PI * getFWHM(rbf->invec, vo,
					R2ANG(rbf->rbfa[dp->ndx].crad),
					rbf_eval, rbf);
}

/* Evaluate FWHM for each incident direction recorded in SIR */
int
main(int argc, char *argv[])
{
	const RBFNODE	*rbf;
	SDData		*sdp;
	FILE		*fp;
	int		ndirs;
	FWHM		*peaka;
	int		i;
						/* set global progname */
	fixargv0(argv[0]);
	if (argc < 2)
		goto userr;

	fp = fopen(argv[1], "rb");		/* load SIR input */
	if (fp == NULL) {
		fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
				progname, argv[1]);
		return(1);
	}
	if (!load_bsdf_rep(fp))
		return(1);
	fclose(fp);
	for (i = 2; i < argc; i++)		/* check/load any XMLs */
		if (SDcacheFile(argv[i]) == NULL)
			return(1);
	ndirs = 0;				/* count input directions */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		if (rbf->nrbf <= 0) {
			ndirs = 0;
			break;
		}
		++ndirs;
	}
	if (!ndirs) {
		fprintf(stderr, "%s: missing/bad RBFs in '%s'\n", progname, argv[1]);
		return(1);
	}
						/* print output header */
	printf("%d incident directions in '%s': %s -> %s\n", ndirs, argv[1],
				input_orient>0 ? "Front" : "Back",
				output_orient>0 ? "Front" : "Back");
	fputs("Incident (theta, phi)\tExiting (theta, phi)\tPeak\tFWHM", stdout);
	for (i = 2; i < argc; i++)
		printf("\t'%s'", argv[i]);
	fputc('\n', stdout);
						/* find SIR peaks */
	peaka = (FWHM *)malloc(sizeof(FWHM)*ndirs);
	if (peaka == NULL) return(1);
	for (i = 0, rbf = dsf_list; i < ndirs; i++, rbf = rbf->next)
		assignFWHM(&peaka[i], rbf);
						/* sort strong to weak */
	qsort(peaka, ndirs, sizeof(FWHM), cmpFWHM);

	for (i = 0; i < ndirs; i++) {		/* report FWHM for each incidence */
		FVECT	vout;
		int	j;
		getOutDir(vout, &peaka[i]);
		printf("%.0f %.0f\t%.0f %.0f", get_theta180(peaka[i].rbs->invec),
				get_phi360(peaka[i].rbs->invec),
				get_theta180(vout), get_phi360(vout));
						/* peak and FWHM from SIR */
		printf("\t%.2e\t%.1f", peaka[i].peakv, peaka[i].width);
						/* FWHM for each XML */
		for (j = 2; j < argc; j++) {
			const SDData	*sd = SDcacheFile(argv[j]);
			double	psa;
			if (SDreportError(
					SDsizeBSDF(&psa, peaka[i].rbs->invec,
						NULL, SDqueryMin, sd),
					stderr))
				return(1);

			printf("\t%.1f", 180./M_PI * getFWHM(peaka[i].rbs->invec,
						vout, sqrt(psa/M_PI),
						bsdf_eval, sd));
			SDfreeCache(sd);
		}
		fputc('\n', stdout);
	}
	/*			we're exiting, anyway...
	SDfreeCache(NULL);
	clear_bsdf_rep();
	*/
	return(0);
userr:
	fprintf(stderr, "Usage: %s bsdf.sir [bsdfrep1.xml ..]\n", progname);
	return(1);
}
