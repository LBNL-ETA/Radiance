#ifndef lint
static const char RCSid[] = "$Id: checkBSDF.c,v 2.2 2021/12/15 02:13:27 greg Exp $";
#endif
/*
 *  checkBSDF.c
 *
 *  Load BSDF XML file and check Helmholtz reciprocity
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "rtio.h"
#include "random.h"
#include "bsdf.h"
#include "bsdf_m.h"
#include "bsdf_t.h"

#define	F_IN_COLOR	0x1
#define	F_ISOTROPIC	0x2
#define F_MATRIX	0x4
#define F_TTREE		0x8

typedef struct {
	double	vmin, vmax;	/* extrema */
	double	vsum;		/* straight sum */
	long	nvals;		/* number of values */
} SimpleStats;

const SimpleStats	SSinit = {FHUGE, -FHUGE, .0, 0};

				/* relative difference formula */
#define rdiff(a,b)	((a)>(b) ? ((a)-(b))/(a) : ((b)-(a))/(b))

/* Figure out BSDF type (and optionally determine if in color) */
const char *
getBSDFtype(const SDData *bsdf, int *flags)
{
	const SDSpectralDF	*df = bsdf->tb;

	if (flags) *flags = 0;
	if (!df) df = bsdf->tf;
	if (!df) df = bsdf->rf;
	if (!df) df = bsdf->rb;
	if (!df) return "Pure_Lambertian";
	if (df->comp[0].func == &SDhandleMtx) {
		const SDMat	*m = (const SDMat *)df->comp[0].dist;
		if (flags) {
			*flags |= F_MATRIX;
			*flags |= F_IN_COLOR*(m->chroma != NULL);
		}
		switch (m->ninc) {
		case 145:
			return "Klems_Full";
		case  73:
			return "Klems_Half";
		case  41:
			return "Klems_Quarter";
		}
		return "Unknown_Matrix";
	}
	if (df->comp[0].func == &SDhandleTre) {
		const SDTre	*t = (const SDTre *)df->comp[0].dist;
		if (flags) {
			*flags |= F_TTREE;
			*flags |= F_IN_COLOR*(t->stc[1] != NULL);
		}
		switch (t->stc[0]->ndim) {
		case 4:
			return "Anisotropic_Tensor_Tree";
		case 3:
			if (flags) *flags |= F_ISOTROPIC;
			return "Isotropic_Tensor_Tree";
		}
		return "Unknown_Tensor_Tree";
	}
	return "Unknown";
}

/* Report details related to one hemisphere distribution */
void
detailComponent(const char *nm, const SDValue *lamb, const SDSpectralDF *df)
{
	printf("%s\t%4.1f %4.1f %4.1f\t\t", nm,
			100.*lamb->cieY*lamb->spec.cx/lamb->spec.cy,
			100.*lamb->cieY,
			100.*lamb->cieY*(1.f - lamb->spec.cx - lamb->spec.cy)/lamb->spec.cy);
	if (df)
		printf("%5.1f%%\t\t%.2f deg\n", 100.*df->maxHemi,
				sqrt(df->minProjSA/M_PI)*(360./M_PI));
	else
		puts("0%\t\t180");
}

/* Add a value to stats */
void
addStat(SimpleStats *ssp, double v)
{
	if (v < ssp->vmin) ssp->vmin = v;
	if (v > ssp->vmax) ssp->vmax = v;
	ssp->vsum += v;
	ssp->nvals++;
}

/* Sample a BSDF hemisphere with callback (quadtree recursion) */
int
qtSampBSDF(double xleft, double ytop, double siz,
		const SDData *bsdf, const int side, const RREAL *v0,
		int (*cf)(const SDData *b, const FVECT v1, const RREAL *v0, void *p),
				void *cdata)
{
	if (siz < 0.124) {			/* make sure we subdivide, first */
		FVECT	vsmp;
		double	sa;
		square2disk(vsmp, xleft + frandom()*siz, ytop + frandom()*siz);
		vsmp[2] = 1. - vsmp[0]*vsmp[0] - vsmp[1]*vsmp[1];
		if (vsmp[2] <= 0) return 0;
		vsmp[2] = side * sqrt(vsmp[2]);
		if (SDreportError( SDsizeBSDF(&sa, vsmp, v0, SDqueryMin, bsdf), stderr))
			return 0;
		if (sa >= M_PI*siz*siz - FTINY)	/* no further division needed */
			return (*cf)(bsdf, vsmp, v0, cdata);
	}
	siz *= .5;				/* 4-branch recursion */
	return(	qtSampBSDF(xleft, ytop, siz, bsdf, side, v0, cf, cdata) &&
		qtSampBSDF(xleft+siz, ytop, siz, bsdf, side, v0, cf, cdata) &&
		qtSampBSDF(xleft, ytop+siz, siz, bsdf, side, v0, cf, cdata) &&
		qtSampBSDF(xleft+siz, ytop+siz, siz, bsdf, side, v0, cf, cdata) );
}

#define sampBSDFhemi(b,s,v0,cf,cd)	qtSampBSDF(0,0,1,b,s,v0,cf,cd)

/* Call-back to compute reciprocity difference */
int
diffRecip(const SDData *bsdf, const FVECT v1, const RREAL *v0, void *p)
{
	SDValue		sdv;
	double		otherY;

	if (SDreportError( SDevalBSDF(&sdv, v0, v1, bsdf), stderr))
		return 0;
	otherY = sdv.cieY;
	if (SDreportError( SDevalBSDF(&sdv, v1, v0, bsdf), stderr))
		return 0;

	addStat((SimpleStats *)p, rdiff(sdv.cieY, otherY));
	return 1;
}

/* Call-back to compute reciprocity over reflected hemisphere */
int
reflHemi(const SDData *bsdf, const FVECT v1, const RREAL *v0, void *p)
{
	return sampBSDFhemi(bsdf, 1 - 2*(v1[2]<0), v1, &diffRecip, p);
}

/* Call-back to compute reciprocity over transmitted hemisphere */
int
transHemi(const SDData *bsdf, const FVECT v1, const RREAL *v0, void *p)
{
	return sampBSDFhemi(bsdf, 1 - 2*(v1[2]>0), v1, &diffRecip, p);
}

/* Report reciprocity errors for the given directions */
void
checkReciprocity(const char *nm, const int side1, const int side2,
			const SDData *bsdf, const int fl)
{
	SimpleStats		myStats = SSinit;
	const SDSpectralDF	*df = bsdf->tf;

	if (side1 == side2) {
		df = (side1 > 0) ? bsdf->rf : bsdf->rb;
		if (!df) goto nothing2do;
	} else if (!bsdf->tf | !bsdf->tb)
		goto nothing2do;

	if (fl & F_MATRIX) {			/* special case for matrix BSDF */
		const SDMat	*m = (const SDMat *)df->comp[0].dist;
		int		i = m->ninc;
		FVECT		vin, vout;
		double		fwdY;
		SDValue		rev;
		while (i--) {
		    int	o = m->nout;
		    if (!mBSDF_incvec(vin, m, i+.5))
			continue;
		    while (o--) {
			if (!mBSDF_outvec(vout, m, o+.5))
				continue;
			fwdY = mBSDF_value(m, o, i);
			if (fwdY <= 1e-4)
				continue;
			if (SDreportError( SDevalBSDF(&rev, vout, vin, bsdf), stderr))
				return;
			if (rev.cieY > 1e-4)
				addStat(&myStats, rdiff(fwdY, rev.cieY));
		    }
		}
	} if (fl & F_ISOTROPIC) {		/* isotropic case */
		const double	stepSize = sqrt(df->minProjSA/M_PI);
		FVECT		vin;
		vin[1] = 0;
		for (vin[0] = 0.5*stepSize; vin[0] < 1; vin[0] += stepSize) {
			vin[2] = side1*sqrt(1. - vin[0]*vin[0]);
			if (!sampBSDFhemi(bsdf, side2, vin, &diffRecip, &myStats))
				return;
		}
	} else if (!sampBSDFhemi(bsdf, side1, NULL,
				(side1==side2) ? &reflHemi : &transHemi, &myStats))
		return;
	if (myStats.nvals) {
		printf("%s\t%5.1f\t%5.1f\t%5.1f\n", nm,
				100.*myStats.vmin,
				100.*myStats.vsum/(double)myStats.nvals,
				100.*myStats.vmax);
		return;
	}
nothing2do:
	printf("%s\t0\t0\t0\n", nm);
}

/* Report on the given BSDF XML file */
int
checkXML(char *fname)
{
	int		flags;
	SDData		myBSDF;
	char		*pth;

	puts("=====================================================");
	printf("File: '%s'\n", fname);
	SDclearBSDF(&myBSDF, fname);
	pth = getpath(fname, getrlibpath(), R_OK);
	if (!pth) {
		fprintf(stderr, "Cannot find file '%s'\n", fname);
		return 0;
	}
	if (SDreportError( SDloadFile(&myBSDF, pth), stderr))
		return 0;
	printf("Manufacturer: '%s'\n", myBSDF.makr);
	printf("BSDF Name: '%s'\n", myBSDF.matn);
	printf("Dimensions (W x H x Thickness): %g x %g x %g cm\n", 100.*myBSDF.dim[0],
				100.*myBSDF.dim[1], 100.*myBSDF.dim[2]);
	printf("Type: %s\n", getBSDFtype(&myBSDF, &flags));
	printf("Color: %d\n", (flags & F_IN_COLOR) != 0);
	printf("Has Geometry: %d\n", (myBSDF.mgf != NULL));
	puts("Component\tLambertian XYZ %\tMax. Dir\tMin. Angle");
	detailComponent("Internal Refl", &myBSDF.rLambFront, myBSDF.rf);
	detailComponent("External Refl", &myBSDF.rLambBack, myBSDF.rb);
	detailComponent("Int->Ext Trans", &myBSDF.tLambFront, myBSDF.tf);
	detailComponent("Ext->Int Trans", &myBSDF.tLambBack, myBSDF.tb);
	puts("Component\tReciprocity Error (min avg max %)");
	checkReciprocity("Front Refl", 1, 1, &myBSDF, flags);
	checkReciprocity("Back Refl", -1, -1, &myBSDF, flags);
	checkReciprocity("Transmission", -1, 1, &myBSDF, flags);
	SDfreeBSDF(&myBSDF);
	return 1;
}

int
main(int argc, char *argv[])
{
	int	i;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s bsdf.xml ..\n", argv[0]);
		return 1;
	}
	for (i = 1; i < argc; i++)
		if (!checkXML(argv[i]))
			return 1;
	return 0;
}
