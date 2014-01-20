#ifndef lint
static const char RCSid[] = "$Id: cmbsdf.c,v 2.1 2014/01/20 21:29:04 greg Exp $";
#endif
/*
 * Load and convert BSDF into color coefficient matrix representation.
 *
 *	G. Ward
 */

#include "standard.h"
#include "cmatrix.h"
#include "paths.h"
#include "bsdf.h"
#include "bsdf_m.h"

/* Convert a BSDF to our matrix representation */
static CMATRIX *
cm_bsdf(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->nout, bsdf->ninc);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const double	dom = mBSDF_incohm(bsdf,c);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			float	f = mBSDF_value(bsdf,c,r);
			COLORV	*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Convert between input and output indices for reciprocity */
static int
recip_out_from_in(const SDMat *bsdf, int in_recip)
{
	FVECT	v;

	if (!mBSDF_incvec(v, bsdf, in_recip+.5))
		return(in_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_outndx(bsdf, v));
}

/* Convert between output and input indices for reciprocity */
static int
recip_in_from_out(const SDMat *bsdf, int out_recip)
{
	FVECT	v;

	if (!mBSDF_outvec(v, bsdf, out_recip+.5))
		return(out_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_incndx(bsdf, v));
}

/* Convert a BSDF to our matrix representation, applying reciprocity */
static CMATRIX *
cm_bsdf_recip(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->ninc, bsdf->nout);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const int	ro = recip_out_from_in(bsdf,c);
		const double	dom = mBSDF_outohm(bsdf,ro);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			const int	ri = recip_in_from_out(bsdf,r);
			float		f = mBSDF_value(bsdf,ri,ro);
			COLORV		*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Load and convert a matrix BSDF from the given XML file */
CMATRIX *
cm_loadBSDF(char *fname, COLOR cLamb)
{
	CMATRIX		*Tmat;
	char		*fpath;
	int		recip;
	SDError		ec;
	SDData		myBSDF;
	SDSpectralDF	*tdf;
	COLOR		bsdfLamb, specCol;
					/* find path to BSDF file */
	fpath = getpath(fname, getrlibpath(), R_OK);
	if (fpath == NULL) {
		sprintf(errmsg, "cannot find BSDF file '%s'", fname);
		error(USER, errmsg);
	}
	SDclearBSDF(&myBSDF, fname);	/* load XML and check type */
	ec = SDloadFile(&myBSDF, fpath);
	if (ec)
		error(USER, transSDError(ec));
	ccy2rgb(&myBSDF.tLamb.spec, myBSDF.tLamb.cieY/PI, bsdfLamb);
	recip = (myBSDF.tb == NULL);
	tdf = recip ? myBSDF.tf : myBSDF.tb;
	if (tdf == NULL) {		/* no non-Lambertian transmission? */
		if (cLamb != NULL)
			copycolor(cLamb, bsdfLamb);
		SDfreeBSDF(&myBSDF);
		return(NULL);
	}
	if (tdf->ncomp != 1 || tdf->comp[0].func != &SDhandleMtx) {
		sprintf(errmsg, "unsupported BSDF '%s'", fpath);
		error(USER, errmsg);
	}
					/* convert BTDF to matrix */
	ccy2rgb(&tdf->comp[0].cspec[0], 1., specCol);
	Tmat = recip ? cm_bsdf_recip(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist)
			: cm_bsdf(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist);
	if (cLamb != NULL)		/* Lambertian is included */
		setcolor(cLamb, .0, .0, .0);
					/* free BSDF and return */
	SDfreeBSDF(&myBSDF);
	return(Tmat);
}
