#ifndef lint
static const char RCSid[] = "$Id: cmbsdf.c,v 2.6 2019/07/31 22:00:08 greg Exp $";
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
cm_bsdf(const COLOR diffBSDF, const SDMat *bsdf)
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
				setcolor(mp, .0f, .0f, .0f);
			} else if (bsdf->chroma != NULL) {
				C_COLOR	cxy;
				c_decodeChroma(&cxy, mBSDF_chroma(bsdf,c,r));
				ccy2rgb(&cxy, f, mp);
			} else
				setcolor(mp, f, f, f);
			addcolor(mp, diffBSDF);
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
cm_bsdf_recip(const COLOR diffBSDF, const SDMat *bsdf)
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
				setcolor(mp, .0f, .0f, .0f);
			} else if (bsdf->chroma != NULL) {
				C_COLOR	cxy;
				c_decodeChroma(&cxy, mBSDF_chroma(bsdf,ri,ro));
				ccy2rgb(&cxy, f, mp);
			} else
				setcolor(mp, f, f, f);
			addcolor(mp, diffBSDF);
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

/* Return a Lambertian (diffuse) matrix */
static CMATRIX *
cm_bsdf_Lamb(const COLOR diffBSDF)
{
	CMATRIX	*cm = cm_alloc(145, 145);	/* this is a hack */
	int	r, c;

	for (c = 0; c < cm->ncols; c++)
		for (r = 0; r < cm->nrows; r++) {
			COLORV	*mp = cm_lval(cm,r,c);
			copycolor(mp, diffBSDF);
		}
	return(cm);
}

/* Load and convert a matrix BTDF from the given XML file */
CMATRIX *
cm_loadBTDF(char *fname)
{
	CMATRIX		*Tmat;
	int		recip;
	SDError		ec;
	SDData		myBSDF;
	SDSpectralDF	*tdf;
	COLOR		diffBSDF;

	SDclearBSDF(&myBSDF, fname);	/* load XML and check type */
	ec = SDloadFile(&myBSDF, fname);
	if (ec)
		error(USER, transSDError(ec));
	ccy2rgb(&myBSDF.tLamb.spec, myBSDF.tLamb.cieY/PI, diffBSDF);
	recip = (myBSDF.tb == NULL);
	tdf = recip ? myBSDF.tf : myBSDF.tb;
	if (tdf == NULL) {		/* no non-Lambertian transmission? */
		SDfreeBSDF(&myBSDF);
		return(cm_bsdf_Lamb(diffBSDF));
	}
	if (tdf->ncomp != 1 || tdf->comp[0].func != &SDhandleMtx) {
		sprintf(errmsg, "unsupported BSDF '%s'", fname);
		error(USER, errmsg);
	}
					/* convert BTDF to matrix */
	Tmat = recip ? cm_bsdf_recip(diffBSDF, (SDMat *)tdf->comp[0].dist)
			: cm_bsdf(diffBSDF, (SDMat *)tdf->comp[0].dist);
					/* free BSDF and return */
	SDfreeBSDF(&myBSDF);
	return(Tmat);
}
