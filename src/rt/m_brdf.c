#ifndef lint
static const char	RCSid[] = "$Id: m_brdf.c,v 2.29 2013/03/20 23:07:51 greg Exp $";
#endif
/*
 *  Shading for materials with arbitrary BRDF's
 */

#include "copyright.h"

#include  "ray.h"
#include  "ambient.h"
#include  "data.h"
#include  "source.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "func.h"

/*
 *	Arguments to this material include the color and specularity.
 *  String arguments include the reflection function and files.
 *  The BRDF is currently used just for the specular component to light
 *  sources.  Reflectance values or data coordinates are functions
 *  of the direction to the light source.  (Data modification functions
 *  are passed the source direction as args 2-4.)
 *	We orient the surface towards the incoming ray, so a single
 *  surface can be used to represent an infinitely thin object.
 *
 *  Arguments for MAT_PFUNC and MAT_MFUNC are:
 *	2+	func	funcfile	transform
 *	0
 *	4+	red	grn	blu	specularity	A5 ..
 *
 *  Arguments for MAT_PDATA and MAT_MDATA are:
 *	4+	func	datafile	funcfile	v0 ..	transform
 *	0
 *	4+	red	grn	blu	specularity	A5 ..
 *
 *  Arguments for MAT_TFUNC are:
 *	2+	func	funcfile	transform
 *	0
 *	6+	red	grn	blu	rspec	trans	tspec	A7 ..
 *
 *  Arguments for MAT_TDATA are:
 *	4+	func	datafile	funcfile	v0 ..	transform
 *	0
 *	6+	red	grn	blu	rspec	trans	tspec	A7 ..
 *
 *  Arguments for the more general MAT_BRTDF are:
 *	10+	rrefl	grefl	brefl
 *		rtrns	gtrns	btrns
 *		rbrtd	gbrtd	bbrtd
 *		funcfile	transform
 *	0
 *	9+	rdf	gdf	bdf
 *		rdb	gdb	bdb
 *		rdt	gdt	bdt	A10 ..
 *
 *	In addition to the normal variables available to functions,
 *  we define the following:
 *		NxP, NyP, NzP -		perturbed surface normal
 *		RdotP -			perturbed ray dot product
 *		CrP, CgP, CbP -		perturbed material color (or pattern)
 */

typedef struct {
	OBJREC  *mp;		/* material pointer */
	RAY  *pr;		/* intersected ray */
	DATARRAY  *dp;		/* data array for PDATA, MDATA or TDATA */
	COLOR  mcolor;		/* material (or pattern) color */
	COLOR  rdiff;		/* diffuse reflection */
	COLOR  tdiff;		/* diffuse transmission */
	double  rspec;		/* specular reflectance (1 for BRDTF) */
	double  trans;		/* transmissivity (.5 for BRDTF) */
	double  tspec;		/* specular transmittance (1 for BRDTF) */
	FVECT  pnorm;		/* perturbed surface normal */
	double  pdot;		/* perturbed dot product */
}  BRDFDAT;		/* BRDF material data */


static int setbrdfunc(BRDFDAT *np);


static void
dirbrdf(		/* compute source contribution */
	COLOR  cval,			/* returned coefficient */
	void  *nnp,		/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BRDFDAT *np = nnp;
	double  ldot;
	double  dtmp;
	COLOR  ctmp;
	FVECT  ldx;
	static double  vldx[5], pt[MAXDIM];
	char	**sa;
	int	i;
#define lddx (vldx+1)

	setcolor(cval, 0.0, 0.0, 0.0);
	
	ldot = DOT(np->pnorm, ldir);

	if (ldot <= FTINY && ldot >= -FTINY)
		return;		/* too close to grazing */

	if (ldot < 0.0 ? np->trans <= FTINY : np->trans >= 1.0-FTINY)
		return;		/* wrong side */

	if (ldot > 0.0) {
		/*
		 *  Compute and add diffuse reflected component to returned
		 *  color.  The diffuse reflected component will always be
		 *  modified by the color of the material.
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega / PI;
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	} else {
		/*
		 *  Diffuse transmitted component.
		 */
		copycolor(ctmp, np->tdiff);
		dtmp = -ldot * omega / PI;
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ldot > 0.0 ? np->rspec <= FTINY : np->tspec <= FTINY)
		return;		/* diffuse only */
					/* set up function */
	setbrdfunc(np);
	sa = np->mp->oargs.sarg;
	errno = 0;
					/* transform light vector */
	multv3(ldx, ldir, funcxf.xfm);
	for (i = 0; i < 3; i++)
		lddx[i] = ldx[i]/funcxf.sca;
	lddx[3] = omega;
					/* compute BRTDF */
	if (np->mp->otype == MAT_BRTDF) {
		if (sa[6][0] == '0')		/* special case */
			colval(ctmp,RED) = 0.0;
		else
			colval(ctmp,RED) = funvalue(sa[6], 4, lddx);
		if (sa[7][0] == '0')
			colval(ctmp,GRN) = 0.0;
		else if (!strcmp(sa[7],sa[6]))
			colval(ctmp,GRN) = colval(ctmp,RED);
		else
			colval(ctmp,GRN) = funvalue(sa[7], 4, lddx);
		if (!strcmp(sa[8],sa[6]))
			colval(ctmp,BLU) = colval(ctmp,RED);
		else if (!strcmp(sa[8],sa[7]))
			colval(ctmp,BLU) = colval(ctmp,GRN);
		else
			colval(ctmp,BLU) = funvalue(sa[8], 4, lddx);
		dtmp = bright(ctmp);
	} else if (np->dp == NULL) {
		dtmp = funvalue(sa[0], 4, lddx);
		setcolor(ctmp, dtmp, dtmp, dtmp);
	} else {
		for (i = 0; i < np->dp->nd; i++)
			pt[i] = funvalue(sa[3+i], 4, lddx);
		vldx[0] = datavalue(np->dp, pt);
		dtmp = funvalue(sa[0], 5, vldx);
		setcolor(ctmp, dtmp, dtmp, dtmp);
	}
	if ((errno == EDOM) | (errno == ERANGE)) {
		objerror(np->mp, WARNING, "compute error");
		return;
	}
	if (dtmp <= FTINY)
		return;
	if (ldot > 0.0) {
		/*
		 *  Compute reflected non-diffuse component.
		 */
		if ((np->mp->otype == MAT_MFUNC) | (np->mp->otype == MAT_MDATA))
			multcolor(ctmp, np->mcolor);
		dtmp = ldot * omega * np->rspec;
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	} else {
		/*
		 *  Compute transmitted non-diffuse component.
		 */
		if ((np->mp->otype == MAT_TFUNC) | (np->mp->otype == MAT_TDATA))
			multcolor(ctmp, np->mcolor);
		dtmp = -ldot * omega * np->tspec;
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
#undef lddx
}


int
m_brdf(			/* color a ray that hit a BRDTfunc material */
	OBJREC  *m,
	RAY  *r
)
{
	int  hitfront = 1;
	BRDFDAT  nd;
	RAY  sr;
	double  mirtest=0, mirdist=0;
	double  transtest, transdist;
	int  hasrefl, hastrans;
	int  hastexture;
	COLOR  ctmp;
	FVECT  vtmp;
	double  d;
	MFUNC  *mf;
	int  i;
						/* check arguments */
	if ((m->oargs.nsargs < 10) | (m->oargs.nfargs < 9))
		objerror(m, USER, "bad # arguments");
	nd.mp = m;
	nd.pr = r;
						/* dummy values */
	nd.rspec = nd.tspec = 1.0;
	nd.trans = 0.5;
						/* diffuse reflectance */
	if (r->rod > 0.0)
		setcolor(nd.rdiff, m->oargs.farg[0],
				m->oargs.farg[1],
				m->oargs.farg[2]);
	else
		setcolor(nd.rdiff, m->oargs.farg[3],
				m->oargs.farg[4],
				m->oargs.farg[5]);
						/* diffuse transmittance */
	setcolor(nd.tdiff, m->oargs.farg[6],
			m->oargs.farg[7],
			m->oargs.farg[8]);
						/* get modifiers */
	raytexture(r, m->omod);
	hastexture = DOT(r->pert,r->pert) > FTINY*FTINY;
	if (hastexture) {			/* perturb normal */
		nd.pdot = raynormal(nd.pnorm, r);
	} else {
		VCOPY(nd.pnorm, r->ron);
		nd.pdot = r->rod;
	}
	if (r->rod < 0.0) {			/* orient perturbed values */
		nd.pdot = -nd.pdot;
		for (i = 0; i < 3; i++) {
			nd.pnorm[i] = -nd.pnorm[i];
			r->pert[i] = -r->pert[i];
		}
		hitfront = 0;
	}
	copycolor(nd.mcolor, r->pcol);		/* get pattern color */
	multcolor(nd.rdiff, nd.mcolor);		/* modify diffuse values */
	multcolor(nd.tdiff, nd.mcolor);
	hasrefl = bright(nd.rdiff) > FTINY;
	hastrans = bright(nd.tdiff) > FTINY;
						/* load cal file */
	nd.dp = NULL;
	mf = getfunc(m, 9, 0x3f, 0);
						/* compute transmitted ray */
	setbrdfunc(&nd);
	errno = 0;
	setcolor(ctmp, evalue(mf->ep[3]),
			evalue(mf->ep[4]),
			evalue(mf->ep[5]));
	if ((errno == EDOM) | (errno == ERANGE))
		objerror(m, WARNING, "compute error");
	else if (rayorigin(&sr, TRANS, r, ctmp) == 0) {
		if (!(r->crtype & SHADOW) && hastexture) {
						/* perturb direction */
			VSUM(sr.rdir, r->rdir, r->pert, -.75);
			if (normalize(sr.rdir) == 0.0) {
				objerror(m, WARNING, "illegal perturbation");
				VCOPY(sr.rdir, r->rdir);
			}
		} else {
			VCOPY(sr.rdir, r->rdir);
		}
		rayvalue(&sr);
		multcolor(sr.rcol, sr.rcoef);
		addcolor(r->rcol, sr.rcol);
		if (!hastexture) {
			transtest = 2.0*bright(sr.rcol);
			transdist = r->rot + sr.rt;
		}
	}
	if (r->crtype & SHADOW)			/* the rest is shadow */
		return(1);
						/* compute reflected ray */
	setbrdfunc(&nd);
	errno = 0;
	setcolor(ctmp, evalue(mf->ep[0]),
			evalue(mf->ep[1]),
			evalue(mf->ep[2]));
	if ((errno == EDOM) | (errno == ERANGE))
		objerror(m, WARNING, "compute error");
	else if (rayorigin(&sr, REFLECTED, r, ctmp) == 0) {
		VSUM(sr.rdir, r->rdir, nd.pnorm, 2.*nd.pdot);
		checknorm(sr.rdir);
		rayvalue(&sr);
		multcolor(sr.rcol, sr.rcoef);
		addcolor(r->rcol, sr.rcol);
		if (!hastexture && r->ro != NULL && isflat(r->ro->otype)) {
			mirtest = 2.0*bright(sr.rcol);
			mirdist = r->rot + sr.rt;
		}
	}
						/* compute ambient */
	if (hasrefl) {
		if (!hitfront)
			flipsurface(r);
		copycolor(ctmp, nd.rdiff);
		multambient(ctmp, r, nd.pnorm);
		addcolor(r->rcol, ctmp);	/* add to returned color */
		if (!hitfront)
			flipsurface(r);
	}
	if (hastrans) {				/* from other side */
		if (hitfront)
			flipsurface(r);
		vtmp[0] = -nd.pnorm[0];
		vtmp[1] = -nd.pnorm[1];
		vtmp[2] = -nd.pnorm[2];
		copycolor(ctmp, nd.tdiff);
		multambient(ctmp, r, vtmp);
		addcolor(r->rcol, ctmp);
		if (hitfront)
			flipsurface(r);
	}
	if (hasrefl | hastrans || m->oargs.sarg[6][0] != '0')
		direct(r, dirbrdf, &nd);	/* add direct component */

	d = bright(r->rcol);			/* set effective distance */
	if (transtest > d)
		r->rt = transdist;
	else if (mirtest > d)
		r->rt = mirdist;

	return(1);
}



int
m_brdf2(			/* color a ray that hit a BRDF material */
	OBJREC  *m,
	RAY  *r
)
{
	BRDFDAT  nd;
	COLOR  ctmp;
	FVECT  vtmp;
	double  dtmp;
						/* always a shadow */
	if (r->crtype & SHADOW)
		return(1);
						/* check arguments */
	if ((m->oargs.nsargs < (hasdata(m->otype)?4:2)) | (m->oargs.nfargs <
			((m->otype==MAT_TFUNC)|(m->otype==MAT_TDATA)?6:4)))
		objerror(m, USER, "bad # arguments");
						/* check for back side */
	if (r->rod < 0.0) {
		if (!backvis && m->otype != MAT_TFUNC
				&& m->otype != MAT_TDATA) {
			raytrans(r);
			return(1);
		}
		raytexture(r, m->omod);
		flipsurface(r);			/* reorient if backvis */
	} else
		raytexture(r, m->omod);

	nd.mp = m;
	nd.pr = r;
						/* get material color */
	setcolor(nd.mcolor, m->oargs.farg[0],
			m->oargs.farg[1],
			m->oargs.farg[2]);
						/* get specular component */
	nd.rspec = m->oargs.farg[3];
						/* compute transmittance */
	if ((m->otype == MAT_TFUNC) | (m->otype == MAT_TDATA)) {
		nd.trans = m->oargs.farg[4]*(1.0 - nd.rspec);
		nd.tspec = nd.trans * m->oargs.farg[5];
		dtmp = nd.trans - nd.tspec;
		setcolor(nd.tdiff, dtmp, dtmp, dtmp);
	} else {
		nd.tspec = nd.trans = 0.0;
		setcolor(nd.tdiff, 0.0, 0.0, 0.0);
	}
						/* compute reflectance */
	dtmp = 1.0 - nd.trans - nd.rspec;
	setcolor(nd.rdiff, dtmp, dtmp, dtmp);
	nd.pdot = raynormal(nd.pnorm, r);	/* perturb normal */
	multcolor(nd.mcolor, r->pcol);		/* modify material color */
	multcolor(nd.rdiff, nd.mcolor);
	multcolor(nd.tdiff, nd.mcolor);
						/* load auxiliary files */
	if (hasdata(m->otype)) {
		nd.dp = getdata(m->oargs.sarg[1]);
		getfunc(m, 2, 0, 0);
	} else {
		nd.dp = NULL;
		getfunc(m, 1, 0, 0);
	}
						/* compute ambient */
	if (nd.trans < 1.0-FTINY) {
		copycolor(ctmp, nd.mcolor);	/* modified by material color */
		scalecolor(ctmp, 1.0-nd.trans);
		multambient(ctmp, r, nd.pnorm);
		addcolor(r->rcol, ctmp);	/* add to returned color */
	}
	if (nd.trans > FTINY) {		/* from other side */
		flipsurface(r);
		vtmp[0] = -nd.pnorm[0];
		vtmp[1] = -nd.pnorm[1];
		vtmp[2] = -nd.pnorm[2];
		copycolor(ctmp, nd.mcolor);
		scalecolor(ctmp, nd.trans);
		multambient(ctmp, r, vtmp);
		addcolor(r->rcol, ctmp);
		flipsurface(r);
	}
						/* add direct component */
	direct(r, dirbrdf, &nd);

	return(1);
}


static int
setbrdfunc(			/* set up brdf function and variables */
	BRDFDAT  *np
)
{
	FVECT  vec;

	if (setfunc(np->mp, np->pr) == 0)
		return(0);	/* it's OK, setfunc says we're done */
				/* else (re)assign special variables */
	multv3(vec, np->pnorm, funcxf.xfm);
	varset("NxP", '=', vec[0]/funcxf.sca);
	varset("NyP", '=', vec[1]/funcxf.sca);
	varset("NzP", '=', vec[2]/funcxf.sca);
	varset("RdotP", '=', np->pdot <= -1.0 ? -1.0 :
			np->pdot >= 1.0 ? 1.0 : np->pdot);
	varset("CrP", '=', colval(np->mcolor,RED));
	varset("CgP", '=', colval(np->mcolor,GRN));
	varset("CbP", '=', colval(np->mcolor,BLU));
	return(1);
}
