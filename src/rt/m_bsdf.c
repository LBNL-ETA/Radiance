#ifndef lint
static const char RCSid[] = "$Id: m_bsdf.c,v 2.16 2011/08/24 04:31:13 greg Exp $";
#endif
/*
 *  Shading for materials with BSDFs taken from XML data files
 */

#include "copyright.h"

#include  "ray.h"
#include  "ambient.h"
#include  "source.h"
#include  "func.h"
#include  "bsdf.h"
#include  "random.h"

/*
 *	Arguments to this material include optional diffuse colors.
 *  String arguments include the BSDF and function files.
 *	A non-zero thickness causes the strange but useful behavior
 *  of translating transmitted rays this distance beneath the surface
 *  (opposite the surface normal) to bypass any intervening geometry.
 *  Translation only affects scattered, non-source-directed samples.
 *  A non-zero thickness has the further side-effect that an unscattered
 *  (view) ray will pass right through our material if it has any
 *  non-diffuse transmission, making the BSDF surface invisible.  This
 *  shows the proxied geometry instead. Thickness has the further
 *  effect of turning off reflection on the hidden side so that rays
 *  heading in the opposite direction pass unimpeded through the BSDF
 *  surface.  A paired surface may be placed on the opposide side of
 *  the detail geometry, less than this thickness away, if a two-way
 *  proxy is desired.  Note that the sign of the thickness is important.
 *  A positive thickness hides geometry behind the BSDF surface and uses
 *  front reflectance and transmission properties.  A negative thickness
 *  hides geometry in front of the surface when rays hit from behind,
 *  and applies only the transmission and backside reflectance properties.
 *  Reflection is ignored on the hidden side, as those rays pass through.
 *	The "up" vector for the BSDF is given by three variables, defined
 *  (along with the thickness) by the named function file, or '.' if none.
 *  Together with the surface normal, this defines the local coordinate
 *  system for the BSDF.
 *	We do not reorient the surface, so if the BSDF has no back-side
 *  reflectance and none is given in the real arguments, a BSDF surface
 *  with zero thickness will appear black when viewed from behind
 *  unless backface visibility is off.
 *	The diffuse arguments are added to components in the BSDF file,
 *  not multiplied.  However, patterns affect this material as a multiplier
 *  on everything except non-diffuse reflection.
 *
 *  Arguments for MAT_BSDF are:
 *	6+	thick	BSDFfile	ux uy uz	funcfile	transform
 *	0
 *	0|3|6|9	rdf	gdf	bdf
 *		rdb	gdb	bdb
 *		rdt	gdt	bdt
 */

/*
 * Note that our reverse ray-tracing process means that the positions
 * of incoming and outgoing vectors may be reversed in our calls
 * to the BSDF library.  This is fine, since the bidirectional nature
 * of the BSDF (that's what the 'B' stands for) means it all works out.
 */

typedef struct {
	OBJREC	*mp;		/* material pointer */
	RAY	*pr;		/* intersected ray */
	FVECT	pnorm;		/* perturbed surface normal */
	FVECT	vray;		/* local outgoing (return) vector */
	double	sr_vpsa[2];	/* sqrt of BSDF projected solid angle extrema */
	RREAL	toloc[3][3];	/* world to local BSDF coords */
	RREAL	fromloc[3][3];	/* local BSDF coords to world */
	double	thick;		/* surface thickness */
	SDData	*sd;		/* loaded BSDF data */
	COLOR	runsamp;	/* BSDF hemispherical reflection */
	COLOR	rdiff;		/* added diffuse reflection */
	COLOR	tunsamp;	/* BSDF hemispherical transmission */
	COLOR	tdiff;		/* added diffuse transmission */
}  BSDFDAT;		/* BSDF material data */

#define	cvt_sdcolor(cv, svp)	ccy2rgb(&(svp)->spec, (svp)->cieY, cv)

/* Jitter ray sample according to projected solid angle and specjitter */
static void
bsdf_jitter(FVECT vres, BSDFDAT *ndp, double sr_psa)
{
	VCOPY(vres, ndp->vray);
	if (specjitter < 1.)
		sr_psa *= specjitter;
	if (sr_psa <= FTINY)
		return;
	vres[0] += sr_psa*(.5 - frandom());
	vres[1] += sr_psa*(.5 - frandom());
	normalize(vres);
}

/* Evaluate BSDF for direct component, returning true if OK to proceed */
static int
direct_bsdf_OK(COLOR cval, FVECT ldir, double omega, BSDFDAT *ndp)
{
	int	nsamp, ok = 0;
	FVECT	vsrc, vsmp, vjit;
	double	tomega;
	double	sf, tsr, sd[2];
	COLOR	csmp;
	SDValue	sv;
	SDError	ec;
	int	i;
					/* transform source direction */
	if (SDmapDir(vsrc, ndp->toloc, ldir) != SDEnone)
		return(0);
					/* assign number of samples */
	ec = SDsizeBSDF(&tomega, ndp->vray, vsrc, SDqueryMin, ndp->sd);
	if (ec)
		goto baderror;
					/* check indirect over-counting */
	if (ndp->thick != 0 && ndp->pr->crtype & (SPECULAR|AMBIENT)
				&& vsrc[2] > 0 ^ ndp->vray[2] > 0) {
		double	dx = vsrc[0] + ndp->vray[0];
		double	dy = vsrc[1] + ndp->vray[1];
		if (dx*dx + dy*dy <= omega+tomega)
			return(0);
	}
	sf = specjitter * ndp->pr->rweight;
	if (25.*tomega <= omega)
		nsamp = 100.*sf + .5;
	else
		nsamp = 4.*sf*omega/tomega + .5;
	nsamp += !nsamp;
	setcolor(cval, .0, .0, .0);	/* sample our source area */
	sf = sqrt(omega);
	tsr = sqrt(tomega);
	for (i = nsamp; i--; ) {
		VCOPY(vsmp, vsrc);	/* jitter query directions */
		if (nsamp > 1) {
			multisamp(sd, 2, (i + frandom())/(double)nsamp);
			vsmp[0] += (sd[0] - .5)*sf;
			vsmp[1] += (sd[1] - .5)*sf;
			if (normalize(vsmp) == 0) {
				--nsamp;
				continue;
			}
		}
		bsdf_jitter(vjit, ndp, tsr);
					/* compute BSDF */
		ec = SDevalBSDF(&sv, vjit, vsmp, ndp->sd);
		if (ec)
			goto baderror;
		if (sv.cieY <= FTINY)	/* worth using? */
			continue;
		cvt_sdcolor(csmp, &sv);
		addcolor(cval, csmp);	/* average it in */
		++ok;
	}
	sf = 1./(double)nsamp;
	scalecolor(cval, sf);
	return(ok);
baderror:
	objerror(ndp->mp, USER, transSDError(ec));
}

/* Compute source contribution for BSDF (reflected & transmitted) */
static void
dir_bsdf(
	COLOR  cval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	COLOR		ctmp;

	setcolor(cval, .0, .0, .0);

	ldot = DOT(np->pnorm, ldir);
	if ((-FTINY <= ldot) & (ldot <= FTINY))
		return;

	if (ldot > 0 && bright(np->rdiff) > FTINY) {
		/*
		 *  Compute added diffuse reflected component.
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ldot < 0 && bright(np->tdiff) > FTINY) {
		/*
		 *  Compute added diffuse transmission.
		 */
		copycolor(ctmp, np->tdiff);
		dtmp = -ldot * omega * (1.0/PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	/*
	 *  Compute scattering coefficient using BSDF.
	 */
	if (!direct_bsdf_OK(ctmp, ldir, omega, np))
		return;
	if (ldot > 0) {		/* pattern only diffuse reflection */
		COLOR	ctmp1, ctmp2;
		dtmp = (np->pr->rod > 0) ? np->sd->rLambFront.cieY
					: np->sd->rLambBack.cieY;
					/* diffuse fraction */
		dtmp /= PI * bright(ctmp);
		copycolor(ctmp2, np->pr->pcol);
		scalecolor(ctmp2, dtmp);
		setcolor(ctmp1, 1.-dtmp, 1.-dtmp, 1.-dtmp);
		addcolor(ctmp1, ctmp2);
		multcolor(ctmp, ctmp1);	/* apply derated pattern */
		dtmp = ldot * omega;
	} else {			/* full pattern on transmission */
		multcolor(ctmp, np->pr->pcol);
		dtmp = -ldot * omega;
	}
	scalecolor(ctmp, dtmp);
	addcolor(cval, ctmp);
}

/* Compute source contribution for BSDF (reflected only) */
static void
dir_brdf(
	COLOR  cval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	COLOR		ctmp, ctmp1, ctmp2;

	setcolor(cval, .0, .0, .0);

	ldot = DOT(np->pnorm, ldir);
	
	if (ldot <= FTINY)
		return;

	if (bright(np->rdiff) > FTINY) {
		/*
		 *  Compute added diffuse reflected component.
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	/*
	 *  Compute reflection coefficient using BSDF.
	 */
	if (!direct_bsdf_OK(ctmp, ldir, omega, np))
		return;
					/* pattern only diffuse reflection */
	dtmp = (np->pr->rod > 0) ? np->sd->rLambFront.cieY
				: np->sd->rLambBack.cieY;
	dtmp /= PI * bright(ctmp);	/* diffuse fraction */
	copycolor(ctmp2, np->pr->pcol);
	scalecolor(ctmp2, dtmp);
	setcolor(ctmp1, 1.-dtmp, 1.-dtmp, 1.-dtmp);
	addcolor(ctmp1, ctmp2);
	multcolor(ctmp, ctmp1);		/* apply derated pattern */
	dtmp = ldot * omega;
	scalecolor(ctmp, dtmp);
	addcolor(cval, ctmp);
}

/* Compute source contribution for BSDF (transmitted only) */
static void
dir_btdf(
	COLOR  cval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	COLOR		ctmp;

	setcolor(cval, .0, .0, .0);

	ldot = DOT(np->pnorm, ldir);

	if (ldot >= -FTINY)
		return;

	if (bright(np->tdiff) > FTINY) {
		/*
		 *  Compute added diffuse transmission.
		 */
		copycolor(ctmp, np->tdiff);
		dtmp = -ldot * omega * (1.0/PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	/*
	 *  Compute scattering coefficient using BSDF.
	 */
	if (!direct_bsdf_OK(ctmp, ldir, omega, np))
		return;
					/* full pattern on transmission */
	multcolor(ctmp, np->pr->pcol);
	dtmp = -ldot * omega;
	scalecolor(ctmp, dtmp);
	addcolor(cval, ctmp);
}

/* Sample separate BSDF component */
static int
sample_sdcomp(BSDFDAT *ndp, SDComponent *dcp, int usepat)
{
	int	nstarget = 1;
	int	nsent;
	SDError	ec;
	SDValue bsv;
	double	xrand;
	FVECT	vsmp;
	RAY	sr;
						/* multiple samples? */
	if (specjitter > 1.5) {
		nstarget = specjitter*ndp->pr->rweight + .5;
		nstarget += !nstarget;
	}
						/* run through our samples */
	for (nsent = 0; nsent < nstarget; nsent++) {
		if (nstarget == 1) {		/* stratify random variable */
			xrand = urand(ilhash(dimlist,ndims)+samplendx);
			if (specjitter < 1.)
				xrand = .5 + specjitter*(xrand-.5);
		} else {
			xrand = (nsent + frandom())/(double)nstarget;
		}
		SDerrorDetail[0] = '\0';	/* sample direction & coef. */
		bsdf_jitter(vsmp, ndp, ndp->sr_vpsa[0]);
		ec = SDsampComponent(&bsv, vsmp, xrand, dcp);
		if (ec)
			objerror(ndp->mp, USER, transSDError(ec));
		if (bsv.cieY <= FTINY)		/* zero component? */
			break;
						/* map vector to world */
		if (SDmapDir(sr.rdir, ndp->fromloc, vsmp) != SDEnone)
			break;
						/* spawn a specular ray */
		if (nstarget > 1)
			bsv.cieY /= (double)nstarget;
		cvt_sdcolor(sr.rcoef, &bsv);	/* use sample color */
		if (usepat)			/* apply pattern? */
			multcolor(sr.rcoef, ndp->pr->pcol);
		if (rayorigin(&sr, SPECULAR, ndp->pr, sr.rcoef) < 0) {
			if (maxdepth > 0)
				break;
			continue;		/* Russian roulette victim */
		}
						/* need to offset origin? */
		if (ndp->thick != 0 && ndp->pr->rod > 0 ^ vsmp[2] > 0)
			VSUM(sr.rorg, sr.rorg, ndp->pr->ron, -ndp->thick);
		rayvalue(&sr);			/* send & evaluate sample */
		multcolor(sr.rcol, sr.rcoef);
		addcolor(ndp->pr->rcol, sr.rcol);
	}
	return(nsent);
}

/* Sample non-diffuse components of BSDF */
static int
sample_sdf(BSDFDAT *ndp, int sflags)
{
	int		n, ntotal = 0;
	SDSpectralDF	*dfp;
	COLORV		*unsc;

	if (sflags == SDsampSpT) {
		unsc = ndp->tunsamp;
		dfp = ndp->sd->tf;
		cvt_sdcolor(unsc, &ndp->sd->tLamb);
	} else /* sflags == SDsampSpR */ {
		unsc = ndp->runsamp;
		if (ndp->pr->rod > 0) {
			dfp = ndp->sd->rf;
			cvt_sdcolor(unsc, &ndp->sd->rLambFront);
		} else {
			dfp = ndp->sd->rb;
			cvt_sdcolor(unsc, &ndp->sd->rLambBack);
		}
	}
	multcolor(unsc, ndp->pr->pcol);
	if (dfp == NULL)			/* no specular component? */
		return(0);
						/* below sampling threshold? */
	if (dfp->maxHemi <= specthresh+FTINY) {
		if (dfp->maxHemi > FTINY) {	/* XXX no color from BSDF */
			FVECT	vjit;
			double	d;
			COLOR	ctmp;
			bsdf_jitter(vjit, ndp, ndp->sr_vpsa[1]);
			d = SDdirectHemi(vjit, sflags, ndp->sd);
			if (sflags == SDsampSpT) {
				copycolor(ctmp, ndp->pr->pcol);
				scalecolor(ctmp, d);
			} else			/* no pattern on reflection */
				setcolor(ctmp, d, d, d);
			addcolor(unsc, ctmp);
		}
		return(0);
	}
						/* else need to sample */
	dimlist[ndims++] = (int)(size_t)ndp->mp;
	ndims++;
	for (n = dfp->ncomp; n--; ) {		/* loop over components */
		dimlist[ndims-1] = n + 9438;
		ntotal += sample_sdcomp(ndp, &dfp->comp[n], sflags==SDsampSpT);
	}
	ndims -= 2;
	return(ntotal);
}

/* Color a ray that hit a BSDF material */
int
m_bsdf(OBJREC *m, RAY *r)
{
	int	hitfront;
	COLOR	ctmp;
	SDError	ec;
	FVECT	upvec, vtmp;
	MFUNC	*mf;
	BSDFDAT	nd;
						/* check arguments */
	if ((m->oargs.nsargs < 6) | (m->oargs.nfargs > 9) |
				(m->oargs.nfargs % 3))
		objerror(m, USER, "bad # arguments");
						/* record surface struck */
	hitfront = (r->rod > 0);
						/* load cal file */
	mf = getfunc(m, 5, 0x1d, 1);
						/* get thickness */
	nd.thick = evalue(mf->ep[0]);
	if ((-FTINY <= nd.thick) & (nd.thick <= FTINY))
		nd.thick = .0;
						/* check shadow */
	if (r->crtype & SHADOW) {
		if (nd.thick != 0)
			raytrans(r);		/* pass-through */
		return(1);			/* or shadow */
	}
						/* check other rays to pass */
	if (nd.thick != 0 && (!(r->crtype & (SPECULAR|AMBIENT)) ||
				nd.thick > 0 ^ hitfront)) {
		raytrans(r);			/* hide our proxy */
		return(1);
	}
						/* get BSDF data */
	nd.sd = loadBSDF(m->oargs.sarg[1]);
						/* diffuse reflectance */
	if (hitfront) {
		if (m->oargs.nfargs < 3)
			setcolor(nd.rdiff, .0, .0, .0);
		else
			setcolor(nd.rdiff, m->oargs.farg[0],
					m->oargs.farg[1],
					m->oargs.farg[2]);
	} else {
		if (m->oargs.nfargs < 6) {	/* check invisible backside */
			if (!backvis && (nd.sd->rb == NULL) &
						(nd.sd->tf == NULL)) {
				SDfreeCache(nd.sd);
				raytrans(r);
				return(1);
			}
			setcolor(nd.rdiff, .0, .0, .0);
		} else
			setcolor(nd.rdiff, m->oargs.farg[3],
					m->oargs.farg[4],
					m->oargs.farg[5]);
	}
						/* diffuse transmittance */
	if (m->oargs.nfargs < 9)
		setcolor(nd.tdiff, .0, .0, .0);
	else
		setcolor(nd.tdiff, m->oargs.farg[6],
				m->oargs.farg[7],
				m->oargs.farg[8]);
	nd.mp = m;
	nd.pr = r;
						/* get modifiers */
	raytexture(r, m->omod);
						/* modify diffuse values */
	multcolor(nd.rdiff, r->pcol);
	multcolor(nd.tdiff, r->pcol);
						/* get up vector */
	upvec[0] = evalue(mf->ep[1]);
	upvec[1] = evalue(mf->ep[2]);
	upvec[2] = evalue(mf->ep[3]);
						/* return to world coords */
	if (mf->f != &unitxf) {
		multv3(upvec, upvec, mf->f->xfm);
		nd.thick *= mf->f->sca;
	}
	raynormal(nd.pnorm, r);
						/* compute local BSDF xform */
	ec = SDcompXform(nd.toloc, nd.pnorm, upvec);
	if (!ec) {
		nd.vray[0] = -r->rdir[0];
		nd.vray[1] = -r->rdir[1];
		nd.vray[2] = -r->rdir[2];
		ec = SDmapDir(nd.vray, nd.toloc, nd.vray);
	}
	if (!ec)
		ec = SDinvXform(nd.fromloc, nd.toloc);
						/* determine BSDF resolution */
	if (!ec)
		ec = SDsizeBSDF(nd.sr_vpsa, nd.vray, NULL,
						SDqueryMin+SDqueryMax, nd.sd);
	if (ec) {
		objerror(m, WARNING, transSDError(ec));
		SDfreeCache(nd.sd);
		return(1);
	}
	nd.sr_vpsa[0] = sqrt(nd.sr_vpsa[0]);
	nd.sr_vpsa[1] = sqrt(nd.sr_vpsa[1]);
	if (!hitfront) {			/* perturb normal towards hit */
		nd.pnorm[0] = -nd.pnorm[0];
		nd.pnorm[1] = -nd.pnorm[1];
		nd.pnorm[2] = -nd.pnorm[2];
	}
						/* sample reflection */
	sample_sdf(&nd, SDsampSpR);
						/* sample transmission */
	sample_sdf(&nd, SDsampSpT);
						/* compute indirect diffuse */
	copycolor(ctmp, nd.rdiff);
	addcolor(ctmp, nd.runsamp);
	if (bright(ctmp) > FTINY) {		/* ambient from reflection */
		if (!hitfront)
			flipsurface(r);
		multambient(ctmp, r, nd.pnorm);
		addcolor(r->rcol, ctmp);
		if (!hitfront)
			flipsurface(r);
	}
	copycolor(ctmp, nd.tdiff);
	addcolor(ctmp, nd.tunsamp);
	if (bright(ctmp) > FTINY) {		/* ambient from other side */
		FVECT  bnorm;
		if (hitfront)
			flipsurface(r);
		bnorm[0] = -nd.pnorm[0];
		bnorm[1] = -nd.pnorm[1];
		bnorm[2] = -nd.pnorm[2];
		if (nd.thick != 0) {		/* proxy with offset? */
			VCOPY(vtmp, r->rop);
			VSUM(r->rop, vtmp, r->ron, -nd.thick);
			multambient(ctmp, r, bnorm);
			VCOPY(r->rop, vtmp);
		} else
			multambient(ctmp, r, bnorm);
		addcolor(r->rcol, ctmp);
		if (hitfront)
			flipsurface(r);
	}
						/* add direct component */
	if ((bright(nd.tdiff) <= FTINY) & (nd.sd->tf == NULL)) {
		direct(r, dir_brdf, &nd);	/* reflection only */
	} else if (nd.thick == 0) {
		direct(r, dir_bsdf, &nd);	/* thin surface scattering */
	} else {
		direct(r, dir_brdf, &nd);	/* reflection first */
		VCOPY(vtmp, r->rop);		/* offset for transmitted */
		VSUM(r->rop, vtmp, r->ron, -nd.thick);
		direct(r, dir_btdf, &nd);	/* separate transmission */
		VCOPY(r->rop, vtmp);
	}
						/* clean up */
	SDfreeCache(nd.sd);
	return(1);
}
