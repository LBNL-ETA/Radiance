#ifndef lint
static const char RCSid[] = "$Id: m_bsdf.c,v 2.4 2011/02/19 23:42:09 greg Exp $";
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
 *	A thickness variable causes the strange but useful behavior
 *  of translating transmitted rays this distance past the surface
 *  intersection in the normal direction to bypass intervening geometry.
 *  This only affects scattered, non-source directed samples.  Thus,
 *  thickness is relevant only if there is a transmitted component.
 *  A positive thickness has the further side-effect that an unscattered
 *  (view) ray will pass right through our material if it has any
 *  non-diffuse transmission, making our BSDF invisible.  This allows the
 *  underlying geometry to become visible.  A matching surface should be
 *  placed on the other side, less than the thickness away, if the backside
 *  reflectance is non-zero.
 *	The "up" vector for the BSDF is given by three variables, defined
 *  (along with the thickness) by the named function file, or '.' if none.
 *  Together with the surface normal, this defines the local coordinate
 *  system for the BSDF.
 *	We do not reorient the surface, so if the BSDF has no back-side
 *  reflectance and none is given in the real arguments, the surface will
 *  appear as black when viewed from behind (unless backvis is false).
 *  The diffuse compnent arguments are added to components in the BSDF file,
 *  not multiplied.  However, patterns affect this material as a multiplier
 *  on everything except non-diffuse reflection.
 *
 *  Arguments for MAT_BSDF are:
 *	6+	thick	BSDFfile	ux uy uz	funcfile	transform
 *	0
 *	0|3|9	rdf	gdf	bdf
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
	double	sr_vpsa;	/* sqrt of BSDF projected solid angle */
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
bsdf_jitter(FVECT vres, BSDFDAT *ndp)
{
	double	sr_psa = ndp->sr_vpsa;

	VCOPY(vres, ndp->vray);
	if (specjitter < 1.)
		sr_psa *= specjitter;
	if (sr_psa <= FTINY)
		return;
	vres[0] += sr_psa*(.5 - frandom());
	vres[1] += sr_psa*(.5 - frandom());
	normalize(vres);
}

/* Compute source contribution for BSDF */
static void
dirbsdf(
	COLOR  cval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	SDError		ec;
	SDValue		sv;
	FVECT		vsrc;
	FVECT		vjit;
	double		ldot;
	double		dtmp;
	COLOR		ctmp;

	setcolor(cval, .0, .0, .0);

	ldot = DOT(np->pnorm, ldir);
	if ((-FTINY <= ldot) & (ldot <= FTINY))
		return;

	if (ldot > .0 && bright(np->rdiff) > FTINY) {
		/*
		 *  Compute added diffuse reflected component.
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ldot < .0 && bright(np->tdiff) > FTINY) {
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
	if (SDmapDir(vsrc, np->toloc, ldir) != SDEnone)
		return;
	bsdf_jitter(vjit, np);
	ec = SDevalBSDF(&sv, vjit, vsrc, np->sd);
	if (ec)
		objerror(np->mp, USER, transSDError(ec));

	if (sv.cieY <= FTINY)		/* not worth using? */
		return;
	cvt_sdcolor(ctmp, &sv);
	if (ldot > .0) {		/* pattern only diffuse reflection */
		COLOR	ctmp1, ctmp2;
		dtmp = (np->pr->rod > .0) ? np->sd->rLambFront.cieY
					: np->sd->rLambBack.cieY;
		dtmp /= PI * sv.cieY;	/* diffuse fraction */
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

/* Sample separate BSDF component */
static int
sample_sdcomp(BSDFDAT *ndp, SDComponent *dcp, int usepat)
{
	int	nstarget = 1;
	int	nsent = 0;
	SDError	ec;
	SDValue bsv;
	double	sthick;
	FVECT	vjit, vsmp;
	RAY	sr;
	int	ntrials;
						/* multiple samples? */
	if (specjitter > 1.5) {
		nstarget = specjitter*ndp->pr->rweight + .5;
		if (nstarget < 1)
			nstarget = 1;
	}
						/* run through our trials */
	for (ntrials = 0; nsent < nstarget && ntrials < 9*nstarget; ntrials++) {
		SDerrorDetail[0] = '\0';
						/* sample direction & coef. */
		bsdf_jitter(vjit, ndp);
		ec = SDsampComponent(&bsv, vsmp, vjit, ntrials ? frandom()
				: urand(ilhash(dimlist,ndims)+samplendx), dcp);
		if (ec)
			objerror(ndp->mp, USER, transSDError(ec));
						/* zero component? */
		if (bsv.cieY <= FTINY)
			break;
						/* map vector to world */
		if (SDmapDir(sr.rdir, ndp->fromloc, vsmp) != SDEnone)
			break;
						/* unintentional penetration? */
		if (DOT(sr.rdir, ndp->pr->ron) > .0 ^ vsmp[2] > .0)
			continue;
						/* spawn a specular ray */
		if (nstarget > 1)
			bsv.cieY /= (double)nstarget;
		cvt_sdcolor(sr.rcoef, &bsv);	/* use color */
		if (usepat)			/* pattern on transmission */
			multcolor(sr.rcoef, ndp->pr->pcol);
		if (rayorigin(&sr, SPECULAR, ndp->pr, sr.rcoef) < 0) {
			if (maxdepth  > 0)
				break;
			++nsent;		/* Russian roulette victim */
			continue;
		}
		if (ndp->thick > FTINY) {	/* need to move origin? */
			sthick = (ndp->pr->rod > .0) ? -ndp->thick : ndp->thick;
			if (sthick < .0 ^ vsmp[2] > .0)
				VSUM(sr.rorg, sr.rorg, ndp->pr->ron, sthick);
		}
		rayvalue(&sr);			/* send & evaluate sample */
		multcolor(sr.rcol, sr.rcoef);
		addcolor(ndp->pr->rcol, sr.rcol);
		++nsent;
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
		if (ndp->pr->rod > .0) {
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
			bsdf_jitter(vjit, ndp);
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
	COLOR	ctmp;
	SDError	ec;
	FVECT	upvec, outVec;
	MFUNC	*mf;
	BSDFDAT	nd;
						/* check arguments */
	if ((m->oargs.nsargs < 6) | (m->oargs.nfargs > 9) |
				(m->oargs.nfargs % 3))
		objerror(m, USER, "bad # arguments");

						/* get BSDF data */
	nd.sd = loadBSDF(m->oargs.sarg[1]);
						/* load cal file */
	mf = getfunc(m, 5, 0x1d, 1);
						/* get thickness */
	nd.thick = evalue(mf->ep[0]);
	if (nd.thick < .0)
		nd.thick = .0;
						/* check shadow */
	if (r->crtype & SHADOW) {
		if ((nd.thick > FTINY) & (nd.sd->tf != NULL))
			raytrans(r);		/* pass-through */
		SDfreeCache(nd.sd);
		return(1);			/* else shadow */
	}
						/* check unscattered ray */
	if (!(r->crtype & (SPECULAR|AMBIENT)) &&
			(nd.thick > FTINY) & (nd.sd->tf != NULL)) {
		SDfreeCache(nd.sd);
		raytrans(r);			/* pass-through */
		return(1);
	}
						/* diffuse reflectance */
	if (r->rod > .0) {
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
	if (bright(r->pcol) <= FTINY) {		/* black pattern?! */
		SDfreeCache(nd.sd);
		return(1);
	}
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
		ec = SDsizeBSDF(&nd.sr_vpsa, nd.vray, SDqueryMin, nd.sd);
	if (!ec)
		nd.sr_vpsa = sqrt(nd.sr_vpsa);
	else {
		objerror(m, WARNING, transSDError(ec));
		SDfreeCache(nd.sd);
		return(1);
	}
	
	if (r->rod < .0) {			/* perturb normal towards hit */
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
	if (bright(ctmp) > FTINY) {		/* ambient from this side */
		if (r->rod < .0)
			flipsurface(r);
		multambient(ctmp, r, nd.pnorm);
		addcolor(r->rcol, ctmp);
		if (r->rod < .0)
			flipsurface(r);
	}
	copycolor(ctmp, nd.tdiff);
	addcolor(ctmp, nd.tunsamp);
	if (bright(ctmp) > FTINY) {		/* ambient from other side */
		FVECT  bnorm;
		if (r->rod > .0)
			flipsurface(r);
		bnorm[0] = -nd.pnorm[0];
		bnorm[1] = -nd.pnorm[1];
		bnorm[2] = -nd.pnorm[2];
		multambient(ctmp, r, bnorm);
		addcolor(r->rcol, ctmp);
		if (r->rod > .0)
			flipsurface(r);
	}
						/* add direct component */
	direct(r, dirbsdf, &nd);
						/* clean up */
	SDfreeCache(nd.sd);
	return(1);
}
