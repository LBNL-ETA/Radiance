#ifndef lint
static const char RCSid[] = "$Id: m_bsdf.c,v 2.1 2011/02/18 00:40:25 greg Exp $";
#endif
/*
 *  Shading for materials with BSDFs taken from XML data files
 */

#include "copyright.h"

#include  "ray.h"
#include  "paths.h"
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

typedef struct {
	OBJREC	*mp;		/* material pointer */
	RAY	*pr;		/* intersected ray */
	FVECT	pnorm;		/* perturbed surface normal */
	FVECT	vinc;		/* local incident vector */
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

/* Convert error from BSDF library */
static char *
cvt_sderr(SDError ec)
{
	if (!SDerrorDetail[0])
		return(strcpy(errmsg, SDerrorEnglish[ec]));
	sprintf(errmsg, "%s: %s", SDerrorEnglish[ec], SDerrorDetail);
	return(errmsg);
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
	BSDFDAT		*np = nnp;
	SDError		ec;
	SDValue		sv;
	FVECT		vout;
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
	if (SDmapDir(vout, np->toloc, ldir) != SDEnone)
		return;
	ec = SDevalBSDF(&sv, vout, np->vinc, np->sd);
	if (ec)
		objerror(np->mp, USER, cvt_sderr(ec));

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
		multcolor(ctmp, ctmp1);	/* apply desaturated pattern */
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
	FVECT	vout;
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
		ec = SDsampComponent(&bsv, vout, ndp->vinc,
				ntrials ? frandom()
					: urand(ilhash(dimlist,ndims)+samplendx),
						dcp);
		if (ec)
			objerror(ndp->mp, USER, cvt_sderr(ec));
						/* zero component? */
		if (bsv.cieY <= FTINY)
			break;
						/* map vector to world */
		if (SDmapDir(sr.rdir, ndp->fromloc, vout) != SDEnone)
			break;
						/* unintentional penetration? */
		if (DOT(sr.rdir, ndp->pr->ron) > .0 ^ vout[2] > .0)
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
						/* need to move origin? */
		sthick = (ndp->pr->rod > .0) ? -ndp->thick : ndp->thick;
		if (sthick < .0 ^ vout[2] > .0)
			VSUM(sr.rorg, sr.rorg, ndp->pr->ron, sthick);

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
		if (dfp->maxHemi > FTINY) {	/* XXX no color from BSDF! */
			double	d = SDdirectHemi(ndp->vinc, sflags, ndp->sd);
			COLOR	ctmp;
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

	SDerrorDetail[0] = '\0';		/* get BSDF data */
	nd.sd = SDgetCache(m->oargs.sarg[1]);
	if (nd.sd == NULL)
		error(SYSTEM, "out of memory in m_bsdf");
	if (!SDisLoaded(nd.sd)) {
		char	*pname = getpath(m->oargs.sarg[1], getrlibpath(), R_OK);
		if (pname == NULL) {
			sprintf(errmsg, "cannot find BSDF file \"%s\"",
						m->oargs.sarg[1]);
			objerror(m, USER, errmsg);
		}
		ec = SDloadFile(nd.sd, pname);
		if (ec)
			objerror(m, USER, cvt_sderr(ec));
		SDretainSet = SDretainAll;
	}
						/* load cal file */
	mf = getfunc(m, 5, 0x1d, 1);
						/* get thickness */
	nd.thick = evalue(mf->ep[0]);
	if (nd.thick < .0)
		nd.thick = .0;
						/* check shadow */
	if (r->crtype & SHADOW) {
		SDfreeCache(nd.sd);
		if (nd.thick > FTINY && nd.sd->tf != NULL &&
				nd.sd->tf->maxHemi > FTINY)
			raytrans(r);		/* pass-through */
		return(1);			/* else shadow */
	}
						/* check unscattered ray */
	if (!(r->crtype & (SPECULAR|AMBIENT)) && nd.thick > FTINY &&
			nd.sd->tf != NULL && nd.sd->tf->maxHemi > FTINY) {
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
			if (!backvis && (nd.sd->rb == NULL ||
						nd.sd->rb->maxHemi <= FTINY) &&
					(nd.sd->tf == NULL ||
						nd.sd->tf->maxHemi <= FTINY)) {
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
		nd.vinc[0] = -r->rdir[0];
		nd.vinc[1] = -r->rdir[1];
		nd.vinc[2] = -r->rdir[2];
		ec = SDmapDir(nd.vinc, nd.toloc, nd.vinc);
	}
	if (!ec)
		ec = SDinvXform(nd.fromloc, nd.toloc);
	if (ec) {
		objerror(m, WARNING, cvt_sderr(ec));
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
