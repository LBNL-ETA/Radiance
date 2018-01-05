#ifndef lint
static const char RCSid[] = "$Id: m_bsdf.c,v 2.45 2018/01/05 02:47:46 greg Exp $";
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
#include  "pmapmat.h"

/*
 *	Arguments to this material include optional diffuse colors.
 *  String arguments include the BSDF and function files.
 *	A non-zero thickness causes the strange but useful behavior
 *  of translating transmitted rays this distance beneath the surface
 *  (opposite the surface normal) to bypass any intervening geometry.
 *  Translation only affects scattered, non-source-directed samples.
 *  A non-zero thickness has the further side-effect that an unscattered
 *  (view) ray will pass right through our material, making the BSDF
 *  surface invisible and showing the proxied geometry instead. Thickness
 *  has the further effect of turning off reflection on the reverse side so
 *  rays heading in the opposite direction pass unimpeded through the BSDF
 *  surface.  A paired surface may be placed on the opposide side of
 *  the detail geometry, less than this thickness away, if a two-way
 *  proxy is desired.  Note that the sign of the thickness is important.
 *  A positive thickness hides geometry behind the BSDF surface and uses
 *  front reflectance and transmission properties.  A negative thickness
 *  hides geometry in front of the surface when rays hit from behind,
 *  and applies only the transmission and backside reflectance properties.
 *  Reflection is ignored on the hidden side, as those rays pass through.
 *	When thickness is set to zero, shadow rays will be blocked unless
 *  a BTDF has a strong "through" component in the source direction.
 *  A separate test prevents over-counting by dropping samples that are
 *  too close to this "through" direction.  BSDFs with such a through direction
 *  will also have a view component, meaning they are somewhat see-through.
 *	The "up" vector for the BSDF is given by three variables, defined
 *  (along with the thickness) by the named function file, or '.' if none.
 *  Together with the surface normal, this defines the local coordinate
 *  system for the BSDF.
 *	We do not reorient the surface, so if the BSDF has no back-side
 *  reflectance and none is given in the real arguments, a BSDF surface
 *  with zero thickness will appear black when viewed from behind
 *  unless backface visibility is on, when it becomes invisible.
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
 * to the BSDF library.  This is usually fine, since the bidirectional nature
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
	COLOR	cthru;		/* "through" component multiplier */
	SDData	*sd;		/* loaded BSDF data */
	COLOR	rdiff;		/* diffuse reflection */
	COLOR	runsamp;	/* BSDF hemispherical reflection */
	COLOR	tdiff;		/* diffuse transmission */
	COLOR	tunsamp;	/* BSDF hemispherical transmission */
}  BSDFDAT;		/* BSDF material data */

#define	cvt_sdcolor(cv, svp)	ccy2rgb(&(svp)->spec, (svp)->cieY, cv)

/* Compute "through" component color */
static void
compute_through(BSDFDAT *ndp)
{
#define NDIR2CHECK	13
	static const float	dir2check[NDIR2CHECK][2] = {
					{0, 0},
					{-0.8, 0},
					{0, 0.8},
					{0, -0.8},
					{0.8, 0},
					{-0.8, 0.8},
					{-0.8, -0.8},
					{0.8, 0.8},
					{0.8, -0.8},
					{-1.6, 0},
					{0, 1.6},
					{0, -1.6},
					{1.6, 0},
				};
	const double	peak_over = 1.5;
	SDSpectralDF	*dfp;
	FVECT		pdir;
	double		tomega, srchrad;
	COLOR		vpeak, vsum;
	int		i;
	SDError		ec;

	setcolor(ndp->cthru, 0, 0, 0);		/* starting assumption */

	if (ndp->pr->rod > 0)
		dfp = (ndp->sd->tf != NULL) ? ndp->sd->tf : ndp->sd->tb;
	else
		dfp = (ndp->sd->tb != NULL) ? ndp->sd->tb : ndp->sd->tf;

	if (dfp == NULL)
		return;				/* no specular transmission */
	if (bright(ndp->pr->pcol) <= FTINY)
		return;				/* pattern is black, here */
	srchrad = sqrt(dfp->minProjSA);		/* else search for peak */
	setcolor(vpeak, 0, 0, 0);
	setcolor(vsum, 0, 0, 0);
	pdir[2] = 0.0;
	for (i = 0; i < NDIR2CHECK; i++) {
		FVECT	tdir;
		SDValue	sv;
		COLOR	vcol;
		tdir[0] = -ndp->vray[0] + dir2check[i][0]*srchrad;
		tdir[1] = -ndp->vray[1] + dir2check[i][1]*srchrad;
		tdir[2] = -ndp->vray[2];
		normalize(tdir);
		ec = SDevalBSDF(&sv, tdir, ndp->vray, ndp->sd);
		if (ec)
			goto baderror;
		cvt_sdcolor(vcol, &sv);
		addcolor(vsum, vcol);
		if (sv.cieY > bright(vpeak)) {
			copycolor(vpeak, vcol);
			VCOPY(pdir, tdir);
		}
	}
	if (pdir[2] == 0.0)
		return;				/* zero neighborhood */
	ec = SDsizeBSDF(&tomega, pdir, ndp->vray, SDqueryMin, ndp->sd);
	if (ec)
		goto baderror;
	if (tomega > 1.5*dfp->minProjSA)
		return;				/* not really a peak? */
	tomega /= fabs(pdir[2]);		/* remove cosine factor */
	if ((bright(vpeak) - ndp->sd->tLamb.cieY*(1./PI))*tomega <= .001)
		return;				/* < 0.1% transmission */
	for (i = 3; i--; )			/* remove peak from average */
		colval(vsum,i) -= colval(vpeak,i);
	if (peak_over*bright(vsum) >= (NDIR2CHECK-1)*bright(vpeak))
		return;				/* not peaky enough */
	copycolor(ndp->cthru, vpeak);		/* else use it */
	scalecolor(ndp->cthru, tomega);
	multcolor(ndp->cthru, ndp->pr->pcol);	/* modify by pattern */
	return;
baderror:
	objerror(ndp->mp, USER, transSDError(ec));
#undef NDIR2CHECK
}

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

/* Get BSDF specular for direct component, returning true if OK to proceed */
static int
direct_specular_OK(COLOR cval, FVECT ldir, double omega, BSDFDAT *ndp)
{
	int	nsamp;
	double	wtot = 0;
	FVECT	vsrc, vsmp, vjit;
	double	tomega, tomega2;
	double	sf, tsr, sd[2];
	COLOR	csmp, cdiff;
	double	diffY;
	SDValue	sv;
	SDError	ec;
	int	i;
					/* in case we fail */
	setcolor(cval,  0, 0, 0);
					/* transform source direction */
	if (SDmapDir(vsrc, ndp->toloc, ldir) != SDEnone)
		return(0);
					/* will discount diffuse portion */
	switch ((vsrc[2] > 0)<<1 | (ndp->vray[2] > 0)) {
	case 3:
		if (ndp->sd->rf == NULL)
			return(0);	/* all diffuse */
		sv = ndp->sd->rLambFront;
		break;
	case 0:
		if (ndp->sd->rb == NULL)
			return(0);	/* all diffuse */
		sv = ndp->sd->rLambBack;
		break;
	default:
		if ((ndp->sd->tf == NULL) & (ndp->sd->tb == NULL))
			return(0);	/* all diffuse */
		sv = ndp->sd->tLamb;
		break;
	}
	if (sv.cieY > FTINY) {
		diffY = sv.cieY *= 1./PI;
		cvt_sdcolor(cdiff, &sv);
	} else {
		diffY = 0;
		setcolor(cdiff,  0, 0, 0);
	}
					/* need projected solid angles */
	omega *= fabs(vsrc[2]);
	ec = SDsizeBSDF(&tomega, ndp->vray, vsrc, SDqueryMin, ndp->sd);
	if (ec)
		goto baderror;
					/* check indirect over-counting */
	if ((vsrc[2] > 0) ^ (ndp->vray[2] > 0) && bright(ndp->cthru) > FTINY) {
		double	dx = vsrc[0] + ndp->vray[0];
		double	dy = vsrc[1] + ndp->vray[1];
		if (dx*dx + dy*dy <= (4./PI)*(omega + tomega +
						2.*sqrt(omega*tomega)))
			return(0);
	}
					/* assign number of samples */
	sf = specjitter * ndp->pr->rweight;
	if (tomega <= 0)
		nsamp = 1;
	else if (25.*tomega <= omega)
		nsamp = 100.*sf + .5;
	else
		nsamp = 4.*sf*omega/tomega + .5;
	nsamp += !nsamp;
	sf = sqrt(omega);		/* sample our source area */
	tsr = sqrt(tomega);
	for (i = nsamp; i--; ) {
		VCOPY(vsmp, vsrc);	/* jitter query directions */
		if (nsamp > 1) {
			multisamp(sd, 2, (i + frandom())/(double)nsamp);
			vsmp[0] += (sd[0] - .5)*sf;
			vsmp[1] += (sd[1] - .5)*sf;
			normalize(vsmp);
		}
		bsdf_jitter(vjit, ndp, tsr);
					/* compute BSDF */
		ec = SDevalBSDF(&sv, vjit, vsmp, ndp->sd);
		if (ec)
			goto baderror;
		if (sv.cieY - diffY <= FTINY)
			continue;	/* no specular part */
					/* check for variable resolution */
		ec = SDsizeBSDF(&tomega2, vjit, vsmp, SDqueryMin, ndp->sd);
		if (ec)
			goto baderror;
		if (tomega2 < .12*tomega)
			continue;	/* not safe to include */
		cvt_sdcolor(csmp, &sv);

		if (sf < 2.5*tsr) {	/* weight by Y for small sources */
			scalecolor(csmp, sv.cieY);
			wtot += sv.cieY;
		} else
			wtot += 1.;
		addcolor(cval, csmp);
	}
	if (wtot <= FTINY)		/* no valid specular samples? */
		return(0);

	sf = 1./wtot;			/* weighted average BSDF */
	scalecolor(cval, sf);
					/* subtract diffuse contribution */
	for (i = 3*(diffY > FTINY); i--; )
		if ((colval(cval,i) -= colval(cdiff,i)) < 0)
			colval(cval,i) = 0;
	return(1);
baderror:
	objerror(ndp->mp, USER, transSDError(ec));
	return(0);			/* gratis return */
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

	setcolor(cval,  0, 0, 0);

	ldot = DOT(np->pnorm, ldir);
	if ((-FTINY <= ldot) & (ldot <= FTINY))
		return;

	if (ldot > 0 && bright(np->rdiff) > FTINY) {
		/*
		 *  Compute diffuse reflected component
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ldot < 0 && bright(np->tdiff) > FTINY) {
		/*
		 *  Compute diffuse transmission
		 */
		copycolor(ctmp, np->tdiff);
		dtmp = -ldot * omega * (1.0/PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular scattering coefficient using BSDF
	 */
	if (!direct_specular_OK(ctmp, ldir, omega, np))
		return;
	if (ldot < 0) {		/* pattern for specular transmission */
		multcolor(ctmp, np->pr->pcol);
		dtmp = -ldot * omega;
	} else
		dtmp = ldot * omega;
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

	setcolor(cval,  0, 0, 0);

	ldot = DOT(np->pnorm, ldir);
	
	if (ldot <= FTINY)
		return;

	if (bright(np->rdiff) > FTINY) {
		/*
		 *  Compute diffuse reflected component
		 */
		copycolor(ctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular reflection coefficient using BSDF
	 */
	if (!direct_specular_OK(ctmp, ldir, omega, np))
		return;
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

	setcolor(cval,  0, 0, 0);

	ldot = DOT(np->pnorm, ldir);

	if (ldot >= -FTINY)
		return;

	if (bright(np->tdiff) > FTINY) {
		/*
		 *  Compute diffuse transmission
		 */
		copycolor(ctmp, np->tdiff);
		dtmp = -ldot * omega * (1.0/PI);
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular scattering coefficient using BSDF
	 */
	if (!direct_specular_OK(ctmp, ldir, omega, np))
		return;
					/* full pattern on transmission */
	multcolor(ctmp, np->pr->pcol);
	dtmp = -ldot * omega;
	scalecolor(ctmp, dtmp);
	addcolor(cval, ctmp);
}

/* Sample separate BSDF component */
static int
sample_sdcomp(BSDFDAT *ndp, SDComponent *dcp, int xmit)
{
	const int	hasthru = (xmit &&
				!(ndp->pr->crtype & (SPECULAR|AMBIENT)) &&
				bright(ndp->cthru) > FTINY);
	int		nstarget = 1;
	int		nsent = 0;
	int		n;
	SDError		ec;
	SDValue		bsv;
	double		xrand;
	FVECT		vsmp, vinc;
	RAY		sr;
						/* multiple samples? */
	if (specjitter > 1.5) {
		nstarget = specjitter*ndp->pr->rweight + .5;
		nstarget += !nstarget;
	}
						/* run through our samples */
	for (n = 0; n < nstarget; n++) {
		if (nstarget == 1) {		/* stratify random variable */
			xrand = urand(ilhash(dimlist,ndims)+samplendx);
			if (specjitter < 1.)
				xrand = .5 + specjitter*(xrand-.5);
		} else {
			xrand = (n + frandom())/(double)nstarget;
		}
		SDerrorDetail[0] = '\0';	/* sample direction & coef. */
		bsdf_jitter(vsmp, ndp, ndp->sr_vpsa[0]);
		VCOPY(vinc, vsmp);		/* to compare after */
		ec = SDsampComponent(&bsv, vsmp, xrand, dcp);
		if (ec)
			objerror(ndp->mp, USER, transSDError(ec));
		if (bsv.cieY <= FTINY)		/* zero component? */
			break;
		if (hasthru) {			/* check for view ray */
			double	dx = vinc[0] + vsmp[0];
			double	dy = vinc[1] + vsmp[1];
			if (dx*dx + dy*dy <= ndp->sr_vpsa[0]*ndp->sr_vpsa[0])
				continue;	/* exclude view sample */
		}
						/* map non-view sample->world */
		if (SDmapDir(sr.rdir, ndp->fromloc, vsmp) != SDEnone)
			break;
						/* spawn a specular ray */
		if (nstarget > 1)
			bsv.cieY /= (double)nstarget;
		cvt_sdcolor(sr.rcoef, &bsv);	/* use sample color */
		if (xmit)			/* apply pattern on transmit */
			multcolor(sr.rcoef, ndp->pr->pcol);
		if (rayorigin(&sr, SPECULAR, ndp->pr, sr.rcoef) < 0) {
			if (maxdepth > 0)
				break;
			continue;		/* Russian roulette victim */
		}
		if (xmit && ndp->thick != 0)	/* need to offset origin? */
			VSUM(sr.rorg, sr.rorg, ndp->pr->ron, -ndp->thick);
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
	int		hasthru = (sflags == SDsampSpT
				    && !(ndp->pr->crtype & (SPECULAR|AMBIENT))
				    && bright(ndp->cthru) > FTINY);
	int		n, ntotal = 0;
	double		b = 0;
	SDSpectralDF	*dfp;
	COLORV		*unsc;

	if (sflags == SDsampSpT) {
		unsc = ndp->tunsamp;
		if (ndp->pr->rod > 0)
			dfp = (ndp->sd->tf != NULL) ? ndp->sd->tf : ndp->sd->tb;
		else
			dfp = (ndp->sd->tb != NULL) ? ndp->sd->tb : ndp->sd->tf;
	} else /* sflags == SDsampSpR */ {
		unsc = ndp->runsamp;
		if (ndp->pr->rod > 0)
			dfp = ndp->sd->rf;
		else
			dfp = ndp->sd->rb;
	}
	setcolor(unsc,  0, 0, 0);
	if (dfp == NULL)			/* no specular component? */
		return(0);

	if (hasthru) {				/* separate view sample? */
		RAY	tr;
		if (rayorigin(&tr, TRANS, ndp->pr, ndp->cthru) == 0) {
			VCOPY(tr.rdir, ndp->pr->rdir);
			rayvalue(&tr);
			multcolor(tr.rcol, tr.rcoef);
			addcolor(ndp->pr->rcol, tr.rcol);
			++ntotal;
			b = bright(ndp->cthru);
		} else
			hasthru = 0;
	}
	if (dfp->maxHemi - b <= FTINY) {	/* have specular to sample? */
		b = 0;
	} else {
		FVECT	vjit;
		bsdf_jitter(vjit, ndp, ndp->sr_vpsa[1]);
		b = SDdirectHemi(vjit, sflags, ndp->sd) - b;
		if (b < 0) b = 0;
	}
	if (b <= specthresh+FTINY) {		/* below sampling threshold? */
		if (b > FTINY) {		/* XXX no color from BSDF */
			if (sflags == SDsampSpT) {
				copycolor(unsc, ndp->pr->pcol);
				scalecolor(unsc, b);
			} else			/* no pattern on reflection */
				setcolor(unsc, b, b, b);
		}
		return(ntotal);
	}
	dimlist[ndims] = (int)(size_t)ndp->mp;	/* else sample specular */
	ndims += 2;
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
	setfunc(m, r);
						/* get thickness */
	nd.thick = evalue(mf->ep[0]);
	if ((-FTINY <= nd.thick) & (nd.thick <= FTINY))
		nd.thick = 0;
						/* check backface visibility */
	if (!hitfront & !backvis) {
		raytrans(r);
		return(1);
	}
						/* check other rays to pass */
	if (nd.thick != 0 && (r->crtype & SHADOW ||
				!(r->crtype & (SPECULAR|AMBIENT)) ||
				(nd.thick > 0) ^ hitfront)) {
		raytrans(r);			/* hide our proxy */
		return(1);
	}
	nd.mp = m;
	nd.pr = r;
						/* get BSDF data */
	nd.sd = loadBSDF(m->oargs.sarg[1]);
						/* early shadow check */
	if (r->crtype & SHADOW && (nd.sd->tf == NULL) & (nd.sd->tb == NULL))
		return(1);
						/* diffuse reflectance */
	if (hitfront) {
		cvt_sdcolor(nd.rdiff, &nd.sd->rLambFront);
		if (m->oargs.nfargs >= 3) {
			setcolor(ctmp, m->oargs.farg[0],
					m->oargs.farg[1],
					m->oargs.farg[2]);
			addcolor(nd.rdiff, ctmp);
		}
	} else {
		cvt_sdcolor(nd.rdiff, &nd.sd->rLambBack);
		if (m->oargs.nfargs >= 6) {
			setcolor(ctmp, m->oargs.farg[3],
					m->oargs.farg[4],
					m->oargs.farg[5]);
			addcolor(nd.rdiff, ctmp);
		}
	}
						/* diffuse transmittance */
	cvt_sdcolor(nd.tdiff, &nd.sd->tLamb);
	if (m->oargs.nfargs >= 9) {
		setcolor(ctmp, m->oargs.farg[6],
				m->oargs.farg[7],
				m->oargs.farg[8]);
		addcolor(nd.tdiff, ctmp);
	}
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
	if (mf->fxp != &unitxf) {
		multv3(upvec, upvec, mf->fxp->xfm);
		nd.thick *= mf->fxp->sca;
	}
	if (r->rox != NULL) {
		multv3(upvec, upvec, r->rox->f.xfm);
		nd.thick *= r->rox->f.sca;
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
	if (ec) {
		objerror(m, WARNING, "Illegal orientation vector");
		return(1);
	}
	compute_through(&nd);			/* compute through component */
	if (r->crtype & SHADOW) {
		RAY	tr;			/* attempt to pass shadow ray */
		if (rayorigin(&tr, TRANS, r, nd.cthru) < 0)
			return(1);		/* blocked */
		VCOPY(tr.rdir, r->rdir);
		rayvalue(&tr);			/* transmit with scaling */
		multcolor(tr.rcol, tr.rcoef);
		copycolor(r->rcol, tr.rcol);
		return(1);			/* we're done */
	}
	ec = SDinvXform(nd.fromloc, nd.toloc);
	if (!ec)				/* determine BSDF resolution */
		ec = SDsizeBSDF(nd.sr_vpsa, nd.vray, NULL,
					SDqueryMin+SDqueryMax, nd.sd);
	if (ec)
		objerror(m, USER, transSDError(ec));

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
			VSUM(r->rop, vtmp, r->ron, nd.thick);
			multambient(ctmp, r, bnorm);
			VCOPY(r->rop, vtmp);
		} else
			multambient(ctmp, r, bnorm);
		addcolor(r->rcol, ctmp);
		if (hitfront)
			flipsurface(r);
	}
						/* add direct component */
	if ((bright(nd.tdiff) <= FTINY) & (nd.sd->tf == NULL) &
					(nd.sd->tb == NULL)) {
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
