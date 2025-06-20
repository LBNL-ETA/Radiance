#ifndef lint
static const char RCSid[] = "$Id: m_bsdf.c,v 2.76 2025/06/20 17:03:01 greg Exp $";
#endif
/*
 *  Shading for materials with BSDFs taken from XML data files
 */

#include "copyright.h"

#include  "ray.h"
#include  "otypes.h"
#include  "ambient.h"
#include  "source.h"
#include  "func.h"
#include  "bsdf.h"
#include  "random.h"
#include  "pmapmat.h"

/*
 *	Arguments to this material include optional diffuse colors.
 *  String arguments include the BSDF and function files.
 *	For the MAT_BSDF type, a non-zero thickness causes the useful behavior
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
 *	For the MAT_ABSDF type, we check for a strong "through" component.
 *  Such a component will cause direct rays to pass through unscattered.
 *  A separate test prevents over-counting by dropping samples that are
 *  too close to this "through" direction.  BSDFs with such a through direction
 *  will also have a view component, meaning they are somewhat see-through.
 *  A MAT_BSDF type with zero thickness behaves the same as a MAT_ABSDF
 *  type with no strong through component.
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
 *  Arguments for MAT_ABSDF are:
 *	5+	BSDFfile	ux uy uz	funcfile	transform
 *	0
 *	0|3|6|9	rdf	gdf	bdf
 *		rdb	gdb	bdb
 *		rdt	gdt	bdt
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
	SCOLOR	cthru;		/* "through" component for MAT_ABSDF */
	SCOLOR	cthru_surr;	/* surround for "through" component */
	SDData	*sd;		/* loaded BSDF data */
	SCOLOR	rdiff;		/* diffuse reflection */
	SCOLOR	runsamp;	/* BSDF hemispherical reflection */
	SCOLOR	tdiff;		/* diffuse transmission */
	SCOLOR	tunsamp;	/* BSDF hemispherical transmission */
}  BSDFDAT;		/* BSDF material data */

#define	cvt_sdcolor(scv, svp)	ccy2scolor(&(svp)->spec, (svp)->cieY, scv)

typedef struct {
	double	vy;		/* brightness (for sorting) */
	FVECT	tdir;		/* through sample direction (normalized) */
	SCOLOR	vcol;		/* BTDF color */
}  PEAKSAMP;		/* BTDF peak sample */

/* Comparison function to put near-peak values in descending order */
static int
cmp_psamp(const void *p1, const void *p2)
{
	double	diff = (*(const PEAKSAMP *)p1).vy - (*(const PEAKSAMP *)p2).vy;
	if (diff > 0) return(-1);
	if (diff < 0) return(1);
	return(0);
}

/* Compute "through" component color for MAT_ABSDF */
static void
compute_through(BSDFDAT *ndp)
{
#define NDIR2CHECK	29
	static const float	dir2check[NDIR2CHECK][2] = {
					{0, 0}, {-0.6, 0}, {0, 0.6},
					{0, -0.6}, {0.6, 0}, {-0.6, 0.6},
					{-0.6, -0.6}, {0.6, 0.6}, {0.6, -0.6},
					{-1.2, 0}, {0, 1.2}, {0, -1.2},
					{1.2, 0}, {-1.2, 1.2}, {-1.2, -1.2},
					{1.2, 1.2}, {1.2, -1.2}, {-1.8, 0},
					{0, 1.8}, {0, -1.8}, {1.8, 0},
					{-1.8, 1.8}, {-1.8, -1.8}, {1.8, 1.8},
					{1.8, -1.8}, {-2.4, 0}, {0, 2.4},
					{0, -2.4}, {2.4, 0},
				};
	PEAKSAMP	psamp[NDIR2CHECK];
	SDSpectralDF	*dfp;
	FVECT		pdir;
	double		tomega, srchrad;
	double		tomsum, tomsurr;
	SCOLOR		vpeak, vsurr, btdiff;
	double		vypeak;
	int		i, ns;
	SDError		ec;

	if (ndp->pr->rod > 0)
		dfp = (ndp->sd->tf != NULL) ? ndp->sd->tf : ndp->sd->tb;
	else
		dfp = (ndp->sd->tb != NULL) ? ndp->sd->tb : ndp->sd->tf;

	if (dfp == NULL)
		return;				/* no specular transmission */
	if (sintens(ndp->pr->pcol) <= FTINY)
		return;				/* pattern is black, here */
	srchrad = sqrt(dfp->minProjSA);		/* else evaluate peak */
	for (i = 0; i < NDIR2CHECK; i++) {
		SDValue	sv;
		psamp[i].tdir[0] = -ndp->vray[0] + dir2check[i][0]*srchrad;
		psamp[i].tdir[1] = -ndp->vray[1] + dir2check[i][1]*srchrad;
		psamp[i].tdir[2] = -ndp->vray[2];
		normalize(psamp[i].tdir);
		ec = SDevalBSDF(&sv, ndp->vray, psamp[i].tdir, ndp->sd);
		if (ec)
			goto baderror;
		cvt_sdcolor(psamp[i].vcol, &sv);
		psamp[i].vy = sv.cieY;
	}
	qsort(psamp, NDIR2CHECK, sizeof(PEAKSAMP), cmp_psamp);
	if (psamp[0].vy <= FTINY)
		return;				/* zero BTDF here */
	scolorblack(vpeak);
	scolorblack(vsurr);
	vypeak = tomsum = tomsurr = 0;		/* combine top unique values */
	ns = 0;
	for (i = 0; i < NDIR2CHECK; i++) {
		if (i && psamp[i].vy == psamp[i-1].vy)
			continue;		/* assume duplicate sample */

		ec = SDsizeBSDF(&tomega, ndp->vray, psamp[i].tdir,
						SDqueryMin, ndp->sd);
		if (ec)
			goto baderror;

		scalescolor(psamp[i].vcol, tomega);
						/* not part of peak? */
		if (tomega > 1.5*dfp->minProjSA ||
					vypeak > 8.*psamp[i].vy*ns) {
			if (!i) return;		/* abort */
			saddscolor(vsurr, psamp[i].vcol);
			tomsurr += tomega;
			continue;
		}
		saddscolor(vpeak, psamp[i].vcol);
		tomsum += tomega;
		vypeak += psamp[i].vy;
		++ns;
	}
	if (tomsurr < 0.2*tomsum)		/* insufficient surround? */
		return;
	scalescolor(vsurr, 1./tomsurr);		/* surround is avg. BTDF */
	if (ndp->vray[2] > 0)			/* get diffuse BTDF */
		cvt_sdcolor(btdiff, &ndp->sd->tLambFront);
	else
		cvt_sdcolor(btdiff, &ndp->sd->tLambBack);
	scalescolor(btdiff, (1./PI));
	for (i = NCSAMP; i--; ) {		/* remove diffuse contrib. */
		if ((vpeak[i] -= tomsum*btdiff[i]) < 0)
			vpeak[i] = 0;
		if ((vsurr[i] -= btdiff[i]) < 0)
			vsurr[i] = 0;
	}
	if (pbright(vpeak) < .0005)		/* < 0.05% specular? */
		return;
	smultscolor(vsurr, ndp->pr->pcol);	/* modify by pattern */
	smultscolor(vpeak, ndp->pr->pcol);
	copyscolor(ndp->cthru_surr, vsurr);
	copyscolor(ndp->cthru, vpeak);
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
direct_specular_OK(SCOLOR scval, FVECT ldir, double omega, BSDFDAT *ndp)
{
	int	nsamp = 1;
	int	scnt = 0;
	FVECT	vsrc;
	double	tomega, tomega2;
	double	tsr, sd[2];
	SCOLOR	csmp, cdiff;
	double	diffY;
	SDValue	sv;
	SDError	ec;
	int	i;
					/* in case we fail */
	scolorblack(scval);
					/* transform source direction */
	if (SDmapDir(vsrc, ndp->toloc, ldir) != SDEnone)
		return(0);
					/* check indirect over-counting */
	if ((vsrc[2] > 0) ^ (ndp->vray[2] > 0) && sintens(ndp->cthru) > FTINY) {
		double		dx = vsrc[0] + ndp->vray[0];
		double		dy = vsrc[1] + ndp->vray[1];
		SDSpectralDF	*dfp = (ndp->pr->rod > 0) ?
			((ndp->sd->tf != NULL) ? ndp->sd->tf : ndp->sd->tb) :
			((ndp->sd->tb != NULL) ? ndp->sd->tb : ndp->sd->tf) ;

		tomega = omega*fabs(vsrc[2]);
		if (dx*dx + dy*dy <= (2.5*4./PI)*(tomega + dfp->minProjSA +
						2.*sqrt(tomega*dfp->minProjSA))) {
			if (sintens(ndp->cthru_surr) <= FTINY)
				return(0);
			copyscolor(scval, ndp->cthru_surr);
			return(1);	/* return non-zero surround BTDF */
		}
	}
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
	case 1:
		if ((ndp->sd->tf == NULL) & (ndp->sd->tb == NULL))
			return(0);	/* all diffuse */
		sv = ndp->sd->tLambFront;
		break;
	case 2:
		if ((ndp->sd->tf == NULL) & (ndp->sd->tb == NULL))
			return(0);	/* all diffuse */
		sv = ndp->sd->tLambBack;
		break;
	}
	if (sv.cieY > FTINY) {
		diffY = sv.cieY *= 1./PI;
		cvt_sdcolor(cdiff, &sv);
	} else {
		diffY = 0;
		scolorblack(cdiff);
	}
	ec = SDsizeBSDF(&tomega, ndp->vray, vsrc, SDqueryMin, ndp->sd);
	if (ec)
		goto baderror;
					/* check if sampling BSDF */
	if ((tsr = sqrt(tomega)) > 0) {
		nsamp = 4.*specjitter*ndp->pr->rweight + .5;
		nsamp += !nsamp;
	}
	for (i = nsamp; i--; ) {	/* jitter to fuzz BSDF cells */
		FVECT	vjit;
		bsdf_jitter(vjit, ndp, tsr);
					/* compute BSDF */
		ec = SDevalBSDF(&sv, vjit, vsrc, ndp->sd);
		if (ec)
			goto baderror;
		if (sv.cieY - diffY <= FTINY) {
			++scnt;		/* still counts as 0 contribution */
			continue;
		}
					/* check for variable resolution */
		ec = SDsizeBSDF(&tomega2, vjit, vsrc, SDqueryMin, ndp->sd);
		if (ec)
			goto baderror;
		if (tomega2 < .12*tomega)
			continue;	/* not safe to include */
		cvt_sdcolor(csmp, &sv);
		saddscolor(scval, csmp);
		++scnt;
	}
	if (!scnt)			/* no valid specular samples? */
		return(0);

	scalescolor(scval, 1./scnt);	/* weighted average BSDF */
					/* subtract diffuse contribution */
	for (i = NCSAMP*(diffY > FTINY); i--; )
		if ((scval[i] -= cdiff[i]) < 0)
			scval[i] = 0;
	return(1);
baderror:
	objerror(ndp->mp, USER, transSDError(ec));
	return(0);			/* gratis return */
}

/* Compute source contribution for BSDF (reflected & transmitted) */
static void
dir_bsdf(
	SCOLOR  scval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	SCOLOR		sctmp;

	scolorblack(scval);

	ldot = DOT(np->pnorm, ldir);
	if ((-FTINY <= ldot) & (ldot <= FTINY))
		return;

	if (ldot > 0 && sintens(np->rdiff) > FTINY) {
		/*
		 *  Compute diffuse reflected component
		 */
		copyscolor(sctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
	if (ldot < 0 && sintens(np->tdiff) > FTINY) {
		/*
		 *  Compute diffuse transmission
		 */
		copyscolor(sctmp, np->tdiff);
		dtmp = -ldot * omega * (1./PI);
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular scattering coefficient using BSDF
	 */
	if (!direct_specular_OK(sctmp, ldir, omega, np))
		return;
	if (ldot < 0) {		/* pattern for specular transmission */
		smultscolor(sctmp, np->pr->pcol);
		dtmp = -ldot * omega;
	} else
		dtmp = ldot * omega;
	scalescolor(sctmp, dtmp);
	saddscolor(scval, sctmp);
}

/* Compute source contribution for BSDF (reflected only) */
static void
dir_brdf(
	SCOLOR  scval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	SCOLOR		sctmp;

	scolorblack(scval);

	ldot = DOT(np->pnorm, ldir);
	
	if (ldot <= FTINY)
		return;

	if (sintens(np->rdiff) > FTINY) {
		/*
		 *  Compute diffuse reflected component
		 */
		copyscolor(sctmp, np->rdiff);
		dtmp = ldot * omega * (1./PI);
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular reflection coefficient using BSDF
	 */
	if (!direct_specular_OK(sctmp, ldir, omega, np))
		return;
	dtmp = ldot * omega;
	scalescolor(sctmp, dtmp);
	saddscolor(scval, sctmp);
}

/* Compute source contribution for BSDF (transmitted only) */
static void
dir_btdf(
	SCOLOR  scval,			/* returned coefficient */
	void  *nnp,			/* material data */
	FVECT  ldir,			/* light source direction */
	double  omega			/* light source size */
)
{
	BSDFDAT		*np = (BSDFDAT *)nnp;
	double		ldot;
	double		dtmp;
	SCOLOR		sctmp;

	scolorblack(scval);

	ldot = DOT(np->pnorm, ldir);

	if (ldot >= -FTINY)
		return;

	if (sintens(np->tdiff) > FTINY) {
		/*
		 *  Compute diffuse transmission
		 */
		copyscolor(sctmp, np->tdiff);
		dtmp = -ldot * omega * (1./PI);
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
	if (ambRayInPmap(np->pr))
		return;		/* specular already in photon map */
	/*
	 *  Compute specular scattering coefficient using BSDF
	 */
	if (!direct_specular_OK(sctmp, ldir, omega, np))
		return;
					/* full pattern on transmission */
	smultscolor(sctmp, np->pr->pcol);
	dtmp = -ldot * omega;
	scalescolor(sctmp, dtmp);
	saddscolor(scval, sctmp);
}

/* Sample separate BSDF component */
static int
sample_sdcomp(BSDFDAT *ndp, SDComponent *dcp, int xmit)
{
	const int	hasthru = (xmit &&
					!(ndp->pr->crtype & (SPECULAR|AMBIENT))
					&& sintens(ndp->cthru) > FTINY);
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
			smultscolor(sr.rcoef, ndp->pr->pcol);
		if (rayorigin(&sr, xmit ? TSPECULAR : RSPECULAR, ndp->pr, sr.rcoef) < 0) {
			if (!n & (nstarget > 1)) {
				n = nstarget;	/* avoid infinitue loop */
				nstarget = nstarget*sr.rweight/(minweight + 1e-20);
				if (n == nstarget) break;
				n = -1;		/* moved target */
			}
			continue;		/* try again */
		}
		if (xmit && ndp->thick != 0)	/* need to offset origin? */
			VSUM(sr.rorg, sr.rorg, ndp->pr->ron, -ndp->thick);
		rayvalue(&sr);			/* send & evaluate sample */
		smultscolor(sr.rcol, sr.rcoef);
		saddscolor(ndp->pr->rcol, sr.rcol);
		++nsent;
	}
	return(nsent);
}

/* Sample non-diffuse components of BSDF */
static int
sample_sdf(BSDFDAT *ndp, int sflags)
{
	int		hasthru = (sflags == SDsampSpT &&
					!(ndp->pr->crtype & (SPECULAR|AMBIENT))
					&& sintens(ndp->cthru) > FTINY);
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
	scolorblack(unsc);
	if (dfp == NULL)			/* no specular component? */
		return(0);

	if (hasthru) {				/* separate view sample? */
		RAY	tr;
		if (rayorigin(&tr, TRANS, ndp->pr, ndp->cthru) == 0) {
			VCOPY(tr.rdir, ndp->pr->rdir);
			rayvalue(&tr);
			smultscolor(tr.rcol, tr.rcoef);
			saddscolor(ndp->pr->rcol, tr.rcol);
			ndp->pr->rxt = ndp->pr->rot + raydistance(&tr);
			++ntotal;
			b = pbright(ndp->cthru);
		} else
			hasthru = 0;
	}
	if (dfp->maxHemi - b <= FTINY) {	/* have specular to sample? */
		b = 0;
	} else {
		FVECT	vjit;
		bsdf_jitter(vjit, ndp, ndp->sr_vpsa[1]);
		b = SDdirectHemi(vjit, sflags, ndp->sd) - b;
		b *= (b > 0);
	}
	if (b <= specthresh+FTINY) {		/* below sampling threshold? */
		if (b > FTINY) {		/* XXX no color from BSDF */
			if (sflags == SDsampSpT) {
				copyscolor(unsc, ndp->pr->pcol);
				scalescolor(unsc, b);
			} else			/* no pattern on reflection */
				setscolor(unsc, b, b, b);
		}
		return(ntotal);
	}
	dimlist[ndims_inc] = (int)(size_t)ndp->mp;	/* else sample specular */
	ndims_inc;
	for (n = dfp->ncomp; n--; ) {		/* loop over components */
		dimlist[ndims-1] = n + 9438;
		ntotal += sample_sdcomp(ndp, &dfp->comp[n], sflags==SDsampSpT);
	}
	dec_ndims; dec_ndims;
	return(ntotal);
}

/* Color a ray that hit a BSDF material */
int
m_bsdf(OBJREC *m, RAY *r)
{
	const int	hasthick = (m->otype == MAT_BSDF);
	const int	hitfront = (r->rod > 0);
	SCOLOR	sctmp;
	SDError	ec;
	FVECT	upvec, vtmp;
	MFUNC	*mf;
	BSDFDAT	nd;
						/* check backface visibility */
	if (!hitfront & !backvis) {
		raytrans(r);
		return(1);
	}
						/* check arguments */
	if ((m->oargs.nsargs < hasthick+5) | (m->oargs.nfargs > 9) |
				(m->oargs.nfargs % 3))
		objerror(m, USER, "bad # arguments");
						/* load cal file */
	mf = hasthick	? getfunc(m, 5, 0x1d, 1)
			: getfunc(m, 4, 0xe, 1) ;
	setfunc(m, r);
	nd.thick = 0;				/* set thickness */
	if (hasthick) {
		nd.thick = evalue(mf->ep[0]);
		if ((-FTINY <= nd.thick) & (nd.thick <= FTINY))
			nd.thick = 0;
	}
						/* check other rays to pass */
	if (nd.thick != 0 && (r->crtype & SHADOW ||
				!(r->crtype & (SPECULAR|AMBIENT)) ||
				(nd.thick > 0) ^ hitfront)) {
		raytrans(r);			/* hide our proxy */
		return(1);
	}
	if (hasthick && r->crtype & SHADOW)	/* early shadow check #1 */
		return(1);
	nd.mp = m;
	nd.pr = r;
						/* get BSDF data */
	nd.sd = loadBSDF(m->oargs.sarg[hasthick]);
						/* early shadow check #2 */
	if (r->crtype & SHADOW && (nd.sd->tf == NULL) & (nd.sd->tb == NULL)) {
		SDfreeCache(nd.sd);
		return(1);
	}
						/* diffuse components */
	if (hitfront) {
		cvt_sdcolor(nd.rdiff, &nd.sd->rLambFront);
		if (m->oargs.nfargs >= 3) {
			setscolor(sctmp, m->oargs.farg[0],
					m->oargs.farg[1],
					m->oargs.farg[2]);
			saddscolor(nd.rdiff, sctmp);
		}
		cvt_sdcolor(nd.tdiff, &nd.sd->tLambFront);
	} else {
		cvt_sdcolor(nd.rdiff, &nd.sd->rLambBack);
		if (m->oargs.nfargs >= 6) {
			setscolor(sctmp, m->oargs.farg[3],
					m->oargs.farg[4],
					m->oargs.farg[5]);
			saddscolor(nd.rdiff, sctmp);
		}
		cvt_sdcolor(nd.tdiff, &nd.sd->tLambBack);
	}
	if (m->oargs.nfargs >= 9) {		/* add diffuse transmittance? */
		setscolor(sctmp, m->oargs.farg[6],
				m->oargs.farg[7],
				m->oargs.farg[8]);
		saddscolor(nd.tdiff, sctmp);
	}
						/* get modifiers */
	raytexture(r, m->omod);
						/* modify diffuse values */
	smultscolor(nd.rdiff, r->pcol);
	smultscolor(nd.tdiff, r->pcol);
						/* get up vector */
	upvec[0] = evalue(mf->ep[hasthick+0]);
	upvec[1] = evalue(mf->ep[hasthick+1]);
	upvec[2] = evalue(mf->ep[hasthick+2]);
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
		SDfreeCache(nd.sd);
		return(1);
	}
	scolorblack(nd.cthru);			/* consider through component */
	scolorblack(nd.cthru_surr);
	if (m->otype == MAT_ABSDF) {
		compute_through(&nd);
		if (r->crtype & SHADOW) {
			RAY	tr;		/* attempt to pass shadow ray */
			SDfreeCache(nd.sd);
			if (rayorigin(&tr, TRANS, r, nd.cthru) < 0)
				return(1);	/* no through component */
			VCOPY(tr.rdir, r->rdir);
			rayvalue(&tr);		/* transmit with scaling */
			smultscolor(tr.rcol, tr.rcoef);
			copyscolor(r->rcol, tr.rcol);
			return(1);		/* we're done */
		}
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
	copyscolor(sctmp, nd.rdiff);
	saddscolor(sctmp, nd.runsamp);
	if (sintens(sctmp) > FTINY) {		/* ambient from reflection */
		multambient(sctmp, r, nd.pnorm);
		saddscolor(r->rcol, sctmp);
	}
	copyscolor(sctmp, nd.tdiff);
	saddscolor(sctmp, nd.tunsamp);
	if (sintens(sctmp) > FTINY) {		/* ambient from other side */
		FVECT  bnorm;
		bnorm[0] = -nd.pnorm[0];
		bnorm[1] = -nd.pnorm[1];
		bnorm[2] = -nd.pnorm[2];
		if (nd.thick != 0) {		/* proxy with offset? */
			VCOPY(vtmp, r->rop);
			VSUM(r->rop, vtmp, r->ron, nd.thick);
			multambient(sctmp, r, bnorm);
			VCOPY(r->rop, vtmp);
		} else
			multambient(sctmp, r, bnorm);
		saddscolor(r->rcol, sctmp);
	}
						/* add direct component */
	if ((nd.sd->tf == NULL) & (nd.sd->tb == NULL) &&
			sintens(nd.tdiff) <= FTINY) {
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
