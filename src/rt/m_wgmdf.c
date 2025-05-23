#ifndef lint
static const char RCSid[] = "$Id: m_wgmdf.c,v 2.10 2025/05/23 17:09:26 greg Exp $";
#endif
/*
 *  Shading function for programmable Ward-Geisler-Moroder-Duer material.
 */

#include "copyright.h"

#include  "ray.h"
#include  "ambient.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "source.h"
#include  "func.h"
#include  "random.h"
#include  "pmapmat.h"

#ifndef  MAXITER
#define  MAXITER	10		/* maximum # specular ray attempts */
#endif
					/* estimate of Fresnel function */
#define  FRESNE(ci)	(exp(-5.85*(ci)) - 0.00202943064)
#define  FRESTHRESH	0.017999	/* minimum specularity for approx. */

/*
 *	This routine implements the anisotropic Gaussian
 *  model described by Ward in a 1992 Siggraph article and updated by
 *  Geisler-Moroder and Duer in a 2010 article in High Performance Graphics.
 *	We do not reorient incoming ray, using side in part to determine
 *  reflectance values.  Most parameters are programmable with their own
 *  modifiers and/or value expressions.
 *
 *  Arguments for MAT_WGMDF are:
 *	13+	rs_mod  rs  rs_urough rs_vrough
 *		ts_mod  ts  ts_urough ts_vrough
 *		td_mod
 *		ux uy uz  funcfile  transform
 *	0
 *	9+	rfdif gfdif bfdif
 *		rbdif gbdif bbdif
 *		rtdif gtdif btdif
 *		A10 ..
 *
 *  Where the rs_urough or rs_vrough expression yields zero, mirror-Fresnel
 *  effects are computed, similar to MAT_PLASTIC and MAT_METAL.  The
 *  rs* expressions should not vary with incident angle, or the material
 *  will not be physically valid.  Similarly, the ts* expressions should
 *  give the same value for coincident direction vectors from either side.
 *	There are independent modifiers for specular reflection,
 *  transmission, and diffuse transmission.  Diffuse reflection
 *  applies the material's main modifier, which doesn't apply to
 *  anything else by default.  However, any of the modifiers may be
 *  ALIASMOD, which will use the main material modifier, or VOIDID,
 *  which will just be white.
 *	Diffuse reflection and transmission colors and patterns add to
 *  the specular components, and are only adjusted with mirror-Fresnel
 *  reflection if specular reflection is greater than FRESHTHRESH.  The
 *  specular transmission is likewise adjusted in such cases.  Specified
 *  values for all components should sum to less than 1, but like other
 *  Radiance materials, this is not enforced, nor is a warning issued.
 */
				/* specularity flags */
#define  SP_REFL	01		/* has reflected specular component */
#define  SP_TRAN	02		/* has transmitted specular */
#define  SP_RPURE	04		/* mirror reflection */
#define  SP_TPURE	010		/* has view component */
#define  SP_FLAT	020		/* flat reflecting surface */
#define  SP_RBLT	040		/* reflection below sample threshold */
#define  SP_TBLT	0100		/* transmission below threshold */

typedef struct {
	char		*nam;		/* modifier name */
	int		hastexture;	/* has a texture? */
	FVECT		pnorm;		/* perturbed normal direction */
	double  	pdot;		/* perturbed dot product */
	SCOLOR		pcol;		/* pattern color */
} MODVAL;		/* modifier-derived values */

typedef struct {
	MODVAL		mo;		/* modifier parameters */
	SCOLOR		scol;		/* modified diffuse color */
} DCOMP;		/* diffuse component parameters */

typedef struct {
	MODVAL		mo;		/* modifier parameters */
	SCOLOR		scol;		/* modified specular color */
	FVECT		u, v;		/* u and v in-plane vectors */
	double  	u_alpha;	/* u roughness */
	double  	v_alpha;	/* v roughness */
} SCOMP;		/* specular component parameters */

typedef struct {
	RAY		*rp;		/* ray pointer */
	OBJREC		*mtp;		/* material pointer */
	MFUNC		*mf;		/* pointer to expression list */
	int		specfl;		/* specularity flags, defined above */
	FVECT		ulocal;		/* u-vector in local coordinates */
	DCOMP		rd, td;		/* diffuse component params */
	SCOMP		rs, ts;		/* specular component params */
	FVECT		prdir;		/* vector in transmitted direction */
} WGMDDAT;		/* WGMD material data */

#define clr_comps(wp)	((wp)->specfl = 0, \
			(wp)->rd.mo.nam = (wp)->td.mo.nam = \
			(wp)->rs.mo.nam = (wp)->ts.mo.nam = "")

/* assign modifier values */
static int
set_modval(MODVAL *mp, OBJECT omod, const RAY *r)
{
	RAY	tr;

	if (!mp->nam[0])
		mp->nam = (omod == OVOID) ? VOIDID : objptr(omod)->oname;
	else if (!strcmp(mp->nam, VOIDID))
		omod = OVOID;
	else if (omod == OVOID)
		return(0);
	tr = *r;			/* independent modifier */
	raytexture(&tr, omod);
	if (DOT(tr.pert,tr.pert) > FTINY*FTINY) {
		mp->pdot = raynormal(mp->pnorm, &tr);
		mp->hastexture = 1;
	} else {
		VCOPY(mp->pnorm, tr.ron);
		mp->pdot = tr.rod;
		mp->hastexture = 0;
	}
	copyscolor(mp->pcol, tr.pcol);
	return(1);
}

/* fill modifier values, using previous setting if found */
static int
fill_modval(MODVAL *mp, const WGMDDAT *wp)
{
	if (mp == &wp->rd.mo) {		/* special case (should be first) */
		set_modval(mp, wp->mtp->omod, wp->rp);
		return(1);
	}				/* use main modifier? */
	if (!strcmp(mp->nam, ALIASMOD) || !strcmp(mp->nam, wp->rd.mo.nam)) {
		*mp = wp->rd.mo;
		return(1);
	}				/* check others */
	if (mp != &wp->td.mo && !strcmp(mp->nam, wp->td.mo.nam)) {
		*mp = wp->td.mo;
		return(1);
	}
	if (mp != &wp->rs.mo && !strcmp(mp->nam, wp->rs.mo.nam)) {
		*mp = wp->rs.mo;
		return(1);
	}
	if (mp != &wp->ts.mo && !strcmp(mp->nam, wp->ts.mo.nam)) {
		*mp = wp->ts.mo;
		return(1);
	}				/* new modifier */
	return(set_modval(mp, lastmod(objndx(wp->mtp), mp->nam), wp->rp));
}

/* set calculation context for given component of MAT_WGMDF */
static int
setWGMDfunc(MODVAL *mp, const WGMDDAT *wp)
{
	static char	lastMod[MAXSTR];
	double		sf;
	FVECT		vec;

	if (setfunc(wp->mtp, wp->rp) == 0 &&
			!strcmp(mp->nam, lastMod))
		return(0);	/* already set */
	strcpy(lastMod, mp->nam);
				/* else (re)assign special variables */
	sf = 1 - 2*(wp->rp->rod < 0);
	varset("RdotP`", '=', mp->pdot*sf);
	multv3(vec, mp->pnorm, funcxf.xfm);
	sf /= funcxf.sca;
	varset("NxP`", '=', vec[0]*sf);
	varset("NyP`", '=', vec[1]*sf);
	varset("NzP`", '=', vec[2]*sf);
	return(1);
}

/* assign indicated diffuse component (do !trans first) */
static void
set_dcomp(WGMDDAT *wp, int trans)
{
	DCOMP		*dp = trans ? &wp->td : &wp->rd;
	const int	offs = trans ? 6 : 3*(wp->rp->rod < 0);

	if (trans) {			/* transmitted diffuse? */
		if (intens(wp->mtp->oargs.farg+offs) <= FTINY) {
			scolorblack(dp->scol);
			return;
		}
		dp->mo.nam = wp->mtp->oargs.sarg[8];
		if (!fill_modval(&dp->mo, wp)) {
			sprintf(errmsg,
			"unknown diffuse transmission modifier '%s'",
					dp->mo.nam);
			objerror(wp->mtp, USER, errmsg);
		}
	} else				/* no priors for main mod */
		fill_modval(&dp->mo, wp);

	setscolor(dp->scol, wp->mtp->oargs.farg[offs],
			wp->mtp->oargs.farg[offs+1],
			wp->mtp->oargs.farg[offs+2]);
	smultscolor(dp->scol, dp->mo.pcol);
}

/* assign indicated specular component */
static void
set_scomp(WGMDDAT *wp, int trans)
{
	SCOMP	*sp = trans ? &wp->ts : &wp->rs;
	EPNODE  **exa = wp->mf->ep + 3*(trans != 0);
	double	coef;
					/* constant zero check */
	if (exa[0]->type == NUM && exa[0]->v.num <= FTINY)
		goto blackout;
					/* need modifier */
	sp->mo.nam = wp->mtp->oargs.sarg[4*(trans != 0)];
	if (!fill_modval(&sp->mo, wp)) {
		sprintf(errmsg, "unknown specular %s modifier '%s'",
			trans ? "transmission" : "reflection", sp->mo.nam);
		objerror(wp->mtp, USER, errmsg);
	}
	if (sintens(sp->mo.pcol) <= FTINY)
		goto blackout;		/* got black pattern */
	setWGMDfunc(&sp->mo, wp);	/* else compute coefficient */
	errno = 0;
	coef = evalue(exa[0]);
	if ((errno == EDOM) | (errno == ERANGE)) {
		objerror(wp->mtp, WARNING, "specular compute error");
		goto blackout;
	}
	if (coef <= FTINY)		/* negligible value? */
		goto blackout;
	copyscolor(sp->scol, sp->mo.pcol);
	scalescolor(sp->scol, coef);
	errno = 0;			/* else get roughness */
	sp->u_alpha = evalue(exa[1]);
	sp->v_alpha = (sp->u_alpha > FTINY) ? evalue(exa[2]) : 0.0;
	if ((errno == EDOM) | (errno == ERANGE)) {
		objerror(wp->mtp, WARNING, "roughness compute error");
		goto blackout;
	}				/* we have something... */
	wp->specfl |= trans ? SP_TRAN : SP_REFL;
	if (sp->v_alpha <= FTINY) {	/* is it pure specular? */
		wp->specfl |= trans ? SP_TPURE : SP_RPURE;
		sp->u_alpha = sp->v_alpha = 0.0;
		return;
	}				/* else get aniso coordinates */
	fcross(sp->v, sp->mo.pnorm, wp->ulocal);
	if (normalize(sp->v) == 0.0) {	/* orientation vector==normal? */
		if (fabs(sp->u_alpha - sp->v_alpha) > 0.001)
			objerror(wp->mtp, WARNING, "bad orientation vector");
		getperpendicular(sp->u, sp->mo.pnorm, 0);	/* punting */
		fcross(sp->v, sp->mo.pnorm, sp->u);
		sp->u_alpha = sp->v_alpha = sqrt( 0.5 *
			(sp->u_alpha*sp->u_alpha + sp->v_alpha*sp->v_alpha) );
	} else
		fcross(sp->u, sp->v, sp->mo.pnorm);
	return;
blackout:
	scolorblack(sp->scol);		/* zero out component */
}

/* sample anisotropic Gaussian specular */
static void
agaussamp(WGMDDAT *wp)
{
	RAY	sr;
	FVECT	h;
	double	rv[2];
	double	d, sinp, cosp;
	int	maxiter, ntrials, nstarget, nstaken;
	int	i;
					/* compute reflection */
	if ((wp->specfl & (SP_REFL|SP_RPURE|SP_RBLT)) == SP_REFL &&
			rayorigin(&sr, RSPECULAR, wp->rp, wp->rs.scol) == 0) {
		SCOLOR	scol;
		nstarget = 1;
		if (specjitter > 1.5) {	/* multiple samples? */
			nstarget = specjitter*wp->rp->rweight + .5;
			if (sr.rweight <= minweight*nstarget)
				nstarget = sr.rweight/minweight;
			if (nstarget > 1) {
				d = 1./nstarget;
				scalescolor(sr.rcoef, d);
				sr.rweight *= d;
			} else
				nstarget = 1;
		}
		scolorblack(scol);
		dimlist[ndims++] = (int)(size_t)wp->mtp;
		maxiter = MAXITER*nstarget;
		for (nstaken = ntrials = 0; (nstaken < nstarget) &
						(ntrials < maxiter); ntrials++) {
			if (ntrials)
				d = frandom();
			else
				d = urand(ilhash(dimlist,ndims)+samplendx);
			multisamp(rv, 2, d);
			d = 2.0*PI * rv[0];
			cosp = tcos(d) * wp->rs.u_alpha;
			sinp = tsin(d) * wp->rs.v_alpha;
			d = 1./sqrt(cosp*cosp + sinp*sinp);
			cosp *= d;
			sinp *= d;
			if ((0. <= specjitter) & (specjitter < 1.))
				rv[1] = 1.0 - specjitter*rv[1];
			d = (rv[1] <= FTINY) ? 1.0 : sqrt( -log(rv[1]) /
					(cosp*cosp/(wp->rs.u_alpha*wp->rs.u_alpha) +
					 sinp*sinp/(wp->rs.v_alpha*wp->rs.v_alpha)) );
			for (i = 0; i < 3; i++)
				h[i] = wp->rs.mo.pnorm[i] +
					d*(cosp*wp->rs.u[i] + sinp*wp->rs.v[i]);
			d = -2.0 * DOT(h, wp->rp->rdir) / (1.0 + d*d);
			VSUM(sr.rdir, wp->rp->rdir, h, d);
						/* sample rejection test */
			d = DOT(sr.rdir, wp->rp->ron);
			if ((d > 0) ^ (wp->rp->rod > 0))
				continue;
			checknorm(sr.rdir);
			if (nstarget > 1) {	/* W-G-M-D adjustment */
				if (nstaken) rayclear(&sr);
				rayvalue(&sr);
				d = 2./(1. + wp->rp->rod/d);
				scalescolor(sr.rcol, d);
				saddscolor(scol, sr.rcol);
			} else {
				rayvalue(&sr);
				smultscolor(sr.rcol, sr.rcoef);
				saddscolor(wp->rp->rcol, sr.rcol);
			}
			++nstaken;
		}
		if (nstarget > 1) {		/* final W-G-M-D weighting */
			smultscolor(scol, sr.rcoef);
			d = (double)nstarget/ntrials;
			scalescolor(scol, d);
			saddscolor(wp->rp->rcol, scol);
		}
		ndims--;
	}
					/* compute transmission */
	if ((wp->specfl & (SP_TRAN|SP_TPURE|SP_TBLT)) == SP_TRAN &&
			rayorigin(&sr, TSPECULAR, wp->rp, wp->ts.scol) == 0) {
		nstarget = 1;
		if (specjitter > 1.5) {	/* multiple samples? */
			nstarget = specjitter*wp->rp->rweight + .5;
			if (sr.rweight <= minweight*nstarget)
				nstarget = sr.rweight/minweight;
			if (nstarget > 1) {
				d = 1./nstarget;
				scalescolor(sr.rcoef, d);
				sr.rweight *= d;
			} else
				nstarget = 1;
		}
		dimlist[ndims++] = (int)(size_t)wp->mtp;
		maxiter = MAXITER*nstarget;
		for (nstaken = ntrials = 0; (nstaken < nstarget) &
						(ntrials < maxiter); ntrials++) {
			if (ntrials)
				d = frandom();
			else
				d = urand(ilhash(dimlist,ndims)+1823+samplendx);
			multisamp(rv, 2, d);
			d = 2.0*PI * rv[0];
			cosp = tcos(d) * wp->ts.u_alpha;
			sinp = tsin(d) * wp->ts.v_alpha;
			d = 1./sqrt(cosp*cosp + sinp*sinp);
			cosp *= d;
			sinp *= d;
			if ((0. <= specjitter) & (specjitter < 1.))
				rv[1] = 1.0 - specjitter*rv[1];
			if (rv[1] <= FTINY)
				d = 1.0;
			else
				d = sqrt(-log(rv[1]) /
					(cosp*cosp/(wp->ts.u_alpha*wp->ts.u_alpha) +
					 sinp*sinp/(wp->ts.v_alpha*wp->ts.v_alpha)));
			for (i = 0; i < 3; i++)
				sr.rdir[i] = wp->prdir[i] +
						d*(cosp*wp->ts.u[i] + sinp*wp->ts.v[i]);
						/* rejection test */
			if ((DOT(sr.rdir,wp->rp->ron) > 0) == (wp->rp->rod > 0))
				continue;
			normalize(sr.rdir);	/* OK, normalize */
			if (nstaken)		/* multi-sampling? */
				rayclear(&sr);
			rayvalue(&sr);
			smultscolor(sr.rcol, sr.rcoef);
			saddscolor(wp->rp->rcol, sr.rcol);
			++nstaken;
		}
		ndims--;
	}
}

/* compute source contribution for MAT_WGMDF */
static void
dirwgmdf(SCOLOR scval, void *uwp, FVECT ldir, double omega)
{
	WGMDDAT		*wp = (WGMDDAT *)uwp;
	const int	hitfront = (wp->rp->rod > 0);
	double		fresadj = 1.;
	double		ldot;
	double		dtmp, dtmp1, dtmp2;
	FVECT		h;
	double		au2, av2;
	SCOLOR		sctmp;

	scolorblack(scval);		/* will add component coefficients */

					/* XXX ignores which side is lit */
	if (wp->specfl & SP_RPURE && pbright(wp->rs.scol) >= FRESTHRESH)
		fresadj = 1. - FRESNE(fabs(DOT(wp->rs.mo.pnorm,ldir)));

	if (sintens(wp->rd.scol) > FTINY &&
			((ldot = DOT(wp->rd.mo.pnorm,ldir)) > 0) == hitfront) {
		/*
		 *  Compute diffuse reflection coefficient for source.
		 */
		copyscolor(sctmp, wp->rd.scol);
		dtmp = fabs(ldot) * omega * (1.0/PI) * fresadj;
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
	if (sintens(wp->td.scol) > FTINY &&
			((ldot = DOT(wp->td.mo.pnorm,ldir)) > 0) ^ hitfront) {
		/*
		 *  Compute diffuse transmission coefficient for source.
		 */
		copyscolor(sctmp, wp->td.scol);
		dtmp = fabs(ldot) * omega * (1.0/PI) * fresadj;
		scalescolor(sctmp, dtmp);
		saddscolor(scval, sctmp);
	}
#if 0	/* XXX not yet implemented */
	if (ambRayInPmap(wp->rp))
		return;		/* specular accounted for in photon map */
#endif
	if ((wp->specfl & (SP_REFL|SP_RPURE)) == SP_REFL &&
			((ldot = DOT(wp->rs.mo.pnorm,ldir)) > 0) == hitfront) {
		/*
		 *  Compute specular reflection coefficient for source using
		 *  anisotropic Gaussian distribution model.
		 */
						/* add source width if flat */
		if (wp->specfl & SP_FLAT)
			au2 = av2 = (1. - dstrsrc) * omega * (0.25/PI);
		else
			au2 = av2 = 0.0;
		au2 += wp->rs.u_alpha*wp->rs.u_alpha;
		av2 += wp->rs.v_alpha*wp->rs.v_alpha;
						/* half vector */
		VSUB(h, ldir, wp->rp->rdir);
						/* ellipse */
		dtmp1 = DOT(wp->rs.u, h);
		dtmp1 *= dtmp1 / au2;
		dtmp2 = DOT(wp->rs.v, h);
		dtmp2 *= dtmp2 / av2;
						/* W-G-M-D model */
		dtmp = DOT(wp->rs.mo.pnorm, h);
		dtmp *= dtmp;
		dtmp1 = (dtmp1 + dtmp2) / dtmp;
		dtmp = exp(-dtmp1) * DOT(h,h) /
				(PI * dtmp*dtmp * sqrt(au2*av2));

		if (dtmp > FTINY) {		/* worth using? */
			copyscolor(sctmp, wp->rs.scol);
			dtmp *= fabs(ldot) * omega;
			scalescolor(sctmp, dtmp);
			saddscolor(scval, sctmp);
		}
	}
	if ((wp->specfl & (SP_TRAN|SP_TPURE)) == SP_TRAN &&
			((ldot = DOT(wp->ts.mo.pnorm,ldir)) > 0) ^ hitfront) {
		/*
		 *  Compute specular transmission coefficient for source.
		 */
						/* roughness + source */
		au2 = av2 = omega * (1.0/PI);
		au2 += wp->ts.u_alpha*wp->ts.u_alpha;
		av2 += wp->ts.v_alpha*wp->ts.v_alpha;
						/* "half vector" */
		VSUB(h, ldir, wp->prdir);
		dtmp = DOT(h,h);
		if (dtmp > FTINY*FTINY) {
			dtmp1 = DOT(h,wp->ts.mo.pnorm);
			dtmp = 1.0 - dtmp1*dtmp1/dtmp;
		}
		if (dtmp > FTINY*FTINY) {
			dtmp1 = DOT(h,wp->ts.u);
			dtmp1 *= dtmp1 / au2;
			dtmp2 = DOT(h,wp->ts.v);
			dtmp2 *= dtmp2 / av2;
			dtmp = (dtmp1 + dtmp2) / dtmp;
			dtmp = exp(-dtmp);
		} else
			dtmp = 1.0;
						/* Gaussian */
		dtmp *= (1.0/PI) * sqrt(-ldot/(wp->ts.mo.pdot*au2*av2));

		if (dtmp > FTINY) {		/* worth using? */
			copyscolor(sctmp, wp->ts.scol);
			dtmp *= omega;
			scalescolor(sctmp, dtmp);
			saddscolor(scval, sctmp);
		}
	}
}

/* color a ray that hit a programmable WGMD material */
int
m_wgmdf(OBJREC *m, RAY *r)
{
	RAY		lr;
	WGMDDAT		wd;
	SCOLOR		sctmp;
	FVECT		anorm;
	int		i;

	if (!backvis & (r->rod < 0.0)) {
		raytrans(r);
		return(1);		/* backside invisible */
	}
	if ((m->oargs.nsargs < 13) | (m->oargs.nfargs < 9))
		objerror(m, USER, "bad number of arguments");

	if (r->crtype & SHADOW && !strcmp(m->oargs.sarg[5], "0"))
		return(1);		/* first shadow test */
	clr_comps(&wd);
	wd.rp = r;
	wd.mtp = m;
	wd.mf = getfunc(m, 12, 0xEEE, 1);
	set_dcomp(&wd, 0);		/* gets main modifier */
	setWGMDfunc(&wd.rd.mo, &wd);	/* get local u vector */
	errno = 0;
	for (i = 0; i < 3; i++)
		wd.ulocal[i] = evalue(wd.mf->ep[6+i]);
	if ((errno == EDOM) | (errno == ERANGE))
		wd.ulocal[0] = wd.ulocal[1] = wd.ulocal[2] = 0.0;
	else if (wd.mf->fxp != &unitxf)
		multv3(wd.ulocal, wd.ulocal, wd.mf->fxp->xfm);

	set_scomp(&wd, 1);		/* sets SP_TPURE */
	if (r->crtype & SHADOW && !(wd.specfl & SP_TPURE))
		return(1);		/* second shadow test */
	set_dcomp(&wd, 1);
	set_scomp(&wd, 0);
	wd.specfl |= SP_FLAT*(!wd.rs.mo.hastexture &&
				r->ro != NULL && isflat(r->ro->otype));
					/* apply Fresnel adjustments? */
	if (wd.specfl & SP_RPURE && pbright(wd.rs.scol) >= FRESTHRESH) {
		const double	fest = FRESNE(fabs(wd.rs.mo.pdot));
		for (i = NCSAMP; i--; )
			wd.rs.scol[i] += fest*(1. - wd.rs.scol[i]);
		scalescolor(wd.rd.scol, 1.-fest);
		scalescolor(wd.ts.scol, 1.-fest);
		scalescolor(wd.td.scol, 1.-fest);
	}
					/* check specular thresholds */
	wd.specfl |= SP_RBLT*((wd.specfl & (SP_REFL|SP_RPURE)) == SP_REFL &&
				specthresh >= pbright(wd.rs.scol)-FTINY);
	wd.specfl |= SP_TBLT*((wd.specfl & (SP_TRAN|SP_TPURE)) == SP_TRAN &&
				specthresh >= pbright(wd.ts.scol)-FTINY);
					/* get through direction */
	if (wd.specfl & SP_TRAN && wd.ts.mo.hastexture &&
			!(r->crtype & (SHADOW|AMBIENT))) {
		for (i = 0; i < 3; i++)	/* perturb */
			wd.prdir[i] = r->rdir[i] - wd.ts.mo.pnorm[i] + r->ron[i];
		if ((DOT(wd.prdir,r->ron) > 0) ^ (r->rod > 0))
			normalize(wd.prdir);	/* OK */
		else			/* too much */
			VCOPY(wd.prdir, r->rdir);
	} else
		VCOPY(wd.prdir, r->rdir);
					/* transmitted view ray? */
	if ((wd.specfl & (SP_TRAN|SP_TPURE|SP_TBLT)) == (SP_TRAN|SP_TPURE) &&
			rayorigin(&lr, TRANS, r, wd.ts.scol) == 0) {
		VCOPY(lr.rdir, wd.prdir);
		rayvalue(&lr);
		smultscolor(lr.rcol, lr.rcoef);
		saddscolor(r->rcol, lr.rcol);
		if (scolor_mean(wd.ts.scol) >= 0.999) {
					/* completely transparent */
			smultscolor(lr.mcol, lr.rcoef);
			copyscolor(r->mcol, lr.mcol);
			r->rmt = r->rot + lr.rmt;
			r->rxt = r->rot + lr.rxt;
		} else if (pbright(wd.ts.scol) >
				pbright(wd.td.scol) + pbright(wd.rd.scol))
			r->rxt = r->rot + raydistance(&lr);
	}
	if (r->crtype & SHADOW)
		return(1);		/* the rest is shadow */
					/* mirror ray? */
	if ((wd.specfl & (SP_REFL|SP_RPURE|SP_RBLT)) == (SP_REFL|SP_RPURE) &&
			rayorigin(&lr, REFLECTED, r, wd.rs.scol) == 0) {
		VSUM(lr.rdir, r->rdir, wd.rs.mo.pnorm, 2.*wd.rs.mo.pdot);
					/* fall back if would penetrate */
		if (wd.rs.mo.hastexture &&
				(DOT(lr.rdir,r->ron) > 0) ^ (r->rod > 0))
			VSUM(lr.rdir, r->rdir, r->ron, 2.*r->rod);
		checknorm(lr.rdir);
		rayvalue(&lr);
		smultscolor(lr.rcol, lr.rcoef);
		copyscolor(r->mcol, lr.rcol);
		saddscolor(r->rcol, lr.rcol);
		r->rmt = r->rot;
		if (wd.specfl & SP_FLAT && r->crtype & AMBIENT)
			r->rmt += raydistance(&lr);
	}
	if (wd.specfl & (SP_REFL|SP_TRAN))	/* specularly scattered rays */
		agaussamp(&wd);		/* checks *BLT flags */

	if (sintens(wd.rd.scol) > FTINY) {	/* ambient from this side */
		if (r->rod > 0) {
			VCOPY(anorm, wd.rd.mo.pnorm);
		} else {
			anorm[0] = -wd.rd.mo.pnorm[0];
			anorm[1] = -wd.rd.mo.pnorm[1];
			anorm[2] = -wd.rd.mo.pnorm[2];
		}
		copyscolor(sctmp, wd.rd.scol);
		if (wd.specfl & SP_RBLT)	/* add in specular as well? */
			saddscolor(sctmp, wd.rs.scol);
		multambient(sctmp, r, anorm);
		saddscolor(r->rcol, sctmp);	/* add to returned color */
	}
	if (sintens(wd.td.scol) > FTINY) {	/* ambient from other side */
		if (r->rod > 0) {
			anorm[0] = -wd.td.mo.pnorm[0];
			anorm[1] = -wd.td.mo.pnorm[1];
			anorm[2] = -wd.td.mo.pnorm[2];
		} else {
			VCOPY(anorm, wd.td.mo.pnorm);
		}
		copyscolor(sctmp, wd.td.scol);
		if (wd.specfl & SP_TBLT)	/* add in specular as well? */
			saddscolor(sctmp, wd.ts.scol)
		multambient(sctmp, r, anorm);
		saddscolor(r->rcol, sctmp);
	}
	direct(r, dirwgmdf, &wd);	/* add direct component last */
	return(1);
}
