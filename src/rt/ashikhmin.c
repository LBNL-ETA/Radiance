#ifndef lint
static const char RCSid[] = "$Id: ashikhmin.c,v 2.6 2015/09/02 18:59:01 greg Exp $";
#endif
/*
 *  Shading functions for Ashikhmin-Shirley anisotropic materials.
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

/*
 *  Ashikhmin-Shirley model
 *
 *  Arguments for MAT_ASHIKHMIN are:
 *  4+ ux uy uz	funcfile [transform...]
 *  0
 *  8  dred dgrn dblu sred sgrn sblu u-power v-power
 */

				/* specularity flags */
#define  SPA_REFL	01		/* has reflected specular component */
#define  SPA_FLAT	02		/* reflecting surface is flat */
#define  SPA_RBLT	010		/* reflection below sample threshold */

typedef struct {
	OBJREC  *mp;		/* material pointer */
	RAY  *rp;		/* ray pointer */
	short  specfl;		/* specularity flags, defined above */
	COLOR  mcolor;		/* color of this material */
	COLOR  scolor;		/* color of specular component */
	FVECT  u, v;		/* u and v vectors orienting anisotropy */
	double  u_power;	/* u power */
	double  v_power;	/* v power */
	FVECT  pnorm;		/* perturbed surface normal */
	double  pdot;		/* perturbed dot product */
}  ASHIKDAT;		/* anisotropic material data */

static void getacoords_as(ASHIKDAT *np);
static void ashiksamp(ASHIKDAT *np);

#undef MAX
#define MAX(a,b)	((a)>(b) ? (a) : (b))


static double
schlick_fres(double dprod)
{
	double	pf = 1. - dprod;

	return(pf*pf*pf*pf*pf);
}


static void
dirashik(		/* compute source contribution */
	COLOR  cval,		/* returned coefficient */
	void  *nnp,		/* material data */
	FVECT  ldir,		/* light source direction */
	double  omega		/* light source size */
)
{
	ASHIKDAT *np = nnp;
	double  ldot;
	double  dtmp, dtmp1, dtmp2;
	FVECT  h;
	COLOR  ctmp;

	setcolor(cval, 0.0, 0.0, 0.0);

	ldot = DOT(np->pnorm, ldir);

	if (ldot < 0.0)
		return;		/* wrong side */

	/*
	 *  Compute and add diffuse reflected component to returned
	 *  color.
	 */
	copycolor(ctmp, np->mcolor);
	dtmp = ldot * omega * (1.0/PI) * (1. - schlick_fres(ldot));
	scalecolor(ctmp, dtmp);		
	addcolor(cval, ctmp);

	if (!(np->specfl & SPA_REFL) || ambRayInPmap(np->rp))
		return;
	/*
	 *  Compute specular reflection coefficient
	 */
					/* half vector */
	VSUB(h, ldir, np->rp->rdir);
	normalize(h);
					/* ellipse */
	dtmp1 = DOT(np->u, h);
	dtmp1 *= dtmp1 * np->u_power;
	dtmp2 = DOT(np->v, h);
	dtmp2 *= dtmp2 * np->v_power;
					/* Ashikhmin-Shirley model*/
	dtmp = DOT(np->pnorm, h);
	dtmp = pow(dtmp, (dtmp1+dtmp2)/(1.-dtmp*dtmp));
	dtmp *= sqrt((np->u_power+1.)*(np->v_power+1.));
	dtmp /= 8.*PI * DOT(ldir,h) * MAX(ldot,np->pdot);
					/* worth using? */
	if (dtmp > FTINY) {
		copycolor(ctmp, np->scolor);
		dtmp *= ldot * omega;
		scalecolor(ctmp, dtmp);
		addcolor(cval, ctmp);
	}
}


int
m_ashikhmin(			/* shade ray that hit something anisotropic */
	OBJREC  *m,
	RAY  *r
)
{
	ASHIKDAT  nd;
	COLOR  ctmp;
	double	fres;
	int  i;
						/* easy shadow test */
	if (r->crtype & SHADOW)
		return(1);

	if (m->oargs.nfargs != 8)
		objerror(m, USER, "bad number of real arguments");
						/* check for back side */
	if (r->rod < 0.0) {
		if (!backvis) {
			raytrans(r);
			return(1);
		}
		raytexture(r, m->omod);
		flipsurface(r);			/* reorient if backvis */
	} else
		raytexture(r, m->omod);
						/* get material color */
	nd.mp = m;
	nd.rp = r;
	setcolor(nd.mcolor, m->oargs.farg[0],
			   m->oargs.farg[1],
			   m->oargs.farg[2]);
	setcolor(nd.scolor, m->oargs.farg[3],
			   m->oargs.farg[4],
			   m->oargs.farg[5]);
						/* get specular power */
	nd.specfl = 0;
	nd.u_power = m->oargs.farg[6];
	nd.v_power = m->oargs.farg[7];

	nd.pdot = raynormal(nd.pnorm, r);	/* perturb normal */
	if (nd.pdot < .001)
		nd.pdot = .001;			/* non-zero for dirashik() */
	multcolor(nd.mcolor, r->pcol);		/* modify diffuse color */

	if (bright(nd.scolor) > FTINY) {	/* adjust specular color */
		nd.specfl |= SPA_REFL;
						/* check threshold */
		if (specthresh >= bright(nd.scolor)-FTINY)
			nd.specfl |= SPA_RBLT;
		fres = schlick_fres(nd.pdot);	/* Schick's Fresnel approx */
		for (i = 0; i < 3; i++)
			colval(nd.scolor,i) += (1.-colval(nd.scolor,i))*fres;
	}
	if (r->ro != NULL && isflat(r->ro->otype))
		nd.specfl |= SPA_FLAT;
						/* set up coordinates */
	getacoords_as(&nd);
						/* specular sampling? */
	if ((nd.specfl & (SPA_REFL|SPA_RBLT)) == SPA_REFL)
		ashiksamp(&nd);
						/* diffuse interreflection */
	if (bright(nd.mcolor) > FTINY) {
		copycolor(ctmp, nd.mcolor);	/* modified by material color */		
		if (nd.specfl & SPA_RBLT)	/* add in specular as well? */
			addcolor(ctmp, nd.scolor);
		multambient(ctmp, r, nd.pnorm);
		addcolor(r->rcol, ctmp);	/* add to returned color */
	}
	direct(r, dirashik, &nd);		/* add direct component */

	return(1);
}


static void
getacoords_as(		/* set up coordinate system */
	ASHIKDAT  *np
)
{
	MFUNC  *mf;
	int  i;

	mf = getfunc(np->mp, 3, 0x7, 1);
	setfunc(np->mp, np->rp);
	errno = 0;
	for (i = 0; i < 3; i++)
		np->u[i] = evalue(mf->ep[i]);
	if ((errno == EDOM) | (errno == ERANGE))
		np->u[0] = np->u[1] = np->u[2] = 0.0;
	if (mf->fxp != &unitxf)
		multv3(np->u, np->u, mf->fxp->xfm);
	fcross(np->v, np->pnorm, np->u);
	if (normalize(np->v) == 0.0) {
		if (fabs(np->u_power - np->v_power) > 0.1)
			objerror(np->mp, WARNING, "bad orientation vector");
		getperpendicular(np->u, np->pnorm, 1);	/* punting */
		fcross(np->v, np->pnorm, np->u);
		np->u_power = np->v_power =
			2./(1./(np->u_power+1e-5) + 1./(np->v_power+1e-5));
	} else
		fcross(np->u, np->v, np->pnorm);
}


static void
ashiksamp(		/* sample anisotropic Ashikhmin-Shirley specular */
	ASHIKDAT  *np
)
{
	RAY  sr;
	FVECT  h;
	double  rv[2], dtmp;
	double  cosph, sinph, costh, sinth;
	int  maxiter, ntrials, nstarget, nstaken;
	int  i;

	if (rayorigin(&sr, SPECULAR, np->rp, np->scolor) < 0)
		return;

	nstarget = 1;
	if (specjitter > 1.5) {			/* multiple samples? */
		nstarget = specjitter*np->rp->rweight + .5;
		if (sr.rweight <= minweight*nstarget)
			nstarget = sr.rweight/minweight;
		if (nstarget > 1) {
			dtmp = 1./nstarget;
			scalecolor(sr.rcoef, dtmp);
			sr.rweight *= dtmp;
		} else
			nstarget = 1;
	}
	dimlist[ndims++] = (int)(size_t)np->mp;
	maxiter = MAXITER*nstarget;
	for (nstaken = ntrials = 0; nstaken < nstarget &&
					ntrials < maxiter; ntrials++) {
		if (ntrials)
			dtmp = frandom();
		else
			dtmp = urand(ilhash(dimlist,ndims)+647+samplendx);
		multisamp(rv, 2, dtmp);
		dtmp = 2.*PI * rv[0];
		cosph = sqrt(np->v_power + 1.) * tcos(dtmp);
		sinph = sqrt(np->u_power + 1.) * tsin(dtmp);
		dtmp = 1./sqrt(cosph*cosph + sinph*sinph);
		cosph *= dtmp;
		sinph *= dtmp;
		costh = pow(rv[1], 1./(np->u_power*cosph*cosph+np->v_power*sinph*sinph+1.));
		if (costh <= FTINY)
			continue;
		sinth = sqrt(1. - costh*costh);
		for (i = 0; i < 3; i++)
			h[i] = cosph*sinth*np->u[i] + sinph*sinth*np->v[i] + costh*np->pnorm[i];

		if (nstaken)
			rayclear(&sr);
		dtmp = -2.*DOT(h, np->rp->rdir);
		VSUM(sr.rdir, np->rp->rdir, h, dtmp);				
						/* sample rejection test */
		if (DOT(sr.rdir, np->rp->ron) <= FTINY)
			continue;
		checknorm(sr.rdir);
		rayvalue(&sr);
		multcolor(sr.rcol, sr.rcoef);
		addcolor(np->rp->rcol, sr.rcol);
		++nstaken;
	}
	ndims--;
}
