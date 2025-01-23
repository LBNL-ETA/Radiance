static const char	RCSid[] = "$Id: ambient.c,v 2.129 2025/01/23 02:15:33 greg Exp $";
/*
 *  ambient.c - routines dealing with ambient (inter-reflected) component.
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include <string.h>

#include  "platform.h"
#include  "ray.h"
#include  "otypes.h"
#include  "otspecial.h"
#include  "resolu.h"
#include  "ambient.h"
#include  "random.h"
#include  "pmapamb.h"

#ifndef  OCTSCALE
#define	 OCTSCALE	1.0	/* ceil((valid rad.)/(cube size)) */
#endif

#ifndef  MAXASET
#define	 MAXASET	4095	/* maximum number of elements in ambient set */
#endif
OBJECT	ambset[MAXASET+1]={0};	/* ambient include/exclude set */

double	maxarad;		/* maximum ambient radius */
double	minarad;		/* minimum ambient radius */

static AMBTREE	atrunk;		/* our ambient trunk node */

static FILE  *ambfp = NULL;	/* ambient file pointer */
static int  nunflshed;		/* number of unflushed ambient values */
static FILE  *ambinp = NULL;	/* input pointer for ambient i/o */

static double  avsum = 0.;		/* computed ambient value sum (log) */
static unsigned int  navsum = 0;	/* number of values in avsum */
static unsigned int  nambvals = 0;	/* total number of indirect values */
static off_t  lastpos = -1;		/* last flush position */

#define	 AMBFLUSH	(BUFSIZ/AMBVALSIZ)

#define  AVSIZE		(sizeof(AMBVAL)-sizeof(SCOLOR)+sizeof(COLORV)*NCSAMP)
#define	 newambval()	(AMBVAL *)malloc(AVSIZE)

#define  tfunc(x0, x, x1)	(((x)-(x0))/((x1)-(x0)))

static void initambfile(int cre8);
static void avsave(AMBVAL *av);
static AMBVAL *avstore(AMBVAL  *aval);
static AMBTREE *newambtree(void);
static void freeambtree(AMBTREE  *atp);

typedef void unloadtf_t(AMBVAL *);
static unloadtf_t avinsert;
static unloadtf_t avfree;
static void unloadatree(AMBTREE  *at, unloadtf_t *f);

static void sortambvals(void);

static int	plugaleak(RAY *r, AMBVAL *ap, FVECT anorm, double ang);
static double	sumambient(SCOLOR acol, RAY *r, FVECT rn, int al,
				AMBTREE *at, FVECT c0, double s);
static int	makeambient(SCOLOR acol, RAY *r, FVECT rn, int al);
static int	extambient(SCOLOR cr, AMBVAL *ap, FVECT pv, FVECT nv,
				FVECT uvw[3]);


void
setambres(				/* set ambient resolution */
	int  ar
)
{
	ambres = ar < 0 ? 0 : ar;		/* may be done already */
						/* set min & max radii */
	if (ar <= 0) {
		minarad = 0;
		maxarad = thescene.cusize*0.2;
	} else {
		minarad = thescene.cusize / ar;
		maxarad = 64.0 * minarad;		/* heuristic */
		if (maxarad > thescene.cusize*0.2)
			maxarad = thescene.cusize*0.2;
	}
	if (minarad <= FTINY)
		minarad = 10.0*FTINY;
	if (maxarad <= minarad)
		maxarad = 64.0 * minarad;
}


void
setambacc(				/* set ambient accuracy */
	double  newa
)
{
	static double	olda;		/* remember previous setting here */
	
	newa *= (newa > 0);
	if (fabs(newa - olda) >= .05*(newa + olda)) {
		ambacc = newa;
		if (ambacc > FTINY && nambvals > 0)
			sortambvals();		/* rebuild tree */
	}
}


void
setambient(void)				/* initialize calculation */
{
	int	exists;
	off_t	flen;
	AMBVAL	amb;
						/* make sure we're fresh */
	ambdone();
						/* init ambient limits */
	setambres(ambres);
	setambacc(ambacc);
	if (ambfile == NULL || !ambfile[0])
		return;
	if (ambacc <= FTINY) {
		sprintf(errmsg, "zero ambient accuracy so \"%s\" not opened",
				ambfile);
		error(WARNING, errmsg);
		return;
	}
	exists = access(ambfile, F_OK) == 0;	/* check existence, first */
	ambfp = fopen(ambfile, "a+");		/* try read/append */
	if (!exists & (ambfp == NULL)) {
		sprintf(errmsg, "cannot create ambient file \"%s\"", ambfile);
		error(SYSTEM, errmsg);
	}
	if (ambfp == NULL) {			/* try opening read-only? */
		if ((ambfp = fopen(ambfile, "r")) == NULL) {
			sprintf(errmsg,
			"cannot open ambient file \"%s\" for reading",
					ambfile);
			error(SYSTEM, errmsg);
		}
		exists = -1;			/* flag read-only */
	} else if (exists)
		rewind(ambfp);	/* XXX not necessary? */

	if (exists) {
		initambfile(0);			/* file already exists */
		lastpos = ftell(ambfp);
		while (readambval(&amb, ambfp))
			avstore(&amb);		/* load what we can */
		if (exists < 0) {		/* read-only? */
			sprintf(errmsg,
				"loaded %u values from read-only ambient file",
					nambvals);
			error(WARNING, errmsg);
			fclose(ambfp);		/* close file so no writes */
			ambfp = NULL;
			return;
		}
						/* align file pointer */
		lastpos += (off_t)nambvals*AMBVALSIZ;
		flen = lseek(fileno(ambfp), 0, SEEK_END);
		if (flen != lastpos) {
			sprintf(errmsg,
			"ignoring last %lu values in ambient file (corrupted)",
				(unsigned long)((flen - lastpos)/AMBVALSIZ));
			error(WARNING, errmsg);
						/* fseek() not needed? */
			fseeko(ambfp, lastpos, SEEK_SET);
			ftruncate(fileno(ambfp), lastpos);
		}
	} else {
		initambfile(1);			/* else start new file */
		fflush(ambfp);
		lastpos = ftell(ambfp);
	}
}


void
ambdone(void)			/* close ambient file and free memory */
{
	if (ambfp != NULL) {		/* close ambient file */
		fclose(ambfp);		/* don't call ambsync() */
		ambfp = NULL;
		lastpos = -1;
		if (ambinp != NULL) {
			fclose(ambinp);
			ambinp = NULL;
		}
	}
					/* free ambient tree */
	unloadatree(&atrunk, avfree);
					/* reset state variables */
	avsum = 0.;
	navsum = 0;
	nambvals = 0;
}


void
ambnotify(			/* record new modifier */
	OBJECT	obj
)
{
	static int  hitlimit = 0;
	OBJREC	 *o;
	char  **amblp;

	if (obj == OVOID) {		/* starting over */
		ambset[0] = 0;
		hitlimit = 0;
		return;
	}
	o = objptr(obj);
	if (hitlimit || !ismodifier(o->otype))
		return;
	for (amblp = amblist; *amblp != NULL; amblp++)
		if (!strcmp(o->oname, *amblp)) {
			if (ambset[0] >= MAXASET) {
				error(WARNING, "too many modifiers in ambient list");
				hitlimit++;
				return;		/* should this be fatal? */
			}
			insertelem(ambset, obj);
			return;
		}
}


void
multambient(		/* compute ambient component & multiply by coef. */
	SCOLOR  aval,
	RAY  *r,
	FVECT  nrm
)
{
	static double  logAvgAbsorp = 1;
	static int  rdepth = 0;			/* ambient recursion */
	SCOLOR	acol, caustic;
	int	i, ok;
	double	d, l;

	/* PMAP: Factor in ambient from photon map, if enabled and ray is
	 * ambient. Return as all ambient components accounted for, else
	 * continue. */
	if (ambPmap(aval, r, rdepth))
		return;

	if (logAvgAbsorp > 0)			/* exclude in -aw to avoid growth */
		logAvgAbsorp = log(1.-AVGREFL);

	/* PMAP: Factor in specular-diffuse ambient (caustics) from photon
	 * map, if enabled and ray is primary, else caustic is zero.  Continue
	 * with RADIANCE ambient calculation */
{/* XXX TEMPORARY */
	COLOR	pmc;
	scolor_color(pmc, aval);
	ambPmapCaustic(pmc, r, rdepth);
	setscolor(caustic, colval(pmc,RED), colval(pmc,GRN), colval(pmc,BLU));
}
	if (ambdiv <= 0)			/* no ambient calculation */
		goto dumbamb;
						/* check number of bounces */
	if (rdepth >= ambounce)
		goto dumbamb;
						/* check ambient list */
	if (ambincl != -1 && r->ro != NULL &&
			ambincl != inset(ambset, r->ro->omod))
		goto dumbamb;

	if (ambacc <= FTINY) {			/* no ambient storage? */
		double	rdot = DOT(nrm,r->ron);
		int	sgn = 1 - 2*(rdot < 0);
		float	dgrad[2], *dgp = NULL;
		FVECT	uvd[2];

		if (sgn*rdot < 0.9999)
			dgp = dgrad;		/* compute rotational grad. */
		copyscolor(acol, aval);
		rdepth++;
		ok = doambient(acol, r, r->rweight*sgn,
				uvd, NULL, NULL, dgp, NULL);
		rdepth--;
		if (!ok)
			goto dumbamb;
		if ((ok > 0) & (dgp != NULL)) {	/* apply texture */
			FVECT	v1;
			VCROSS(v1, r->ron, nrm);
			d = 1.0;
			for (i = 3; i--; )
				d += sgn*v1[i] * (dgp[0]*uvd[0][i] + dgp[1]*uvd[1][i]);
			if (d >= 0.05)
				scalescolor(acol, d);
		}
		copyscolor(aval, acol);

		/* PMAP: add in caustic */
		saddscolor(aval, caustic);
		return;
	}
						/* interpolate ambient value */
	scolorblack(acol);
	d = sumambient(acol, r, nrm, rdepth,
			&atrunk, thescene.cuorg, thescene.cusize);
			
	if (d > FTINY) {
		scalescolor(acol, 1.0/d);
		smultscolor(aval, acol);

		/* PMAP: add in caustic */
		saddscolor(aval, caustic);
		return;
	}
	
	rdepth++;				/* need to cache new value */
	ok = makeambient(acol, r, nrm, rdepth-1);
	rdepth--;
	
	if (ok) {
		smultscolor(aval, acol);	/* computed new value */

		/* PMAP: add in caustic */
		saddscolor(aval, caustic);
		return;
	}
	
dumbamb:					/* return global value */
	if ((ambvwt <= 0) | (navsum == 0)) {
		smultcolor(aval, ambval);
		
		/* PMAP: add in caustic */
		saddscolor(aval, caustic);
		return;
	}
	
	l = bright(ambval);			/* average in computations */	
	if (l > FTINY) {
		d = (log(l)*(double)ambvwt + avsum + logAvgAbsorp*navsum) /
				(double)(ambvwt + navsum);
		d = exp(d) / l;
		scalescolor(aval, d);
		smultcolor(aval, ambval);	/* apply color of ambval */
	} else {
		d = exp( avsum/(double)navsum + logAvgAbsorp );
		scalescolor(aval, d);		/* neutral color */
	}
}


/* Plug a potential leak where ambient cache value is occluded */
static int
plugaleak(RAY *r, AMBVAL *ap, FVECT anorm, double ang)
{
	const double	cost70sq = 0.1169778;	/* cos(70deg)^2 */
	RAY		rtst;
	FVECT		vdif;
	double		normdot, ndotd, nadotd;
	double		a, b, c, t[2];

	ang += 2.*PI*(ang < 0);			/* check direction flags */
	if ( !(ap->corral>>(int)(ang*(16./PI)) & 1) )
		return(0);
	/*
	 * Generate test ray, targeting 20 degrees above sample point plane
	 * along surface normal from cache position.  This should be high
	 * enough to miss local geometry we don't really care about.
	 */
	VSUB(vdif, ap->pos, r->rop);
	normdot = DOT(anorm, r->ron);
	ndotd = DOT(vdif, r->ron);
	nadotd = DOT(vdif, anorm);
	a = normdot*normdot - cost70sq;
	b = 2.0*(normdot*ndotd - nadotd*cost70sq);
	c = ndotd*ndotd - DOT(vdif,vdif)*cost70sq;
	if (quadratic(t, a, b, c) != 2)
		return(1);			/* should rarely happen */
	if (t[1] <= FTINY)
		return(0);			/* should fail behind test */
	rayorigin(&rtst, SHADOW, r, NULL);
	VSUM(rtst.rdir, vdif, anorm, t[1]);	/* further dist. > plane */
	rtst.rmax = normalize(rtst.rdir);	/* short ray test */
	while (localhit(&rtst, &thescene)) {	/* check for occluder */
		OBJREC	*m = findmaterial(rtst.ro);
		if (m != NULL && !istransp(m) && !isBSDFproxy(m) &&
				(rtst.clipset == NULL ||
					!inset(rtst.clipset, rtst.ro->omod)))
			return(1);		/* plug light leak */
		VCOPY(rtst.rorg, rtst.rop);	/* skip invisible surface */
		rtst.rmax -= rtst.rot;
		rayclear(&rtst);
	}
	return(0);				/* seems we're OK */
}


static double
sumambient(		/* get interpolated ambient value */
	SCOLOR  acol,
	RAY  *r,
	FVECT  rn,
	int  al,
	AMBTREE	 *at,
	FVECT  c0,
	double	s
)
{			/* initial limit is 10 degrees plus ambacc radians */
	const double	minangle = 10.0 * PI/180.;
	const int	sgn = 1 - 2*(DOT(r->ron,rn) < 0);
	double		maxangle = minangle + ambacc;
	double		wsum = 0.0;
	FVECT		ck0;
	int		i, j;
	AMBVAL		*av;

	if (at->kid != NULL) {		/* sum children first */				
		s *= 0.5;
		for (i = 0; i < 8; i++) {
			for (j = 0; j < 3; j++) {
				ck0[j] = c0[j];
				if (1<<j & i)
					ck0[j] += s;
				if (r->rop[j] < ck0[j] - OCTSCALE*s)
					break;
				if (r->rop[j] > ck0[j] + (1.0+OCTSCALE)*s)
					break;
			}
			if (j == 3)
				wsum += sumambient(acol, r, rn, al,
							at->kid+i, ck0, s);
		}
					/* good enough? */
		if ((wsum >= 0.05) & (s*ambacc > minarad))
			return(wsum);
	}
					/* adjust maximum angle */
	if (at->alist != NULL && (at->alist->lvl <= al) & (r->rweight < 0.6))
		maxangle = (maxangle - PI/2.)*pow(r->rweight,0.13) + PI/2.;
					/* sum this node */
	for (av = at->alist; av != NULL; av = av->next) {
		double	u, v, d, delta_r2, delta_t2;
		SCOLOR	sct;
		FVECT	uvw[3];
		/*
		 *  Ambient level test
		 */
		if (av->lvl > al ||	/* list sorted, so this works */
				(av->lvl == al) & (av->weight < 0.9*r->rweight))
			break;
		/*
		 *  Direction test using unperturbed normal
		 */
		decodedir(uvw[2], av->ndir);
		d = sgn * DOT(uvw[2], r->ron);
		if (d <= 0.0)		/* >= 90 degrees */
			continue;
		delta_r2 = 2.0 - 2.0*d;	/* approx. radians^2 */
		if (delta_r2 >= maxangle*maxangle)
			continue;
		/*
		 *  Modified ray behind test
		 */
		VSUB(ck0, r->rop, av->pos);
		d = DOT(ck0, uvw[2]);
		if (d < -minarad*ambacc)
			continue;
		d /= av->rad[0];
		delta_t2 = d*d;
		if (delta_t2 >= ambacc*ambacc)
			continue;
		/*
		 *  Elliptical radii test based on Hessian
		 */
		decodedir(uvw[0], av->udir);
		VCROSS(uvw[1], uvw[2], uvw[0]);
		d = (u = DOT(ck0, uvw[0])) / av->rad[0];
		delta_t2 += d*d;
		d = (v = DOT(ck0, uvw[1])) / av->rad[1];
		delta_t2 += d*d;
		if (delta_t2 >= ambacc*ambacc)
			continue;
		/*
		 *  Test for potential light leak
		 */
		if (av->corral && plugaleak(r, av, uvw[2], atan2a(v,u)))
			continue;
		/*
		 *  Extrapolate value and compute final weight (hat function)
		 */
		if (!extambient(sct, av, r->rop, rn, uvw))
			continue;
		d = tfunc(maxangle, sqrt(delta_r2), 0.0) *
			tfunc(ambacc, sqrt(delta_t2), 0.0);
		scalescolor(sct, d);
		saddscolor(acol, sct);
		wsum += d;
	}
	return(wsum);
}


static int
makeambient(		/* make a new ambient value for storage */
	SCOLOR  acol,
	RAY  *r,
	FVECT  rn,
	int  al
)
{
	int	sgn = 1 - 2*(DOT(r->ron,rn) < 0);
	AMBVAL	amb;
	FVECT	uvw[3];
	int	i;

	amb.weight = 1.0;			/* compute weight */
	for (i = al; i-- > 0; )
		amb.weight *= AVGREFL;
	if (r->rweight < 0.1*amb.weight)	/* heuristic override */
		amb.weight = 1.25*r->rweight;
	setscolor(acol, AVGREFL, AVGREFL, AVGREFL);
						/* compute ambient */
	i = doambient(acol, r, amb.weight*sgn,
			uvw, amb.rad, amb.gpos, amb.gdir, &amb.corral);
	scalescolor(acol, 1./AVGREFL);		/* undo assumed reflectance */
	if (i <= 0 || amb.rad[0] <= FTINY)	/* no Hessian or zero radius */
		return(i);
	uvw[2][0] = sgn*r->ron[0];		/* orient unperturbed normal */
	uvw[2][1] = sgn*r->ron[1];
	uvw[2][2] = sgn*r->ron[2];
						/* store value */
	VCOPY(amb.pos, r->rop);
	amb.ndir = encodedir(uvw[2]);
	amb.udir = encodedir(uvw[0]);
	amb.lvl = al;
	copyscolor(amb.val, acol);
	avsave(&amb);				/* insert and save to file */
	if (DOT(uvw[2],rn) < 0.9999)		/* texture? */
		extambient(acol, &amb, r->rop, rn, uvw);
	return(1);
}


static int
extambient(		/* extrapolate value at pv, nv */
	SCOLOR  scr,
	AMBVAL	 *ap,
	FVECT  pv,
	FVECT  nv,
	FVECT  uvw[3]
)
{
	const double	min_d = 0.05;
	const double	max_d = 20.;
	static FVECT	my_uvw[3];
	FVECT		v1;
	int		i;
	double		d = 1.0;	/* zeroeth order */

	if (uvw == NULL) {		/* need local coordinates? */
		decodedir(my_uvw[2], ap->ndir);
		decodedir(my_uvw[0], ap->udir);
		VCROSS(my_uvw[1], my_uvw[2], my_uvw[0]);
		uvw = my_uvw;
	}
	for (i = 3; i--; )		/* gradient due to translation */
		d += (pv[i] - ap->pos[i]) *
			(ap->gpos[0]*uvw[0][i] + ap->gpos[1]*uvw[1][i]);

	VCROSS(v1, uvw[2], nv);		/* gradient due to rotation */
	for (i = 3; i--; )
		d += v1[i] * (ap->gdir[0]*uvw[0][i] + ap->gdir[1]*uvw[1][i]);
	
	if (d < min_d)			/* clamp min/max scaling */
		d = min_d;
	else if (d > max_d)
		d = max_d;
	copyscolor(scr, ap->val);
	scalescolor(scr, d);
	return(d > min_d);
}


static void
avinsert(				/* insert ambient value in our tree */
	AMBVAL *av
)
{
	AMBTREE  *at;
	AMBVAL  *ap;
	AMBVAL  avh;
	FVECT  ck0;
	double	s;
	int  branch;
	int  i;

	if (av->rad[0] <= FTINY)
		error(CONSISTENCY, "zero ambient radius in avinsert");
	at = &atrunk;
	VCOPY(ck0, thescene.cuorg);
	s = thescene.cusize;
	while (s*(OCTSCALE/2) > av->rad[1]*ambacc) {
		if (at->kid == NULL)
			if ((at->kid = newambtree()) == NULL)
				error(SYSTEM, "out of memory in avinsert");
		s *= 0.5;
		branch = 0;
		for (i = 0; i < 3; i++)
			if (av->pos[i] > ck0[i] + s) {
				ck0[i] += s;
				branch |= 1 << i;
			}
		at = at->kid + branch;
	}
	avh.next = at->alist;		/* order by increasing level */
	for (ap = &avh; ap->next != NULL; ap = ap->next)
		if ( ap->next->lvl > av->lvl ||
				(ap->next->lvl == av->lvl) &
				(ap->next->weight <= av->weight) )
			break;
	av->next = ap->next;
	ap->next = (AMBVAL*)av;
	at->alist = avh.next;
}


static void
initambfile(		/* initialize ambient file */
	int  cre8
)
{
	extern char  *progname, *octname;
	static char  *mybuf = NULL;
	int  ntries = 3;

	if (!AMBFLUSH)
		error(INTERNAL, "BUFSIZ too small in initambfile");
	SET_FILE_BINARY(ambfp);
	if (mybuf == NULL)
		mybuf = (char *)bmalloc(BUFSIZ);
	setbuf(ambfp, mybuf);
	nunflshed = 0;
retry:
	if (cre8) {			/* new file */
		newheader("RADIANCE", ambfp);
		fprintf(ambfp, "%s -av %g %g %g -aw %d -ab %d -aa %g ",
				progname, colval(ambval,RED),
				colval(ambval,GRN), colval(ambval,BLU),
				ambvwt, ambounce, ambacc);
		fprintf(ambfp, "-ad %d -as %d -ar %d ",
				ambdiv, ambssamp, ambres);
		fprintf(ambfp, "-dr %d -ds %g -dt %g -dc %g ", directrelay,
				srcsizerat, shadthresh, shadcert);
		fprintf(ambfp, "-ss %g -st %g -lr %d -lw %g ", specjitter,
				specthresh, maxdepth, minweight);
		fprintf(ambfp, "-cw %g %g -cs %d ", WLPART[3], WLPART[0], NCSAMP);
		if (octname != NULL)
			fputs(octname, ambfp);
		fputc('\n', ambfp);	/* end of command line, not header! */
		fprintf(ambfp, "SOFTWARE= %s\n", VersionID);
		fputnow(ambfp);
		AMB_CNDX = CNDX;	/* use current spectral sampling */
		AMB_WLPART = WLPART;
		fputwlsplit(WLPART, ambfp);
		fputncomp(NCSAMP, ambfp);
		fputformat(AMBFMT, ambfp);
		fputc('\n', ambfp);
		putambmagic(ambfp);
	} else if (getheader(ambfp, amb_headline, NULL) < 0 || !hasambmagic(ambfp)) {
		if (--ntries > 0 && ftell(ambfp) == 0) {
			clearerr(ambfp);
			sleep(2);
			goto retry;
		}
		error(USER, "bad/incompatible ambient file");
	}
	if ((AMB_CNDX != CNDX) | (AMB_WLPART != WLPART)) {
		if (setspectrsamp(AMB_CNDX, AMB_WLPART) < 0)
			error(USER, "bad wavelength sampling in ambient file");
		if (AMB_CNDX[3] == CNDX[3] && FABSEQ(AMB_WLPART[0],WLPART[0]) &&
					FABSEQ(AMB_WLPART[3],WLPART[3])) {
			AMB_CNDX = CNDX;
			AMB_WLPART = WLPART;		/* just the same */
		} else
			error(WARNING, "different ambient file wavelength sampling");
	}
}


static void
avsave(				/* insert and save an ambient value */
	AMBVAL	*av
)
{
	avstore(av);
	if (ambfp == NULL)
		return;
	if (writambval(av, ambfp) < 0)
		goto writerr;
	if (++nunflshed >= AMBFLUSH)
		if (ambsync() == EOF)
			goto writerr;
	return;
writerr:
	error(SYSTEM, "error writing to ambient file");
}


static AMBVAL *
avstore(				/* allocate memory and save aval */
	AMBVAL  *aval
)
{
	AMBVAL  *av;
	double	d;

	if ((av = newambval()) == NULL)
		error(SYSTEM, "out of memory in avstore");
	memcpy(av, aval, AVSIZE);	/* AVSIZE <= sizeof(AMBVAL) */
	av->next = NULL;
	nambvals++;
	d = pbright(av->val);
	if (d > FTINY) {		/* add to log sum for averaging */
		avsum += log(d);
		navsum++;
	}
	avinsert(av);			/* insert in our cache tree */
	return(av);
}


#define ATALLOCSZ	512		/* #/8 trees to allocate at once */

static AMBTREE  *atfreelist = NULL;	/* free ambient tree structures */


static AMBTREE *
newambtree(void)				/* allocate 8 ambient tree structs */
{
	AMBTREE  *atp, *upperlim;

	if (atfreelist == NULL) {	/* get more nodes */
		atfreelist = (AMBTREE *)malloc(ATALLOCSZ*8*sizeof(AMBTREE));
		if (atfreelist == NULL)
			return(NULL);
					/* link new free list */
		upperlim = atfreelist + 8*(ATALLOCSZ-1);
		for (atp = atfreelist; atp < upperlim; atp += 8)
			atp->kid = atp + 8;
		atp->kid = NULL;
	}
	atp = atfreelist;
	atfreelist = atp->kid;
	memset(atp, 0, 8*sizeof(AMBTREE));
	return(atp);
}


static void
freeambtree(			/* free 8 ambient tree structs */
	AMBTREE  *atp
)
{
	atp->kid = atfreelist;
	atfreelist = atp;
}


static void
unloadatree(			/* unload an ambient value tree */
	AMBTREE  *at,
	unloadtf_t *f
)
{
	AMBVAL  *av;
	int  i;
					/* transfer values at this node */
	for (av = at->alist; av != NULL; av = at->alist) {
		at->alist = av->next;
		av->next = NULL;
		(*f)(av);
	}
	if (at->kid == NULL)
		return;
	for (i = 0; i < 8; i++)		/* transfer and free children */
		unloadatree(at->kid+i, f);
	freeambtree(at->kid);
	at->kid = NULL;
}


static void
avfree(AMBVAL *av)
{
	free(av);
}


static void
sortambvals(void)			/* resort ambient values */
{
	AMBTREE  oldatrunk = atrunk;

	atrunk.alist = NULL;
	atrunk.kid = NULL;
	unloadatree(&oldatrunk, avinsert);
}


int
ambsync(void)			/* synchronize ambient file */
{
	off_t	newpos;
	int	n;
	AMBVAL	avs;

	if (ambfp == NULL)	/* no ambient file? */
		return(0);

	if (nunflshed > 0) {	/* append new values? */
		if (fflush(ambfp) < 0)
			return(EOF);
	} else if (fseeko(ambfp, 0, SEEK_END) < 0)
		goto seekerr;

	if ((newpos = ftello(ambfp)) < 0)
		goto seekerr;
				/* how many others added? */
	n = (newpos - lastpos)/AMBVALSIZ - nunflshed;
	nunflshed = 0;
	if (n <= 0) {		/* no one helping this time? */
		lastpos = newpos;
		return(0);
	}
	if (ambinp == NULL) {	/* else need to open for input? */
		ambinp = fopen(ambfile, "r");
		if (ambinp == NULL) {
			sprintf(errmsg, "cannot reopen ambient file \"%s\"",
					ambfile);
			error(SYSTEM, errmsg);
		}
		SET_FILE_BINARY(ambinp);
	}
				/* read from last endpoint */
	if (fseeko(ambinp, lastpos, SEEK_SET) < 0)
		goto seekerr;
	while (n-- > 0) {	/* load new contributed values */
		if (!readambval(&avs, ambinp)) {
			sprintf(errmsg, "ambient file \"%s\" corrupted",
					ambfile);
			error(WARNING, errmsg);
			break;
		}
		avstore(&avs);
	}
	lastpos = newpos;	/* update endpoint */
	return(0);
seekerr:
	error(SYSTEM, "seek failed in ambsync");
	return(EOF);	/* pro forma return */
}
