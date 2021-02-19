static const char	RCSid[] = "$Id: ambient.c,v 2.110 2021/02/19 22:05:46 greg Exp $";
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

extern char  *shm_boundary;	/* memory sharing boundary */

#ifndef  MAXASET
#define	 MAXASET	4095	/* maximum number of elements in ambient set */
#endif
OBJECT	ambset[MAXASET+1]={0};	/* ambient include/exclude set */

double	maxarad;		/* maximum ambient radius */
double	minarad;		/* minimum ambient radius */

static AMBTREE	atrunk;		/* our ambient trunk node */

static FILE  *ambfp = NULL;	/* ambient file pointer */
static int  nunflshed = 0;	/* number of unflushed ambient values */

#ifndef SORT_THRESH
#ifdef SMLMEM
#define SORT_THRESH	((16L<<20)/sizeof(AMBVAL))
#else
#define SORT_THRESH	((64L<<20)/sizeof(AMBVAL))
#endif
#endif
#ifndef SORT_INTVL
#define SORT_INTVL	(SORT_THRESH<<1)
#endif
#ifndef MAX_SORT_INTVL
#define MAX_SORT_INTVL	(SORT_INTVL<<6)
#endif


static double  avsum = 0.;		/* computed ambient value sum (log) */
static unsigned int  navsum = 0;	/* number of values in avsum */
static unsigned int  nambvals = 0;	/* total number of indirect values */
static unsigned int  nambshare = 0;	/* number of values from file */
static unsigned long  ambclock = 0;	/* ambient access clock */
static unsigned long  lastsort = 0;	/* time of last value sort */
static long  sortintvl = SORT_INTVL;	/* time until next sort */
static FILE  *ambinp = NULL;		/* auxiliary file for input */
static long  lastpos = -1;		/* last flush position */

#define MAXACLOCK	(1L<<30)	/* clock turnover value */
	/*
	 * Track access times unless we are sharing ambient values
	 * through memory on a multiprocessor, when we want to avoid
	 * claiming our own memory (copy on write).  Go ahead anyway
	 * if more than two thirds of our values are unshared.
	 * Compile with -Dtracktime=0 to turn this code off.
	 */
#ifndef tracktime
#define tracktime	(shm_boundary == NULL || nambvals > 3*nambshare)
#endif

#define	 AMBFLUSH	(BUFSIZ/AMBVALSIZ)

#define	 newambval()	(AMBVAL *)malloc(sizeof(AMBVAL))

#define  tfunc(lwr, x, upr)	(((x)-(lwr))/((upr)-(lwr)))

static void initambfile(int creat);
static void avsave(AMBVAL *av);
static AMBVAL *avstore(AMBVAL  *aval);
static AMBTREE *newambtree(void);
static void freeambtree(AMBTREE  *atp);

typedef void unloadtf_t(AMBVAL *);
static unloadtf_t avinsert;
static unloadtf_t av2list;
static unloadtf_t avfree;
static void unloadatree(AMBTREE  *at, unloadtf_t *f);

static int aposcmp(const void *avp1, const void *avp2);
static int avlmemi(AMBVAL *avaddr);
static void sortambvals(int always);

static int	plugaleak(RAY *r, AMBVAL *ap, FVECT anorm, double ang);
static double	sumambient(COLOR acol, RAY *r, FVECT rn, int al,
				AMBTREE *at, FVECT c0, double s);
static int	makeambient(COLOR acol, RAY *r, FVECT rn, int al);
static int	extambient(COLOR cr, AMBVAL *ap, FVECT pv, FVECT nv, 
				FVECT uvw[3]);

#ifdef  F_SETLKW
static void aflock(int  typ);
#endif


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
			sortambvals(1);		/* rebuild tree */
	}
}


void
setambient(void)				/* initialize calculation */
{
	int	readonly = 0;
	long	flen;
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
						/* open ambient file */
	if ((ambfp = fopen(ambfile, "r+")) == NULL)
		readonly = (ambfp = fopen(ambfile, "r")) != NULL;
	if (ambfp != NULL) {
		initambfile(0);			/* file exists */
		lastpos = ftell(ambfp);
		while (readambval(&amb, ambfp))
			avstore(&amb);
		nambshare = nambvals;		/* share loaded values */
		if (readonly) {
			sprintf(errmsg,
				"loaded %u values from read-only ambient file",
					nambvals);
			error(WARNING, errmsg);
			fclose(ambfp);		/* close file so no writes */
			ambfp = NULL;
			return;			/* avoid ambsync() */
		}
						/* align file pointer */
		lastpos += (long)nambvals*AMBVALSIZ;
		flen = lseek(fileno(ambfp), (off_t)0, SEEK_END);
		if (flen != lastpos) {
			sprintf(errmsg,
			"ignoring last %ld values in ambient file (corrupted)",
					(flen - lastpos)/AMBVALSIZ);
			error(WARNING, errmsg);
			fseek(ambfp, lastpos, SEEK_SET);
			ftruncate(fileno(ambfp), (off_t)lastpos);
		}
	} else if ((ambfp = fopen(ambfile, "w+")) != NULL) {
		initambfile(1);			/* else create new file */
		fflush(ambfp);
		lastpos = ftell(ambfp);
	} else {
		sprintf(errmsg, "cannot open ambient file \"%s\"", ambfile);
		error(SYSTEM, errmsg);
	}
#ifdef  F_SETLKW
	aflock(F_UNLCK);			/* release file */
#endif
}


void
ambdone(void)			/* close ambient file and free memory */
{
	if (ambfp != NULL) {		/* close ambient file */
		ambsync();
		fclose(ambfp);
		ambfp = NULL;
		if (ambinp != NULL) {	
			fclose(ambinp);
			ambinp = NULL;
		}
		lastpos = -1;
	}
					/* free ambient tree */
	unloadatree(&atrunk, avfree);
					/* reset state variables */
	avsum = 0.;
	navsum = 0;
	nambvals = 0;
	nambshare = 0;
	ambclock = 0;
	lastsort = 0;
	sortintvl = SORT_INTVL;
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
	COLOR  aval,
	RAY  *r,
	FVECT  nrm
)
{
	static int  rdepth = 0;			/* ambient recursion */
	COLOR	acol, caustic;
	int	i, ok;
	double	d, l;

	/* PMAP: Factor in ambient from photon map, if enabled and ray is
	 * ambient. Return as all ambient components accounted for, else
	 * continue. */
	if (ambPmap(aval, r, rdepth))
		return;

	/* PMAP: Factor in specular-diffuse ambient (caustics) from photon
	 * map, if enabled and ray is primary, else caustic is zero.  Continue
	 * with RADIANCE ambient calculation */
	copycolor(caustic, aval);
	ambPmapCaustic(caustic, r, rdepth);
	
	if (ambdiv <= 0)			/* no ambient calculation */
		goto dumbamb;
						/* check number of bounces */
	if (rdepth >= ambounce)
		goto dumbamb;
						/* check ambient list */
	if (ambincl != -1 && r->ro != NULL &&
			ambincl != inset(ambset, r->ro->omod))
		goto dumbamb;

	if (ambacc <= FTINY) {			/* no ambient storage */
		FVECT	uvd[2];
		float	dgrad[2], *dgp = NULL;

		if (nrm != r->ron && DOT(nrm,r->ron) < 0.9999)
			dgp = dgrad;		/* compute rotational grad. */
		copycolor(acol, aval);
		rdepth++;
		ok = doambient(acol, r, r->rweight,
				uvd, NULL, NULL, dgp, NULL);
		rdepth--;
		if (!ok)
			goto dumbamb;
		if ((ok > 0) & (dgp != NULL)) {	/* apply texture */
			FVECT	v1;
			VCROSS(v1, r->ron, nrm);
			d = 1.0;
			for (i = 3; i--; )
				d += v1[i] * (dgp[0]*uvd[0][i] + dgp[1]*uvd[1][i]);
			if (d >= 0.05)
				scalecolor(acol, d);
		}
		copycolor(aval, acol);

		/* PMAP: add in caustic */
		addcolor(aval, caustic);
		return;
	}

	if (tracktime)				/* sort to minimize thrashing */
		sortambvals(0);
						/* interpolate ambient value */
	setcolor(acol, 0.0, 0.0, 0.0);
	d = sumambient(acol, r, nrm, rdepth,
			&atrunk, thescene.cuorg, thescene.cusize);
			
	if (d > FTINY) {
		d = 1.0/d;
		scalecolor(acol, d);
		multcolor(aval, acol);

		/* PMAP: add in caustic */
		addcolor(aval, caustic);
		return;
	}
	
	rdepth++;				/* need to cache new value */
	ok = makeambient(acol, r, nrm, rdepth-1);
	rdepth--;
	
	if (ok) {
		multcolor(aval, acol);		/* computed new value */

		/* PMAP: add in caustic */
		addcolor(aval, caustic);
		return;
	}
	
dumbamb:					/* return global value */
	if ((ambvwt <= 0) | (navsum == 0)) {
		multcolor(aval, ambval);
		
		/* PMAP: add in caustic */
		addcolor(aval, caustic);
		return;
	}
	
	l = bright(ambval);			/* average in computations */	
	if (l > FTINY) {
		d = (log(l)*(double)ambvwt + avsum) /
				(double)(ambvwt + navsum);
		d = exp(d) / l;
		scalecolor(aval, d);
		multcolor(aval, ambval);	/* apply color of ambval */
	} else {
		d = exp( avsum / (double)navsum );
		scalecolor(aval, d);		/* neutral color */
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
		if (m != NULL && !istransp(m->otype) && !isBSDFproxy(m) &&
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
	COLOR  acol,
	RAY  *r,
	FVECT  rn,
	int  al,
	AMBTREE	 *at,
	FVECT  c0,
	double	s
)
{			/* initial limit is 10 degrees plus ambacc radians */
	const double	minangle = 10.0 * PI/180.;
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
		if (wsum >= 0.05 && s > minarad*10.0)
			return(wsum);
	}
					/* adjust maximum angle */
	if (at->alist != NULL && (at->alist->lvl <= al) & (r->rweight < 0.6))
		maxangle = (maxangle - PI/2.)*pow(r->rweight,0.13) + PI/2.;
					/* sum this node */
	for (av = at->alist; av != NULL; av = av->next) {
		double	u, v, d, delta_r2, delta_t2;
		COLOR	ct;
		FVECT	uvw[3];
					/* record access */
		if (tracktime)
			av->latick = ambclock;
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
		d = DOT(uvw[2], r->ron);
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
		if (d < -minarad*ambacc-.001)
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
		if (!extambient(ct, av, r->rop, rn, uvw))
			continue;
		d = tfunc(maxangle, sqrt(delta_r2), 0.0) *
			tfunc(ambacc, sqrt(delta_t2), 0.0);
		scalecolor(ct, d);
		addcolor(acol, ct);
		wsum += d;
	}
	return(wsum);
}


static int
makeambient(		/* make a new ambient value for storage */
	COLOR  acol,
	RAY  *r,
	FVECT  rn,
	int  al
)
{
	AMBVAL	amb;
	FVECT	uvw[3];
	int	i;

	amb.weight = 1.0;			/* compute weight */
	for (i = al; i-- > 0; )
		amb.weight *= AVGREFL;
	if (r->rweight < 0.1*amb.weight)	/* heuristic override */
		amb.weight = 1.25*r->rweight;
	setcolor(acol, AVGREFL, AVGREFL, AVGREFL);
						/* compute ambient */
	i = doambient(acol, r, amb.weight,
			uvw, amb.rad, amb.gpos, amb.gdir, &amb.corral);
	scalecolor(acol, 1./AVGREFL);		/* undo assumed reflectance */
	if (i <= 0 || amb.rad[0] <= FTINY)	/* no Hessian or zero radius */
		return(i);
						/* store value */
	VCOPY(amb.pos, r->rop);
	amb.ndir = encodedir(r->ron);
	amb.udir = encodedir(uvw[0]);
	amb.lvl = al;
	copycolor(amb.val, acol);
						/* insert into tree */
	avsave(&amb);				/* and save to file */
	if (rn != r->ron) {			/* texture */
		VCOPY(uvw[2], r->ron);
		extambient(acol, &amb, r->rop, rn, uvw);
	}
	return(1);
}


static int
extambient(		/* extrapolate value at pv, nv */
	COLOR  cr,
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
	copycolor(cr, ap->val);
	scalecolor(cr, d);
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

#ifdef	F_SETLKW
	aflock(cre8 ? F_WRLCK : F_RDLCK);
#endif
	SET_FILE_BINARY(ambfp);
	if (mybuf == NULL)
		mybuf = (char *)bmalloc(BUFSIZ+8);
	setbuf(ambfp, mybuf);
	if (cre8) {			/* new file */
		newheader("RADIANCE", ambfp);
		fprintf(ambfp, "%s -av %g %g %g -aw %d -ab %d -aa %g ",
				progname, colval(ambval,RED),
				colval(ambval,GRN), colval(ambval,BLU),
				ambvwt, ambounce, ambacc);
		fprintf(ambfp, "-ad %d -as %d -ar %d ",
				ambdiv, ambssamp, ambres);
		if (octname != NULL)
			fputs(octname, ambfp);
		fputc('\n', ambfp);
		fprintf(ambfp, "SOFTWARE= %s\n", VersionID);
		fputnow(ambfp);
		fputformat(AMBFMT, ambfp);
		fputc('\n', ambfp);
		putambmagic(ambfp);
	} else if (checkheader(ambfp, AMBFMT, NULL) < 0 || !hasambmagic(ambfp))
		error(USER, "bad ambient file");
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
	*av = *aval;
	av->latick = ambclock;
	av->next = NULL;
	nambvals++;
	d = bright(av->val);
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


static struct avl {
	AMBVAL	*p;
	unsigned long	t;
}	*avlist1;			/* ambient value list with ticks */
static AMBVAL	**avlist2;		/* memory positions for sorting */
static int	i_avlist;		/* index for lists */

static int alatcmp(const void *av1, const void *av2);

static void
avfree(AMBVAL *av)
{
	free(av);
}

static void
av2list(
	AMBVAL *av
)
{
#ifdef DEBUG
	if (i_avlist >= nambvals)
		error(CONSISTENCY, "too many ambient values in av2list1");
#endif
	avlist1[i_avlist].p = avlist2[i_avlist] = (AMBVAL*)av;
	avlist1[i_avlist++].t = av->latick;
}


static int
alatcmp(			/* compare ambient values for MRA */
	const void *av1,
	const void *av2
)
{
	long  lc = ((struct avl *)av2)->t - ((struct avl *)av1)->t;
	return(lc<0 ? -1 : lc>0 ? 1 : 0);
}


/* GW NOTE 2002/10/3:
 * I used to compare AMBVAL pointers, but found that this was the
 * cause of a serious consistency error with gcc, since the optimizer
 * uses some dangerous trick in pointer subtraction that
 * assumes pointers differ by exact struct size increments.
 */
static int
aposcmp(			/* compare ambient value positions */
	const void	*avp1,
	const void	*avp2
)
{
	long	diff = *(char * const *)avp1 - *(char * const *)avp2;
	if (diff < 0)
		return(-1);
	return(diff > 0);
}


static int
avlmemi(				/* find list position from address */
	AMBVAL	*avaddr
)
{
	AMBVAL  **avlpp;

	avlpp = (AMBVAL **)bsearch(&avaddr, avlist2,
			nambvals, sizeof(AMBVAL *), aposcmp);
	if (avlpp == NULL)
		error(CONSISTENCY, "address not found in avlmemi");
	return(avlpp - avlist2);
}


static void
sortambvals(			/* resort ambient values */
	int	always
)
{
	AMBTREE  oldatrunk;
	AMBVAL	tav, *tap, *pnext;
	int	i, j;
					/* see if it's time yet */
	if (!always && (ambclock++ < lastsort+sortintvl ||
			nambvals < SORT_THRESH))
		return;
	/*
	 * The idea here is to minimize memory thrashing
	 * in VM systems by improving reference locality.
	 * We do this by periodically sorting our stored ambient
	 * values in memory in order of most recently to least
	 * recently accessed.  This ordering was chosen so that new
	 * ambient values (which tend to be less important) go into
	 * higher memory with the infrequently accessed values.
	 *	Since we expect our values to need sorting less
	 * frequently as the process continues, we double our
	 * waiting interval after each call.
	 * 	This routine is also called by setambacc() with
	 * the "always" parameter set to 1 so that the ambient
	 * tree will be rebuilt with the new accuracy parameter.
	 */
	if (tracktime) {		/* allocate pointer arrays to sort */
		avlist2 = (AMBVAL **)malloc(nambvals*sizeof(AMBVAL *));
		avlist1 = (struct avl *)malloc(nambvals*sizeof(struct avl));
	} else {
		avlist2 = NULL;
		avlist1 = NULL;
	}
	if (avlist1 == NULL) {		/* no time tracking -- rebuild tree? */
		if (avlist2 != NULL)
			free(avlist2);
		if (always) {		/* rebuild without sorting */
			oldatrunk = atrunk;
			atrunk.alist = NULL;
			atrunk.kid = NULL;
			unloadatree(&oldatrunk, avinsert);
		}
	} else {			/* sort memory by last access time */
		/*
		 * Sorting memory is tricky because it isn't contiguous.
		 * We have to sort an array of pointers by MRA and also
		 * by memory position.  We then copy values in "loops"
		 * to minimize memory hits.  Nevertheless, we will visit
		 * everyone at least twice, and this is an expensive process
		 * when we're thrashing, which is when we need to do it.
		 */
#ifdef DEBUG
		sprintf(errmsg, "sorting %u ambient values at ambclock=%lu...",
				nambvals, ambclock);
		eputs(errmsg);
#endif
		i_avlist = 0;
		unloadatree(&atrunk, av2list);	/* empty current tree */
#ifdef DEBUG
		if (i_avlist < nambvals)
			error(CONSISTENCY, "missing ambient values in sortambvals");
#endif
		qsort(avlist1, nambvals, sizeof(struct avl), alatcmp);
		qsort(avlist2, nambvals, sizeof(AMBVAL *), aposcmp);
		for (i = 0; i < nambvals; i++) {
			if (avlist1[i].p == NULL)
				continue;
			tap = avlist2[i];
			tav = *tap;
			for (j = i; (pnext = avlist1[j].p) != tap;
					j = avlmemi(pnext)) {
				*(avlist2[j]) = *pnext;
				avinsert(avlist2[j]);
				avlist1[j].p = NULL;
			}
			*(avlist2[j]) = tav;
			avinsert(avlist2[j]);
			avlist1[j].p = NULL;
		}
		free(avlist1);
		free(avlist2);
						/* compute new sort interval */
		sortintvl = ambclock - lastsort;
		if (sortintvl >= MAX_SORT_INTVL/2)
			sortintvl = MAX_SORT_INTVL;
		else
			sortintvl <<= 1;	/* wait twice as long next */
#ifdef DEBUG
		eputs("done\n");
#endif
	}
	if (ambclock >= MAXACLOCK)
		ambclock = MAXACLOCK/2;
	lastsort = ambclock;
}


#ifdef	F_SETLKW

static void
aflock(			/* lock/unlock ambient file */
	int  typ
)
{
	static struct flock  fls;	/* static so initialized to zeroes */

	if (typ == fls.l_type)		/* already called? */
		return;

	fls.l_type = typ;
	do
		if (fcntl(fileno(ambfp), F_SETLKW, &fls) != -1)
			return;
	while (errno == EINTR);
	
	error(SYSTEM, "cannot (un)lock ambient file");
}


int
ambsync(void)			/* synchronize ambient file */
{
	long  flen;
	AMBVAL	avs;
	int  n;

	if (ambfp == NULL)	/* no ambient file? */
		return(0);
				/* gain appropriate access */
	aflock(nunflshed ? F_WRLCK : F_RDLCK);
				/* see if file has grown */
	if ((flen = lseek(fileno(ambfp), (off_t)0, SEEK_END)) < 0)
		goto seekerr;
	if ((n = flen - lastpos) > 0) {		/* file has grown */
		if (ambinp == NULL) {		/* get new file pointer */
			ambinp = fopen(ambfile, "rb");
			if (ambinp == NULL)
				error(SYSTEM, "fopen failed in ambsync");
		}
		if (fseek(ambinp, lastpos, SEEK_SET) < 0)
			goto seekerr;
		while (n >= AMBVALSIZ) {	/* load contributed values */
			if (!readambval(&avs, ambinp)) {
				sprintf(errmsg,
			"ambient file \"%s\" corrupted near character %ld",
						ambfile, flen - n);
				error(WARNING, errmsg);
				break;
			}
			avstore(&avs);
			n -= AMBVALSIZ;
		}
		lastpos = flen - n;		/* check alignment */
		if (n && lseek(fileno(ambfp), (off_t)lastpos, SEEK_SET) < 0)
			goto seekerr;
	}
	n = fflush(ambfp);			/* calls write() at last */
	lastpos += (long)nunflshed*AMBVALSIZ;
	aflock(F_UNLCK);			/* release file */
	nunflshed = 0;
	return(n);
seekerr:
	error(SYSTEM, "seek failed in ambsync");
	return(EOF);	/* pro forma return */
}

#else	/* ! F_SETLKW */

int
ambsync(void)			/* flush ambient file */
{
	if (ambfp == NULL)
		return(0);
	nunflshed = 0;
	return(fflush(ambfp));
}

#endif	/* ! F_SETLKW */
