#ifndef lint
static const char RCSid[] = "$Id: source.c,v 2.44 2004/03/30 16:13:01 schorsch Exp $";
#endif
/*
 *  source.c - routines dealing with illumination sources.
 *
 *  External symbols declared in source.h
 */

#include  "ray.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "source.h"
#include  "random.h"

extern double  ssampdist;		/* scatter sampling distance */

#ifndef MAXSSAMP
#define MAXSSAMP	16		/* maximum samples per ray */
#endif

/*
 * Structures used by direct()
 */

typedef struct {
	int  sno;		/* source number */
	FVECT  dir;		/* source direction */
	COLOR  coef;		/* material coefficient */
	COLOR  val;		/* contribution */
}  CONTRIB;		/* direct contribution */

typedef struct {
	int  sndx;		/* source index (to CONTRIB array) */
	float  brt;		/* brightness (for comparison) */
}  CNTPTR;		/* contribution pointer */

static CONTRIB  *srccnt;		/* source contributions in direct() */
static CNTPTR  *cntord;			/* source ordering in direct() */
static int  maxcntr = 0;		/* size of contribution arrays */

static int cntcmp(const void *p1, const void *p2);


extern OBJREC *			/* find an object's actual material */
findmaterial(register OBJREC *o)
{
	while (!ismaterial(o->otype)) {
		if (ismixture(o->otype))
			return(NULL);   /* reject mixed materials */
		if (o->otype == MOD_ALIAS && o->oargs.nsargs) {
			OBJECT  aobj;
			OBJREC  *ao;
			aobj = lastmod(objndx(o), o->oargs.sarg[0]);
			if (aobj < 0)
				objerror(o, USER, "bad reference");
			ao = objptr(aobj);
			if (ismaterial(ao->otype))
				return(ao);
		}
		if (o->omod == OVOID)
			return(NULL);
		o = objptr(o->omod);
	}
	return(o);
}


extern void
marksources(void)			/* find and mark source objects */
{
	int  foundsource = 0;
	int  i;
	register OBJREC  *o, *m;
	register int  ns;
					/* initialize dispatch table */
	initstypes();
					/* find direct sources */
	for (i = 0; i < nsceneobjs; i++) {
	
		o = objptr(i);

		if (!issurface(o->otype) || o->omod == OVOID)
			continue;
					/* find material */
		m = findmaterial(objptr(o->omod));
		if (m == NULL || !islight(m->otype))
			continue;	/* not source modifier */
	
		if (m->oargs.nfargs != (m->otype == MAT_GLOW ? 4 :
				m->otype == MAT_SPOT ? 7 : 3))
			objerror(m, USER, "bad # arguments");

		if (m->otype == MAT_GLOW &&
				o->otype != OBJ_SOURCE &&
				m->oargs.farg[3] <= FTINY)
			continue;			/* don't bother */
		if (m->oargs.farg[0] <= FTINY && m->oargs.farg[1] <= FTINY &&
				m->oargs.farg[2] <= FTINY)
			continue;			/* don't bother */

		if (sfun[o->otype].of == NULL ||
				sfun[o->otype].of->setsrc == NULL)
			objerror(o, USER, "illegal material");

		if ((ns = newsource()) < 0)
			goto memerr;

		setsource(&source[ns], o);

		if (m->otype == MAT_GLOW) {
			source[ns].sflags |= SPROX;
			source[ns].sl.prox = m->oargs.farg[3];
			if (source[ns].sflags & SDISTANT)
				source[ns].sflags |= SSKIP;
		} else if (m->otype == MAT_SPOT) {
			source[ns].sflags |= SSPOT;
			if ((source[ns].sl.s = makespot(m)) == NULL)
				goto memerr;
			if (source[ns].sflags & SFLAT &&
				!checkspot(source[ns].sl.s,source[ns].snorm)) {
				objerror(o, WARNING,
					"invalid spotlight direction");
				source[ns].sflags |= SSKIP;
			}
		}
#if  SHADCACHE
		source[ns].obscache = NULL;
#endif
		if (!(source[ns].sflags & SSKIP))
			foundsource++;
	}
	if (!foundsource) {
		error(WARNING, "no light sources found");
		return;
	}
	markvirtuals();			/* find and add virtual sources */
				/* allocate our contribution arrays */
	maxcntr = nsources + MAXSPART;	/* start with this many */
	srccnt = (CONTRIB *)malloc(maxcntr*sizeof(CONTRIB));
	cntord = (CNTPTR *)malloc(maxcntr*sizeof(CNTPTR));
	if ((srccnt == NULL) | (cntord == NULL))
		goto memerr;
	return;
memerr:
	error(SYSTEM, "out of memory in marksources");
}


extern void
freesources(void)			/* free all source structures */
{
	if (nsources > 0) {
#if SHADCACHE
		while (nsources--)
			freeobscache(&source[nsources]);
#endif
		free((void *)source);
		source = NULL;
		nsources = 0;
	}
	if (maxcntr <= 0)
		return;
	free((void *)srccnt);
	srccnt = NULL;
	free((void *)cntord);
	cntord = NULL;
	maxcntr = 0;
}


extern int
srcray(				/* send a ray to a source, return domega */
	register RAY  *sr,		/* returned source ray */
	RAY  *r,			/* ray which hit object */
	SRCINDEX  *si			/* source sample index */
)
{
    double  d;				/* distance to source */
    register SRCREC  *srcp;

    rayorigin(sr, r, SHADOW, 1.0);		/* ignore limits */

    while ((d = nextssamp(sr, si)) != 0.0) {
	sr->rsrc = si->sn;			/* remember source */
	srcp = source + si->sn;
	if (srcp->sflags & SDISTANT) {
		if (srcp->sflags & SSPOT && spotout(sr, srcp->sl.s))
			continue;
		return(1);		/* sample OK */
	}
				/* local source */
						/* check proximity */
	if (srcp->sflags & SPROX && d > srcp->sl.prox)
		continue;
						/* check angle */
	if (srcp->sflags & SSPOT) {
		if (spotout(sr, srcp->sl.s))
			continue;
					/* adjust solid angle */
		si->dom *= d*d;
		d += srcp->sl.s->flen;
		si->dom /= d*d;
	}
	return(1);			/* sample OK */
    }
    return(0);			/* no more samples */
}


extern void
srcvalue(			/* punch ray to source and compute value */
	register RAY  *r
)
{
	register SRCREC  *sp;

	sp = &source[r->rsrc];
	if (sp->sflags & SVIRTUAL) {	/* virtual source */
					/* check intersection */
		if (!(*ofun[sp->so->otype].funp)(sp->so, r))
			return;
		if (!rayshade(r, r->ro->omod))	/* compute contribution */
			goto nomat;
		rayparticipate(r);
		return;
	}
					/* compute intersection */
	if (sp->sflags & SDISTANT ? sourcehit(r) :
			(*ofun[sp->so->otype].funp)(sp->so, r)) {
		if (sp->sa.success >= 0)
			sp->sa.success++;
		if (!rayshade(r, r->ro->omod))	/* compute contribution */
			goto nomat;
		rayparticipate(r);
		return;
	}
					/* we missed our mark! */
	if (sp->sa.success < 0)
		return;			/* bitched already */
	sp->sa.success -= AIMREQT;
	if (sp->sa.success >= 0)
		return;			/* leniency */
	sprintf(errmsg, "aiming failure for light source \"%s\"",
			sp->so->oname);
	error(WARNING, errmsg);		/* issue warning */
	return;
nomat:
	objerror(r->ro, USER, "material not found");
}


extern int
sourcehit(			/* check to see if ray hit distant source */
	register RAY  *r
)
{
	int  first, last;
	register int  i;

	if (r->rsrc >= 0) {		/* check only one if aimed */
		first = last = r->rsrc;
	} else {			/* otherwise check all */
		first = 0; last = nsources-1;
	}
	for (i = first; i <= last; i++)
		if ((source[i].sflags & (SDISTANT|SVIRTUAL)) == SDISTANT)
			/*
			 * Check to see if ray is within
			 * solid angle of source.
			 */
			if (2.0*PI * (1.0 - DOT(source[i].sloc,r->rdir))
					<= source[i].ss2) {
				r->ro = source[i].so;
				if (!(source[i].sflags & SSKIP))
					break;
			}

	if (r->ro != NULL) {
		r->robj = objndx(r->ro);
		for (i = 0; i < 3; i++)
			r->ron[i] = -r->rdir[i];
		r->rod = 1.0;
		r->pert[0] = r->pert[1] = r->pert[2] = 0.0;
		r->uv[0] = r->uv[1] = 0.0;
		r->rox = NULL;
		return(1);
	}
	return(0);
}


static int
cntcmp(				/* contribution compare (descending) */
	const void *p1,
	const void *p2
)
{
	register const CNTPTR  *sc1 = (const CNTPTR *)p1;
	register const CNTPTR  *sc2 = (const CNTPTR *)p2;

	if (sc1->brt > sc2->brt)
		return(-1);
	if (sc1->brt < sc2->brt)
		return(1);
	return(0);
}


extern void
direct(					/* add direct component */
	RAY  *r,			/* ray that hit surface */
	srcdirf_t *f,			/* direct component coefficient function */
	void  *p			/* data for f */
)
{
	register int  sn;
	register CONTRIB  *scp;
	SRCINDEX  si;
	int  nshadcheck, ncnts;
	int  nhits;
	double  prob, ourthresh, hwt;
	RAY  sr;
			/* NOTE: srccnt and cntord global so no recursion */
	if (nsources <= 0)
		return;		/* no sources?! */
						/* potential contributions */
	initsrcindex(&si);
	for (sn = 0; srcray(&sr, r, &si); sn++) {
		if (sn >= maxcntr) {
			maxcntr = sn + MAXSPART;
			srccnt = (CONTRIB *)realloc((void *)srccnt,
					maxcntr*sizeof(CONTRIB));
			cntord = (CNTPTR *)realloc((void *)cntord,
					maxcntr*sizeof(CNTPTR));
			if ((srccnt == NULL) | (cntord == NULL))
				error(SYSTEM, "out of memory in direct");
		}
		cntord[sn].sndx = sn;
		scp = srccnt + sn;
		scp->sno = sr.rsrc;
						/* compute coefficient */
		(*f)(scp->coef, p, sr.rdir, si.dom);
		cntord[sn].brt = bright(scp->coef);
		if (cntord[sn].brt <= 0.0)
			continue;
#if SHADCACHE
						/* check shadow cache */
		if (si.np == 1 && srcblocked(&sr)) {
			cntord[sn].brt = 0.0;
			continue;
		}
#endif
		VCOPY(scp->dir, sr.rdir);
						/* compute potential */
		sr.revf = srcvalue;
		rayvalue(&sr);
		copycolor(scp->val, sr.rcol);
		multcolor(scp->val, scp->coef);
		cntord[sn].brt = bright(scp->val);
	}
						/* sort contributions */
	qsort(cntord, sn, sizeof(CNTPTR), cntcmp);
	{					/* find last */
		register int  l, m;

		ncnts = l = sn;
		sn = 0;
		while ((m = (sn + ncnts) >> 1) != l) {
			if (cntord[m].brt > 0.0)
				sn = m;
			else
				ncnts = m;
			l = m;
		}
	}
	if (ncnts == 0)
		return;		/* no contributions! */
                                                /* accumulate tail */
        for (sn = ncnts-1; sn > 0; sn--)
                cntord[sn-1].brt += cntord[sn].brt;
						/* compute number to check */
	nshadcheck = pow((double)ncnts, shadcert) + .5;
						/* modify threshold */
	ourthresh = shadthresh / r->rweight;
						/* test for shadows */
	for (nhits = 0, hwt = 0.0, sn = 0; sn < ncnts;
			hwt += (double)source[scp->sno].nhits /
				(double)source[scp->sno].ntests,
			sn++) {
						/* check threshold */
		if ((sn+nshadcheck>=ncnts ? cntord[sn].brt :
				cntord[sn].brt-cntord[sn+nshadcheck].brt)
				< ourthresh*bright(r->rcol))
			break;
		scp = srccnt + cntord[sn].sndx;
						/* test for hit */
		rayorigin(&sr, r, SHADOW, 1.0);
		VCOPY(sr.rdir, scp->dir);
		sr.rsrc = scp->sno;
						/* keep statistics */
		if (source[scp->sno].ntests++ > 0xfffffff0) {
			source[scp->sno].ntests >>= 1;
			source[scp->sno].nhits >>= 1;
		}
		if (localhit(&sr, &thescene) &&
				( sr.ro != source[scp->sno].so ||
				source[scp->sno].sflags & SFOLLOW )) {
						/* follow entire path */
			raycont(&sr);
			rayparticipate(&sr);
			if (trace != NULL)
				(*trace)(&sr);	/* trace execution */
			if (bright(sr.rcol) <= FTINY) {
#if SHADCACHE
				if ((scp <= srccnt || scp[-1].sno != scp->sno)
						&& (scp >= srccnt+ncnts-1 ||
						    scp[1].sno != scp->sno))
					srcblocker(&sr);
#endif
				continue;	/* missed! */
			}
			copycolor(scp->val, sr.rcol);
			multcolor(scp->val, scp->coef);
		}
						/* add contribution if hit */
		addcolor(r->rcol, scp->val);
		nhits++;
		source[scp->sno].nhits++;
	}
					/* source hit rate */
	if (hwt > FTINY)
		hwt = (double)nhits / hwt;
	else
		hwt = 0.5;
#ifdef DEBUG
	sprintf(errmsg, "%d tested, %d untested, %f conditional hit rate\n",
			sn, ncnts-sn, hwt);
	eputs(errmsg);
#endif
					/* add in untested sources */
	for ( ; sn < ncnts; sn++) {
		scp = srccnt + cntord[sn].sndx;
		prob = hwt * (double)source[scp->sno].nhits /
				(double)source[scp->sno].ntests;
		if (prob > 1.0)
			prob = 1.0;
		scalecolor(scp->val, prob);
		addcolor(r->rcol, scp->val);
	}
}


extern void
srcscatter(			/* compute source scattering into ray */
	register RAY  *r
)
{
	int  oldsampndx;
	int  nsamps;
	RAY  sr;
	SRCINDEX  si;
	double  t, d;
	double  re, ge, be;
	COLOR  cvext;
	int  i, j;

	if (r->slights == NULL || r->slights[0] == 0
			|| r->gecc >= 1.-FTINY || r->rot >= FHUGE)
		return;
	if (ssampdist <= FTINY || (nsamps = r->rot/ssampdist + .5) < 1)
		nsamps = 1;
#if MAXSSAMP
	else if (nsamps > MAXSSAMP)
		nsamps = MAXSSAMP;
#endif
	oldsampndx = samplendx;
	samplendx = random()&0x7fff;		/* randomize */
	for (i = r->slights[0]; i > 0; i--) {	/* for each source */
		for (j = 0; j < nsamps; j++) {	/* for each sample position */
			samplendx++;
			t = r->rot * (j+frandom())/nsamps;
							/* extinction */
			re = t*colval(r->cext,RED);
			ge = t*colval(r->cext,GRN);
			be = t*colval(r->cext,BLU);
			setcolor(cvext,	re > 92. ? 0. : exp(-re),
					ge > 92. ? 0. : exp(-ge),
					be > 92. ? 0. : exp(-be));
			if (intens(cvext) <= FTINY)
				break;			/* too far away */
			sr.rorg[0] = r->rorg[0] + r->rdir[0]*t;
			sr.rorg[1] = r->rorg[1] + r->rdir[1]*t;
			sr.rorg[2] = r->rorg[2] + r->rdir[2]*t;
			sr.rmax = 0.;
			initsrcindex(&si);	/* sample ray to this source */
			si.sn = r->slights[i];
			nopart(&si, &sr);
			if (!srcray(&sr, NULL, &si) ||
					sr.rsrc != r->slights[i])
				continue;		/* no path */
			copycolor(sr.cext, r->cext);
			copycolor(sr.albedo, r->albedo);
			sr.gecc = r->gecc;
			sr.slights = r->slights;
			rayvalue(&sr);			/* eval. source ray */
			if (bright(sr.rcol) <= FTINY)
				continue;
			if (r->gecc <= FTINY)		/* compute P(theta) */
				d = 1.;
			else {
				d = DOT(r->rdir, sr.rdir);
				d = 1. + r->gecc*r->gecc - 2.*r->gecc*d;
				d = (1. - r->gecc*r->gecc) / (d*sqrt(d));
			}
							/* other factors */
			d *= si.dom * r->rot / (4.*PI*nsamps);
			multcolor(sr.rcol, r->cext);
			multcolor(sr.rcol, r->albedo);
			scalecolor(sr.rcol, d);
			multcolor(sr.rcol, cvext);
			addcolor(r->rcol, sr.rcol);	/* add it in */
		}
	}
	samplendx = oldsampndx;
}


/****************************************************************
 * The following macros were separated from the m_light() routine
 * because they are very nasty and difficult to understand.
 */

/* illumblock *
 *
 * We cannot allow an illum to pass to another illum, because that
 * would almost certainly constitute overcounting.
 * However, we do allow an illum to pass to another illum
 * that is actually going to relay to a virtual light source.
 * We also prevent an illum from passing to a glow; this provides a
 * convenient mechanism for defining detailed light source
 * geometry behind (or inside) an effective radiator.
 */

static int
weaksrcmat(OBJECT obj)		/* identify material */
{
	OBJREC *o = findmaterial(objptr(obj));
	
	if (o == NULL) return 0;
	return((o->otype==MAT_ILLUM)|(o->otype==MAT_GLOW));
}

#define  illumblock(m, r)	(!(source[r->rsrc].sflags&SVIRTUAL) && \
				r->rod > 0.0 && \
				weaksrcmat(source[r->rsrc].so->omod))

/* wrongsource *
 *
 * This source is the wrong source (ie. overcounted) if we are
 * aimed to a different source than the one we hit and the one
 * we hit is not an illum that should be passed.
 */

#define  wrongsource(m, r)	(r->rsrc>=0 && source[r->rsrc].so!=r->ro && \
				(m->otype!=MAT_ILLUM || illumblock(m,r)))

/* distglow *
 *
 * A distant glow is an object that sometimes acts as a light source,
 * but is too far away from the test point to be one in this case.
 * (Glows with negative radii should NEVER participate in illumination.)
 */

#define  distglow(m, r, d)	(m->otype==MAT_GLOW && \
				m->oargs.farg[3] >= -FTINY && \
				d > m->oargs.farg[3])

/* badcomponent *
 *
 * We must avoid counting light sources in the ambient calculation,
 * since the direct component is handled separately.  Therefore, any
 * ambient ray which hits an active light source must be discarded.
 * The same is true for stray specular samples, since the specular
 * contribution from light sources is calculated separately.
 */

#define  badcomponent(m, r)	(r->crtype&(AMBIENT|SPECULAR) && \
				!(r->crtype&SHADOW || r->rod < 0.0 || \
		/* not 100% correct */	distglow(m, r, r->rot)))

/* passillum *
 *
 * An illum passes to another material type when we didn't hit it
 * on purpose (as part of a direct calculation), or it is relaying
 * a virtual light source.
 */

#define  passillum(m, r)	(m->otype==MAT_ILLUM && \
				(r->rsrc<0 || source[r->rsrc].so!=r->ro || \
				source[r->rsrc].sflags&SVIRTUAL))

/* srcignore *
 *
 * The -dv flag is normally on for sources to be visible.
 */

#define  srcignore(m, r)	!(directvis || r->crtype&SHADOW || \
				distglow(m, r, raydist(r,PRIMARY)))


extern int
m_light(				/* ray hit a light source */
	register OBJREC  *m,
	register RAY  *r
)
{
						/* check for over-counting */
	if (badcomponent(m, r))
		return(1);
	if (wrongsource(m, r))
		return(1);
						/* check for passed illum */
	if (passillum(m, r)) {
		if (m->oargs.nsargs && strcmp(m->oargs.sarg[0], VOIDID))
			return(rayshade(r,lastmod(objndx(m),m->oargs.sarg[0])));
		raytrans(r);
		return(1);
	}
					/* otherwise treat as source */
						/* check for behind */
	if (r->rod < 0.0)
		return(1);
						/* check for invisibility */
	if (srcignore(m, r))
		return(1);
						/* check for outside spot */
	if (m->otype==MAT_SPOT && spotout(r, makespot(m)))
		return(1);
						/* get distribution pattern */
	raytexture(r, m->omod);
						/* get source color */
	setcolor(r->rcol, m->oargs.farg[0],
			  m->oargs.farg[1],
			  m->oargs.farg[2]);
						/* modify value */
	multcolor(r->rcol, r->pcol);
	return(1);
}
