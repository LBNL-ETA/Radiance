#ifndef lint
static const char	RCSid[] = "$Id: ambient.c,v 2.61 2005/05/31 18:01:09 greg Exp $";
#endif
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
#include  "resolu.h"
#include  "ambient.h"
#include  "random.h"

#ifndef  OCTSCALE
#define	 OCTSCALE	1.0	/* ceil((valid rad.)/(cube size)) */
#endif

extern char  *shm_boundary;	/* memory sharing boundary */

#ifndef  MAXASET
#define	 MAXASET	2047	/* maximum number of elements in ambient set */
#endif
OBJECT	ambset[MAXASET+1]={0};	/* ambient include/exclude set */

double	maxarad;		/* maximum ambient radius */
double	minarad;		/* minimum ambient radius */

static AMBTREE	atrunk;		/* our ambient trunk node */

static FILE  *ambfp = NULL;	/* ambient file pointer */
static int  nunflshed = 0;	/* number of unflushed ambient values */

#ifndef SORT_THRESH
#ifdef SMLMEM
#define SORT_THRESH	((3L<<20)/sizeof(AMBVAL))
#else
#define SORT_THRESH	((9L<<20)/sizeof(AMBVAL))
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
#define  freeav(av)	free((void *)av);

static void initambfile(int creat);
static void avsave(AMBVAL *av);
static AMBVAL *avstore(AMBVAL  *aval);
static AMBTREE *newambtree(void);
static void freeambtree(AMBTREE  *atp);

typedef void unloadtf_t(void *);
static unloadtf_t avinsert;
static unloadtf_t av2list;
static void unloadatree(AMBTREE  *at, unloadtf_t *f);

static int aposcmp(const void *avp1, const void *avp2);
static int avlmemi(AMBVAL *avaddr);
static void sortambvals(int always);

#ifdef  F_SETLKW
static void aflock(int  typ);
#endif


extern void
setambres(				/* set ambient resolution */
	int  ar
)
{
	ambres = ar < 0 ? 0 : ar;		/* may be done already */
						/* set min & max radii */
	if (ar <= 0) {
		minarad = 0;
		maxarad = thescene.cusize / 2.0;
	} else {
		minarad = thescene.cusize / ar;
		maxarad = 64 * minarad;			/* heuristic */
		if (maxarad > thescene.cusize / 2.0)
			maxarad = thescene.cusize / 2.0;
	}
	if (minarad <= FTINY)
		minarad = 10*FTINY;
	if (maxarad <= minarad)
		maxarad = 64 * minarad;
}


extern void
setambacc(				/* set ambient accuracy */
	double  newa
)
{
	double  ambdiff;

	if (newa < 0.0)
		newa = 0.0;
	ambdiff = fabs(newa - ambacc);
	if (ambdiff >= .01 && (ambacc = newa) > FTINY && nambvals > 0)
		sortambvals(1);			/* rebuild tree */
}


extern void
setambient(void)				/* initialize calculation */
{
	int	readonly = 0;
	long  pos, flen;
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
		pos = ftell(ambfp);
		while (readambval(&amb, ambfp))
			avinsert(avstore(&amb));
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
		pos += (long)nambvals*AMBVALSIZ;
		flen = lseek(fileno(ambfp), (off_t)0, SEEK_END);
		if (flen != pos) {
			sprintf(errmsg,
			"ignoring last %ld values in ambient file (corrupted)",
					(flen - pos)/AMBVALSIZ);
			error(WARNING, errmsg);
			fseek(ambfp, pos, 0);
#ifndef _WIN32 /* XXX we need a replacement for that one */
			ftruncate(fileno(ambfp), (off_t)pos);
#endif
		}
	} else if ((ambfp = fopen(ambfile, "w+")) != NULL) {
		initambfile(1);			/* else create new file */
	} else {
		sprintf(errmsg, "cannot open ambient file \"%s\"", ambfile);
		error(SYSTEM, errmsg);
	}
	nunflshed++;	/* lie */
	ambsync();
}


extern void
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
	unloadatree(&atrunk, free);
					/* reset state variables */
	avsum = 0.;
	navsum = 0;
	nambvals = 0;
	nambshare = 0;
	ambclock = 0;
	lastsort = 0;
	sortintvl = SORT_INTVL;
}


extern void
ambnotify(			/* record new modifier */
	OBJECT	obj
)
{
	static int  hitlimit = 0;
	register OBJREC	 *o;
	register char  **amblp;

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


extern void
multambient(		/* compute ambient component & multiply by coef. */
	COLOR  aval,
	register RAY  *r,
	FVECT  nrm
)
{
	static int  rdepth = 0;			/* ambient recursion */
	COLOR	acol;
	double	d, l;

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
		copycolor(acol, aval);
		rdepth++;
		d = doambient(acol, r, r->rweight, NULL, NULL);
		rdepth--;
		if (d <= FTINY)
			goto dumbamb;
		copycolor(aval, acol);
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
		return;
	}
	rdepth++;				/* need to cache new value */
	d = makeambient(acol, r, nrm, rdepth-1);
	rdepth--;
	if (d > FTINY) {
		multcolor(aval, acol);		/* got new value */
		return;
	}
dumbamb:					/* return global value */
	if ((ambvwt <= 0) | (navsum == 0)) {
		multcolor(aval, ambval);
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


extern double
sumambient(	/* get interpolated ambient value */
	COLOR  acol,
	register RAY  *r,
	FVECT  rn,
	int  al,
	AMBTREE	 *at,
	FVECT  c0,
	double	s
)
{
	double	d, e1, e2, wt, wsum;
	COLOR  ct;
	FVECT  ck0;
	int  i;
	register int  j;
	register AMBVAL	 *av;

	wsum = 0.0;
					/* do this node */
	for (av = at->alist; av != NULL; av = av->next) {
		double	rn_dot = -2.0;
		if (tracktime)
			av->latick = ambclock;
		/*
		 *  Ambient level test.
		 */
		if (av->lvl > al)	/* list sorted, so this works */
			break;
		if (av->weight < 0.9*r->rweight)
			continue;
		/*
		 *  Ambient radius test.
		 */
		d = av->pos[0] - r->rop[0];
		e1 = d * d;
		d = av->pos[1] - r->rop[1];
		e1 += d * d;
		d = av->pos[2] - r->rop[2];
		e1 += d * d;
		e1 /= av->rad * av->rad;
		if (e1 > ambacc*ambacc*1.21)
			continue;
		/*
		 *  Direction test using closest normal.
		 */
		d = DOT(av->dir, r->ron);
		if (rn != r->ron) {
			rn_dot = DOT(av->dir, rn);
			if (rn_dot > 1.0-FTINY)
				rn_dot = 1.0-FTINY;
			if (rn_dot >= d-FTINY) {
				d = rn_dot;
				rn_dot = -2.0;
			}
		}
		e2 = (1.0 - d) * r->rweight;
		if (e2 < 0.0) e2 = 0.0;
		if (e1 + e2 > ambacc*ambacc*1.21)
			continue;
		/*
		 *  Ray behind test.
		 */
		d = 0.0;
		for (j = 0; j < 3; j++)
			d += (r->rop[j] - av->pos[j]) *
					(av->dir[j] + r->ron[j]);
		if (d*0.5 < -minarad*ambacc-.001)
			continue;
		/*
		 *  Jittering final test reduces image artifacts.
		 */
		e1 = sqrt(e1);
		e2 = sqrt(e2);
		wt = e1 + e2;
		if (wt > ambacc*(.9+.2*urand(9015+samplendx)))
			continue;
		/*
		 *  Recompute directional error using perturbed normal
		 */
		if (rn_dot > 0.0) {
			e2 = sqrt((1.0 - rn_dot)*r->rweight);
			wt = e1 + e2;
		}
		if (wt <= 1e-3)
			wt = 1e3;
		else
			wt = 1.0 / wt;
		wsum += wt;
		extambient(ct, av, r->rop, rn);
		scalecolor(ct, wt);
		addcolor(acol, ct);
	}
	if (at->kid == NULL)
		return(wsum);
					/* do children */
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
	return(wsum);
}


extern double
makeambient(		/* make a new ambient value for storage */
	COLOR  acol,
	RAY  *r,
	FVECT  rn,
	int  al
)
{
	AMBVAL	amb;
	FVECT	gp, gd;
	int	i;

	amb.weight = 1.0;			/* compute weight */
	for (i = al; i-- > 0; )
		amb.weight *= AVGREFL;
	if (r->rweight < 0.1*amb.weight)	/* heuristic override */
		amb.weight = r->rweight;
	setcolor(acol, AVGREFL, AVGREFL, AVGREFL);
						/* compute ambient */
	amb.rad = doambient(acol, r, amb.weight, gp, gd);
	if (amb.rad <= FTINY) {
		setcolor(acol, 0.0, 0.0, 0.0);
		return(0.0);
	}
	scalecolor(acol, 1./AVGREFL);		/* undo assumed reflectance */
						/* store value */
	VCOPY(amb.pos, r->rop);
	VCOPY(amb.dir, r->ron);
	amb.lvl = al;
	copycolor(amb.val, acol);
	VCOPY(amb.gpos, gp);
	VCOPY(amb.gdir, gd);
						/* insert into tree */
	avsave(&amb);				/* and save to file */
	if (rn != r->ron)
		extambient(acol, &amb, r->rop, rn);	/* texture */
	return(amb.rad);
}


extern void
extambient(		/* extrapolate value at pv, nv */
	COLOR  cr,
	register AMBVAL	 *ap,
	FVECT  pv,
	FVECT  nv
)
{
	FVECT  v1;
	register int  i;
	double	d;

	d = 1.0;			/* zeroeth order */
					/* gradient due to translation */
	for (i = 0; i < 3; i++)
		d += ap->gpos[i]*(pv[i]-ap->pos[i]);
					/* gradient due to rotation */
	VCROSS(v1, ap->dir, nv);
	d += DOT(ap->gdir, v1);
	if (d <= 0.0) {
		setcolor(cr, 0.0, 0.0, 0.0);
		return;
	}
	copycolor(cr, ap->val);
	scalecolor(cr, d);
}


static void
initambfile(		/* initialize ambient file */
	int  creat
)
{
	extern char  *progname, *octname;
	static char  *mybuf = NULL;

#ifdef	F_SETLKW
	aflock(creat ? F_WRLCK : F_RDLCK);
#endif
	SET_FILE_BINARY(ambfp);
	if (mybuf == NULL)
		mybuf = (char *)bmalloc(BUFSIZ+8);
	setbuf(ambfp, mybuf);
	if (creat) {			/* new file */
		newheader("RADIANCE", ambfp);
		fprintf(ambfp, "%s -av %g %g %g -aw %d -ab %d -aa %g ",
				progname, colval(ambval,RED),
				colval(ambval,GRN), colval(ambval,BLU),
				ambvwt, ambounce, ambacc);
		fprintf(ambfp, "-ad %d -as %d -ar %d ",
				ambdiv, ambssamp, ambres);
		if (octname != NULL)
			printargs(1, &octname, ambfp);
		else
			fputc('\n', ambfp);
		fprintf(ambfp, "SOFTWARE= %s\n", VersionID);
		fputnow(ambfp);
		fputformat(AMBFMT, ambfp);
		putc('\n', ambfp);
		putambmagic(ambfp);
	} else if (checkheader(ambfp, AMBFMT, NULL) < 0 || !hasambmagic(ambfp))
		error(USER, "bad ambient file");
}


static void
avsave(				/* insert and save an ambient value */
	AMBVAL	*av
)
{
	avinsert(avstore(av));
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
avstore(				/* allocate memory and store aval */
	register AMBVAL  *aval
)
{
	register AMBVAL  *av;
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
	return(av);
}


#define ATALLOCSZ	512		/* #/8 trees to allocate at once */

static AMBTREE  *atfreelist = NULL;	/* free ambient tree structures */


static AMBTREE *
newambtree(void)				/* allocate 8 ambient tree structs */
{
	register AMBTREE  *atp, *upperlim;

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
	memset((char *)atp, '\0', 8*sizeof(AMBTREE));
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
avinsert(				/* insert ambient value in our tree */
	void *av
)
{
	register AMBTREE  *at;
	register AMBVAL  *ap;
	AMBVAL  avh;
	FVECT  ck0;
	double	s;
	int  branch;
	register int  i;

	if (((AMBVAL*)av)->rad <= FTINY)
		error(CONSISTENCY, "zero ambient radius in avinsert");
	at = &atrunk;
	VCOPY(ck0, thescene.cuorg);
	s = thescene.cusize;
	while (s*(OCTSCALE/2) > ((AMBVAL*)av)->rad*ambacc) {
		if (at->kid == NULL)
			if ((at->kid = newambtree()) == NULL)
				error(SYSTEM, "out of memory in avinsert");
		s *= 0.5;
		branch = 0;
		for (i = 0; i < 3; i++)
			if (((AMBVAL*)av)->pos[i] > ck0[i] + s) {
				ck0[i] += s;
				branch |= 1 << i;
			}
		at = at->kid + branch;
	}
	avh.next = at->alist;		/* order by increasing level */
	for (ap = &avh; ap->next != NULL; ap = ap->next)
		if (ap->next->lvl >= ((AMBVAL*)av)->lvl)
			break;
	((AMBVAL*)av)->next = ap->next;
	ap->next = (AMBVAL*)av;
	at->alist = avh.next;
}


static void
unloadatree(			/* unload an ambient value tree */
	register AMBTREE  *at,
	unloadtf_t *f
)
{
	register AMBVAL  *av;
	register int  i;
					/* transfer values at this node */
	for (av = at->alist; av != NULL; av = at->alist) {
		at->alist = av->next;
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
av2list(
	void *av
)
{
#ifdef DEBUG
	if (i_avlist >= nambvals)
		error(CONSISTENCY, "too many ambient values in av2list1");
#endif
	avlist1[i_avlist].p = avlist2[i_avlist] = (AMBVAL*)av;
	avlist1[i_avlist++].t = ((AMBVAL*)av)->latick;
}


static int
alatcmp(			/* compare ambient values for MRA */
	const void *av1,
	const void *av2
)
{
	register long  lc = ((struct avl *)av2)->t - ((struct avl *)av1)->t;
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
	register long	diff = *(char * const *)avp1 - *(char * const *)avp2;
	if (diff < 0)
		return(-1);
	return(diff > 0);
}

#if 1
static int
avlmemi(				/* find list position from address */
	AMBVAL	*avaddr
)
{
	register AMBVAL  **avlpp;

	avlpp = (AMBVAL **)bsearch((char *)&avaddr, (char *)avlist2,
			nambvals, sizeof(AMBVAL *), aposcmp);
	if (avlpp == NULL)
		error(CONSISTENCY, "address not found in avlmemi");
	return(avlpp - avlist2);
}
#else
#define avlmemi(avaddr)	((AMBVAL **)bsearch((char *)&avaddr,(char *)avlist2, \
				nambvals,sizeof(AMBVAL *),aposcmp) - avlist2)
#endif


static void
sortambvals(			/* resort ambient values */
	int	always
)
{
	AMBTREE  oldatrunk;
	AMBVAL	tav, *tap, *pnext;
	register int	i, j;
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
			free((void *)avlist2);
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
		qsort((char *)avlist1, nambvals, sizeof(struct avl), alatcmp);
		qsort((char *)avlist2, nambvals, sizeof(AMBVAL *), aposcmp);
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
		free((void *)avlist1);
		free((void *)avlist2);
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

	fls.l_type = typ;
	if (fcntl(fileno(ambfp), F_SETLKW, &fls) < 0)
		error(SYSTEM, "cannot (un)lock ambient file");
}


extern int
ambsync(void)			/* synchronize ambient file */
{
	long  flen;
	AMBVAL	avs;
	register int  n;

	if (nunflshed == 0)
		return(0);
	if (lastpos < 0)	/* initializing (locked in initambfile) */
		goto syncend;
				/* gain exclusive access */
	aflock(F_WRLCK);
				/* see if file has grown */
	if ((flen = lseek(fileno(ambfp), (off_t)0, SEEK_END)) < 0)
		goto seekerr;
	if ( (n = flen - lastpos) ) {		/* file has grown */
		if (ambinp == NULL) {		/* use duplicate filedes */
			ambinp = fdopen(dup(fileno(ambfp)), "r");
			if (ambinp == NULL)
				error(SYSTEM, "fdopen failed in ambsync");
		}
		if (fseek(ambinp, lastpos, 0) < 0)
			goto seekerr;
		while (n >= AMBVALSIZ) {	/* load contributed values */
			if (!readambval(&avs, ambinp)) {
				sprintf(errmsg,
			"ambient file \"%s\" corrupted near character %ld",
						ambfile, flen - n);
				error(WARNING, errmsg);
				break;
			}
			avinsert(avstore(&avs));
			n -= AMBVALSIZ;
		}
		/*** seek always as safety measure
		if (n) ***/			/* alignment */
			if (lseek(fileno(ambfp), (off_t)(flen-n), SEEK_SET) < 0)
				goto seekerr;
	}
#ifdef  DEBUG
	if (ambfp->_ptr - ambfp->_base != nunflshed*AMBVALSIZ) {
		sprintf(errmsg, "ambient file buffer at %d rather than %d",
				ambfp->_ptr - ambfp->_base,
				nunflshed*AMBVALSIZ);
		error(CONSISTENCY, errmsg);
	}
#endif
syncend:
	n = fflush(ambfp);			/* calls write() at last */
	if ((lastpos = lseek(fileno(ambfp), (off_t)0, SEEK_CUR)) < 0)
		goto seekerr;
	aflock(F_UNLCK);			/* release file */
	nunflshed = 0;
	return(n);
seekerr:
	error(SYSTEM, "seek failed in ambsync");
	return -1; /* pro forma return */
}

#else

extern int
ambsync(void)			/* flush ambient file */
{
	if (nunflshed == 0)
		return(0);
	nunflshed = 0;
	return(fflush(ambfp));
}

#endif
