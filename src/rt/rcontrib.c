#ifndef lint
static const char RCSid[] = "$Id: rcontrib.c,v 2.11 2012/06/19 01:27:13 greg Exp $";
#endif
/*
 * Accumulate ray contributions for a set of materials
 * Initialization and calculation routines
 */

#include "rcontrib.h"
#include "source.h"
#include "otypes.h"
#include "platform.h"

char	*shm_boundary = NULL;		/* boundary of shared memory */

CUBE	thescene;			/* our scene */
OBJECT	nsceneobjs;			/* number of objects in our scene */

int	dimlist[MAXDIM];		/* sampling dimensions */
int	ndims = 0;			/* number of sampling dimensions */
int	samplendx = 0;			/* index for this sample */

static void	trace_contrib(RAY *r);	/* our trace callback */
void	(*trace)() = trace_contrib;

int	do_irrad = 0;			/* compute irradiance? */

int	rand_samp = 1;			/* pure Monte Carlo sampling? */

double	dstrsrc = 0.9;			/* square source distribution */
double	shadthresh = .03;		/* shadow threshold */
double	shadcert = .75;			/* shadow certainty */
int	directrelay = 3;		/* number of source relays */
int	vspretest = 512;		/* virtual source pretest density */
int	directvis = 1;			/* sources visible? */
double	srcsizerat = .2;		/* maximum ratio source size/dist. */

COLOR	cextinction = BLKCOLOR;		/* global extinction coefficient */
COLOR	salbedo = BLKCOLOR;		/* global scattering albedo */
double	seccg = 0.;			/* global scattering eccentricity */
double	ssampdist = 0.;			/* scatter sampling distance */

double	specthresh = .15;		/* specular sampling threshold */
double	specjitter = 1.;		/* specular sampling jitter */

int	backvis = 1;			/* back face visibility */

int	maxdepth = -10;			/* maximum recursion depth */
double	minweight = 2e-3;		/* minimum ray weight */

char	*ambfile = NULL;		/* ambient file name */
COLOR	ambval = BLKCOLOR;		/* ambient value */
int	ambvwt = 0;			/* initial weight for ambient value */
double	ambacc = 0;			/* ambient accuracy */
int	ambres = 256;			/* ambient resolution */
int	ambdiv = 350;			/* ambient divisions */
int	ambssamp = 0;			/* ambient super-samples */
int	ambounce = 1;			/* ambient bounces */
char	*amblist[AMBLLEN+1];		/* ambient include/exclude list */
int	ambincl = -1;			/* include == 1, exclude == 0 */

int	account;			/* current accumulation count */
RNUMBER	raysleft;			/* number of rays left to trace */
long	waitflush;			/* how long until next flush */

RNUMBER	lastray = 0;			/* last ray number sent */
RNUMBER	lastdone = 0;			/* last ray output */

static void mcfree(void *p) { epfree((*(MODCONT *)p).binv); free(p); }

LUTAB	modconttab = LU_SINIT(NULL,mcfree);	/* modifier lookup table */

/************************** INITIALIZATION ROUTINES ***********************/

char *
formstr(				/* return format identifier */
	int  f
)
{
	switch (f) {
	case 'a': return("ascii");
	case 'f': return("float");
	case 'd': return("double");
	case 'c': return(COLRFMT);
	}
	return("unknown");
}


/* Add modifier to our list to track */
MODCONT *
addmodifier(char *modn, char *outf, char *binv, int bincnt)
{
	LUENT	*lep = lu_find(&modconttab,modn);
	MODCONT	*mp;
	EPNODE	*ebinv;
	int	i;
	
	if (lep->data != NULL) {
		sprintf(errmsg, "duplicate modifier '%s'", modn);
		error(USER, errmsg);
	}
	if (nmods >= MAXMODLIST)
		error(INTERNAL, "too many modifiers");
	modname[nmods++] = modn;	/* XXX assumes static string */
	lep->key = modn;		/* XXX assumes static string */
	if (binv == NULL)
		binv = "0";		/* use single bin if unspecified */
	ebinv = eparse(binv);
	if (ebinv->type == NUM) {	/* check value if constant */
		bincnt = (int)(evalue(ebinv) + 1.5);
		if (bincnt != 1) {
			sprintf(errmsg, "illegal non-zero constant for bin (%s)",
					binv);
			error(USER, errmsg);
		}
	} else if (bincnt <= 0) {
		sprintf(errmsg,
			"unspecified or illegal bin count for modifier '%s'",
				modn);
		error(USER, errmsg);
	}
					/* initialize results holder */
	mp = (MODCONT *)malloc(sizeof(MODCONT)+sizeof(DCOLOR)*(bincnt-1));
	if (mp == NULL)
		error(SYSTEM, "out of memory in addmodifier");
	mp->outspec = outf;		/* XXX assumes static string */
	mp->modname = modn;		/* XXX assumes static string */
	mp->binv = ebinv;
	mp->nbins = bincnt;
	memset(mp->cbin, 0, sizeof(DCOLOR)*bincnt);
					/* allocate output streams */
	for (i = bincnt; i-- > 0; )
		getostream(mp->outspec, mp->modname, i, 1);
	lep->data = (char *)mp;
	return(mp);
}


/* Add modifiers from a file list */
void
addmodfile(char *fname, char *outf, char *binv, int bincnt)
{
	char	*mname[MAXMODLIST];
	int	i;
					/* find the file & store strings */
	if (wordfile(mname, getpath(fname, getrlibpath(), R_OK)) < 0) {
		sprintf(errmsg, "cannot find modifier file '%s'", fname);
		error(SYSTEM, errmsg);
	}
	for (i = 0; mname[i]; i++)	/* add each one */
		addmodifier(mname[i], outf, binv, bincnt);
}


void
quit(			/* quit program */
	int  code
)
{
	if (nchild > 0)		/* close children if any */
		end_children(code != 0);
	exit(code);
}


/* Initialize our process(es) */
static void
rcinit()
{
	int	i;

	if (nproc > MAXPROCESS)
		sprintf(errmsg, "too many processes requested -- reducing to %d",
				nproc = MAXPROCESS);
	if (nproc > 1) {
		preload_objs();		/* preload auxiliary data */
					/* set shared memory boundary */
		shm_boundary = strcpy((char *)malloc(16), "SHM_BOUNDARY");
	}
	if (yres > 0) {			/* set up flushing & ray counts */
		if (xres > 0)
			raysleft = (RNUMBER)xres*yres;
		else
			raysleft = yres;
	} else
		raysleft = 0;
	if ((account = accumulate) > 1)
		raysleft *= accumulate;
	waitflush = (yres > 0) & (xres > 1) ? 0 : xres;
					/* tracing to sources as well */
	for (i = 0; i < nsources; i++)
		source[i].sflags |= SFOLLOW;

	if (nproc > 1 && in_rchild())	/* forked child? */
		return;			/* return to main processing loop */

	if (recover) {			/* recover previous output? */
		if (accumulate <= 0)
			reload_output();
		else
			recover_output();
	}
	if (nproc == 1)			/* single process? */
		return;
					/* else run appropriate controller */
	if (accumulate <= 0)
		feeder_loop();
	else
		parental_loop();
	quit(0);			/* parent musn't return! */
}

/************************** MAIN CALCULATION PROCESS ***********************/

/* Our trace call to sum contributions */
static void
trace_contrib(RAY *r)
{
	MODCONT	*mp;
	int	bn;
	RREAL	contr[3];

	if (r->ro == NULL || r->ro->omod == OVOID)
		return;

	mp = (MODCONT *)lu_find(&modconttab,objptr(r->ro->omod)->oname)->data;

	if (mp == NULL)				/* not in our list? */
		return;

	worldfunc(RCCONTEXT, r);		/* get bin number */
	bn = (int)(evalue(mp->binv) + .5);
	if ((bn < 0) | (bn >= mp->nbins)) {
		error(WARNING, "bad bin number (ignored)");
		return;
	}
	raycontrib(contr, r, PRIMARY);
	if (contrib)
		multcolor(contr, r->rcol);
	addcolor(mp->cbin[bn], contr);
}


/* Evaluate irradiance contributions */
static void
eval_irrad(FVECT org, FVECT dir)
{
	RAY	thisray;

	VSUM(thisray.rorg, org, dir, 1.1e-4);
	thisray.rdir[0] = -dir[0];
	thisray.rdir[1] = -dir[1];
	thisray.rdir[2] = -dir[2];
	thisray.rmax = 0.0;
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	thisray.rot = 1e-5;		/* pretend we hit surface */
	thisray.rod = 1.0;
	VSUM(thisray.rop, org, dir, 1e-4);
	samplendx++;			/* compute result */
	(*ofun[Lamb.otype].funp)(&Lamb, &thisray);
}


/* Evaluate radiance contributions */
static void
eval_rad(FVECT org, FVECT dir, double dmax)
{
	RAY	thisray;
					/* set up ray */
	VCOPY(thisray.rorg, org);
	VCOPY(thisray.rdir, dir);
	thisray.rmax = dmax;
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	samplendx++;			/* call ray evaluation */
	rayvalue(&thisray);
}


/* Accumulate and/or output ray contributions (child or only process) */
static void
done_contrib()
{
	MODCONT	*mp;
	int	i;

	if (account <= 0 || --account)
		return;			/* not time yet */

	for (i = 0; i < nmods; i++) {	/* output records & clear */
		mp = (MODCONT *)lu_find(&modconttab,modname[i])->data;
		mod_output(mp);
		memset(mp->cbin, 0, sizeof(DCOLOR)*mp->nbins);
	}
	end_record();			/* end lines & flush if time */

	account = accumulate;		/* reset accumulation counter */
}


/* Principal calculation loop (called by main) */
void
rcontrib()
{
	static int	ignore_warning_given = 0;
	FVECT		orig, direc;
	double		d;
					/* initialize (& fork more of us) */
	rcinit();
					/* load rays from stdin & process */
#ifdef getc_unlocked
	flockfile(stdin);		/* avoid mutex overhead */
#endif
	while (getvec(orig) == 0 && getvec(direc) == 0) {
		d = normalize(direc);
		if (nchild != -1 && (d == 0.0) & (accumulate != 1)) {
			if (!ignore_warning_given++)
				error(WARNING,
				"dummy ray(s) ignored during accumulation\n");
			continue;
		}
		if (lastray+1 < lastray)
			lastray = lastdone = 0;
		++lastray;
		if (d == 0.0) {				/* zero ==> flush */
			if ((yres <= 0) | (xres <= 0))
				waitflush = 1;		/* flush right after */
			account = 1;
		} else if (imm_irrad) {			/* else compute */
			eval_irrad(orig, direc);
		} else {
			eval_rad(orig, direc, lim_dist ? d : 0.0);
		}
		done_contrib();		/* accumulate/output */
		++lastdone;
		if (raysleft && !--raysleft)
			break;		/* preemptive EOI */
	}
	if (nchild != -1 && (accumulate <= 0) | (account < accumulate)) {
		if (account < accumulate) {
			error(WARNING, "partial accumulation in final record");
			accumulate -= account;
		}
		account = 1;		/* output accumulated totals */
		done_contrib();
	}
	lu_done(&ofiletab);		/* close output files */
	if (raysleft)
		error(USER, "unexpected EOF on input");
}
