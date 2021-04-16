#ifndef lint
static const char	RCSid[] = "$Id: rtrace.c,v 2.102 2021/04/16 16:31:35 greg Exp $";
#endif
/*
 *  rtrace.c - program and variables for individual ray tracing.
 */

#include "copyright.h"

/*
 *  Input is in the form:
 *
 *	xorg	yorg	zorg	xdir	ydir	zdir
 *
 *  The direction need not be normalized.  Output is flexible.
 *  If the direction vector is (0,0,0), then the output is flushed.
 *  All values default to ascii representation of real
 *  numbers.  Binary representations can be selected
 *  with '-ff' for float or '-fd' for double.  By default,
 *  radiance is computed.  The '-i' or '-I' options indicate that
 *  irradiance values are desired.
 */

#include  <time.h>

#include  "platform.h"
#include  "ray.h"
#include  "ambient.h"
#include  "source.h"
#include  "otypes.h"
#include  "otspecial.h"
#include  "resolu.h"
#include  "random.h"

extern int  inform;			/* input format */
extern int  outform;			/* output format */
extern char  *outvals;			/* output values */

extern int  imm_irrad;			/* compute immediate irradiance? */
extern int  lim_dist;			/* limit distance? */

extern char  *tralist[];		/* list of modifers to trace (or no) */
extern int  traincl;			/* include == 1, exclude == 0 */

extern int  hresolu;			/* horizontal resolution */
extern int  vresolu;			/* vertical resolution */

int  castonly = 0;			/* only doing ray-casting? */

#ifndef  MAXTSET
#define	 MAXTSET	8191		/* maximum number in trace set */
#endif
OBJECT	traset[MAXTSET+1]={0};		/* trace include/exclude set */

static int  Tflag = 0;			/* source tracing enabled? */

static RAY  thisray;			/* for our convenience */

static FILE  *inpfp = NULL;		/* input stream pointer */

static FVECT	*inp_queue = NULL;	/* ray input queue if flushing */
static int	inp_qpos = 0;		/* next ray to return */
static int	inp_qend = 0;		/* number of rays in this work group */

typedef void putf_t(RREAL *v, int n);
static putf_t puta, putd, putf, putrgbe;

typedef void oputf_t(RAY *r);
static oputf_t  oputo, oputd, oputv, oputV, oputl, oputL, oputc, oputp,
		oputr, oputR, oputx, oputX, oputn, oputN, oputs,
		oputw, oputW, oputm, oputM, oputtilde;

extern void tranotify(OBJECT obj);
static int is_fifo(FILE *fp);
static void bogusray(void);
static void raycast(RAY *r);
static void rayirrad(RAY *r);
static void rtcompute(FVECT org, FVECT dir, double dmax);
static int printvals(RAY *r);
static int getvec(FVECT vec, int fmt, FILE *fp);
static int iszerovec(const FVECT vec);
static double nextray(FVECT org, FVECT dir);
static void tabin(RAY *r);
static void ourtrace(RAY *r);

static oputf_t *ray_out[32], *every_out[32];
static putf_t *putreal;


void
quit(			/* quit program */
	int  code
)
{
	if (ray_pnprocs > 0)	/* close children if any */
		ray_pclose(0);		
	else if (ray_pnprocs < 0)
		_exit(code);	/* avoid flush() in child */
#ifndef  NON_POSIX
	else {
		headclean();	/* delete header file */
		pfclean();	/* clean up persist files */
	}
#endif
	exit(code);
}


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


static void
trace_sources(void)			/* trace rays to light sources, also */
{
	int	sn;

	for (sn = 0; sn < nsources; sn++)
		source[sn].sflags |= SFOLLOW;
}


void
rtrace(				/* trace rays from file */
	char  *fname,
	int  nproc
)
{
	unsigned long  vcount = (hresolu > 1) ? (unsigned long)hresolu*vresolu
					      : (unsigned long)vresolu;
	long  nextflush = (!vresolu | (hresolu <= 1)) * hresolu;
	int  something2flush = 0;
	FILE  *fp;
	double	d;
	FVECT  orig, direc;
					/* set up input */
	if (fname == NULL)
		inpfp = stdin;
	else if ((inpfp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open input file \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
#ifdef getc_unlocked
	flockfile(inpfp);		/* avoid lock/unlock overhead */
	flockfile(stdout);
#endif
	if (inform != 'a')
		SET_FILE_BINARY(inpfp);
					/* set up output */
	if (castonly || every_out[0] != NULL)
		nproc = 1;		/* don't bother multiprocessing */
	if (Tflag && every_out[0] != NULL)
		trace_sources();	/* asking to trace light sources */
	if ((nextflush > 0) & (nproc > nextflush)) {
		error(WARNING, "reducing number of processes to match flush interval");
		nproc = nextflush;
	}
	switch (outform) {
	case 'a': putreal = puta; break;
	case 'f': putreal = putf; break;
	case 'd': putreal = putd; break;
	case 'c': 
		if (outvals[1] || !strchr("vrx", outvals[0]))
			error(USER, "color format only with -ov, -or, -ox");
		putreal = putrgbe; break;
	default:
		error(CONSISTENCY, "botched output format");
	}
	if (nproc > 1) {		/* start multiprocessing */
		ray_popen(nproc);
		ray_fifo_out = printvals;
	}
	if (hresolu > 0) {
		if (vresolu > 0)
			fprtresolu(hresolu, vresolu, stdout);
		else
			fflush(stdout);
	}
					/* process input rays */
	while ((d = nextray(orig, direc)) >= 0.0) {
		if (d == 0.0) {				/* flush request? */
			if (something2flush) {
				if (ray_pnprocs > 1 && ray_fifo_flush() < 0)
					error(USER, "child(ren) died");
				bogusray();
				fflush(stdout);
				nextflush = (!vresolu | (hresolu <= 1)) * hresolu;
				something2flush = 0;
			} else
				bogusray();
		} else {				/* compute and print */
			rtcompute(orig, direc, lim_dist ? d : 0.0);
			if (!--nextflush) {		/* flush if time */
				if (ray_pnprocs > 1 && ray_fifo_flush() < 0)
					error(USER, "child(ren) died");
				fflush(stdout);
				nextflush = hresolu;
			} else
				something2flush = 1;
		}
		if (ferror(stdout))
			error(SYSTEM, "write error");
		if (vcount && !--vcount)		/* check for end */
			break;
	}
	if (ray_pnprocs > 1) {				/* clean up children */
		if (ray_fifo_flush() < 0)
			error(USER, "unable to complete processing");
		ray_pclose(0);
	}
	if (vcount)
		error(WARNING, "unexpected EOF on input");
	if (fflush(stdout) < 0)
		error(SYSTEM, "write error");
	if (fname != NULL) {
		fclose(inpfp);
		inpfp = NULL;
	}
	nextray(NULL, NULL);
}


int
setrtoutput(void)			/* set up output tables, return #comp */
{
	char  *vs = outvals;
	oputf_t **table = ray_out;
	int  ncomp = 0;

	if (!*vs)
		error(USER, "empty output specification");

	castonly = 1;			/* sets castonly as side-effect */
	do
		switch (*vs) {
		case 'T':				/* trace sources */
			Tflag++;
			/* fall through */
		case 't':				/* trace */
			if (!vs[1]) break;
			*table = NULL;
			table = every_out;
			trace = ourtrace;
			castonly = 0;
			break;
		case 'o':				/* origin */
			*table++ = oputo;
			ncomp += 3;
			break;
		case 'd':				/* direction */
			*table++ = oputd;
			ncomp += 3;
			break;
		case 'r':				/* reflected contrib. */
			*table++ = oputr;
			ncomp += 3;
			castonly = 0;
			break;
		case 'R':				/* reflected distance */
			*table++ = oputR;
			ncomp++;
			castonly = 0;
			break;
		case 'x':				/* xmit contrib. */
			*table++ = oputx;
			ncomp += 3;
			castonly = 0;
			break;
		case 'X':				/* xmit distance */
			*table++ = oputX;
			ncomp++;
			castonly = 0;
			break;
		case 'v':				/* value */
			*table++ = oputv;
			ncomp += 3;
			castonly = 0;
			break;
		case 'V':				/* contribution */
			*table++ = oputV;
			ncomp += 3;
			castonly = 0;
			if (ambounce > 0 && (ambacc > FTINY || ambssamp > 0))
				error(WARNING,
					"-otV accuracy depends on -aa 0 -as 0");
			break;
		case 'l':				/* effective distance */
			*table++ = oputl;
			ncomp++;
			castonly = 0;
			break;
		case 'c':				/* local coordinates */
			*table++ = oputc;
			ncomp += 2;
			break;
		case 'L':				/* single ray length */
			*table++ = oputL;
			ncomp++;
			break;
		case 'p':				/* point */
			*table++ = oputp;
			ncomp += 3;
			break;
		case 'n':				/* perturbed normal */
			*table++ = oputn;
			ncomp += 3;
			castonly = 0;
			break;
		case 'N':				/* unperturbed normal */
			*table++ = oputN;
			ncomp += 3;
			break;
		case 's':				/* surface */
			*table++ = oputs;
			ncomp++;
			break;
		case 'w':				/* weight */
			*table++ = oputw;
			ncomp++;
			break;
		case 'W':				/* coefficient */
			*table++ = oputW;
			ncomp += 3;
			castonly = 0;
			if (ambounce > 0 && (ambacc > FTINY) | (ambssamp > 0))
				error(WARNING,
					"-otW accuracy depends on -aa 0 -as 0");
			break;
		case 'm':				/* modifier */
			*table++ = oputm;
			ncomp++;
			break;
		case 'M':				/* material */
			*table++ = oputM;
			ncomp++;
			break;
		case '~':				/* tilde */
			*table++ = oputtilde;
			break;
		default:
			sprintf(errmsg, "unrecognized output option '%c'", *vs);
			error(USER, errmsg);
		}
	while (*++vs);

	*table = NULL;
	if (*every_out != NULL)
		ncomp = 0;
							/* compatibility */
	if ((do_irrad | imm_irrad) && castonly)
		error(USER, "-I+ and -i+ options require some value output");
	for (table = ray_out; *table != NULL; table++) {
		if ((*table == oputV) | (*table == oputW))
			error(WARNING, "-oVW options require trace mode");
		if ((do_irrad | imm_irrad) &&
				(*table == oputr) | (*table == oputR) |
				(*table == oputx) | (*table == oputX))
			error(WARNING, "-orRxX options incompatible with -I+ and -i+");
	}
	return(ncomp);
}


static void
bogusray(void)			/* print out empty record */
{
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	printvals(&thisray);
}


static void
raycast(			/* compute first ray intersection only */
	RAY *r
)
{
	if (!localhit(r, &thescene)) {
		if (r->ro == &Aftplane) {	/* clipped */
			r->ro = NULL;
			r->rot = FHUGE;
		} else
			sourcehit(r);
	}
}


static void
rayirrad(			/* compute irradiance rather than radiance */
	RAY *r
)
{
	void	(*old_revf)(RAY *) = r->revf;
					/* pretend we hit surface */
	r->rxt = r->rot = 1e-5;
	VSUM(r->rop, r->rorg, r->rdir, r->rot);
	r->ron[0] = -r->rdir[0];
	r->ron[1] = -r->rdir[1];
	r->ron[2] = -r->rdir[2];
	r->rod = 1.0;
					/* compute result */
	r->revf = raytrace;
	(*ofun[Lamb.otype].funp)(&Lamb, r);
	r->revf = old_revf;
}


static void
rtcompute(			/* compute and print ray value(s) */
	FVECT  org,
	FVECT  dir,
	double	dmax
)
{
					/* set up ray */
	if (imm_irrad) {
		VSUM(thisray.rorg, org, dir, 1.1e-4);
		thisray.rdir[0] = -dir[0];
		thisray.rdir[1] = -dir[1];
		thisray.rdir[2] = -dir[2];
		thisray.rmax = 0.0;
	} else {
		VCOPY(thisray.rorg, org);
		VCOPY(thisray.rdir, dir);
		thisray.rmax = dmax;
	}
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	if (imm_irrad)
		thisray.revf = rayirrad;
	else if (castonly)
		thisray.revf = raycast;
	if (ray_pnprocs > 1) {		/* multiprocessing FIFO? */
		if (ray_fifo_in(&thisray) < 0)
			error(USER, "lost children");
		return;
	}
	samplendx++;			/* else do it ourselves */
	rayvalue(&thisray);
	printvals(&thisray);
}


static int
printvals(			/* print requested ray values */
	RAY  *r
)
{
	oputf_t **tp;

	if (ray_out[0] == NULL)
		return(0);
	for (tp = ray_out; *tp != NULL; tp++)
		(**tp)(r);
	if (outform == 'a')
		putchar('\n');
	return(1);
}


static int
is_fifo(		/* check if file pointer connected to pipe */
	FILE *fp
)
{
#ifdef S_ISFIFO
	struct stat  sbuf;

	if (fstat(fileno(fp), &sbuf) < 0)
		error(SYSTEM, "fstat() failed on input stream");
	return(S_ISFIFO(sbuf.st_mode));
#else
	return (fp == stdin);		/* just a guess, really */
#endif
}


static int
getvec(			/* get a vector from fp */
	FVECT  vec,
	int  fmt,
	FILE  *fp
)
{
	static float  vf[3];
	static double  vd[3];
	char  buf[32];
	int  i;

	switch (fmt) {
	case 'a':					/* ascii */
		for (i = 0; i < 3; i++) {
			if (fgetword(buf, sizeof(buf), fp) == NULL ||
					!isflt(buf))
				return(-1);
			vec[i] = atof(buf);
		}
		break;
	case 'f':					/* binary float */
		if (getbinary(vf, sizeof(float), 3, fp) != 3)
			return(-1);
		VCOPY(vec, vf);
		break;
	case 'd':					/* binary double */
		if (getbinary(vd, sizeof(double), 3, fp) != 3)
			return(-1);
		VCOPY(vec, vd);
		break;
	default:
		error(CONSISTENCY, "botched input format");
	}
	return(0);
}


static int
iszerovec(const FVECT vec)
{
	return (vec[0] == 0.0) & (vec[1] == 0.0) & (vec[2] == 0.0);
}


static double
nextray(		/* return next ray in work group (-1.0 if EOF) */
	FVECT org,
	FVECT dir
)
{
	const size_t	qlength = !vresolu * hresolu;

	if ((org == NULL) | (dir == NULL)) {
		if (inp_queue != NULL)	/* asking to free queue */
			free(inp_queue);
		inp_queue = NULL;
		inp_qpos = inp_qend = 0;
		return(-1.);
	}
	if (!inp_qend) {		/* initialize FIFO queue */
		int	rsiz = 6*20;	/* conservative ascii ray size */
		if (inform == 'f') rsiz = 6*sizeof(float);
		else if (inform == 'd') rsiz = 6*sizeof(double);
					/* check against pipe limit */
		if (qlength*rsiz > 512 && is_fifo(inpfp))
			inp_queue = (FVECT *)malloc(sizeof(FVECT)*2*qlength);
		inp_qend = -(inp_queue == NULL);	/* flag for no queue */
	}
	if (inp_qend < 0) {		/* not queuing? */
		if (getvec(org, inform, inpfp) < 0 ||
				getvec(dir, inform, inpfp) < 0)
			return(-1.);
		return normalize(dir);
	}
	if (inp_qpos >= inp_qend) {	/* need to refill input queue? */
		for (inp_qend = 0; inp_qend < qlength; inp_qend++) {
			if (getvec(inp_queue[2*inp_qend], inform, inpfp) < 0
					|| getvec(inp_queue[2*inp_qend+1],
							inform, inpfp) < 0)
				break;		/* hit EOF */
			if (iszerovec(inp_queue[2*inp_qend+1])) {
				++inp_qend;	/* flush request */
				break;
			}
		}
		inp_qpos = 0;
	}
	if (inp_qpos >= inp_qend)	/* unexpected EOF? */
		return(-1.);
	VCOPY(org, inp_queue[2*inp_qpos]);
	VCOPY(dir, inp_queue[2*inp_qpos+1]);
	++inp_qpos;
	return normalize(dir);
}


void
tranotify(			/* record new modifier */
	OBJECT	obj
)
{
	static int  hitlimit = 0;
	OBJREC	 *o = objptr(obj);
	char  **tralp;

	if (obj == OVOID) {		/* starting over */
		traset[0] = 0;
		hitlimit = 0;
		return;
	}
	if (hitlimit || !ismodifier(o->otype))
		return;
	for (tralp = tralist; *tralp != NULL; tralp++)
		if (!strcmp(o->oname, *tralp)) {
			if (traset[0] >= MAXTSET) {
				error(WARNING, "too many modifiers in trace list");
				hitlimit++;
				return;		/* should this be fatal? */
			}
			insertelem(traset, obj);
			return;
		}
}


static void
ourtrace(				/* print ray values */
	RAY  *r
)
{
	oputf_t **tp;

	if (every_out[0] == NULL)
		return;
	if (r->ro == NULL) {
		if (traincl == 1)
			return;
	} else if (traincl != -1 && traincl != inset(traset, r->ro->omod))
		return;
	tabin(r);
	for (tp = every_out; *tp != NULL; tp++)
		(**tp)(r);
	if (outform == 'a')
		putchar('\n');
}


static void
tabin(				/* tab in appropriate amount */
	RAY  *r
)
{
	const RAY  *rp;

	for (rp = r->parent; rp != NULL; rp = rp->parent)
		putchar('\t');
}


static void
oputo(				/* print origin */
	RAY  *r
)
{
	(*putreal)(r->rorg, 3);
}


static void
oputd(				/* print direction */
	RAY  *r
)
{
	(*putreal)(r->rdir, 3);
}


static void
oputr(				/* print mirrored contribution */
	RAY  *r
)
{
	RREAL	cval[3];

	cval[0] = colval(r->mcol,RED);
	cval[1] = colval(r->mcol,GRN);
	cval[2] = colval(r->mcol,BLU);
	(*putreal)(cval, 3);
}



static void
oputR(				/* print mirrored distance */
	RAY  *r
)
{
	(*putreal)(&r->rmt, 1);
}


static void
oputx(				/* print unmirrored contribution */
	RAY  *r
)
{
	RREAL	cval[3];

	cval[0] = colval(r->rcol,RED) - colval(r->mcol,RED);
	cval[1] = colval(r->rcol,GRN) - colval(r->mcol,GRN);
	cval[2] = colval(r->rcol,BLU) - colval(r->mcol,BLU);
	(*putreal)(cval, 3);
}


static void
oputX(				/* print unmirrored distance */
	RAY  *r
)
{
	(*putreal)(&r->rxt, 1);
}


static void
oputv(				/* print value */
	RAY  *r
)
{
	RREAL	cval[3];

	cval[0] = colval(r->rcol,RED);
	cval[1] = colval(r->rcol,GRN);
	cval[2] = colval(r->rcol,BLU);
	(*putreal)(cval, 3);
}


static void
oputV(				/* print value contribution */
	RAY *r
)
{
	RREAL	contr[3];

	raycontrib(contr, r, PRIMARY);
	multcolor(contr, r->rcol);
	(*putreal)(contr, 3);
}


static void
oputl(				/* print effective distance */
	RAY  *r
)
{
	RREAL	d = raydistance(r);

	(*putreal)(&d, 1);
}


static void
oputL(				/* print single ray length */
	RAY  *r
)
{
	(*putreal)(&r->rot, 1);
}


static void
oputc(				/* print local coordinates */
	RAY  *r
)
{
	(*putreal)(r->uv, 2);
}


static RREAL	vdummy[3] = {0.0, 0.0, 0.0};


static void
oputp(				/* print intersection point */
	RAY  *r
)
{
	(*putreal)(r->rop, 3);	/* set to ray origin if distant or no hit */
}


static void
oputN(				/* print unperturbed normal */
	RAY  *r
)
{
	if (r->ro == NULL) {	/* zero vector if clipped or no hit */
		(*putreal)(vdummy, 3);
		return;
	}
	if (r->rflips & 1) {	/* undo any flippin' flips */
		FVECT	unrm;
		unrm[0] = -r->ron[0];
		unrm[1] = -r->ron[1];
		unrm[2] = -r->ron[2];
		(*putreal)(unrm, 3);
	} else
		(*putreal)(r->ron, 3);
}


static void
oputn(				/* print perturbed normal */
	RAY  *r
)
{
	FVECT  pnorm;

	if (r->ro == NULL) {	/* clipped or no hit */
		(*putreal)(vdummy, 3);
		return;
	}
	raynormal(pnorm, r);
	(*putreal)(pnorm, 3);
}


static void
oputs(				/* print name */
	RAY  *r
)
{
	if (r->ro != NULL)
		fputs(r->ro->oname, stdout);
	else
		putchar('*');
	putchar('\t');
}


static void
oputw(				/* print weight */
	RAY  *r
)
{
	RREAL	rwt = r->rweight;
	
	(*putreal)(&rwt, 1);
}


static void
oputW(				/* print coefficient */
	RAY  *r
)
{
	RREAL	contr[3];
				/* shadow ray not on source? */
	if (r->rsrc >= 0 && source[r->rsrc].so != r->ro)
		setcolor(contr, 0.0, 0.0, 0.0);
	else
		raycontrib(contr, r, PRIMARY);

	(*putreal)(contr, 3);
}


static void
oputm(				/* print modifier */
	RAY  *r
)
{
	if (r->ro != NULL)
		if (r->ro->omod != OVOID)
			fputs(objptr(r->ro->omod)->oname, stdout);
		else
			fputs(VOIDID, stdout);
	else
		putchar('*');
	putchar('\t');
}


static void
oputM(				/* print material */
	RAY  *r
)
{
	OBJREC	*mat;

	if (r->ro != NULL) {
		if ((mat = findmaterial(r->ro)) != NULL)
			fputs(mat->oname, stdout);
		else
			fputs(VOIDID, stdout);
	} else
		putchar('*');
	putchar('\t');
}


static void
oputtilde(			/* output tilde (spacer) */
	RAY  *r
)
{
	fputs("~\t", stdout);
}


static void
puta(				/* print ascii value(s) */
	RREAL *v, int n
)
{
	if (n == 3) {
		printf("%e\t%e\t%e\t", v[0], v[1], v[2]);
		return;
	}
	while (n--)
		printf("%e\t", *v++);
}


static void
putd(RREAL *v, int n)		/* output binary double(s) */
{
#ifdef	SMLFLT
	double	da[3];
	int	i;

	if (n > 3)
		error(INTERNAL, "code error in putd()");
	for (i = n; i--; )
		da[i] = v[i];
	putbinary(da, sizeof(double), n, stdout);
#else
	putbinary(v, sizeof(RREAL), n, stdout);
#endif
}


static void
putf(RREAL *v, int n)		/* output binary float(s) */
{
#ifndef	SMLFLT
	float	fa[3];
	int	i;

	if (n > 3)
		error(INTERNAL, "code error in putf()");
	for (i = n; i--; )
		fa[i] = v[i];
	putbinary(fa, sizeof(float), n, stdout);
#else
	putbinary(v, sizeof(RREAL), n, stdout);
#endif
}


static void
putrgbe(RREAL *v, int n)	/* output RGBE color */
{
	COLR  cout;

	if (n != 3)
		error(INTERNAL, "putrgbe() not called with 3 components");
	setcolr(cout, v[0], v[1], v[2]);
	putbinary(cout, sizeof(cout), 1, stdout);
}
