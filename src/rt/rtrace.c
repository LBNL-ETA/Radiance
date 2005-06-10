#ifndef lint
static const char	RCSid[] = "$Id: rtrace.c,v 2.45 2005/06/10 16:49:42 greg Exp $";
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
#include  "resolu.h"

CUBE  thescene;				/* our scene */
OBJECT	nsceneobjs;			/* number of objects in our scene */

int  dimlist[MAXDIM];			/* sampling dimensions */
int  ndims = 0;				/* number of sampling dimensions */
int  samplendx = 0;			/* index for this sample */

int  imm_irrad = 0;			/* compute immediate irradiance? */
int  lim_dist = 0;			/* limit distance? */

int  inform = 'a';			/* input format */
int  outform = 'a';			/* output format */
char  *outvals = "v";			/* output specification */

int  do_irrad = 0;			/* compute irradiance? */

void  (*trace)() = NULL;		/* trace call */

char  *tralist[128];			/* list of modifers to trace (or no) */
int  traincl = -1;			/* include == 1, exclude == 0 */
#define	 MAXTSET	511		/* maximum number in trace set */
OBJECT	traset[MAXTSET+1]={0};		/* trace include/exclude set */

int  hresolu = 0;			/* horizontal (scan) size */
int  vresolu = 0;			/* vertical resolution */

double  dstrsrc = 0.0;			/* square source distribution */
double  shadthresh = .03;		/* shadow threshold */
double  shadcert = .75;			/* shadow certainty */
int  directrelay = 2;			/* number of source relays */
int  vspretest = 512;			/* virtual source pretest density */
int  directvis = 1;			/* sources visible? */
double  srcsizerat = .2;		/* maximum ratio source size/dist. */

COLOR  cextinction = BLKCOLOR;		/* global extinction coefficient */
COLOR  salbedo = BLKCOLOR;		/* global scattering albedo */
double  seccg = 0.;			/* global scattering eccentricity */
double  ssampdist = 0.;			/* scatter sampling distance */

double  specthresh = .15;		/* specular sampling threshold */
double  specjitter = 1.;		/* specular sampling jitter */

int  backvis = 1;			/* back face visibility */

int  maxdepth = -10;			/* maximum recursion depth */
double  minweight = 2e-3;		/* minimum ray weight */

char  *ambfile = NULL;			/* ambient file name */
COLOR  ambval = BLKCOLOR;		/* ambient value */
int  ambvwt = 0;			/* initial weight for ambient value */
double  ambacc = 0.15;			/* ambient accuracy */
int  ambres = 256;			/* ambient resolution */
int  ambdiv = 1024;			/* ambient divisions */
int  ambssamp = 512;			/* ambient super-samples */
int  ambounce = 0;			/* ambient bounces */
char  *amblist[AMBLLEN];		/* ambient include/exclude list */
int  ambincl = -1;			/* include == 1, exclude == 0 */

static int  castonly = 0;

static RAY  thisray;			/* for our convenience */

typedef void putf_t(double v);
static putf_t puta, putd, putf;

typedef void oputf_t(RAY *r);
static oputf_t  oputo, oputd, oputv, oputl, oputL, oputc, oputp,
		oputn, oputN, oputs, oputw, oputW, oputm, oputM, oputtilde;

static void setoutput(char *vs);
static void tranotify(OBJECT obj);
static void bogusray(void);
static void rad(FVECT  org, FVECT  dir, double	dmax);
static void irrad(FVECT  org, FVECT  dir);
static void printvals(RAY  *r);
static int getvec(FVECT  vec, int  fmt, FILE  *fp);
static void tabin(RAY  *r);
static void ourtrace(RAY  *r);

static oputf_t *ray_out[16], *every_out[16];
static putf_t *putreal;

void  (*addobjnotify[])() = {ambnotify, tranotify, NULL};


void
quit(			/* quit program */
	int  code
)
{
#ifndef  NON_POSIX /* XXX we don't clean up elsewhere? */
	headclean();		/* delete header file */
	pfclean();		/* clean up persist files */
#endif
	exit(code);
}


extern char *
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


extern void
rtrace(				/* trace rays from file */
	char  *fname
)
{
	unsigned long  vcount = (hresolu > 1) ? (unsigned long)hresolu*vresolu
					      : vresolu;
	long  nextflush = hresolu;
	FILE  *fp;
	double	d;
	FVECT  orig, direc;
					/* set up input */
	if (fname == NULL)
		fp = stdin;
	else if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open input file \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
	if (inform != 'a')
		SET_FILE_BINARY(fp);
					/* set up output */
	setoutput(outvals);
	switch (outform) {
	case 'a': putreal = puta; break;
	case 'f': putreal = putf; break;
	case 'd': putreal = putd; break;
	case 'c': 
		if (strcmp(outvals, "v"))
			error(USER, "color format with value output only");
		break;
	default:
		error(CONSISTENCY, "botched output format");
	}
	if (hresolu > 0) {
		if (vresolu > 0)
			fprtresolu(hresolu, vresolu, stdout);
		fflush(stdout);
	}
					/* process file */
	while (getvec(orig, inform, fp) == 0 &&
			getvec(direc, inform, fp) == 0) {

		d = normalize(direc);
		if (d == 0.0) {				/* zero ==> flush */
			bogusray();
			if (--nextflush <= 0 || !vcount) {
				fflush(stdout);
				nextflush = hresolu;
			}
		} else {
			samplendx++;
							/* compute and print */
			if (imm_irrad)
				irrad(orig, direc);
			else
				rad(orig, direc, lim_dist ? d : 0.0);
							/* flush if time */
			if (!--nextflush) {
				fflush(stdout);
				nextflush = hresolu;
			}
		}
		if (ferror(stdout))
			error(SYSTEM, "write error");
		if (vcount && !--vcount)		/* check for end */
			break;
	}
	fflush(stdout);
	if (vcount)
		error(USER, "unexpected EOF on input");
	if (fname != NULL)
		fclose(fp);
}


static void
trace_sources(void)			/* trace rays to light sources, also */
{
	int	sn;
	
	for (sn = 0; sn < nsources; sn++)
		source[sn].sflags |= SFOLLOW;
}


static void
setoutput(				/* set up output tables */
	register char  *vs
)
{
	register oputf_t **table = ray_out;

	castonly = 1;
	while (*vs)
		switch (*vs++) {
		case 'T':				/* trace sources */
			if (!*vs) break;
			trace_sources();
			/* fall through */
		case 't':				/* trace */
			if (!*vs) break;
			*table = NULL;
			table = every_out;
			trace = ourtrace;
			castonly = 0;
			break;
		case 'o':				/* origin */
			*table++ = oputo;
			break;
		case 'd':				/* direction */
			*table++ = oputd;
			break;
		case 'v':				/* value */
			*table++ = oputv;
			castonly = 0;
			break;
		case 'l':				/* effective distance */
			*table++ = oputl;
			castonly = 0;
			break;
		case 'c':				/* local coordinates */
			*table++ = oputc;
			break;
		case 'L':				/* single ray length */
			*table++ = oputL;
			break;
		case 'p':				/* point */
			*table++ = oputp;
			break;
		case 'n':				/* perturbed normal */
			*table++ = oputn;
			castonly = 0;
			break;
		case 'N':				/* unperturbed normal */
			*table++ = oputN;
			break;
		case 's':				/* surface */
			*table++ = oputs;
			break;
		case 'w':				/* weight */
			*table++ = oputw;
			break;
		case 'W':				/* coefficient */
			*table++ = oputW;
			if (ambounce > 0 && (ambacc > FTINY || ambssamp > 0))
				error(WARNING,
					"-otW accuracy depends on -aa 0 -as 0");
			break;
		case 'm':				/* modifier */
			*table++ = oputm;
			break;
		case 'M':				/* material */
			*table++ = oputM;
			break;
		case '~':				/* tilde */
			*table++ = oputtilde;
			break;
		}
	*table = NULL;
}


static void
bogusray(void)			/* print out empty record */
{
	thisray.rorg[0] = thisray.rorg[1] = thisray.rorg[2] =
	thisray.rdir[0] = thisray.rdir[1] = thisray.rdir[2] = 0.0;
	thisray.rmax = 0.0;
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	printvals(&thisray);
}


static void
rad(		/* compute and print ray value(s) */
	FVECT  org,
	FVECT  dir,
	double	dmax
)
{
	VCOPY(thisray.rorg, org);
	VCOPY(thisray.rdir, dir);
	thisray.rmax = dmax;
	rayorigin(&thisray, PRIMARY, NULL, NULL);
	if (castonly) {
		if (!localhit(&thisray, &thescene)) {
			if (thisray.ro == &Aftplane) {	/* clipped */
				thisray.ro = NULL;
				thisray.rot = FHUGE;
			} else
				sourcehit(&thisray);
		}
	} else
		rayvalue(&thisray);
	printvals(&thisray);
}


static void
irrad(			/* compute immediate irradiance value */
	FVECT  org,
	FVECT  dir
)
{
	register int  i;

	for (i = 0; i < 3; i++) {
		thisray.rorg[i] = org[i] + dir[i];
		thisray.rdir[i] = -dir[i];
	}
	thisray.rmax = 0.0;
	rayorigin(&thisray, PRIMARY, NULL, NULL);
					/* pretend we hit surface */
	thisray.rot = 1.0-1e-4;
	thisray.rod = 1.0;
	VCOPY(thisray.ron, dir);
	for (i = 0; i < 3; i++)		/* fudge factor */
		thisray.rop[i] = org[i] + 1e-4*dir[i];
					/* compute and print */
	(*ofun[Lamb.otype].funp)(&Lamb, &thisray);
	printvals(&thisray);
}


static void
printvals(			/* print requested ray values */
	RAY  *r
)
{
	register oputf_t **tp;

	if (ray_out[0] == NULL)
		return;
	for (tp = ray_out; *tp != NULL; tp++)
		(**tp)(r);
	if (outform == 'a')
		putchar('\n');
}


static int
getvec(		/* get a vector from fp */
	register FVECT  vec,
	int  fmt,
	FILE  *fp
)
{
	static float  vf[3];
	static double  vd[3];
	char  buf[32];
	register int  i;

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
		if (fread((char *)vf, sizeof(float), 3, fp) != 3)
			return(-1);
		vec[0] = vf[0]; vec[1] = vf[1]; vec[2] = vf[2];
		break;
	case 'd':					/* binary double */
		if (fread((char *)vd, sizeof(double), 3, fp) != 3)
			return(-1);
		vec[0] = vd[0]; vec[1] = vd[1]; vec[2] = vd[2];
		break;
	default:
		error(CONSISTENCY, "botched input format");
	}
	return(0);
}


static void
tranotify(			/* record new modifier */
	OBJECT	obj
)
{
	static int  hitlimit = 0;
	register OBJREC	 *o = objptr(obj);
	register char  **tralp;

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
	register oputf_t **tp;

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
	(*putreal)(r->rorg[0]);
	(*putreal)(r->rorg[1]);
	(*putreal)(r->rorg[2]);
}


static void
oputd(				/* print direction */
	RAY  *r
)
{
	(*putreal)(r->rdir[0]);
	(*putreal)(r->rdir[1]);
	(*putreal)(r->rdir[2]);
}


static void
oputv(				/* print value */
	RAY  *r
)
{
	if (outform == 'c') {
		COLR  cout;
		setcolr(cout,	colval(r->rcol,RED),
				colval(r->rcol,GRN),
				colval(r->rcol,BLU));
		fwrite((char *)cout, sizeof(cout), 1, stdout);
		return;
	}
	(*putreal)(colval(r->rcol,RED));
	(*putreal)(colval(r->rcol,GRN));
	(*putreal)(colval(r->rcol,BLU));
}


static void
oputl(				/* print effective distance */
	RAY  *r
)
{
	(*putreal)(r->rt);
}


static void
oputL(				/* print single ray length */
	RAY  *r
)
{
	(*putreal)(r->rot);
}


static void
oputc(				/* print local coordinates */
	RAY  *r
)
{
	(*putreal)(r->uv[0]);
	(*putreal)(r->uv[1]);
}


static void
oputp(				/* print point */
	RAY  *r
)
{
	if (r->rot < FHUGE) {
		(*putreal)(r->rop[0]);
		(*putreal)(r->rop[1]);
		(*putreal)(r->rop[2]);
	} else {
		(*putreal)(0.0);
		(*putreal)(0.0);
		(*putreal)(0.0);
	}
}


static void
oputN(				/* print unperturbed normal */
	RAY  *r
)
{
	if (r->rot < FHUGE) {
		(*putreal)(r->ron[0]);
		(*putreal)(r->ron[1]);
		(*putreal)(r->ron[2]);
	} else {
		(*putreal)(0.0);
		(*putreal)(0.0);
		(*putreal)(0.0);
	}
}


static void
oputn(				/* print perturbed normal */
	RAY  *r
)
{
	FVECT  pnorm;

	if (r->rot >= FHUGE) {
		(*putreal)(0.0);
		(*putreal)(0.0);
		(*putreal)(0.0);
		return;
	}
	raynormal(pnorm, r);
	(*putreal)(pnorm[0]);
	(*putreal)(pnorm[1]);
	(*putreal)(pnorm[2]);
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
	(*putreal)(r->rweight);
}


static void
oputW(				/* print contribution */
	RAY  *r
)
{
	COLOR	contr;

	raycontrib(contr, r, PRIMARY);
	(*putreal)(colval(contr,RED));
	(*putreal)(colval(contr,GRN));
	(*putreal)(colval(contr,BLU));
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
puta(				/* print ascii value */
	double  v
)
{
	printf("%e\t", v);
}


static void
putd(v)				/* print binary double */
double  v;
{
	fwrite((char *)&v, sizeof(v), 1, stdout);
}


static void
putf(v)				/* print binary float */
double  v;
{
	float f = v;

	fwrite((char *)&f, sizeof(f), 1, stdout);
}
