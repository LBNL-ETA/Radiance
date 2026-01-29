#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Calculate flux transfer matrix or matrices using RcontribSimulManager class
 *
 * Much of this is cribbed from util/rfluxmtx.c
 * Implementation here improves efficiency and overcomes receiver limits
 */

#include "copyright.h"

#include <ctype.h>
#include <signal.h>
#include <time.h>
#include "RcontribSimulManager.h"
#include "bsdf.h"
#include "bsdf_m.h"
#include "func.h"
#include "random.h"
#include "triangulate.h"
#include "view.h"

int		verby = 0;		/* verbose mode? */

const char	*reinhfn = "reinhartb.cal";
const char	*shirchiufn = "disk2square.cal";
const char	*kfullfn = "klems_full.cal";
const char	*khalffn = "klems_half.cal";
const char	*kquarterfn = "klems_quarter.cal";
const char	*ciefn = "cieskyscan.cal";

const char	*sigerr[NSIG];		/* signal error messages */

int	inpfmt = 'a';			/* input format */
int	outfmt = 'f';			/* output format */

int	report_intvl = 0;		/* reporting interval (seconds) */

RcontribSimulManager	myRCmanager;	// global rcontrib simulation manager

#define	PARAMSNAME	"rfluxmtx"	/* string indicating parameters */

				/* surface type IDs */
#define ST_NONE		0
#define ST_POLY		1
#define ST_RING		2
#define ST_SOURCE	3
				/* surface structure */
struct SURF {
	SURF		*next;		/* next surface in list */
	void		*priv;		/* private data (malloc'ed) */
	char		sname[32];	/* surface name */
	FVECT		snrm;		/* surface normal */
	double		area;		/* surface area / proj. solid angle */
	short		styp;		/* surface type */
	short		nfargs;		/* number of real arguments */
	double		farg[1];	/* real values (extends struct) */
};
			/* triangle */
struct PTRI {
	double	afrac;			/* fraction of total area */
	short	vndx[3];		/* vertex indices */
};
			/* triangulated polygon */
struct POLYTRIS {
	FVECT	uva[2];			/* tangent axes */
	int	ntris;			/* number of triangles */
	PTRI	tri[1];			/* triangle array (extends struct) */
};
			/* sender/receiver parameters */
struct PARAMS {
	char		sign;		/* '-' for axis reversal */
	char		hemis[31];	/* hemispherical sampling spec. */
	int		hsiz;		/* hemisphere basis size */
	int		nsurfs;		/* number of surfaces */
	SURF		*slist;		/* list of surfaces */
	FVECT		vup;		/* up vector (zero if unset) */
	FVECT		nrm;		/* average normal direction */
	FVECT		udir, vdir;	/* tangent axes */
	char		outfn[256];	/* output file name (receiver) */
	int		(*sample_basis)(PARAMS *p, int b, FVECT orig_dir[]);
};

static PARAMS		curparams;
static char		curmod[MAXSTR];
static char		newparams[1024];

typedef int	SURFSAMP(FVECT, SURF *, double);

static SURFSAMP	ssamp_bad, ssamp_poly, ssamp_ring;

static SURFSAMP	*orig_in_surf[4] = {
		ssamp_bad, ssamp_poly, ssamp_ring, ssamp_bad
	};

VIEW		ourview = STDVIEW;	// view parameters (if set)
double		dstrpix = 0;		// pixel jitter
double		dblur = 0;		// depth-of-field

/* Clear parameter set */
static void
clear_params(PARAMS *p, bool reset_only = true)
{
	while (p->slist != NULL) {
		SURF	*sdel = p->slist;
		p->slist = sdel->next;
		if (sdel->priv != NULL)
			free(sdel->priv);
		free(sdel);
	}
	if (reset_only) {
		p->nsurfs = 0;
		p->slist = NULL;
		memset(p->nrm, 0, sizeof(FVECT));
		memset(p->vup, 0, sizeof(FVECT));
	} else
		memset(p, 0, sizeof(PARAMS));
}

/* Get surface type from name */
static int
surf_type(const char *otype)
{
	if (!strcmp(otype, "polygon"))
		return(ST_POLY);
	if (!strcmp(otype, "ring"))
		return(ST_RING);
	if (!strcmp(otype, "source"))
		return(ST_SOURCE);
	return(ST_NONE);
}

/* Add arguments to oconv command */
static char *
oconv_command(int ac, char *av[])
{
	static char	oconvbuf[4096] = "!oconv -f ";
	char		*cp = oconvbuf + 10;
	char		*recv = *av++;
	
	if (ac-- <= 0)
		return(NULL);
	if (erract[WARNING].pf == NULL) {	/* warnings off? */
		strcpy(cp, "-w ");
		cp += 3;
	}
	while (ac-- > 0) {	/* copy each argument */
		int	len = strlen(*av);
		if (cp+len+4 >= oconvbuf+sizeof(oconvbuf))
			goto overrun;
		if (matchany(*av, SPECIALS)) {
			*cp++ = QUOTCHAR;
			strcpy(cp, *av++);
			cp += len;
			*cp++ = QUOTCHAR;
		} else {
			strcpy(cp, *av++);
			cp += len;
		}
		*cp++ = ' ';
	}
				/* receiver goes last */
	if (matchany(recv, SPECIALS)) {
		*cp++ = QUOTCHAR;
		while (*recv) {
			if (cp >= oconvbuf+(sizeof(oconvbuf)-3))
				goto overrun;
			*cp++ = *recv++;
		}
		*cp++ = QUOTCHAR;
		*cp = '\0';
	} else
		strcpy(cp, recv);
	return(oconvbuf);
overrun:
	error(USER, "too many file arguments!");
	return(NULL);	/* pro forma return */
}

/* Get normalized direction vector from string specification */
static int
get_direction(FVECT dv, const char *s)
{
	int	sign = 1;
	int	axis = 0;

	memset(dv, 0, sizeof(FVECT));
nextchar:
	switch (*s) {
	case '+':
		++s;
		goto nextchar;
	case '-':
		sign = -sign;
		++s;
		goto nextchar;
	case 'z':
	case 'Z':
		++axis;
	case 'y':
	case 'Y':
		++axis;
	case 'x':
	case 'X':
		dv[axis] = sign;
		return(!s[1] | isspace(s[1]));
	default:
		break;
	}
#ifdef SMLFLT
	if (sscanf(s, "%f,%f,%f", &dv[0], &dv[1], &dv[2]) != 3)
#else
	if (sscanf(s, "%lf,%lf,%lf", &dv[0], &dv[1], &dv[2]) != 3)
#endif
		return(0);
	dv[0] *= (RREAL)sign;
	return(normalize(dv) > 0);
}

/* Parse program parameters (directives) */
static int
parse_params(PARAMS *p, char *pargs)
{
	char	*cp = pargs;
	int	nparams = 0;
	int	quot;
	int	i;

	for ( ; ; ) {
		switch (*cp++) {
		case 'h':
			if (*cp++ != '=')
				break;
			if ((*cp == '+') | (*cp == '-'))
				p->sign = *cp++;
			else
				p->sign = '+';
			p->hsiz = 0;
			i = 0;
			while (*cp && !isspace(*cp)) {
				if (isdigit(*cp))
					p->hsiz = 10*p->hsiz + *cp - '0';
				p->hemis[i++] = *cp++;
			}
			if (!i)
				break;
			p->hemis[i] = '\0';
			p->hsiz += !p->hsiz;
			++nparams;
			continue;
		case 'u':
			if (*cp++ != '=')
				break;
			if (!get_direction(p->vup, cp))
				break;
			while (*cp && !isspace(*cp++))
				;
			++nparams;
			continue;
		case 'o':
			if (*cp++ != '=')
				break;
			quot = 0;
			if ((*cp == '"') | (*cp == '\''))
				quot = *cp++;
			i = 0;
			while (*cp && (quot ? (*cp != quot) : !isspace(*cp))) {
				i++; cp++;
			}
			if (!i)
				break;
			if (!*cp) {
				if (quot)
					break;
				cp[1] = '\0';
			}
			*cp = '\0';
			strlcpy(p->outfn, cp-i, sizeof(p->outfn));
			*cp++ = quot ? quot : ' ';
			++nparams;
			continue;
		case ' ':
		case '\t':
		case '\r':
		case '\n':
			continue;
		case '\0':
			return(nparams);
		default:
			break;
		}
		break;
	}
	sprintf(errmsg, "bad parameter string:%s", pargs);
	error(USER, errmsg);
	return(-1);	/* pro forma return */
}

/* Add receiver modifier and associated parameters */
static void
finish_receiver()
{
	bool		uniform = false;
	const char	*calfn = NULL;
	char		binv[64] = "0";
	char		params[128] = "";
	const char	*binf = NULL;
	const char	*nbins = NULL;
					/* check arguments */
	if (!curmod[0])
		error(USER, "missing receiver surface");
	if (!curparams.outfn[0])
		error(USER, "missing output file for receiver");
	if (!curparams.hemis[0])
		error(USER, "missing hemisphere sampling type");
	if (normalize(curparams.nrm) == 0)
		error(USER, "undefined normal for hemisphere sampling");
	if (normalize(curparams.vup) == 0) {
		if (fabs(curparams.nrm[2]) < .7)
			curparams.vup[2] = 1;
		else
			curparams.vup[1] = 1;
	}
					/* determine sample type/bin */
	if ((tolower(curparams.hemis[0]) == 'u') | (curparams.hemis[0] == '1')) {
		if (curparams.slist->styp != ST_SOURCE)
			sprintf(binv, "if(-Dx*%g-Dy*%g-Dz*%g,0,-1)",
				curparams.nrm[0], curparams.nrm[1], curparams.nrm[2]);
		uniform = true;		/* uniform sampling -- one bin */
	} else if (tolower(curparams.hemis[0]) == 's' &&
				tolower(curparams.hemis[1]) == 'c') {
					/* assign parameters */
		if (curparams.hsiz <= 1)
			error(USER, "missing size for Shirley-Chiu sampling");
		calfn = shirchiufn; shirchiufn = NULL;
		sprintf(params, "SCdim=%d,rNx=%g,rNy=%g,rNz=%g,Ux=%g,Uy=%g,Uz=%g,RHS=%c1",
				curparams.hsiz,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2],
			curparams.sign);
		strcpy(binv, "scbin");
		nbins = "SCdim*SCdim";
	} else if ((tolower(curparams.hemis[0]) == 'r') |
			(tolower(curparams.hemis[0]) == 't')) {
		calfn = reinhfn; reinhfn = NULL;
		sprintf(params, "MF=%d,rNx=%g,rNy=%g,rNz=%g,Ux=%g,Uy=%g,Uz=%g,RHS=%c1",
				curparams.hsiz,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2],
			curparams.sign);
		strcpy(binv, "rbin");
		nbins = "Nrbins";
	} else if (tolower(curparams.hemis[0]) == 'k' &&
			!curparams.hemis[1] |
			(tolower(curparams.hemis[1]) == 'f') |
			(curparams.hemis[1] == '1')) {
		calfn = kfullfn; kfullfn = NULL;
		binf = "kbin";
		nbins = "Nkbins";
	} else if (tolower(curparams.hemis[0]) == 'k' &&
			(tolower(curparams.hemis[1]) == 'h') |
			(curparams.hemis[1] == '2')) {
		calfn = khalffn; khalffn = NULL;
		binf = "khbin";
		nbins = "Nkhbins";
	} else if (tolower(curparams.hemis[0]) == 'k' &&
			(tolower(curparams.hemis[1]) == 'q') |
			(curparams.hemis[1] == '4')) {
		calfn = kquarterfn; kquarterfn = NULL;
		binf = "kqbin";
		nbins = "Nkqbins";
	} else if (!strcasecmp(curparams.hemis, "cie")) {
		calfn = ciefn; ciefn = NULL;
		sprintf(params, "rNx=%g,rNy=%g,rNz=%g,Ux=%g,Uy=%g,Uz=%g,RHS=%c1",
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2],
			curparams.sign);
		strcpy(binv, "cbin");
		nbins = "Ncbins";
	} else {
		sprintf(errmsg, "unrecognized hemisphere sampling: h=%s",
				curparams.hemis);
		error(USER, errmsg);
	}
	if (tolower(curparams.hemis[0]) == 'k') {
		sprintf(params, "RHS=%c1", curparams.sign);
	}
	if (!uniform)
		for (SURF *sp = curparams.slist; sp != NULL; sp = sp->next)
			if (sp->styp == ST_SOURCE && fabs(sp->area - PI) > 1e-3) {
				sprintf(errmsg, "source '%s' must be 180-degrees",
						sp->sname);
				error(USER, errmsg);
			}
	if (binf != NULL) {		/* complete bin function? */
		sprintf(binv, "%s(%g,%g,%g,%g,%g,%g)", binf,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2]);
	}
	if (calfn != NULL)		/* load cal file if needed */
		loadfunc(const_cast<char *>(calfn));
	if (params[0])			/* set parameters if any */
		set_eparams(params);
					// add modifier to rcontrib object
	if (!myRCmanager.AddModifier(curmod, curparams.outfn, params, binv,
			nbins==NULL ? 1 : int(eval(const_cast<char *>(nbins)) + .5)))
		error(INTERNAL, "AddModifier() call failed");

	clear_params(&curparams, true);		/* clear for next receiver */
	curmod[0] = '\0';
}

/* Make randomly oriented tangent plane axes for given normal direction */
static void
make_axes(FVECT uva[2], const FVECT nrm)
{
	int	i;

	if (!getperpendicular(uva[0], nrm, 1))
		error(USER, "bad surface normal in make_axes()");
	VCROSS(uva[1], nrm, uva[0]);
}

/* Illegal sender surfaces end up here */
static int
ssamp_bad(FVECT orig, SURF *sp, double x)
{
	sprintf(errmsg, "illegal sender surface '%s'", sp->sname);
	error(INTERNAL, errmsg);
	return(0);
}

/* Generate origin on ring surface from uniform random variable */
static int
ssamp_ring(FVECT orig, SURF *sp, double x)
{
	FVECT	*uva = (FVECT *)sp->priv;
	double	samp2[2];
	double	uv[2];
	int	i;

	if (uva == NULL) {		/* need tangent axes */
		uva = (FVECT *)malloc(sizeof(FVECT)*2);
		if (uva == NULL) {
			error(SYSTEM, "out of memory in ssamp_ring()");
			return(0);
		}
		make_axes(uva, sp->snrm);
		sp->priv = uva;
	}
	SDmultiSamp(samp2, 2, x);
	samp2[0] = sqrt(samp2[0]*sp->area*(1./PI) + sp->farg[6]*sp->farg[6]);
	samp2[1] *= 2.*PI;
	uv[0] = samp2[0]*tcos(samp2[1]);
	uv[1] = samp2[0]*tsin(samp2[1]);
	for (i = 3; i--; )
		orig[i] = sp->farg[i] + uv[0]*uva[0][i] + uv[1]*uva[1][i];
	return(1);
}

/* Add triangle to polygon's list (call-back function) */
static int
add_triangle(const Vert2_list *tp, int a, int b, int c)
{
	POLYTRIS	*ptp = (POLYTRIS *)tp->p;
	PTRI		*trip = ptp->tri + ptp->ntris++;

	trip->vndx[0] = a;
	trip->vndx[1] = b;
	trip->vndx[2] = c;
	return(1);
}

/* Generate origin on polygon surface from uniform random variable */
static int
ssamp_poly(FVECT orig, SURF *sp, double x)
{
	POLYTRIS	*ptp = (POLYTRIS *)sp->priv;
	double		samp2[2];
	double		*v0, *v1, *v2;
	int		i;

	if (ptp == NULL) {		/* need to triangulate */
		ptp = (POLYTRIS *)malloc(sizeof(POLYTRIS) +
				sizeof(PTRI)*(sp->nfargs/3 - 3));
		if (ptp == NULL)
			goto memerr;
		if (sp->nfargs == 3) {	/* simple case */
			ptp->ntris = 1;
			ptp->tri[0].vndx[0] = 0;
			ptp->tri[0].vndx[1] = 1;
			ptp->tri[0].vndx[2] = 2;
			ptp->tri[0].afrac = 1;
		} else {
			Vert2_list	*v2l = polyAlloc(sp->nfargs/3);
			if (v2l == NULL)
				goto memerr;
			make_axes(ptp->uva, sp->snrm);
			for (i = v2l->nv; i--; ) {
				v2l->v[i].mX = DOT(sp->farg+3*i, ptp->uva[0]);
				v2l->v[i].mY = DOT(sp->farg+3*i, ptp->uva[1]);
			}
			ptp->ntris = 0;
			v2l->p = ptp;
			if (!polyTriangulate(v2l, add_triangle)) {
				sprintf(errmsg, "cannot triangulate polygon '%s'",
						sp->sname);
				error(USER, errmsg);
				return(0);
			}
			for (i = ptp->ntris; i--; ) {
				int	a = ptp->tri[i].vndx[0];
				int	b = ptp->tri[i].vndx[1];
				int	c = ptp->tri[i].vndx[2];
				ptp->tri[i].afrac =
					(v2l->v[a].mX*v2l->v[b].mY -
					 v2l->v[b].mX*v2l->v[a].mY +
					 v2l->v[b].mX*v2l->v[c].mY -
					 v2l->v[c].mX*v2l->v[b].mY +
					 v2l->v[c].mX*v2l->v[a].mY -
					 v2l->v[a].mX*v2l->v[c].mY) /
						(2.*sp->area);
			}
			polyFree(v2l);
		}
		sp->priv = ptp;
	}
					/* pick triangle by partial area */
	for (i = 0; i < ptp->ntris-1 && x > ptp->tri[i].afrac; i++)
		x -= ptp->tri[i].afrac;
	SDmultiSamp(samp2, 2, x/ptp->tri[i].afrac);
	samp2[0] *= samp2[1] = sqrt(samp2[1]);
	samp2[1] = 1. - samp2[1];
	v0 = sp->farg + 3*ptp->tri[i].vndx[0];
	v1 = sp->farg + 3*ptp->tri[i].vndx[1];
	v2 = sp->farg + 3*ptp->tri[i].vndx[2];
	for (i = 3; i--; )
		orig[i] = v0[i] + samp2[0]*(v1[i] - v0[i])
				+ samp2[1]*(v2[i] - v0[i]) ;
	return(1);
memerr:
	error(SYSTEM, "out of memory in ssamp_poly");
	return(0);
}

/* Compute sample origin based on projected areas of sender subsurfaces */
static int
sample_origin(PARAMS *p, FVECT orig, const FVECT rdir, double x)
{
	static double	*projsa;
	static int	nall;
	double		tarea = 0;
	int		i;
	SURF		*sp;
					/* special case for lone surface */
	if (p->nsurfs == 1) {
		sp = p->slist;
		if (DOT(sp->snrm, rdir) >= FTINY) {
			sprintf(errmsg, "sample behind sender '%s'", sp->sname);
			error(INTERNAL, errmsg);
			return(0);
		}
		return((*orig_in_surf[sp->styp])(orig, sp, x));
	}
	if (p->nsurfs > nall) {		/* (re)allocate surface area cache */
		if (projsa) free(projsa);
		projsa = (double *)malloc(sizeof(double)*p->nsurfs);
		if (projsa == NULL)
			error(SYSTEM, "out of memory in sample_origin()");
		nall = p->nsurfs;
	}
					/* compute projected areas */
	for (i = 0, sp = p->slist; sp != NULL; i++, sp = sp->next) {
		projsa[i] = -DOT(sp->snrm, rdir) * sp->area;
		tarea += projsa[i] *= (double)(projsa[i] > 0);
	}
	if (tarea < FTINY*FTINY) {	/* wrong side of sender? */
		error(INTERNAL, "sample behind all sender elements!");
		return(0);
	}
	tarea *= x;			/* get surface from list */
	for (i = 0, sp = p->slist; tarea > projsa[i]; sp = sp->next)
		tarea -= projsa[i++];
	return((*orig_in_surf[sp->styp])(orig, sp, tarea/projsa[i]));
}

/* Uniform sample generator */
static int
sample_uniform(PARAMS *p, int b, FVECT orig_dir[])
{
	int	n = myRCmanager.accum;
	double	samp3[3];
	FVECT	duvw;
	int	i;

	if (orig_dir == NULL)		/* just requesting number of bins? */
		return(1);

	while (n--) {			/* stratified hemisphere sampling */
		SDmultiSamp(samp3, 3, (n+frandom())/myRCmanager.accum);
		square2disk(duvw, samp3[1], samp3[2]);
		duvw[2] = -sqrt(1. - duvw[0]*duvw[0] - duvw[1]*duvw[1]);
		for (i = 3; i--; )
			orig_dir[1][i] = duvw[0]*p->udir[i] +
						duvw[1]*p->vdir[i] +
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		orig_dir += 2;
	}
	return(1);
}

/* Shirly-Chiu sample generator */
static int
sample_shirchiu(PARAMS *p, int b, FVECT orig_dir[])
{
	int	n = myRCmanager.accum;
	double	samp3[3];
	FVECT	duvw;
	int	i;

	if (orig_dir == NULL)			/* just requesting number of bins? */
		return(p->hsiz*p->hsiz);

	while (n--) {			/* stratified sampling */
		SDmultiSamp(samp3, 3, (n+frandom())/myRCmanager.accum);
		square2disk(duvw, (b/p->hsiz + samp3[1])/p->hsiz,
				(b%p->hsiz + samp3[2])/p->hsiz);
		duvw[2] = sqrt(1. - duvw[0]*duvw[0] - duvw[1]*duvw[1]);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		orig_dir += 2;
	}
	return(1);
}

/* Reinhart/Tregenza sample generator */
static int
sample_reinhart(PARAMS *p, int b, FVECT orig_dir[])
{
#define T_NALT	7
	static const int	tnaz[T_NALT] = {30, 30, 24, 24, 18, 12, 6};
	const int		RowMax = T_NALT*p->hsiz + 1;
	const double		RAH = (.5*PI)/(RowMax-.5);
#define rnaz(r)			(r >= RowMax-1 ? 1 : p->hsiz*tnaz[r/p->hsiz])
	int			n = myRCmanager.accum;
	int			row, col;
	double			samp3[3];
	double			alt, azi;
	double			duvw[3];
	int			i;

	if (orig_dir == NULL) {		/* just requesting number of bins? */
		n = 0;
		for (row = RowMax; row--; ) n += rnaz(row);
		return(n);
	}
	row = 0;			/* identify row & column */
	col = b;
	while (col >= rnaz(row)) {
		col -= rnaz(row);
		++row;
	}
	while (n--) {			/* stratified sampling */
		SDmultiSamp(samp3, 3, (n+frandom())/myRCmanager.accum);
		if (row >= RowMax-1)	/* avoid crowding at zenith */
			samp3[1] *= samp3[1];
		alt = (row+samp3[1])*RAH;
		azi = (2.*PI)*(col+samp3[2]-.5)/rnaz(row);
		duvw[2] = cos(alt);	/* measured from horizon */
		duvw[0] = tsin(azi)*duvw[2];
		duvw[1] = -tcos(azi)*duvw[2];
		duvw[2] = sqrt(1. - duvw[2]*duvw[2]);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		orig_dir += 2;
	}
	return(1);
#undef rnaz
#undef T_NALT
}

/* Klems sample generator */
static int
sample_klems(PARAMS *p, int b, FVECT orig_dir[])
{
	static const char	bname[4][20] = {
					"LBNL/Klems Full",
					"LBNL/Klems Half",
					"INTERNAL ERROR",
					"LBNL/Klems Quarter"
				};
	static ANGLE_BASIS	*kbasis[4];
	const int		bi = p->hemis[1] - '1';
	int			n = myRCmanager.accum;
	double			samp2[2];
	double			duvw[3];
	int			i;

	if (!kbasis[bi]) {		/* need to get basis, first */
		for (i = 4; i--; )
			if (!strcasecmp(abase_list[i].name, bname[bi])) {
				kbasis[bi] = &abase_list[i];
				break;
			}
		if (i < 0) {
			sprintf(errmsg, "unknown hemisphere basis '%s'",
					bname[bi]);
			error(USER, errmsg);
			return(0);
		}
	}
	if (orig_dir == NULL)		/* just requesting number of bins? */
		return(kbasis[bi]->nangles);

	while (n--) {			/* stratified sampling */
		SDmultiSamp(samp2, 2, (n+frandom())/myRCmanager.accum);
		if (!fo_getvec(duvw, b+samp2[1], kbasis[bi]))
			return(0);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp2[0]))
			return(0);
		orig_dir += 2;
	}
	return(1);
}

/* Prepare hemisphere basis sampler that will send rays to rcontrib */
static int
prepare_sampler(PARAMS *p)
{
	if (p->slist == NULL) {		/* missing sample surface! */
		error(USER, "no sender surface");
		return(-1);
	}
					/* misplaced output file spec. */
	if (p->outfn[0]) {
		sprintf(errmsg, "ignoring output file in sender ('%s')",
				p->outfn);
		error(WARNING, errmsg);
	}
					/* check/set basis hemisphere */
	if (!p->hemis[0]) {
		error(USER, "missing sender sampling type");
		return(-1);
	}
	if (normalize(p->nrm) == 0) {
		error(USER, "undefined normal for sender sampling");
		return(-1);
	}
	if (normalize(p->vup) == 0) {
		if (fabs(p->nrm[2]) < .7)
			p->vup[2] = 1;
		else
			p->vup[1] = 1;
	}
	fcross(p->udir, p->vup, p->nrm);
	if (normalize(p->udir) == 0) {
		error(USER, "up vector coincides with sender normal");
		return(-1);
	}
	fcross(p->vdir, p->nrm, p->udir);
	if (p->sign == '-') {		/* left-handed coordinate system? */
		p->udir[0] *= -1.;
		p->udir[1] *= -1.;
		p->udir[2] *= -1.;
	}
	if ((tolower(p->hemis[0]) == 'u') | (p->hemis[0] == '1'))
		p->sample_basis = sample_uniform;
	else if (tolower(p->hemis[0]) == 's' &&
				tolower(p->hemis[1]) == 'c')
		p->sample_basis = sample_shirchiu;
	else if ((tolower(p->hemis[0]) == 'r') |
			(tolower(p->hemis[0]) == 't'))
		p->sample_basis = sample_reinhart;
	else if (tolower(p->hemis[0]) == 'k') {
		switch (p->hemis[1]) {
		case '1':
		case '2':
		case '4':
			break;
		case 'f':
		case 'F':
		case '\0':
			p->hemis[1] = '1';
			break;
		case 'h':
		case 'H':
			p->hemis[1] = '2';
			break;
		case 'q':
		case 'Q':
			p->hemis[1] = '4';
			break;
		default:
			goto unrecognized;
		}
		p->hemis[2] = '\0';
		p->sample_basis = sample_klems;
	} else
		goto unrecognized;
					/* return number of bins */
	return((*p->sample_basis)(p,0,NULL));
unrecognized:
	sprintf(errmsg, "unrecognized sender sampling: h=%s", p->hemis);
	error(USER, errmsg);
	return(-1);
}

/* Compute normal and area for polygon */
static int
finish_polygon(SURF *p)
{
	const int	nv = p->nfargs / 3;
	FVECT		e1, e2, vc;
	int		i;

	memset(p->snrm, 0, sizeof(FVECT));
	VSUB(e1, p->farg+3, p->farg);
	for (i = 2; i < nv; i++) {
		VSUB(e2, p->farg+3*i, p->farg); 
		VCROSS(vc, e1, e2);
		p->snrm[0] += vc[0];
		p->snrm[1] += vc[1];
		p->snrm[2] += vc[2];
		VCOPY(e1, e2);
	}
	p->area = normalize(p->snrm)*0.5;
	return(p->area > FTINY*FTINY);
}

/* Add a surface to our current parameters */
static void
add_surface(int st, const char *oname, FILE *fp)
{
	SURF	*snew;
	int	n;
					/* get floating-point arguments */
	if (!fscanf(fp, "%d", &n)) return;
	while (n-- > 0) fscanf(fp, "%*s");
	if (!fscanf(fp, "%d", &n)) return;
	while (n-- > 0) fscanf(fp, "%*d");
	if (!fscanf(fp, "%d", &n) || n <= 0) return;
	snew = (SURF *)malloc(sizeof(SURF) + sizeof(double)*(n-1));
	if (snew == NULL)
		error(SYSTEM, "out of memory in add_surface()");
	strncpy(snew->sname, oname, sizeof(snew->sname)-1);
	snew->sname[sizeof(snew->sname)-1] = '\0';
	snew->styp = st;
	snew->priv = NULL;
	snew->nfargs = n;
	for (n = 0; n < snew->nfargs; n++)
		if (fscanf(fp, "%lf", &snew->farg[n]) != 1) {
			sprintf(errmsg, "error reading arguments for '%s'",
					oname);
			error(USER, errmsg);
		}
	switch (st) {
	case ST_RING:
		if (snew->nfargs != 8)
			goto badcount;
		VCOPY(snew->snrm, snew->farg+3);
		if (normalize(snew->snrm) == 0)
			goto badnorm;
		if (snew->farg[7] < snew->farg[6]) {
			double	t = snew->farg[7];
			snew->farg[7] = snew->farg[6];
			snew->farg[6] = t;
		}
		snew->area = PI*(snew->farg[7]*snew->farg[7] -
					snew->farg[6]*snew->farg[6]);
		break;
	case ST_POLY:
		if (snew->nfargs < 9 || snew->nfargs % 3)
			goto badcount;
		finish_polygon(snew);
		break;
	case ST_SOURCE:
		if (snew->nfargs != 4)
			goto badcount;
		for (n = 3; n--; )	/* need to reverse "normal" */
			snew->snrm[n] = -snew->farg[n];
		if (normalize(snew->snrm) == 0)
			goto badnorm;
		snew->area = sin((PI/180./2.)*snew->farg[3]);
		snew->area *= PI*snew->area;
		break;
	}
	if (snew->area <= FTINY*FTINY) {
		sprintf(errmsg, "zero area for surface '%s'", oname);
		error(WARNING, errmsg);
		free(snew);
		return;
	}
	VSUM(curparams.nrm, curparams.nrm, snew->snrm, snew->area);
	snew->next = curparams.slist;
	curparams.slist = snew;
	curparams.nsurfs++;
	return;
badcount:
	sprintf(errmsg, "bad argument count for surface element '%s'", oname);
	error(USER, errmsg);
badnorm:
	sprintf(errmsg, "bad orientation for surface element '%s'", oname);
	error(USER, errmsg);
}

/* Parse a receiver object (look for modifiers to add) */
static int
add_recv_object(FILE *fp)
{
	int		st;
	char		thismod[128], otype[32], oname[128];
	int		n;

	if (fscanf(fp, "%s %s %s", thismod, otype, oname) != 3)
		return(0);		/* must have hit EOF! */
	if (!strcmp(otype, "alias")) {
		fscanf(fp, "%*s");	/* skip alias */
		return(0);
	}
					/* is it a new receiver? */
	if ((st = surf_type(otype)) != ST_NONE) {
		if (strcmp(thismod, curmod)) {
			if (curmod[0])	/* output last receiver? */
				finish_receiver();
			parse_params(&curparams, newparams);
			newparams[0] = '\0';
			strcpy(curmod, thismod);
		}
		add_surface(st, oname, fp);	/* read & store surface */
		return(1);
	}
					/* else skip arguments */
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*s");
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*d");
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*f");
	return(0);
}

/* Parse a sender object */
static int
add_send_object(FILE *fp)
{
	int		st;
	char		thismod[128], otype[32], oname[128];
	int		n;

	if (fscanf(fp, "%s %s %s", thismod, otype, oname) != 3)
		return(0);		/* must have hit EOF! */
	if (!strcmp(otype, "alias")) {
		fscanf(fp, "%*s");	/* skip alias */
		return(0);
	}
					/* is it a new surface? */
	if ((st = surf_type(otype)) != ST_NONE) {
		if (st == ST_SOURCE) {
			error(USER, "cannot use source as a sender");
			return(-1);
		}
		if (strcmp(thismod, curmod)) {
			if (curmod[0])
				error(WARNING, "multiple modifiers in sender");
			strcpy(curmod, thismod);
		}
		parse_params(&curparams, newparams);
		newparams[0] = '\0';
		add_surface(st, oname, fp);	/* read & store surface */
		return(0);
	}
					/* else skip arguments */
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*s");
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*d");
	if (!fscanf(fp, "%d", &n)) return(0);
	while (n-- > 0) fscanf(fp, "%*f");
	return(0);
}

/* Load a Radiance scene using the given callback function for objects */
static int
load_scene(const char *inspec, int (*ocb)(FILE *))
{
	int	rv = 0;
	char	inpbuf[1024];
	FILE	*fp;
	int	c;

	if (*inspec == '!')
		fp = popen(inspec+1, "r");
	else
		fp = fopen(inspec, "r");
	if (fp == NULL) {
		sprintf(errmsg, "cannot load '%s'", inspec);
		error(SYSTEM, errmsg);
		return(-1);
	}
	while ((c = getc(fp)) != EOF) {	/* load sender/receiver data */
		if (isspace(c))		/* skip leading white space */
			continue;
		if (c == '!') {		/* read from a new command */
			inpbuf[0] = c;
			if (fgetline(inpbuf+1, sizeof(inpbuf)-1, fp) != NULL) {
				if ((c = load_scene(inpbuf, ocb)) < 0)
					return(c);
				rv += c;
			}
			continue;
		}
		if (c == '#') {		/* parameters/comment */
			if ((c = getc(fp)) == EOF)
				break;
			if (c == '@' && fscanf(fp, "%s", inpbuf) == 1 &&
					(!strcmp(inpbuf, PARAMSNAME) ||
					 !strcmp(inpbuf, progname))) {
				if (fgets(inpbuf, sizeof(inpbuf), fp) != NULL)
					strlcat(newparams, inpbuf, sizeof(newparams));
				continue;
			}
			while (c != '\n' && (c = getc(fp)) != EOF)
				;	/* ...skipping comment */
			continue;
		}
		ungetc(c, fp);		/* else check object for receiver */
		c = (*ocb)(fp);
		if (c < 0)
			return(c);
		rv += c;
	}
					/* close our input stream */
	c = (*inspec == '!') ? pclose(fp) : fclose(fp);
	if (c != 0) {
		sprintf(errmsg, "error loading '%s'", inspec);
		error(SYSTEM, errmsg);
		return(-1);
	}
	return(rv);
}

void
wputs(				/* warning output function */
	const char	*s
)
{
	if (erract[WARNING].pf == NULL)
		return;
	int  lasterrno = errno;
	eputs(s);
	errno = lasterrno;
}

void
eputs(				/* put string to stderr */
	const char  *s
)
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n') {
		fflush(stderr);
		midline = 0;
	}
}

/* set input/output format */
static void
setformat(const char *fmt)
{
	switch (fmt[0]) {
	case 'f':
	case 'd':
		SET_FILE_BINARY(stdin);
		/* fall through */
	case 'a':
		inpfmt = fmt[0];
		break;
	case 'c':
		if (fmt[1])
			goto fmterr;
		outfmt = fmt[0];	// special case for output-only
		return;
	default:
		goto fmterr;
	}
	switch (fmt[1]) {
	case '\0':
		if (inpfmt == 'a')
			goto fmterr;
		outfmt = inpfmt;
		return;
	case 'f':
	case 'd':
	case 'c':
		outfmt = fmt[1];
		break;
	default:
		goto fmterr;
	}
	if (!fmt[2])
		return;
fmterr:
	sprintf(errmsg, "unsupported i/o format: -f%s", fmt);
	error(USER, errmsg);
}

static inline double
pixjitter()
{
	return(0.5 + dstrpix*(frandom()-0.5));
}

// Compute a set of view rays for the given pixel accumulator
bool
viewRayBundle(FVECT orig_dir[], int x, int y)
{
	for (int n = 0; n < myRCmanager.accum; orig_dir += 2, n++) {
		const double	d = viewray(orig_dir[0], orig_dir[1], &ourview,
						(x+pixjitter())/myRCmanager.xres,
						(y+pixjitter())/myRCmanager.yres);
					// off-view sample?
		if (d < -FTINY || !jitteraperture(orig_dir[0], orig_dir[1],
							&ourview, dblur))
			memset(orig_dir+1, 0, sizeof(FVECT));
		else			// else record distance to aft plane
			for (int i = 3*(d > FTINY); i--; )
				orig_dir[1][i] *= d;
					// accumulating identical samples?
		if (!n && (dstrpix <= 0.05) & (dblur <= FTINY)) {
			while (++n < myRCmanager.accum)
				memcpy(orig_dir[2*n], orig_dir[0], sizeof(FVECT)*2);
			break;
		}
	}
	return true;
}

// skip specified number of bytes, return false if EOF
static bool
skipBytes(size_t n2skip)
{
	while (n2skip--)
		if (getchar() == EOF)
			return false;
	return true;
}

// skip specified number of whitespace-separated words, return false if EOF
static bool
skipWords(int n2skip)
{
	int	c;

	while (n2skip--) {
		do {
			c = getchar();
		} while (isspace(c));
		do {
			if (c == EOF) return false;
			c = getchar();
		} while (!isspace(c));
	}
	return true;
}

// Skip a set of input rays
bool
skipRayBundle()
{
	switch (inpfmt) {
	case 'd':
		return skipBytes(sizeof(double)*6*myRCmanager.accum);
	case 'f':
		return skipBytes(sizeof(float)*6*myRCmanager.accum);
	case 'a':
		return skipWords(6*myRCmanager.accum);
	}
	error(INTERNAL, "unsupported input format");
	return false;
}

// Load a set of rays for accumulation (do not normalize directions)
bool
getRayBundle(FVECT orig_dir[])
{
						// read directly if possible
	if (inpfmt == "_fd"[sizeof(RREAL)/sizeof(float)])
		return (getbinary(orig_dir[0], sizeof(FVECT)*2,
				myRCmanager.accum, stdin) == myRCmanager.accum);

	for (int n = 0; n < myRCmanager.accum; orig_dir += 2, n++)
		switch (inpfmt) {
#ifdef SMLFLT
		case 'd': { double	dvin[6];
			if (getbinary(dvin, sizeof(dvin), 1, stdin) != 1)
				return false;
			for (int i = 6; i--; ) orig_dir[0][i] = dvin[i];
			} break;
#else
		case 'f': { float	fvin[6];
			if (getbinary(fvin, sizeof(fvin), 1, stdin) != 1)
				return false;
			for (int i = 6; i--; ) orig_dir[0][i] = fvin[i];
			} break;
#endif
		case 'a':
			if (scanf(FVFORMAT, &orig_dir[0][0], &orig_dir[0][1],
					&orig_dir[0][2]) != 3)
				return false;
			if (scanf(FVFORMAT, &orig_dir[1][0], &orig_dir[1][1],
					&orig_dir[1][2]) != 3)
				return false;
			break;
		default:
			error(INTERNAL, "unsupported input format");
		}

	return true;
}

/* Set default options */
static void
default_options()
{
	rand_samp = 1;
	dstrsrc = 0.9;
	directrelay = 3;
	vspretest = 512;
	srcsizerat = .2;
	specthresh = .02;
	specjitter = 1.;
	maxdepth = -10;
	minweight = 2e-3;
	ambres = 256;
	ambdiv = 350;
	ambounce = 1;
}

/* Set overriding options */
static void
override_options()
{
	shadthresh = 0;
	ambssamp = 0;
	ambacc = 0;
}

void
onsig(				/* fatal signal */
	int  signo
)
{
	static int  gotsig = 0;

	if (gotsig++)			/* two signals and we're gone! */
		_exit(signo);

#ifdef SIGALRM
	alarm(600);			/* allow 10 minutes to clean up */
	signal(SIGALRM, SIG_DFL);	/* make certain we do die */
#endif
	eputs("signal - ");
	eputs(sigerr[signo]);
	eputs("\n");
	quit(3);
}

void
sigdie(			/* set fatal signal */
	int  signo,
	const char  *msg
)
{
	if (signal(signo, onsig) == SIG_IGN)
		signal(signo, SIG_IGN);
	sigerr[signo] = msg;
}

// report progress
void
report_progress(bool force = false)
{
	static time_t	last_report=0, tstart=time(0);
	time_t		tnow;

	if (!force & (report_intvl <= 0))
		return;

	tnow = time(0);
	if (!force & (tnow < last_report + report_intvl))
		return;

	sprintf(errmsg, "%.2f%% done after %.3f hours\n",
			100. * myRCmanager.GetRowFinished() /
				(double)myRCmanager.GetRowMax(),
			(tnow - tstart)*(1./3600.));
	eputs(errmsg);
	last_report = tnow;
}

/* Run rfluxmtx equivalent without leaning on r[x]contrib */
int
main(int argc, char *argv[])
{
#define	 check(ol,al)		if (argv[a][ol] || \
				badarg(argc-a-1,argv+a+1,al)) \
				goto userr
#define	 check_bool(olen,var)		switch (argv[a][olen]) { \
				case '\0': var = !var; break; \
				case '+': case '1': var = true; break; \
				case '-': case '0': var = false; break; \
				default: goto userr; }
	bool		force_open = false;
	bool		recover = false;
	bool		gotView = false;
	double		pixaspect = 1.;
	int		nproc = 1;
	int		nout = 0;
	char		*outfn = NULL;
	FVECT		*rayarr = NULL;
	const char	*sendfn = "";
	double		binjitter = 0;
	PARAMS		sendparams;
	int		rval;
	int		a, i;
					/* set global progname */
	fixargv0(argv[0]);
					// just asking version?
	if (argc == 2 && !strcmp(argv[1], "-version")) {
		puts(VersionID);
		quit(0);
	}
	initfunc();			/* initialize calcomp routines first */
	calcontext(RCCONTEXT);
					/* set rcontrib defaults */
	default_options();
					/* get command-line options */
	for (a = 1; a < argc; a++) {
					/* check for argument expansion */
		while ((rval = expandarg(&argc, &argv, a)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", argv[a]);
			error(SYSTEM, errmsg);
		}
		if (argv[a][0] != '-' || !argv[a][1])
			break;			// break from options
		rval = getrenderopt(argc-a, argv+a);
		if (rval >= 0) {
			a += rval;
			continue;
		}
		rval = getviewopt(&ourview, argc-a, argv+a);
		if (rval >= 0) {
			a += rval;
			gotView = true;
			continue;
		}
		switch (argv[a][1]) {	/* !! Keep consistent !! */
		case 'b':			/* bin jitter? */
			if (argv[a][2] != 'j')
				goto userr;
			check(3,"f");
			binjitter = atof(argv[++a]);
			break;
		case 'v':			// view file or verbose
			if (strchr("+-10", argv[a][2]) != NULL) {
				check_bool(2,verby);
				break;
			}
			if (argv[a][2] != 'f')
				goto userr;
			check(3,"s");
			rval = viewfile(argv[++a], &ourview, NULL);
			if (rval < 0) {
				sprintf(errmsg, "cannot open view file '%s'", argv[a]);
				error(SYSTEM, errmsg);
			} else if (!rval) {
				sprintf(errmsg, "bad view file '%s'", argv[a]);
				error(USER, errmsg);
			}
			gotView = true;
			break;
		case 'p':			// -pj, -pd, or -pa
			switch (argv[a][2]) {
			case 'j':				/* jitter */
				check(3,"f");
				dstrpix = atof(argv[++a]);
				break;
			case 'd':				/* aperture */
				check(3,"f");
				dblur = atof(argv[++a]);
				break;
			case 'a':				/* aspect */
				check(3,"f");
				pixaspect = atof(argv[++a]);
				break;
			default:
				goto userr;
			}
			break;
		case 'n':			/* number of processes */
			check(2,"i");
			nproc = atoi(argv[++a]);
			if (nproc < 0 && (nproc += RadSimulManager::GetNCores()) <= 0)
				nproc = 1;
			break;
		case 'V':			/* output contributions? */
			rval = myRCmanager.HasFlag(RCcontrib);
			check_bool(2,rval);
			myRCmanager.SetFlag(RCcontrib, rval);
			break;
		case 'x':			/* x resolution */
			check(2,"i");
			myRCmanager.xres = atoi(argv[++a]);
			break;
		case 'y':			/* y resolution */
			check(2,"i");
			myRCmanager.yres = atoi(argv[++a]);
			break;
		case 'w':			/* warnings on/off */
			rval = (erract[WARNING].pf != NULL);
			check_bool(2,rval);
			if (rval) erract[WARNING].pf = wputs;
			else erract[WARNING].pf = NULL;
			break;
		case 'l':			/* limit distance */
			if (argv[a][2] != 'd')
				goto userr;
			rval = myRCmanager.HasFlag(RTlimDist);
			check_bool(3,rval);
			myRCmanager.SetFlag(RTlimDist, rval);
			break;
		case 'I':			/* immed. irradiance */
			rval = myRCmanager.HasFlag(RTimmIrrad);
			check_bool(2,rval);
			myRCmanager.SetFlag(RTimmIrrad, rval);
			break;
		case 'f':			/* format or force overwrite */
			if (argv[a][2] == 'o') {
				check_bool(3,force_open);
				break;
			}
			setformat(argv[a]+2);
			break;
		case 'r':			// recover flag
			check_bool(2,recover);
			break;
		case 'o':			/* output file */
			check(2,"s");
			outfn = argv[++a];
			break;
		case 'c':			/* sample count */
			check(2,"i");
			myRCmanager.accum = atoi(argv[++a]);
			break;
		case 't':			/* reporting interval */
			check(2,"i");
			report_intvl = atoi(argv[++a]);
			break;
		default:		/* anything else is verbotten */
			goto userr;
		}
	}
	if (a > argc-1)
		goto userr;

	override_options();		/* override critical options */

	if (!gotView) {
		sendfn = argv[a++];
		if (a > argc-1)
			goto userr;
	} else if (argv[a][0] == '-')
		error(USER, "view specification incompatible with pass-through mode");
					/* assign sender & receiver inputs */
	if (gotView) {			// picture output?
		const char *	err = setview(&ourview);
		if (err != NULL)
			error(USER, err);
		if (myRCmanager.HasFlag(RTimmIrrad))
			error(USER, "-I+ not allowed in view generation mode");
		normaspect(viewaspect(&ourview), &pixaspect,
				&myRCmanager.xres, &myRCmanager.yres);
		if ((myRCmanager.xres <= 0) | (myRCmanager.yres <= 0))
			error(USER, "missing or illegal -x and -y resolution");
		myRCmanager.SetFlag(RTlimDist, ourview.vaft > FTINY);
		if (myRCmanager.accum <= 0)
			myRCmanager.accum = 1;
	} else if (sendfn[0] == '-') {	// pass-through mode?
		if (sendfn[1]) goto userr;
		sendfn = NULL;
		if (myRCmanager.accum <= 0)
			myRCmanager.accum = 1;
	} else {			// else surface sampling
		if (do_irrad | myRCmanager.HasFlag(RTimmIrrad))
			error(USER, "-i+, -I+ not allowed in surface-sampling mode");
		myRCmanager.SetFlag(RTlimDist, false);
		if (load_scene(sendfn, add_send_object) < 0)
			quit(1);
		sendparams = curparams;	// save sender params & clear current
		memset(&curparams, 0, sizeof(PARAMS));
		curmod[0] = '\0';
		myRCmanager.yres = prepare_sampler(&sendparams);
		if (myRCmanager.yres <= 0)
			quit(1);
		myRCmanager.xres = 0;
		if (myRCmanager.accum <= 1)
			myRCmanager.accum = 10000;
	}
	if (outfn != NULL)		// saved from -o option above?
		strlcpy(curparams.outfn, outfn, sizeof(curparams.outfn));
					// get ready to rock...
	if (setspectrsamp(CNDX, WLPART) < 0)
		error(USER, "unsupported spectral sampling");
					/* set up signal handling */
	sigdie(SIGINT, "Interrupt");
#ifdef SIGHUP
	sigdie(SIGHUP, "Hangup");
#endif
	sigdie(SIGTERM, "Terminate");
#ifdef SIGPIPE
	sigdie(SIGPIPE, "Broken pipe");
#endif
#ifdef SIGALRM
	sigdie(SIGALRM, "Alarm clock");
#endif
#ifdef	SIGXCPU
	sigdie(SIGXCPU, "CPU limit exceeded");
	sigdie(SIGXFSZ, "File size exceeded");
#endif
#ifdef	NICE
	nice(NICE);			/* lower priority */
#endif
	rayarr = new FVECT [myRCmanager.accum*2];
					// load octree
	if (!myRCmanager.LoadOctree(oconv_command(argc-a, argv+a)))
		quit(1);
					// add to header
	myRCmanager.AddHeader(argc, argv);
	{
		char	buf[128] = "SOFTWARE= ";
		strcpy(buf+10, VersionID);
		myRCmanager.AddHeader(buf);
		if (gotView) {
			sprintf(buf, "%s%s", VIEWSTR, viewopt(&ourview));
			if (strlen(buf) > VIEWSTRL+3)
				myRCmanager.AddHeader(buf);
		}
		if ((pixaspect > .005) & ((pixaspect < .995) |
					  (pixaspect > 1.005))) {
			sprintf(buf, "%s%f", ASPECTSTR, pixaspect);
			myRCmanager.AddHeader(buf);
		}
	}
					// set output format
	myRCmanager.SetDataFormat(outfmt);
					/* assign receiver modifiers */
	if (load_scene(argv[a], add_recv_object) < 0)
		quit(1);
	finish_receiver();		// makes final AddModifier() call
	if (binjitter > FTINY)		// global bin jitter?
		varset(const_cast<char *>("JTR"), '=', binjitter);
					// prepare output files
	if (recover) {
		if (force_open) {
			error(WARNING, "-r+ mode overrides -fo+");
			force_open = false;
		}
		myRCmanager.outOp = RCOrecover;
	} else if (force_open)
		myRCmanager.outOp = RCOforce;
	else
		myRCmanager.outOp = RCOnew;
					// rval = # rows recovered
	rval = myRCmanager.PrepOutput();
	if (rval < 0)			// PrepOutput() failure?
		error(USER, "issue loading or creating output");
					// in case output is complete
	if (rval >= myRCmanager.GetRowMax()) {
		error(WARNING, "nothing left to compute");
		quit(0);
	}
	if (verby) {			// get output file count?
		if (rval > 0) {
			sprintf(errmsg, "recovered %d of %d rows\n",
					rval, myRCmanager.GetRowMax());
			eputs(errmsg);
		}
		for (const RcontribOutput *op = myRCmanager.GetOutput();
				op != NULL; op = op->Next())
			++nout;
	}
	if (nproc > 1) {		// set #processes
		if (verby) {
			sprintf(errmsg, "starting %d subprocesses\n", nproc);
			eputs(errmsg);
		}
		myRCmanager.SetThreadCount(nproc);
	}
	if (gotView) {			// picture generation mode?
		i = myRCmanager.GetRowCount();
		if (verby) {
			sprintf(errmsg, "%s %d %dx%d pictures\n",
					i ? "completing" : "computing",
					nout, myRCmanager.xres, myRCmanager.yres);
			if (myRCmanager.accum > 1)
				sprintf(errmsg+strlen(errmsg)-1, " with %d samples/pixel\n",
						myRCmanager.accum);
			eputs(errmsg);
		}
		i /= myRCmanager.xres;
		int	xstart = myRCmanager.GetRowCount() - i*myRCmanager.xres;
		i = myRCmanager.yres - i;
		while (i--) {		// compute pixel rows from top down
			for (int x = xstart; x < myRCmanager.xres; x++) {
				report_progress();
				if (!viewRayBundle(rayarr, x, i))
					quit(1);
				if (myRCmanager.ComputeRecord(rayarr) != myRCmanager.accum)
					error(USER, "failed call to ComputeRecord()");
			}
			xstart = 0;	// all rows after first start at x=0
		}
	} else if (sendfn == NULL) {	// pass-through mode?
#ifdef getc_unlocked
		flockfile(stdin);
#endif
					// skip completed rows
		for (i = 0; i < myRCmanager.GetRowCount(); i++)
			if (!skipRayBundle()) {
				sprintf(errmsg, "read error from stdin at row %d", i);
				error(SYSTEM, errmsg);
			}
		if (verby) {
			sprintf(errmsg, "computing %d%s rows in %d matrices\n",
					myRCmanager.GetRowMax()-i, i ? " remaining" : "", nout);
			if (myRCmanager.accum > 1)
				sprintf(errmsg+strlen(errmsg)-1, " with %d samples/row\n",
						myRCmanager.accum);
			eputs(errmsg);
		}
		for ( ; i < myRCmanager.GetRowMax(); i++) {
			report_progress();
			if (!getRayBundle(rayarr)) {
				sprintf(errmsg, "read error from stdin at row %d of %d",
						i, myRCmanager.GetRowMax());
				error(USER, errmsg);
			}
			if (myRCmanager.ComputeRecord(rayarr) != myRCmanager.accum)
				error(USER, "failed call to ComputeRecord()");
		}
	} else {			// else surface-sampling mode
		i = myRCmanager.GetRowCount();
		if (verby) {
			sprintf(errmsg, "sampling %d%s directions in %d matrices with %d samples/direction\n",
					myRCmanager.yres-i, i ? " remaining" : "",
					nout, myRCmanager.accum);
			if (sendparams.nsurfs > 1)
				sprintf(errmsg+strlen(errmsg)-1, " (%d surface elements)\n", sendparams.nsurfs);
			eputs(errmsg);
		}
		for ( ; i < myRCmanager.yres; i++) {
			report_progress();
			if (!(*sendparams.sample_basis)(&sendparams, i, rayarr))
				quit(1);
			if (myRCmanager.ComputeRecord(rayarr) != myRCmanager.accum)
				error(USER, "failed call to ComputeRecord()");
		}
		clear_params(&sendparams);
	}
	delete [] rayarr;
	myRCmanager.FlushQueue();
	report_progress((report_intvl > 0) | verby);
	quit(0);			/* waits on any children */
userr:
	if (a < argc && argv[a][0] == '-')
		fprintf(stderr, "%s: unsupported/misplaced option '%s'\n", progname, argv[a]);
	fprintf(stderr, "Usage: %s [-W][-bj frac] [rcontrib options] { sender.rad | view | - } receiver.rad [-i system.oct] [system.rad ..]\n",
				progname);
	quit(1);
}

/* Exit program */
void
quit(
	int  code
)
{
	if (!code)
		myRCmanager.ClearModifiers();

	exit(code);
}
