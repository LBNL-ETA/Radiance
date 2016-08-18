#ifndef lint
static const char RCSid[] = "$Id: rfluxmtx.c,v 2.38 2016/08/18 00:52:48 greg Exp $";
#endif
/*
 * Calculate flux transfer matrix or matrices using rcontrib
 */

#include "copyright.h"

#include <ctype.h>
#include <stdlib.h>
#include "rtio.h"
#include "rtmath.h"
#include "paths.h"
#include "bsdf.h"
#include "bsdf_m.h"
#include "random.h"
#include "triangulate.h"
#include "platform.h"

#define MAXRCARG	512

char		*progname;		/* global argv[0] */

int		verbose = 0;		/* verbose mode (< 0 no warnings) */

char		*rcarg[MAXRCARG+1] = {"rcontrib", "-fo+"};
int		nrcargs = 2;

const char	overflowerr[] = "%s: too many arguments!\n";

#define	CHECKARGC(n)	if (nrcargs >= MAXRCARG-(n)) \
	{ fprintf(stderr, overflowerr, progname); exit(1); }

int		sampcnt = 0;		/* sample count (0==unset) */

char		*reinhfn = "reinhartb.cal";
char		*shirchiufn = "disk2square.cal";
char		*kfullfn = "klems_full.cal";
char		*khalffn = "klems_half.cal";
char		*kquarterfn = "klems_quarter.cal";

					/* string indicating parameters */
const char	PARAMSTART[] = "@rfluxmtx";

				/* surface type IDs */
#define ST_NONE		0
#define ST_POLY		1
#define ST_RING		2
#define ST_SOURCE	3

typedef struct surf_s {
	struct surf_s	*next;		/* next surface in list */
	void		*priv;		/* private data (malloc'ed) */
	char		sname[32];	/* surface name */
	FVECT		snrm;		/* surface normal */
	double		area;		/* surface area / proj. solid angle */
	short		styp;		/* surface type */
	short		nfargs;		/* number of real arguments */
	double		farg[1];	/* real values (extends struct) */
} SURF;				/* surface structure */

typedef struct {
	FVECT	uva[2];			/* tangent axes */
	int	ntris;			/* number of triangles */
	struct ptri {
		float	afrac;			/* fraction of total area */
		short	vndx[3];		/* vertex indices */
	}	tri[1];			/* triangle array (extends struct) */
} POLYTRIS;			/* triangulated polygon */

typedef struct param_s {
	char		sign;		/* '-' for axis reversal */
	char		hemis[31];	/* hemispherical sampling spec. */
	int		hsiz;		/* hemisphere basis size */
	int		nsurfs;		/* number of surfaces */
	SURF		*slist;		/* list of surfaces */
	FVECT		vup;		/* up vector (zero if unset) */
	FVECT		nrm;		/* average normal direction */
	FVECT		udir, vdir;	/* tangent axes */
	char		*outfn;		/* output file name (receiver) */
	int		(*sample_basis)(struct param_s *p, int, FILE *);
} PARAMS;			/* sender/receiver parameters */

PARAMS		curparams;
char		curmod[128];
char		newparams[1024];

typedef int	SURFSAMP(FVECT, SURF *, double);

static SURFSAMP	ssamp_bad, ssamp_poly, ssamp_ring;

SURFSAMP	*orig_in_surf[4] = {
		ssamp_bad, ssamp_poly, ssamp_ring, ssamp_bad
	};

/* Clear parameter set */
static void
clear_params(PARAMS *p, int reset_only)
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
		memset(p->nrm, 0, sizeof(FVECT));
		memset(p->vup, 0, sizeof(FVECT));
		p->outfn = NULL;
		return;
	}
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
	static char	oconvbuf[2048] = "!oconv -f ";
	char		*cp = oconvbuf + 10;
	char		*recv = *av++;
	
	if (ac-- <= 0)
		return(NULL);
	while (ac-- > 0) {
		strcpy(cp, *av++);
		while (*cp) cp++;
		*cp++ = ' ';
		if (cp >= oconvbuf+(sizeof(oconvbuf)-32))
			goto overrun;
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
	fputs(progname, stderr);
	fputs(": too many file arguments!\n", stderr);
	exit(1);
}

/* Open a pipe to/from a command given as an argument list */
static FILE *
popen_arglist(char *av[], char *mode)
{
	char	cmd[10240];

	if (!convert_commandline(cmd, sizeof(cmd), av)) {
		fputs(progname, stderr);
		fputs(": command line too long in popen_arglist()\n", stderr);
		return(NULL);
	}
	if (verbose > 0)
		fprintf(stderr, "%s: opening pipe %s: %s\n",
				progname, (*mode=='w') ? "to" : "from", cmd);
	return(popen(cmd, mode));
}

#if defined(_WIN32) || defined(_WIN64)
/* Execute system command (Windows version) */
static int
my_exec(char *av[])
{
	char	cmd[10240];

	if (!convert_commandline(cmd, sizeof(cmd), av)) {
		fputs(progname, stderr);
		fputs(": command line too long in my_exec()\n", stderr);
		return(1);
	}
	if (verbose > 0)
		fprintf(stderr, "%s: running: %s\n", progname, cmd);
	return(system(cmd));
}
#else
/* Execute system command in our stead (Unix version) */
static int
my_exec(char *av[])
{
	char	*compath;

	if ((compath = getpath((char *)av[0], getenv("PATH"), X_OK)) == NULL) {
		fprintf(stderr, "%s: cannot locate %s\n", progname, av[0]);
		return(1);
	}
	if (verbose > 0) {
		char	cmd[4096];
		if (!convert_commandline(cmd, sizeof(cmd), av))
			strcpy(cmd, "COMMAND TOO LONG TO SHOW");
		fprintf(stderr, "%s: running: %s\n", progname, cmd);
	}
	execv(compath, av);	/* successful call never returns */
	perror(compath);
	return(1);
}
#endif

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
			p->outfn = savqstr(cp-i);
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
	fprintf(stderr, "%s: bad parameter string: %s", progname, pargs);
	exit(1);
	return(-1);	/* pro forma return */
}

/* Add receiver arguments (directives) corresponding to the current modifier */
static void
finish_receiver(void)
{
	char	sbuf[256];
	int	uniform = 0;
	char	*calfn = NULL;
	char	*params = NULL;
	char	*binv = NULL;
	char	*binf = NULL;
	char	*nbins = NULL;

	if (!curmod[0]) {
		fputs(progname, stderr);
		fputs(": missing receiver surface!\n", stderr);
		exit(1);
	}
	if (curparams.outfn != NULL) {	/* add output file spec. */
		CHECKARGC(2);
		rcarg[nrcargs++] = "-o";
		rcarg[nrcargs++] = curparams.outfn;
	}
					/* check arguments */
	if (!curparams.hemis[0]) {
		fputs(progname, stderr);
		fputs(": missing hemisphere sampling type!\n", stderr);
		exit(1);
	}
	if (normalize(curparams.nrm) == 0) {
		fputs(progname, stderr);
		fputs(": undefined normal for hemisphere sampling\n", stderr);
		exit(1);
	}
	if (normalize(curparams.vup) == 0) {
		if (fabs(curparams.nrm[2]) < .7)
			curparams.vup[2] = 1;
		else
			curparams.vup[1] = 1;
	}
					/* determine sample type/bin */
	if (tolower(curparams.hemis[0]) == 'u' | curparams.hemis[0] == '1') {
		sprintf(sbuf, "if(-Dx*%g-Dy*%g-Dz*%g,0,-1)",
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2]);
		binv = savqstr(sbuf);
		nbins = "1";		/* uniform sampling -- one bin */
		uniform = 1;
	} else if (tolower(curparams.hemis[0]) == 's' &&
				tolower(curparams.hemis[1]) == 'c') {
					/* assign parameters */
		if (curparams.hsiz <= 1) {
			fputs(progname, stderr);
			fputs(": missing size for Shirley-Chiu sampling!\n", stderr);
			exit(1);
		}
		calfn = shirchiufn; shirchiufn = NULL;
		sprintf(sbuf, "SCdim=%d,rNx=%g,rNy=%g,rNz=%g,Ux=%g,Uy=%g,Uz=%g,RHS=%c1",
				curparams.hsiz,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2],
			curparams.sign);
		params = savqstr(sbuf);
		binv = "scbin";
		nbins = "SCdim*SCdim";
	} else if ((tolower(curparams.hemis[0]) == 'r') |
			(tolower(curparams.hemis[0]) == 't')) {
		calfn = reinhfn; reinhfn = NULL;
		sprintf(sbuf, "MF=%d,rNx=%g,rNy=%g,rNz=%g,Ux=%g,Uy=%g,Uz=%g,RHS=%c1",
				curparams.hsiz,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2],
			curparams.sign);
		params = savqstr(sbuf);
		binv = "rbin";
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
	} else {
		fprintf(stderr, "%s: unrecognized hemisphere sampling: h=%s\n",
				progname, curparams.hemis);
		exit(1);
	}
	if (tolower(curparams.hemis[0]) == 'k') {
		sprintf(sbuf, "RHS=%c1", curparams.sign);
		params = savqstr(sbuf);
	}
	if (!uniform & (curparams.slist->styp == ST_SOURCE)) {
		SURF	*sp;
		for (sp = curparams.slist; sp != NULL; sp = sp->next)
			if (fabs(sp->area - PI) > 1e-3) {
				fprintf(stderr, "%s: source '%s' must be 180-degrees\n",
						progname, sp->sname);
				exit(1);
			}
	}
	if (calfn != NULL) {		/* add cal file if needed */
		CHECKARGC(2);
		rcarg[nrcargs++] = "-f";
		rcarg[nrcargs++] = calfn;
	}
	if (params != NULL) {		/* parameters _after_ cal file */
		CHECKARGC(2);
		rcarg[nrcargs++] = "-p";
		rcarg[nrcargs++] = params;
	}
	if (nbins != NULL) {		/* add #bins if set */
		CHECKARGC(2);
		rcarg[nrcargs++] = "-bn";
		rcarg[nrcargs++] = nbins;
	}
	if (binv != NULL) {
		CHECKARGC(2);		/* assign bin variable */
		rcarg[nrcargs++] = "-b";
		rcarg[nrcargs++] = binv;
	} else if (binf != NULL) {
		CHECKARGC(2);		/* assign bin function */
		rcarg[nrcargs++] = "-b";
		sprintf(sbuf, "%s(%g,%g,%g,%g,%g,%g)", binf,
			curparams.nrm[0], curparams.nrm[1], curparams.nrm[2],
			curparams.vup[0], curparams.vup[1], curparams.vup[2]);
		rcarg[nrcargs++] = savqstr(sbuf);
	}
	CHECKARGC(2);				/* modifier argument goes last */
	rcarg[nrcargs++] = "-m";
	rcarg[nrcargs++] = savqstr(curmod);
}

/* Make randomly oriented tangent plane axes for given normal direction */
static void
make_axes(FVECT uva[2], const FVECT nrm)
{
	int	i;

	if (!getperpendicular(uva[0], nrm, 1)) {
		fputs(progname, stderr);
		fputs(": bad surface normal in make_axes!\n", stderr);
		exit(1);
	}
	fcross(uva[1], nrm, uva[0]);
}

/* Illegal sender surfaces end up here */
static int
ssamp_bad(FVECT orig, SURF *sp, double x)
{
	fprintf(stderr, "%s: illegal sender surface '%s'\n",
			progname, sp->sname);
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
			fputs(progname, stderr);
			fputs(": out of memory in ssamp_ring!\n", stderr);
			return(0);
		}
		make_axes(uva, sp->snrm);
		sp->priv = (void *)uva;
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
	struct ptri	*trip = ptp->tri + ptp->ntris++;

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
				sizeof(struct ptri)*(sp->nfargs/3 - 3));
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
			v2l->p = (void *)ptp;
			if (!polyTriangulate(v2l, add_triangle)) {
				fprintf(stderr,
					"%s: cannot triangulate polygon '%s'\n",
						progname, sp->sname);
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
		sp->priv = (void *)ptp;
	}
					/* pick triangle by partial area */
	for (i = 0; i < ptp->ntris && x > ptp->tri[i].afrac; i++)
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
	fputs(progname, stderr);
	fputs(": out of memory in ssamp_poly!\n", stderr);
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
			fprintf(stderr,
				"%s: internal - sample behind sender '%s'\n",
					progname, sp->sname);
			return(0);
		}
		return((*orig_in_surf[sp->styp])(orig, sp, x));
	}
	if (p->nsurfs > nall) {		/* (re)allocate surface area cache */
		if (projsa) free(projsa);
		projsa = (double *)malloc(sizeof(double)*p->nsurfs);
		if (projsa == NULL) {
			fputs(progname, stderr);
			fputs(": out of memory in sample_origin!\n", stderr);
			exit(1);
		}
		nall = p->nsurfs;
	}
					/* compute projected areas */
	for (i = 0, sp = p->slist; sp != NULL; i++, sp = sp->next) {
		projsa[i] = -DOT(sp->snrm, rdir) * sp->area;
		tarea += projsa[i] *= (double)(projsa[i] > FTINY);
	}
	if (tarea <= FTINY) {		/* wrong side of sender? */
		fputs(progname, stderr);
		fputs(": internal - sample behind all sender elements!\n",
				stderr);
		return(0);
	}
	tarea *= x;			/* get surface from list */
	for (i = 0, sp = p->slist; tarea > projsa[i]; sp = sp->next)
		tarea -= projsa[i++];
	return((*orig_in_surf[sp->styp])(orig, sp, tarea/projsa[i]));
}

/* Uniform sample generator */
static int
sample_uniform(PARAMS *p, int b, FILE *fp)
{
	int	n = sampcnt;
	double	samp3[3];
	double	duvw[3];
	FVECT	orig_dir[2];
	int	i;

	if (fp == NULL)			/* just requesting number of bins? */
		return(1);

	while (n--) {			/* stratified hemisphere sampling */
		SDmultiSamp(samp3, 3, (n+frandom())/sampcnt);
		SDsquare2disk(duvw, samp3[1], samp3[2]);
		duvw[2] = -sqrt(1. - duvw[0]*duvw[0] - duvw[1]*duvw[1]);
		for (i = 3; i--; )
			orig_dir[1][i] = duvw[0]*p->udir[i] +
						duvw[1]*p->vdir[i] +
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		if (putbinary(orig_dir, sizeof(FVECT), 2, fp) != 2)
			return(0);
	}
	return(1);
}

/* Shirly-Chiu sample generator */
static int
sample_shirchiu(PARAMS *p, int b, FILE *fp)
{
	int	n = sampcnt;
	double	samp3[3];
	double	duvw[3];
	FVECT	orig_dir[2];
	int	i;

	if (fp == NULL)			/* just requesting number of bins? */
		return(p->hsiz*p->hsiz);

	while (n--) {			/* stratified sampling */
		SDmultiSamp(samp3, 3, (n+frandom())/sampcnt);
		SDsquare2disk(duvw, (b/p->hsiz + samp3[1])/curparams.hsiz,
				(b%p->hsiz + samp3[2])/curparams.hsiz);
		duvw[2] = sqrt(1. - duvw[0]*duvw[0] - duvw[1]*duvw[1]);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		if (putbinary(orig_dir, sizeof(FVECT), 2, fp) != 2)
			return(0);
	}
	return(1);
}

/* Reinhart/Tregenza sample generator */
static int
sample_reinhart(PARAMS *p, int b, FILE *fp)
{
#define T_NALT	7
	static const int	tnaz[T_NALT] = {30, 30, 24, 24, 18, 12, 6};
	const int		RowMax = T_NALT*p->hsiz + 1;
	const double		RAH = (.5*PI)/(RowMax-.5);
#define rnaz(r)			(r >= RowMax-1 ? 1 : p->hsiz*tnaz[r/p->hsiz])
	int			n = sampcnt;
	int			row, col;
	double			samp3[3];
	double			alt, azi;
	double			duvw[3];
	FVECT			orig_dir[2];
	int			i;

	if (fp == NULL) {		/* just requesting number of bins? */
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
		SDmultiSamp(samp3, 3, (n+frandom())/sampcnt);
		alt = (row+samp3[1])*RAH;
		azi = (2.*PI)*(col+samp3[2]-.5)/rnaz(row);
		duvw[2] = cos(alt);	/* measured from horizon */
		duvw[0] = tsin(azi)*duvw[2];
		duvw[1] = tcos(azi)*duvw[2];
		duvw[2] = sqrt(1. - duvw[2]*duvw[2]);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp3[0]))
			return(0);
		if (putbinary(orig_dir, sizeof(FVECT), 2, fp) != 2)
			return(0);
	}
	return(1);
#undef rnaz
#undef T_NALT
}

/* Klems sample generator */
static int
sample_klems(PARAMS *p, int b, FILE *fp)
{
	static const char	bname[4][20] = {
					"LBNL/Klems Full",
					"LBNL/Klems Half",
					"INTERNAL ERROR",
					"LBNL/Klems Quarter"
				};
	static ANGLE_BASIS	*kbasis[4];
	const int		bi = p->hemis[1] - '1';
	int			n = sampcnt;
	double			samp2[2];
	double			duvw[3];
	FVECT			orig_dir[2];
	int			i;

	if (!kbasis[bi]) {		/* need to get basis, first */
		for (i = 4; i--; )
			if (!strcasecmp(abase_list[i].name, bname[bi])) {
				kbasis[bi] = &abase_list[i];
				break;
			}
		if (i < 0) {
			fprintf(stderr, "%s: unknown hemisphere basis '%s'\n",
					progname, bname[bi]);
			return(0);
		}
	}
	if (fp == NULL)			/* just requesting number of bins? */
		return(kbasis[bi]->nangles);

	while (n--) {			/* stratified sampling */
		SDmultiSamp(samp2, 2, (n+frandom())/sampcnt);
		if (!fo_getvec(duvw, b+samp2[1], kbasis[bi]))
			return(0);
		for (i = 3; i--; )
			orig_dir[1][i] = -duvw[0]*p->udir[i] -
						duvw[1]*p->vdir[i] -
						duvw[2]*p->nrm[i] ;
		if (!sample_origin(p, orig_dir[0], orig_dir[1], samp2[0]))
			return(0);
		if (putbinary(orig_dir, sizeof(FVECT), 2, fp) != 2)
			return(0);
	}
	return(1);
}

/* Prepare hemisphere basis sampler that will send rays to rcontrib */
static int
prepare_sampler(void)
{
	if (curparams.slist == NULL) {	/* missing sample surface! */
		fputs(progname, stderr);
		fputs(": no sender surface!\n", stderr);
		return(-1);
	}
					/* misplaced output file spec. */
	if ((curparams.outfn != NULL) & (verbose >= 0))
		fprintf(stderr, "%s: warning - ignoring output file in sender ('%s')\n",
				progname, curparams.outfn);
					/* check/set basis hemisphere */
	if (!curparams.hemis[0]) {
		fputs(progname, stderr);
		fputs(": missing sender sampling type!\n", stderr);
		return(-1);
	}
	if (normalize(curparams.nrm) == 0) {
		fputs(progname, stderr);
		fputs(": undefined normal for sender sampling\n", stderr);
		return(-1);
	}
	if (normalize(curparams.vup) == 0) {
		if (fabs(curparams.nrm[2]) < .7)
			curparams.vup[2] = 1;
		else
			curparams.vup[1] = 1;
	}
	fcross(curparams.udir, curparams.vup, curparams.nrm);
	if (normalize(curparams.udir) == 0) {
		fputs(progname, stderr);
		fputs(": up vector coincides with sender normal\n", stderr);
		return(-1);
	}
	fcross(curparams.vdir, curparams.nrm, curparams.udir);
	if (curparams.sign == '-') {	/* left-handed coordinate system? */
		curparams.udir[0] *= -1.;
		curparams.udir[1] *= -1.;
		curparams.udir[2] *= -1.;
	}
	if (tolower(curparams.hemis[0]) == 'u' | curparams.hemis[0] == '1')
		curparams.sample_basis = sample_uniform;
	else if (tolower(curparams.hemis[0]) == 's' &&
				tolower(curparams.hemis[1]) == 'c')
		curparams.sample_basis = sample_shirchiu;
	else if ((tolower(curparams.hemis[0]) == 'r') |
			(tolower(curparams.hemis[0]) == 't'))
		curparams.sample_basis = sample_reinhart;
	else if (tolower(curparams.hemis[0]) == 'k') {
		switch (curparams.hemis[1]) {
		case '1':
		case '2':
		case '4':
			break;
		case 'f':
		case 'F':
		case '\0':
			curparams.hemis[1] = '1';
			break;
		case 'h':
		case 'H':
			curparams.hemis[1] = '2';
			break;
		case 'q':
		case 'Q':
			curparams.hemis[1] = '4';
			break;
		default:
			goto unrecognized;
		}
		curparams.hemis[2] = '\0';
		curparams.sample_basis = sample_klems;
	} else
		goto unrecognized;
					/* return number of bins */
	return((*curparams.sample_basis)(&curparams,0,NULL));
unrecognized:
	fprintf(stderr, "%s: unrecognized sender sampling: h=%s\n",
			progname, curparams.hemis);
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
	return(p->area > FTINY);
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
	if (snew == NULL) {
		fputs(progname, stderr);
		fputs(": out of memory in add_surface!\n", stderr);
		exit(1);
	}
	strncpy(snew->sname, oname, sizeof(snew->sname)-1);
	snew->sname[sizeof(snew->sname)-1] = '\0';
	snew->styp = st;
	snew->priv = NULL;
	snew->nfargs = n;
	for (n = 0; n < snew->nfargs; n++)
		if (fscanf(fp, "%lf", &snew->farg[n]) != 1) {
			fprintf(stderr, "%s: error reading arguments for '%s'\n",
					progname, oname);
			exit(1);
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
	if ((snew->area <= FTINY) & (verbose >= 0)) {
		fprintf(stderr, "%s: warning - zero area for surface '%s'\n",
				progname, oname);
		free(snew);
		return;
	}
	VSUM(curparams.nrm, curparams.nrm, snew->snrm, snew->area);
	snew->next = curparams.slist;
	curparams.slist = snew;
	curparams.nsurfs++;
	return;
badcount:
	fprintf(stderr, "%s: bad argument count for surface element '%s'\n",
			progname, oname);
	exit(1);
badnorm:
	fprintf(stderr, "%s: bad orientation for surface element '%s'\n",
			progname, oname);
	exit(1);
}

/* Parse a receiver object (look for modifiers to add to rcontrib command) */
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
		if (curparams.slist != NULL && (st == ST_SOURCE) ^
				(curparams.slist->styp == ST_SOURCE)) {
			fputs(progname, stderr);
			fputs(": cannot mix source/non-source receivers!\n", stderr);
			return(-1);
		}
		if (strcmp(thismod, curmod)) {
			if (curmod[0]) {	/* output last receiver? */
				finish_receiver();
				clear_params(&curparams, 1);
			}
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
			fputs(progname, stderr);
			fputs(": cannot use source as a sender!\n", stderr);
			return(-1);
		}
		if (strcmp(thismod, curmod)) {
			if (curmod[0]) {
				fputs(progname, stderr);
				fputs(": warning - multiple modifiers in sender\n",
						stderr);
			}
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
		fprintf(stderr, "%s: cannot load '%s'\n", progname, inspec);
		return(-1);
	}
	while ((c = getc(fp)) != EOF) {	/* load receiver data */
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
			if ((c = getc(fp)) == EOF || ungetc(c,fp) == EOF)
				break;
			if (!isspace(c) && fscanf(fp, "%s", inpbuf) == 1 &&
					!strcmp(inpbuf, PARAMSTART)) {
				if (fgets(inpbuf, sizeof(inpbuf), fp) != NULL)
					strcat(newparams, inpbuf);
				continue;
			}
			while ((c = getc(fp)) != EOF && c != '\n')
				;	/* else skipping comment */
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
		fprintf(stderr, "%s: error loading '%s'\n", progname, inspec);
		return(-1);
	}
	return(rv);
}

/* Get command arguments and run program */
int
main(int argc, char *argv[])
{
	char	fmtopt[6] = "-faa";	/* default output is ASCII */
	char	*xrs=NULL, *yrs=NULL, *ldopt=NULL;
	char	*iropt = NULL;
	char	*sendfn;
	char	sampcntbuf[32], nsbinbuf[32];
	FILE	*rcfp;
	int	nsbins;
	int	a, i;
					/* screen rcontrib options */
	progname = argv[0];
	for (a = 1; a < argc-2; a++) {
		int	na;
					/* check for argument expansion */
		while ((na = expandarg(&argc, &argv, a)) > 0)
			;
		if (na < 0) {
			fprintf(stderr, "%s: cannot expand '%s'\n",
					progname, argv[a]);
			return(1);
		}
		if (argv[a][0] != '-' || !argv[a][1])
			break;
		na = 1;	
		switch (argv[a][1]) {	/* !! Keep consistent !! */
		case 'v':		/* verbose mode */
			verbose = 1;
			na = 0;
			continue;
		case 'f':		/* special case for -fo, -ff, etc. */
			switch (argv[a][2]) {
			case '\0':		/* cal file */
				goto userr;
			case 'o':		/* force output */
				goto userr;
			case 'a':		/* output format */
			case 'f':
			case 'd':
			case 'c':
				if (!(fmtopt[3] = argv[a][3]))
					fmtopt[3] = argv[a][2];
				fmtopt[2] = argv[a][2];
				na = 0;
				continue;	/* will pass later */
			default:
				goto userr;
			}
			break;
		case 'x':		/* x-resolution */
			xrs = argv[++a];
			na = 0;
			continue;
		case 'y':		/* y-resolution */
			yrs = argv[++a];
			na = 0;
			continue;
		case 'c':		/* number of samples */
			sampcnt = atoi(argv[++a]);
			if (sampcnt <= 0)
				goto userr;
			na = 0;		/* we re-add this later */
			continue;
		case 'I':		/* only for pass-through mode */
		case 'i':
			iropt = argv[a];
			na = 0;
			continue;
		case 'w':		/* options without arguments */
			if (argv[a][2] != '+') verbose = -1;
		case 'V':
		case 'u':
		case 'h':
		case 'r':
			break;
		case 'n':		/* options with 1 argument */
		case 's':
		case 'o':
			na = 2;
			break;
		case 'b':		/* special case */
			if (argv[a][2] != 'v') goto userr;
			break;
		case 'l':		/* special case */
			if (argv[a][2] == 'd') {
				ldopt = argv[a];
				na = 0;
				continue;
			}
			na = 2;
			break;
		case 'd':		/* special case */
			if (argv[a][2] != 'v') na = 2;
			break;
		case 'a':		/* special case */
			if (argv[a][2] == 'p') {
				na = 2;	/* photon map [+ bandwidth(s)] */
				if (a < argc-3 && atoi(argv[a+1]) > 0)
					na += 1 + (a < argc-4 && atoi(argv[a+2]) > 0);
			} else
				na = (argv[a][2] == 'v') ? 4 : 2;
			break;
		case 'm':		/* special case */
			if (!argv[a][2]) goto userr;
			na = (argv[a][2] == 'e') | (argv[a][2] == 'a') ? 4 : 2;
			break;
		default:		/* anything else is verbotten */
			goto userr;
		}
		if (na <= 0) continue;
		CHECKARGC(na);		/* pass on option */
		rcarg[nrcargs++] = argv[a];
		while (--na)		/* + arguments if any */
			rcarg[nrcargs++] = argv[++a];
	}
	if (a > argc-2)
		goto userr;		/* check at end of options */
	sendfn = argv[a++];		/* assign sender & receiver inputs */
	if (sendfn[0] == '-') {		/* user wants pass-through mode? */
		if (sendfn[1]) goto userr;
		sendfn = NULL;
		if (iropt) {
			CHECKARGC(1);
			rcarg[nrcargs++] = iropt;
		}
		if (xrs) {
			CHECKARGC(2);
			rcarg[nrcargs++] = "-x";
			rcarg[nrcargs++] = xrs;
		}
		if (yrs) {
			CHECKARGC(2);
			rcarg[nrcargs++] = "-y";
			rcarg[nrcargs++] = yrs;
		}
		if (ldopt) {
			CHECKARGC(1);
			rcarg[nrcargs++] = ldopt;
		}
		if (sampcnt <= 0) sampcnt = 1;
	} else {			/* else in sampling mode */
		if (iropt) {
			fputs(progname, stderr);
			fputs(": -i, -I supported for pass-through only\n", stderr);
			return(1);
		}
		fmtopt[2] = (sizeof(RREAL)==sizeof(double)) ? 'd' : 'f';
		if (sampcnt <= 0) sampcnt = 10000;
	}
	sprintf(sampcntbuf, "%d", sampcnt);
	CHECKARGC(3);			/* add our format & sample count */
	rcarg[nrcargs++] = fmtopt;
	rcarg[nrcargs++] = "-c";
	rcarg[nrcargs++] = sampcntbuf;
					/* add receiver arguments to rcontrib */
	if (load_scene(argv[a], add_recv_object) < 0)
		return(1);
	finish_receiver();
	if (sendfn == NULL) {		/* pass-through mode? */
		CHECKARGC(1);		/* add octree */
		rcarg[nrcargs++] = oconv_command(argc-a, argv+a);
		rcarg[nrcargs] = NULL;
		return(my_exec(rcarg));	/* rcontrib does everything */
	}
	clear_params(&curparams, 0);	/* else load sender surface & params */
	curmod[0] = '\0';
	if (load_scene(sendfn, add_send_object) < 0)
		return(1);
	if ((nsbins = prepare_sampler()) <= 0)
		return(1);
	CHECKARGC(3);			/* add row count and octree */
	rcarg[nrcargs++] = "-y";
	sprintf(nsbinbuf, "%d", nsbins);
	rcarg[nrcargs++] = nsbinbuf;
	rcarg[nrcargs++] = oconv_command(argc-a, argv+a);
	rcarg[nrcargs] = NULL;
					/* open pipe to rcontrib process */
	if ((rcfp = popen_arglist(rcarg, "w")) == NULL)
		return(1);
	SET_FILE_BINARY(rcfp);
#ifdef getc_unlocked
	flockfile(rcfp);
#endif
	if (verbose > 0) {
		fprintf(stderr, "%s: sampling %d directions", progname, nsbins);
		if (curparams.nsurfs > 1)
			fprintf(stderr, " (%d elements)\n", curparams.nsurfs);
		else
			fputc('\n', stderr);
	}
	for (i = 0; i < nsbins; i++)	/* send rcontrib ray samples */
		if (!(*curparams.sample_basis)(&curparams, i, rcfp))
			return(1);
	return(pclose(rcfp) < 0);	/* all finished! */
userr:
	if (a < argc-2)
		fprintf(stderr, "%s: unsupported option '%s'", progname, argv[a]);
	fprintf(stderr, "Usage: %s [-v][rcontrib options] sender.rad receiver.rad [-i system.oct] [system.rad ..]\n",
				progname);
	return(1);
}
