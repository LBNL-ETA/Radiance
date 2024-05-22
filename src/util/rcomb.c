#ifndef lint
static const char RCSid[] = "$Id: rcomb.c,v 2.14 2024/05/22 15:38:04 greg Exp $";
#endif
/*
 * General component matrix combiner, operating on a row at a time.
 */

#include <signal.h>
#include <math.h>
#include "platform.h"
#include "rtprocess.h"
#include "rtio.h"
#include "rmatrix.h"
#include "calcomp.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define MAXCOMP		MAXCSAMP	/* #components we support */

/* Unary matrix operation(s) */
typedef struct {
	double		cmat[MAXCOMP*MAXCOMP];	/* component transformation */
	double		sca[MAXCOMP];		/* scalar coefficients */
	const char	*csym;			/* symbolic coefficients */
	short		clen;			/* number of coefficients */
	short		nsf;			/* number of scalars */
} RUNARYOP;

/* Input matrix */
typedef struct {
	const char	*inspec;		/* input specification */
	RUNARYOP	preop;			/* transform operation */
	RMATRIX		imx;			/* input matrix header info */
	RMATRIX		*rmp;			/* active single-row matrix */
	FILE		*infp;			/* open input stream */
} ROPMAT;

ROPMAT		*mop = NULL;			/* allocated input array */
int		nall = 0;			/* number allocated */
int		nmats = 0;			/* number of actual inputs */

RMATRIX		*mcat = NULL;			/* final concatenation */
int		mcat_last = 0;			/* goes after trailing ops? */

int		in_nrows;			/* number of input rows (or 0) */
#define in_ncols	(mop[0].rmp->ncols)	/* number of input columns */
#define in_ncomp	(mop[0].rmp->ncomp)	/* input #components */

extern int	nowarn;				/* turn off warnings? */

int		cur_row;			/* current input/output row */
int		cur_col;			/* current input/output column */
int		cur_chan;			/* if we're looping channels */

SUBPROC		*cproc = NULL;			/* child process array */
int		nchildren = 0;			/* # of child processes */
int		inchild = -1;			/* our child ID (-1: parent) */
int		pgid = -1;			/* process group ID */
int		nr_out = 0;			/* # of rows output by kids */

static int	checksymbolic(ROPMAT *rop);

static void
on_sigio(int dummy)
{
	nr_out++;			/* happens when child outputs row */
}

static int
split_input(ROPMAT *rop)
{
	if (rop->rmp == &rop->imx && !(rop->rmp = rmx_copy(&rop->imx))) {
		fputs("Out of memory in split_input()\n", stderr);
		return(0);
	}
	rmx_reset(rop->rmp);
	return(1);
}

/* Check/set transform based on a reference input file */
static int
checkreffile(ROPMAT *rop)
{
	static const char	*curRF = NULL;
	static RMATRIX		refm;
	const int		nc = rop->imx.ncomp;
	int			i;

	if (!curRF || strcmp(rop->preop.csym, curRF)) {
		FILE	*fp = fopen(rop->preop.csym, "rb");
		if (!rmx_load_header(&refm, fp)) {
			fprintf(stderr, "%s: cannot read info header\n",
					rop->preop.csym);
			curRF = NULL;
			if (fp) fclose(fp);
			return(0);
		}
		fclose(fp);
		curRF = rop->preop.csym;
	}
	if (refm.ncomp == 3) {
		rop->preop.csym = (refm.dtype == DTxyze) ? "XYZ" : "RGB";
		return(checksymbolic(rop));
	}
	if (refm.ncomp == 2) {
		fprintf(stderr, "%s: cannot convert to 2 components\n",
				curRF);
		return(0);
	}
	if (refm.ncomp == 1) {
		rop->preop.csym = "Y";		/* XXX big assumption */
		return(checksymbolic(rop));
	}
	if (refm.ncomp == nc &&
			!memcmp(refm.wlpart, rop->imx.wlpart, sizeof(refm.wlpart)))
		return(1);			/* nothing to do */

	if ((nc <= 3) | (nc > MAXCSAMP) | (refm.ncomp > MAXCSAMP)) {
		fprintf(stderr, "%s: cannot resample from %d to %d components\n",
				curRF, nc, refm.ncomp);
		return(0);
	}
	if (!split_input(rop))			/* get our own struct */
		return(0);
	rop->preop.clen = refm.ncomp * nc;	/* compute spec to ref */

	for (i = 0; i < nc; i++) {
		SCOLOR	scstim, scresp;
		int	j;
		memset(scstim, 0, sizeof(COLORV)*nc);
		scstim[i] = 1.f;
		convertscolor(scresp, refm.ncomp, refm.wlpart[0], refm.wlpart[3],
				scstim, nc, rop->imx.wlpart[0], rop->imx.wlpart[3]);
		for (j = refm.ncomp; j-- > 0; )
			rop->preop.cmat[j*nc + i] = scresp[j];
	}
						/* remember new spectral params */
	memcpy(rop->rmp->wlpart, refm.wlpart, sizeof(rop->rmp->wlpart));
	rop->rmp->ncomp = refm.ncomp;
	return(1);
}

/* Compute conversion row from spectrum to one channel of RGB */
static void
rgbrow(ROPMAT *rop, int r, int p)
{
	const int	nc = rop->imx.ncomp;
	const float *	wlp = rop->imx.wlpart;
	int		i;

	for (i = nc; i--; ) {
		int	nmEnd = wlp[0] + (wlp[3] - wlp[0])*i/nc;
		int	nmStart = wlp[0] + (wlp[3] - wlp[0])*(i+1)/nc;
		COLOR	crgb;
		spec_rgb(crgb, nmStart, nmEnd);
		rop->preop.cmat[r*nc+i] = crgb[p];
	}
}

/* Compute conversion row from spectrum to one channel of XYZ */
static void
xyzrow(ROPMAT *rop, int r, int p)
{
	const int	nc = rop->imx.ncomp;
	const float *	wlp = rop->imx.wlpart;
	int		i;

	for (i = nc; i--; ) {
		int	nmEnd = wlp[0] + (wlp[3] - wlp[0])*i/nc;
		int	nmStart = wlp[0] + (wlp[3] - wlp[0])*(i+1)/nc;
		COLOR	cxyz;
		spec_cie(cxyz, nmStart, nmEnd);
		rop->preop.cmat[r*nc+i] = cxyz[p];
	}
}

/* Use the spectral sensitivity function to compute matrix coefficients */
static void
sensrow(ROPMAT *rop, int r, double (*sf)(SCOLOR sc, int ncs, const float wlpt[4]))
{
	const int	nc = rop->imx.ncomp;
	int		i;

	for (i = nc; i--; ) {
		SCOLOR	sclr;
		memset(sclr, 0, sizeof(COLORV)*nc);
		sclr[i] = 1.f;
		rop->preop.cmat[r*nc+i] = (*sf)(sclr, nc, rop->imx.wlpart);
	}
}

/* Check/set symbolic transform */
static int
checksymbolic(ROPMAT *rop)
{
	const int	nc = rop->imx.ncomp;
	const int	dt = rop->imx.dtype;
	double		cf = 1;
	int		i, j;
					/* check suffix => reference file */
	if (strchr(rop->preop.csym, '.') > rop->preop.csym)
		return(checkreffile(rop));

	if (nc < 3) {
		fprintf(stderr, "%s: -c '%s' requires at least 3 components\n",
				rop->inspec, rop->preop.csym);
		return(0);
	}
	rop->preop.clen = strlen(rop->preop.csym) * nc;
	if (rop->preop.clen > MAXCOMP*MAXCOMP) {
		fprintf(stderr, "%s: -c '%s' results in too many components\n",
				rop->inspec, rop->preop.csym);
		return(0);
	}
	for (j = 0; rop->preop.csym[j]; j++) {
		int	comp = 0;
		switch (rop->preop.csym[j]) {
		case 'B':
		case 'b':
			++comp;
			/* fall through */
		case 'G':
		case 'g':
			++comp;
			/* fall through */
		case 'R':
		case 'r':
			if (rop->preop.csym[j] <= 'Z')
				cf = 1./WHTEFFICACY;
			if (dt == DTxyze) {
				for (i = 3; i--; )
					rop->preop.cmat[j*nc+i] = cf*xyz2rgbmat[comp][i];
			} else if (nc == 3)
				rop->preop.cmat[j*nc+comp] = 1.;
			else
				rgbrow(rop, j, comp);
			break;
		case 'Z':
		case 'z':
			++comp;
			/* fall through */
		case 'Y':
		case 'y':
			++comp;
			/* fall through */
		case 'X':
		case 'x':
			if ((rop->preop.csym[j] <= 'Z') & (dt != DTxyze))
				cf = WHTEFFICACY;
			if (dt == DTxyze) {
				rop->preop.cmat[j*nc+comp] = 1.;
			} else if (nc == 3) {
				for (i = 3; i--; )
					rop->preop.cmat[j*nc+i] =
							rgb2xyzmat[comp][i];
			} else if (comp == CIEY)
				sensrow(rop, j, scolor2photopic);
			else
				xyzrow(rop, j, comp);

			for (i = nc*(cf != 1); i--; )
				rop->preop.cmat[j*nc+i] *= cf;
			break;
		case 'S':		/* scotopic (il)luminance */
			cf = WHTSCOTOPIC;
			/* fall through */
		case 's':
			sensrow(rop, j, scolor2scotopic);
			for (i = nc*(cf != 1); i--; )
				rop->preop.cmat[j*nc+i] *= cf;
			break;
		case 'M':		/* melanopic (il)luminance */
			cf = WHTMELANOPIC;
			/* fall through */
		case 'm':
			sensrow(rop, j, scolor2melanopic);
			for (i = nc*(cf != 1); i--; )
				rop->preop.cmat[j*nc+i] *= cf;
			break;
		case 'A':		/* average component */
		case 'a':
			for (i = nc; i--; )
				rop->preop.cmat[j*nc+i] = 1./(double)nc;
			break;
		default:
			fprintf(stderr, "%s: -c '%c' unsupported\n",
				rop->inspec, rop->preop.csym[j]);
			return(0);
		}
	}
	if (!split_input(rop))		/* get our own struct */
		return(0);
	memcpy(rop->rmp->wlpart, WLPART, sizeof(rop->rmp->wlpart));
	rop->rmp->ncomp = rop->preop.clen / nc;
					/* decide on output type */
	if (!strcasecmp(rop->preop.csym, "XYZ")) {
		if (dt <= DTspec)
			rop->rmp->dtype = DTxyze;
	} else if (!strcasecmp(rop->preop.csym, "RGB")) {
		if (dt <= DTspec)
			rop->rmp->dtype = DTrgbe;
	} else if (rop->rmp->dtype == DTspec)
		rop->rmp->dtype = DTfloat;
	return(1);
}

static int
get_component_xfm(ROPMAT *rop)
{
	int	i, j;

	if (rop->rmp != &rop->imx) {		/* reset destination matrix */
		rmx_free(rop->rmp);
		rop->rmp = &rop->imx;
	}
	if (rop->preop.csym &&			/* symbolic transform? */
			!checksymbolic(rop))
		return(0);
						/* undo exposure? */
	if (fabs(1. - bright(rop->rmp->cexp)) > .025) {
		if (rop->rmp->ncomp == 1)
			rop->rmp->cexp[RED] = rop->rmp->cexp[GRN] =
					rop->rmp->cexp[BLU] = bright(rop->rmp->cexp);
		if (rop->preop.nsf <= 0) {
			rop->preop.nsf = i = rop->rmp->ncomp;
			while (i--)
				rop->preop.sca[i] = 1.;
		}
		if (rop->preop.nsf == 1) {
			if (rop->rmp->ncomp == 3) {
				rop->preop.sca[2] = rop->preop.sca[1] =
						rop->preop.sca[0];
				rop->preop.nsf = 3;
			} else
				rop->preop.sca[0] /= bright(rop->rmp->cexp);
		}
		if (rop->preop.nsf == 3) {
			opcolor(rop->preop.sca, /=, rop->rmp->cexp);
		} else if (rop->preop.nsf > 3) {	/* punt */
			double	mult = 1./bright(rop->rmp->cexp);
			for (i = rop->preop.nsf; i--; )
				rop->preop.sca[i] *= mult;
		}
		setcolor(rop->rmp->cexp, 1., 1., 1.);
	}
	if (rop->preop.clen > 0) {		/* use component transform? */
		if (rop->preop.clen % rop->imx.ncomp) {
			fprintf(stderr, "%s: -c must have N x %d coefficients\n",
					rop->inspec, rop->imx.ncomp);
			return(0);
		}
		if (rop->preop.nsf > 0) {	/* scale transform, instead */
			if (rop->preop.nsf == 1) {
				for (i = rop->preop.clen; i--; )
					rop->preop.cmat[i] *= rop->preop.sca[0];
			} else if (rop->preop.nsf*rop->imx.ncomp != rop->preop.clen) {
				fprintf(stderr, "%s: -s must have one or %d factors\n",
						rop->inspec,
						rop->preop.clen/rop->imx.ncomp);
				return(0);
			} else {
				for (i = rop->preop.nsf; i--; )
					for (j = rop->imx.ncomp; j--; )
						rop->preop.cmat[i*rop->imx.ncomp+j]
								*= rop->preop.sca[i];
			}
		}
		rop->preop.nsf = 0;		/* now folded in */
		if (!split_input(rop))		/* get our own struct */
			return(0);
		rop->rmp->ncomp = rop->preop.clen / rop->imx.ncomp;
		if ((rop->rmp->ncomp > 3) & (rop->rmp->dtype <= DTspec)) {
			rop->rmp->dtype = DTfloat;	/* probably not actual spectrum */
			memcpy(rop->rmp->wlpart, WLPART, sizeof(rop->rmp->wlpart));
		}
	} else if (rop->preop.nsf > 0) {	/* else use scalar(s)? */
		if (rop->preop.nsf == 1) {
			for (i = rop->rmp->ncomp; --i; )
				rop->preop.sca[i] = rop->preop.sca[0];
			rop->preop.nsf = rop->rmp->ncomp;
		} else if (rop->preop.nsf != rop->rmp->ncomp) {
			fprintf(stderr, "%s: -s must have one or %d factors\n",
					rop->inspec, rop->rmp->ncomp);
			return(0);
		}
	}
	return(1);
}

static int
apply_op(RMATRIX *dst, const RMATRIX *src, const RUNARYOP *ro)
{
	if (ro->clen > 0) {
		RMATRIX	*res = rmx_transform(src, dst->ncomp, ro->cmat);
		if (!res) {
			fputs("Error in call to rmx_transform()\n", stderr);
			return(0);
		}
		if (!rmx_transfer_data(dst, res, 0))
			return(0);
		rmx_free(res);
	} else if (dst != src)
		memcpy(dst->mtx, src->mtx, rmx_array_size(dst));
	if (ro->nsf == dst->ncomp)
		rmx_scale(dst, ro->sca);
	return(1);
}

static int
open_input(ROPMAT *rop)
{
	int	outtype;

	if (!rop || !rop->inspec || !rop->inspec[0])
		return(0);
	if (rop->inspec == stdin_name)
		rop->infp = stdin;
	else if (rop->inspec[0] == '!')
		rop->infp = popen(rop->inspec+1, "r");
	else
		rop->infp = fopen(rop->inspec, "rb");

	if (!rmx_load_header(&rop->imx, rop->infp)) {
		fprintf(stderr, "Bad header from: %s\n", rop->inspec);
		return(0);
	}
	return(get_component_xfm(rop));
}

/* Return nominal wavelength associated with input component (return nm) */
static double
l_wavelength(char *nam)
{
	double	comp = argument(1);

	if ((comp < -.5) | (comp >= in_ncomp+.5)) {
		errno = EDOM;
		return(.0);
	}
	if (comp < .5)				/* asking for #components? */
		return(in_ncomp);

	if (in_ncomp == 3) {			/* special case for RGB */
		const int	w0 = (int)(comp - .5);
		return(mop[0].rmp->wlpart[w0] +
				(comp-.5)*(mop[0].rmp->wlpart[w0+1] -
					mop[0].rmp->wlpart[w0]));
	}
	return(mop[0].rmp->wlpart[0] +		/* general case, even div. */
		(comp-.5)/(double)in_ncomp *
			(mop[0].rmp->wlpart[3] - mop[0].rmp->wlpart[0]));
}

/* Return ith input with optional channel selector */
static double
l_chanin(char *nam)
{
	double	inp = argument(1);
	int	mi, chan;

	if ((mi = (int)(inp-.5)) < 0 || mi >= nmats) {
		errno = EDOM;
		return(.0);
	}
	if (inp < .5)			/* asking for #inputs? */
		return(nmats);

	if (nargum() >= 2) {
		double	cval = argument(2);
		if (cval < .5 || (chan = (int)(cval-.5)) >= in_ncomp) {
			errno = EDOM;
			return(.0);
		}
	} else
		chan = cur_chan;

	return(mop[mi].rmp->mtx[cur_col*in_ncomp + chan]);
}

static int
initialize(RMATRIX *imp)
{
	int	i;
					/* XXX struct is zeroed coming in */
	setcolor(imp->cexp, 1.f, 1.f, 1.f);
	for (i = 0; i < nmats; i++) {	/* open each input */
		int	restype;
		if (!open_input(&mop[i]))
			return(0);
		restype = mop[i].rmp->dtype;
		if (!imp->dtype || (restype = rmx_newtype(restype, imp->dtype)) > 0)
			imp->dtype = restype;
		else
			fprintf(stderr, "%s: warning - data type mismatch\n",
					mop[i].inspec);
		if (!i) {
			imp->ncols = mop[0].rmp->ncols;
			imp->ncomp = mop[0].rmp->ncomp;
			memcpy(imp->wlpart, mop[0].rmp->wlpart, sizeof(imp->wlpart));
		} else if ((mop[i].rmp->ncols != imp->ncols) |
				(mop[i].rmp->ncomp != imp->ncomp) |
				((in_nrows > 0) & (mop[i].rmp->nrows > 0) &
					(mop[i].rmp->nrows != in_nrows))) {
			fprintf(stderr, "%s: mismatch in size or #components\n",
					mop[i].inspec);
			return(0);
		}			/* XXX should check wlpart? */
		if (in_nrows <= 0)
			in_nrows = imp->nrows = mop[i].rmp->nrows;
	}				/* set up .cal environment */
	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_OUTCHAN|E_INCHAN);
	varset("PI", ':', M_PI);
	varset("nfiles", ':', nmats);
	varset("nrows", ':', in_nrows);
	varset("ncols", ':', in_ncols);
	varset("ncomp", ':', in_ncomp);
	varset("R", ':', 1.);
	varset("G", ':', 2.);
	varset("B", ':', 3.);
	funset("wl", 1, ':', l_wavelength);
	funset("ci", 1, '=', l_chanin);
	scompile("ri(i)=ci(i,R);gi(i)=ci(i,G);bi(i)=ci(i,B)", NULL, 0);
	return(1);
}

static void
output_headinfo(FILE *fp)
{
	int	i;

	for (i = 0; i < nmats; i++) {
		const char	*cp = mop[i].imx.info;
		fputs(mop[i].inspec, fp);
		fputs(":\n", fp);
		if (!cp) continue;
		while (*cp) {
			if (*cp == '\n') {
				cp++;		/* avoid inadvertant terminus */
				continue;
			}
			fputc('\t', fp);	/* indent this input's info */
			do
				putc(*cp, fp);
			while (*cp++ != '\n');
		}
	}
}

static int
spawned_children(int np)
{
	int	i, rv;

#if defined(_WIN32) || defined(_WIN64)
	if (np > 1) {
		fputs("Warning: only one process under Windows\n", stderr);
		np = 1;
	} else
#endif
	if ((in_nrows > 0) & (np > in_nrows))
		np = in_nrows;
				/* we'll be doing a row at a time */
	for (i = 0; i < nmats; i++) {
		mop[i].imx.nrows = 1;
		if (!rmx_prepare(&mop[i].imx))
			goto memerror;
		if (mop[i].rmp != &mop[i].imx) {
			mop[i].rmp->nrows = 1;
			if (!rmx_prepare(mop[i].rmp))
				goto memerror;
		}
	}
				/* prep output row buffer */
	if (mcat || mop[nmats].preop.clen > 0) {
		if (!split_input(&mop[nmats]))	/* need separate buffer */
			return(0);
		if (mop[nmats].preop.clen > 0)
			mop[nmats].rmp->ncomp = mop[nmats].preop.clen /
						mop[nmats].imx.ncomp;
		mop[nmats].rmp->nrows = 1;
		if (!mcat | !mcat_last && !rmx_prepare(mop[nmats].rmp))
			goto memerror;
	}
	mop[nmats].imx.nrows = 1;
	if (!rmx_prepare(&mop[nmats].imx))
		goto memerror;
	if (np <= 1) {		/* single process return point */
#ifdef getc_unlocked
		for (i = 0; i < nmats; i++)
			flockfile(mop[i].infp);
		flockfile(stdout);
#endif
		return(0);
	}
	pgid = setpgrp();	/* set process group ID */
	signal(SIGIO, on_sigio);
	fflush(stdout);		/* flush header & spawn children */
	cproc = (SUBPROC *)malloc(sizeof(SUBPROC)*np);
	if (!cproc)
		goto memerror;
	nchildren = np;
	for (i = 0; i < np; i++) {
		cproc[i].flags = PF_FILT_OUT;
		cproc[i].w = dup(1);
		cproc[i].r = 0;
		cproc[i].pid = -1;
		rv = open_process(&cproc[i], NULL);
		if (rv <= 0) break;
	}
	if (rv > 0)
		return(1);	/* parent return value */
	if (rv < 0) {
		perror("fork");
		exit(1);
	}
	inchild = i;		/* our child index */
	while (i-- > 0)		/* don't share siblings' pipes */
		close(cproc[i].w);
	fpurge(stdin);		/* discard previous matrix input */
#ifdef getc_unlocked
	flockfile(stdin);
#endif
	for (i = 0; i < nmats; i++) {
		if (mop[i].infp != stdin)
			fclose(mop[i].infp);	/* ! pclose() */
		mop[i].infp = stdin;
		mop[i].imx.dtype = DTdouble;
	}
	return(0);		/* child return */
memerror:
	fputs("Out of memory in spawned_children()\n", stderr);
	exit(1);
}

static int
parent_loop()
{
	FILE	**outfp = (FILE **)malloc(nchildren*sizeof(FILE *));
	int	i;

	if (!outfp) goto memerror;
	for (i = 0; i < nchildren; i++) {
		outfp[i] = fdopen(cproc[i].w, "w");
		if (!outfp[i]) goto memerror;
#ifdef getc_unlocked
		flockfile(outfp[i]);
#endif
	}
#ifdef getc_unlocked
	for (i = 0; i < nmats; i++)
		flockfile(mop[i].infp);
#endif
	for (cur_row = 0; (in_nrows <= 0) | (cur_row < in_nrows); cur_row++) {
	    FILE	*ofp = outfp[cur_row % nchildren];
	    for (i = 0; i < nmats; i++)
		if (!rmx_load_row(mop[i].imx.mtx, &mop[i].imx, mop[i].infp)) {
			if (cur_row > in_nrows)	/* unknown #input rows? */
				break;
			fprintf(stderr, "%s: read error at row %d\n",
					mop[i].inspec, cur_row);
			return(0);
		}
	    if (i < nmats)
		break;
	    for (i = 0; i < nmats; i++)
	    	if (!rmx_write_data(mop[i].imx.mtx, mop[i].imx.ncomp,
	    				mop[i].imx.ncols, DTdouble, ofp))
	    		return(0);
	    if (fflush(ofp) == EOF)
	    	return(0);
	}
	for (i = 0; i < nchildren; i++)
		fclose(outfp[i]);
	free(outfp);
	i = close_processes(cproc, nchildren);
	free(cproc); cproc = NULL;
	if (i < 0) {
		fputs("Warning: missing child in parent_loop()\n", stderr);
		return(1);
	}
	if (i > 0) {
		fprintf(stderr, "Child exited with status %d\n", i);
		return(0);
	}
	return(1);
memerror:
	fputs("Out of memory in parent_loop()\n", stderr);
	exit(1);
}

static int
combine_input()
{
	const int	row0 = (inchild >= 0)*inchild;
	const int	rstep = nchildren + !nchildren;
	ROPMAT		*res = &mop[nmats];
	int		set_r, set_c;
	RMATRIX		*tmp = NULL;
	int		co_set;
	sigset_t	iomask;
	int		i;

	if (mcat && mcat_last &&
			!(tmp = rmx_alloc(1, res->imx.ncols, res->rmp->ncomp)))
		goto memerror;
					/* figure out what the user set */
	co_set = fundefined("co");
	if (!co_set)
		co_set = -vardefined("co");
	if (!co_set & (in_ncomp == 3) && vardefined("ro") &&
			vardefined("go") && vardefined("bo")) {
		scompile("co(p)=select(p,ro,go,bo)", NULL, 0);
		co_set = 1;
	}
	if (co_set) {			/* set if user wants, didn't set */
		set_r = varlookup("r") != NULL && !vardefined("r");
		set_c = varlookup("c") != NULL && !vardefined("c");
	} else				/* save a little time */
		set_r = set_c = 0;

	sigemptyset(&iomask);		/* read/process row-by-row */
	sigaddset(&iomask, SIGIO);
	for (cur_row = row0; (in_nrows <= 0) | (cur_row < in_nrows); cur_row += rstep) {
	    RMATRIX	*mres = NULL;
	    if (inchild >= 0) sigprocmask(SIG_BLOCK, &iomask, NULL);
	    for (i = 0; i < nmats; i++)
		if (!rmx_load_row(mop[i].imx.mtx, &mop[i].imx, mop[i].infp)) {
			if (cur_row > in_nrows)	/* unknown #input rows? */
				break;
			fprintf(stderr, "%s: read error at row %d\n",
					mop[i].inspec, cur_row);
			return(0);
		}
	    if (inchild >= 0) sigprocmask(SIG_UNBLOCK, &iomask, NULL);
	    if (i < nmats)
	    	break;
	    for (i = 0; i < nmats; i++)
		if (!apply_op(mop[i].rmp, &mop[i].imx, &mop[i].preop))
			return(0);
	    if (set_r) varset("r", '=', cur_row);
	    for (cur_col = 0; cur_col < in_ncols; cur_col++) {
	    	if (set_c) varset("c", '=', cur_col);
		for (cur_chan = 0; cur_chan < in_ncomp; cur_chan++) {
		    const int	ndx = cur_col*in_ncomp + cur_chan;
		    eclock++;
		    if (!co_set) {	/* just summing elements? */
		    	res->imx.mtx[ndx] = 0;
		    	for (i = nmats; i--; )
		    		res->imx.mtx[ndx] += mop[i].rmp->mtx[ndx];
		    } else if (co_set > 0) {
		    	double	dchan = cur_chan+1;
		    	res->imx.mtx[ndx] = funvalue("co", 1, &dchan);
		    } else
			res->imx.mtx[ndx] = varvalue("co");
		}
	    }				/* final conversions */
	    if (!mcat) {
	    	if (!apply_op(res->rmp, &res->imx, &res->preop))
			return(0);
	    } else if (mcat_last) {
	    	if (!apply_op(tmp, &res->imx, &res->preop))
	    		return(0);
	    	mres = rmx_multiply(tmp, mcat);
	    	if (!mres)
	    		goto multerror;
		if (!rmx_transfer_data(res->rmp, mres, 0))
			return(0);
	    } else /* mcat && !mcat_last */ {
	    	mres = rmx_multiply(&res->imx, mcat);
	    	if (!mres)
	    		goto multerror;
		if (!apply_op(res->rmp, mres, &res->preop))
			return(0);
	    }
	    rmx_free(mres); mres = NULL;
	    if (inchild >= 0) {		/* children share stdout */
	    	while (nr_out < cur_row)
		    pause();		/* wait for our turn */
		sigprocmask(SIG_BLOCK, &iomask, NULL);
	    }
	    if (!rmx_write_data(res->rmp->mtx, res->rmp->ncomp,
	    			res->rmp->ncols, res->rmp->dtype, stdout))
	    	return(0);
	    if (inchild >= 0) {		/* flush and notify group */
	    	if (fflush(stdout) == EOF)
		    return(0);
		sigprocmask(SIG_UNBLOCK, &iomask, NULL);
		killpg(pgid, SIGIO);	/* increments everyone's nr_out */
	    }
	}
	return(inchild >= 0 || fflush(stdout) != EOF);
memerror:
	fputs("Out of buffer space in combine_input()\n", stderr);
	return(0);
multerror:
	fputs("Unexpected matrix multiply error in combine_input()\n", stderr);
	return(0);
}

static int
get_factors(double da[], int n, char *av[])
{
	int	ac;

	for (ac = 0; ac < n && isflt(av[ac]); ac++)
		da[ac] = atof(av[ac]);
	return(ac);
}

static void
resize_inparr(int n2alloc)
{
	int	i;

	if (n2alloc == nall)
		return;
	for (i = nall; i > n2alloc; i--) {
		rmx_reset(&mop[i].imx);
		if (mop[i].rmp != &mop[i].imx)
			rmx_free(mop[i].rmp);
	}
	mop = (ROPMAT *)realloc(mop, n2alloc*sizeof(ROPMAT));
	if (mop == NULL) {
		fputs("Out of memory in resize_inparr()\n", stderr);
		exit(1);
	}
	if (n2alloc > nall)
		memset(mop+nall, 0, (n2alloc-nall)*sizeof(ROPMAT));
	nall = n2alloc;
}

/* Load one or more matrices and operate on them, sending results to stdout */
int
main(int argc, char *argv[])
{

	int		outfmt = DTfromHeader;
	const char	*defCsym = NULL;
	int		echoheader = 1;
	int		stdin_used = 0;
	int		nproc = 1;
	const char	*mcat_spec = NULL;
	int		n2comp = 0;
	uby8		comp_ndx[128];
	int		i;
					/* get starting input array */
	mop = (ROPMAT *)calloc(nall=2, sizeof(ROPMAT));
					/* get options and arguments */
	for (i = 1; i < argc; i++)
		if (argv[i][0] != '-' || !argv[i][1]) {
			if (argv[i][0] == '-') {
				if (stdin_used++) goto stdin_error;
				mop[nmats].inspec = stdin_name;
			} else
				mop[nmats].inspec = argv[i];
			if (!mop[nmats].preop.csym)
				mop[nmats].preop.csym = defCsym;
			if (++nmats >= nall)
				resize_inparr(nmats + (nmats>>2) + 2);
		} else {
			int	n = argc-1 - i;
			switch (argv[i][1]) {	/* get option */
			case 'w':
				nowarn = !nowarn;
				break;
			case 'h':
				echoheader = !echoheader;
				break;
			case 'n':
				nproc = atoi(argv[++i]);
				if (nproc <= 0)
					goto userr;
				break;
			case 'e':
				if (!n) goto userr;
				comp_ndx[n2comp++] = i++;
				break;
			case 'f':
				switch (argv[i][2]) {
				case '\0':
					if (!n) goto userr;
					comp_ndx[n2comp++] = i++;
					break;
				case 'd':
					outfmt = DTdouble;
					break;
				case 'f':
					outfmt = DTfloat;
					break;
				case 'a':
					outfmt = DTascii;
					break;
				case 'c':
					outfmt = DTrgbe;
					break;
				default:
					goto userr;
				}
				break;
			case 's':
				if (n > MAXCOMP) n = MAXCOMP;
				i += mop[nmats].preop.nsf =
					get_factors(mop[nmats].preop.sca,
							n, argv+i+1);
				if (mop[nmats].preop.nsf <= 0) {
					fprintf(stderr, "%s: -s missing arguments\n",
							argv[0]);
					goto userr;
				}
				break;
			case 'C':
				if (!n || isflt(argv[i+1]))
					goto userr;
				defCsym = mop[nmats].preop.csym = argv[++i];
				mop[nmats].preop.clen = 0;
				mcat_last = 0;
				break;
			case 'c':
				if (n && !isflt(argv[i+1])) {
					mop[nmats].preop.csym = argv[++i];
					mop[nmats].preop.clen = 0;
					break;
				}
				if (n > MAXCOMP*MAXCOMP) n = MAXCOMP*MAXCOMP;
				i += mop[nmats].preop.clen =
					get_factors(mop[nmats].preop.cmat,
							n, argv+i+1);
				if (mop[nmats].preop.clen <= 0) {
					fprintf(stderr, "%s: -c missing arguments\n",
							argv[0]);
					goto userr;
				}
				mop[nmats].preop.csym = NULL;
				mcat_last = 0;
				break;
			case 'm':
				if (!n) goto userr;
				if (argv[++i][0] == '-' && !argv[i][1]) {
					if (stdin_used++) goto stdin_error;
					mcat_spec = stdin_name;
				} else
					mcat_spec = argv[i];
				mcat_last = 1;
				break;
			default:
				fprintf(stderr, "%s: unknown option '%s'\n",
						argv[0], argv[i]);
				goto userr;
			}
		}
	if (!nmats) {
		fprintf(stderr, "%s: need at least one input matrix\n", argv[0]);
		goto userr;
	}
	resize_inparr(nmats+1);		/* extra matrix at end for result */
	mop[nmats].inspec = "trailing_ops";
					/* load final concatenation matrix */
	if (mcat_spec && !(mcat = rmx_load(mcat_spec, RMPnone))) {
		fprintf(stderr, "%s: error loading concatenation matrix: %s\n",
				argv[0], mcat_spec);
		return(1);
	}
					/* get/check inputs, set constants */
	if (!initialize(&mop[nmats].imx))
		return(1);

	for (i = 0; i < n2comp; i++)	/* user .cal files and expressions */
		if (argv[comp_ndx[i]][1] == 'f') {
			char	*fpath = getpath(argv[comp_ndx[i]+1],
							getrlibpath(), 0);
			if (fpath == NULL) {
				fprintf(stderr, "%s: cannot find file '%s'\n",
						argv[0], argv[comp_ndx[i]+1]);
				return(1);
			}
			fcompile(fpath);
		} else /* (argv[comp_ndx[i]][1] == 'e') */
			scompile(argv[comp_ndx[i]+1], NULL, 0);

					/* get trailing color transform */
	if (!get_component_xfm(&mop[nmats]))
		return(1);
					/* adjust output dimensions and #components */
	if (mcat) {
		if (mop[nmats].imx.ncols != mcat->nrows) {
			fprintf(stderr,
			"%s: number of input columns does not match number of rows in '%s'\n",
					argv[0], mcat_spec);
			return(1);
		}
		if (mcat->ncomp != (mcat_last ? mop[nmats].rmp->ncomp : mop[nmats].imx.ncomp)) {
			fprintf(stderr,
			"%s: number of components does not match those in '%s'\n",
					argv[0], mcat_spec);
			return(1);
		}
		if (!split_input(&mop[nmats]))
			return(1);
		mop[nmats].rmp->ncols = mcat->ncols;
	}
	newheader("RADIANCE", stdout);	/* write output header */
	if (echoheader)
		output_headinfo(stdout);
	printargs(argc, argv, stdout);
	fputnow(stdout);
	mop[nmats].rmp->dtype = rmx_write_header(mop[nmats].rmp, outfmt, stdout);
	if (!mop[nmats].rmp->dtype) {
		fprintf(stderr, "%s: unsupported output format\n", argv[0]);
		return(1);
	}
	doptimize(1);			/* optimize definitions */
	if (spawned_children(nproc))	/* running in parent process? */
		return(parent_loop() ? 0 : 1);
					/* process & write rows */
	return(combine_input() ? 0 : 1);
stdin_error:
	fprintf(stderr, "%s: %s used for more than one input\n",
			argv[0], stdin_name);
	return(1);
userr:
	fprintf(stderr,
	"Usage: %s [-h][-f{adfc}][-n nproc][-e expr][-f file][-s sf .. | -c ce ..] m1 .. -m mcat > mres\n",
			argv[0]);
	return(1);
}
