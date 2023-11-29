#ifndef lint
static const char RCSid[] = "$Id: rmtxop.c,v 2.24 2023/11/29 18:56:28 greg Exp $";
#endif
/*
 * General component matrix operations.
 */

#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include "rtio.h"
#include "resolu.h"
#include "rmatrix.h"
#include "platform.h"

#define MAXCOMP		MAXCSAMP	/* #components we support */

/* Unary matrix operation(s) */
typedef struct {
	double		sca[MAXCOMP];		/* scalar coefficients */
	double		cmat[MAXCOMP*MAXCOMP];	/* component transformation */
	short		nsf;			/* number of scalars */
	short		clen;			/* number of coefficients */
	char		csym[11];		/* symbolic coefficients */
	char		transpose;		/* do transpose? */
} RUNARYOP;

/* Matrix input source and requested operation(s) */
typedef struct {
	const char	*inspec;		/* input specification */
	RMPref		rmp;			/* matrix preference */
	RUNARYOP	preop;			/* unary operation(s) */
	RMATRIX		*mtx;			/* original matrix if loaded */
	int		binop;			/* binary op with next (or 0) */
} ROPMAT;

int	verbose = 0;			/* verbose reporting? */

/* Load matrix */
static int
loadmatrix(ROPMAT *rop)
{
	if (rop->mtx != NULL)		/* already loaded? */
		return(0);

	rop->mtx = rmx_load(rop->inspec, rop->rmp);

	return(!rop->mtx ? -1 : 1);
}

/* Compute conversion row from spectrum to one channel of RGB */
static void
rgbrow(ROPMAT *rop, int r, int p)
{
	const int	nc = rop->mtx->ncomp;
	const float *	wlp = rop->mtx->wlpart;
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
	const int	nc = rop->mtx->ncomp;
	const float *	wlp = rop->mtx->wlpart;
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
	const int	nc = rop->mtx->ncomp;
	int		i;

	for (i = nc; i--; ) {
		SCOLOR	sclr;
		scolorblack(sclr);
		sclr[i] = 1.f;
		rop->preop.cmat[r*nc+i] = (*sf)(sclr, nc, rop->mtx->wlpart);
	}
}

/* Check/set symbolic transform */
static int
checksymbolic(ROPMAT *rop)
{
	const int	nc = rop->mtx->ncomp;
	const int	dt = rop->mtx->dtype;
	int		i, j;

	if (nc < 3) {
		fprintf(stderr, "%s: -c '%s' requires at least 3 components\n",
				rop->inspec, rop->preop.csym);
		return(-1);
	}
	rop->preop.clen = strlen(rop->preop.csym) * nc;
	if (rop->preop.clen > MAXCOMP*MAXCOMP) {
		fprintf(stderr, "%s: -c '%s' results in too many components\n",
				rop->inspec, rop->preop.csym);
		return(-1);
	}
	for (j = 0; rop->preop.csym[j]; j++) {
		int	comp = 0;
		switch (rop->preop.csym[j]) {
		case 'B':
			++comp;
			/* fall through */
		case 'G':
			++comp;
			/* fall through */
		case 'R':
			if (dt == DTxyze) {
				for (i = 3; i--; )
					rop->preop.cmat[j*nc+i] = 1./WHTEFFICACY *
							xyz2rgbmat[comp][i];
			} else if (nc == 3)
				rop->preop.cmat[j*nc+comp] = 1.;
			else
				rgbrow(rop, j, comp);
			break;
		case 'Z':
			++comp;
			/* fall through */
		case 'Y':
			++comp;
			/* fall through */
		case 'X':
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

			for (i = nc*(dt != DTxyze); i--; )
				rop->preop.cmat[j*nc+i] *= WHTEFFICACY;
			break;
		case 'S':		/* scotopic (il)luminance */
			sensrow(rop, j, scolor2scotopic);
			for (i = nc; i--; )
				rop->preop.cmat[j*nc+i] *= WHTSCOTOPIC;
			break;
		case 'M':		/* melanopic (il)luminance */
			sensrow(rop, j, scolor2melanopic);
			for (i = nc; i--; )
				rop->preop.cmat[j*nc+i] *= WHTMELANOPIC;
			break;
		case 'A':		/* average component */
			for (i = nc; i--; )
				rop->preop->cmat[j*nc+i] = 1./(double)nc;
			break;
		default:
			fprintf(stderr, "%s: -c '%c' unsupported\n",
				rop->inspec, rop->preop.csym[j]);
			return(-1);
		}
	}
					/* return recommended output type */
	if (!strcmp(rop->preop.csym, "XYZ")) {
		if (dt <= DTspec)
			return(DTxyze);
	} else if (!strcmp(rop->preop.csym, "RGB")) {
		if (dt <= DTspec)
			return(DTrgbe);
	}
	if ((nc > 3) & (dt <= DTspec))
		return(DTfloat);	/* probably not actual spectrum */
	return(0);
}

/* Get matrix and perform unary operations */
static RMATRIX *
loadop(ROPMAT *rop)
{
	int	outtype = 0;
	RMATRIX	*mres;
	int	i, j;

	if (loadmatrix(rop) < 0)		/* make sure we're loaded */
		return(NULL);

	if (rop->preop.csym[0] &&		/* symbolic transform? */
			(outtype = checksymbolic(rop)) < 0)
		goto failure;
	if (rop->preop.clen > 0) {		/* apply component transform? */
		if (rop->preop.clen % rop->mtx->ncomp) {
			fprintf(stderr, "%s: -c must have N x %d coefficients\n",
					rop->inspec, rop->mtx->ncomp);
			goto failure;
		}
		if (rop->preop.nsf > 0) {	/* scale transform, first */
			if (rop->preop.nsf == 1) {
				for (i = rop->preop.clen; i--; )
					rop->preop.cmat[i] *= rop->preop.sca[0];
 			} else if (rop->preop.nsf != rop->mtx->ncomp) {
				fprintf(stderr, "%s: -s must have one or %d factors\n",
						rop->inspec, rop->mtx->ncomp);
				goto failure;
			} else {
				for (j = rop->preop.clen/rop->preop.nsf; j--; )
					for (i = rop->preop.nsf; i--; )
						rop->preop.cmat[j*rop->preop.nsf+i] *=
								rop->preop.sca[i];
			}
		}
		mres = rmx_transform(rop->mtx, rop->preop.clen/rop->mtx->ncomp,
					rop->preop.cmat);
		if (mres == NULL) {
			fprintf(stderr, "%s: matrix transform failed\n",
						rop->inspec);
			goto failure;
		}
		if (verbose)
			fprintf(stderr, "%s: applied %d x %d transform%s\n",
					rop->inspec, mres->ncomp,
					rop->mtx->ncomp,
					rop->preop.nsf ? " (* scalar)" : "");
		rop->preop.nsf = 0;
		if ((mres->ncomp > 3) & (mres->dtype <= DTspec))
			outtype = DTfloat;	/* probably not actual spectrum */
		rmx_free(rop->mtx);
		rop->mtx = mres;
	}
	if (rop->preop.nsf > 0) {		/* apply scalar(s)? */
		if (rop->preop.nsf == 1) {
			for (i = rop->mtx->ncomp; --i; )
				rop->preop.sca[i] = rop->preop.sca[0];
		} else if (rop->preop.nsf != rop->mtx->ncomp) {
			fprintf(stderr, "%s: -s must have one or %d factors\n",
					rop->inspec, rop->mtx->ncomp);
			goto failure;
		}
		if (!rmx_scale(rop->mtx, rop->preop.sca)) {
			fputs(rop->inspec, stderr);
			fputs(": scalar operation failed\n", stderr);
			goto failure;
		}
		if (verbose) {
			fputs(rop->inspec, stderr);
			fputs(": applied scalar (", stderr);
			for (i = 0; i < rop->preop.nsf; i++)
				fprintf(stderr, " %f", rop->preop.sca[i]);
			fputs(" )\n", stderr);
		}
	}
	if (rop->preop.transpose) {		/* transpose matrix? */
		mres = rmx_transpose(rop->mtx);
		if (mres == NULL) {
			fputs(rop->inspec, stderr);
			fputs(": transpose failed\n", stderr);
			goto failure;
		}
		if (verbose) {
			fputs(rop->inspec, stderr);
			fputs(": transposed rows and columns\n", stderr);
		}
		rmx_free(rop->mtx);
		rop->mtx = mres;
	}
	mres = rop->mtx;
	rop->mtx = NULL;
	if (outtype)
		mres->dtype = outtype;
	return(mres);
failure:
	rmx_free(rop->mtx);
	return(rop->mtx = NULL);
}

/* Execute binary operation, free matrix arguments and return new result */
static RMATRIX *
binaryop(const char *inspec, RMATRIX *mleft, int op, RMATRIX *mright)
{
	RMATRIX	*mres = NULL;
	int	i;

	if ((mleft == NULL) | (mright == NULL))
		return(NULL);
	switch (op) {
	case '.':			/* concatenate */
		if (mleft->ncomp != mright->ncomp) {
			fputs(inspec, stderr);
			fputs(": # components do not match\n", stderr);
		} else if (mleft->ncols != mright->nrows) {
			fputs(inspec, stderr);
			fputs(": mismatched dimensions\n",
					stderr);
		} else
			mres = rmx_multiply(mleft, mright);
		rmx_free(mleft);
		rmx_free(mright);
		if (mres == NULL) {
			fputs(inspec, stderr);
			fputs(": concatenation failed\n", stderr);
			return(NULL);
		}
		if (verbose) {
			fputs(inspec, stderr);
			fputs(": concatenated matrix\n", stderr);
		}
		break;
	case '+':
		if (!rmx_sum(mleft, mright, NULL)) {
			fputs(inspec, stderr);
			fputs(": matrix sum failed\n", stderr);
			rmx_free(mleft);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose) {
			fputs(inspec, stderr);
			fputs(": added in matrix\n", stderr);
		}
		rmx_free(mright);
		mres = mleft;
		break;
	case '*':
	case '/': {
		const char *	tnam = (op == '/') ?
					"division" : "multiplication";
		errno = 0;
		if (!rmx_elemult(mleft, mright, (op == '/'))) {
			fprintf(stderr, "%s: element-wise %s failed\n",
					inspec, tnam);
			rmx_free(mleft);
			rmx_free(mright);
			return(NULL);
		}
		if (errno)
			fprintf(stderr,
				"%s: warning - error during element-wise %s\n",
					inspec, tnam);
		else if (verbose)
			fprintf(stderr, "%s: element-wise %s\n", inspec, tnam);
		rmx_free(mright);
		mres = mleft;
		} break;
	default:
		fprintf(stderr, "%s: unknown operation '%c'\n", inspec, op);
		rmx_free(mleft);
		rmx_free(mright);
		return(NULL);
	}
	return(mres);
}

/* Perform matrix operations from left to right */
static RMATRIX *
op_left2right(ROPMAT *mop)
{
	RMATRIX	*mleft = loadop(mop);

	while (mop->binop) {
		if (mleft == NULL)
			break;
		mleft = binaryop(mop[1].inspec,
				mleft, mop->binop, loadop(mop+1));
		mop++;
	}
	return(mleft);
}

/* Perform matrix operations from right to left */
static RMATRIX *
op_right2left(ROPMAT *mop)
{
	RMATRIX	*mright;
	int	rpos = 0;
					/* find end of list */
	while (mop[rpos].binop)
		if (mop[rpos++].binop != '.') {
			fputs(
		"Right-to-left evaluation only for matrix multiplication!\n",
					stderr);
			return(NULL);
		}
	mright = loadop(mop+rpos);
	while (rpos-- > 0) {
		if (mright == NULL)
			break;
		mright = binaryop(mop[rpos+1].inspec,
				loadop(mop+rpos), mop[rpos].binop, mright);
	}
	return(mright);
}

#define t_nrows(mop)	((mop)->preop.transpose ? (mop)->mtx->ncols \
						: (mop)->mtx->nrows)
#define t_ncols(mop)	((mop)->preop.transpose ? (mop)->mtx->nrows \
						: (mop)->mtx->ncols)

/* Should we prefer concatenating from rightmost matrix towards left? */
static int
prefer_right2left(ROPMAT *mop)
{
	int	mri = 0;

	while (mop[mri].binop)		/* find rightmost matrix */
		if (mop[mri++].binop != '.')
			return(0);	/* pre-empt reversal for other ops */

	if (mri <= 1)
		return(0);		/* won't matter */

	if (loadmatrix(mop+mri) < 0)	/* load rightmost cat */
		return(1);		/* fail will bail in a moment */

	if (t_ncols(mop+mri) == 1)
		return(1);		/* definitely better R->L */

	if (t_ncols(mop+mri) >= t_nrows(mop+mri))
		return(0);		/* ...probably worse */

	if (loadmatrix(mop) < 0)	/* load leftmost */
		return(0);		/* fail will bail in a moment */

	return(t_ncols(mop+mri) < t_nrows(mop));
}

static int
get_factors(double da[], int n, char *av[])
{
	int	ac;

	for (ac = 0; ac < n && isflt(av[ac]); ac++)
		da[ac] = atof(av[ac]);
	return(ac);
}

static ROPMAT *
grow_moparray(ROPMAT *mop, int n2alloc)
{
	int	nmats = 0;

	while (mop[nmats++].binop)
		;
	mop = (ROPMAT *)realloc(mop, n2alloc*sizeof(ROPMAT));
	if (mop == NULL) {
		fputs("Out of memory in grow_moparray()\n", stderr);
		exit(1);
	}
	if (n2alloc > nmats)
		memset(mop+nmats, 0, (n2alloc-nmats)*sizeof(ROPMAT));
	return(mop);
}

/* Load one or more matrices and operate on them, sending results to stdout */
int
main(int argc, char *argv[])
{
	int	outfmt = DTfromHeader;
	int	nall = 2;
	ROPMAT	*mop = (ROPMAT *)calloc(nall, sizeof(ROPMAT));
	int	nmats = 0;
	RMATRIX	*mres = NULL;
	int	stdin_used = 0;
	int	i;
					/* get options and arguments */
	for (i = 1; i < argc; i++) {
		if (argv[i][0] && !argv[i][1] &&
				strchr(".+*/", argv[i][0]) != NULL) {
			if (!nmats || mop[nmats-1].binop) {
				fprintf(stderr,
			"%s: missing matrix argument before '%c' operation\n",
						argv[0], argv[i][0]);
				return(1);
			}
			mop[nmats-1].binop = argv[i][0];
		} else if (argv[i][0] != '-' || !argv[i][1]) {
			if (argv[i][0] == '-') {
				if (stdin_used++) {
					fprintf(stderr,
			"%s: standard input used for more than one matrix\n",
						argv[0]);
					return(1);
				}
				mop[nmats].inspec = stdin_name;
			} else
				mop[nmats].inspec = argv[i];
			if (nmats > 0 && !mop[nmats-1].binop)
				mop[nmats-1].binop = '.';
			nmats++;
		} else {
			int	n = argc-1 - i;
			switch (argv[i][1]) {	/* get option */
			case 'v':
				verbose++;
				break;
			case 'f':
				switch (argv[i][2]) {
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
			case 't':
				mop[nmats].preop.transpose = 1;
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
			case 'c':
				if (n && isupper(argv[i+1][0])) {
					strlcpy(mop[nmats].preop.csym,
						argv[++i],
						sizeof(mop[0].preop.csym));
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
				mop[nmats].preop.csym[0] = '\0';
				break;
			case 'r':
				if (argv[i][2] == 'f')
					mop[nmats].rmp = RMPreflF;
				else if (argv[i][2] == 'b')
					mop[nmats].rmp = RMPreflB;
				else
					goto userr;
				break;
			default:
				fprintf(stderr, "%s: unknown operation '%s'\n",
						argv[0], argv[i]);
				goto userr;
			}
		}
		if (nmats >= nall)
			mop = grow_moparray(mop, nall += 2);
	}
	if (mop[0].inspec == NULL)	/* nothing to do? */
		goto userr;
	if (mop[nmats-1].binop) {
		fprintf(stderr,
			"%s: missing matrix argument after '%c' operation\n",
				argv[0], mop[nmats-1].binop);
		return(1);
	}
					/* favor quicker concatenation */
	mop[nmats].mtx = prefer_right2left(mop) ? op_right2left(mop)
						: op_left2right(mop);
	if (mop[nmats].mtx == NULL)
		return(1);
					/* apply trailing unary operations */
	mop[nmats].inspec = "trailing_ops";
	mres = loadop(mop+nmats);
	if (mres == NULL)
		return(1);
	if (outfmt == DTfromHeader)	/* check data type */
		outfmt = mres->dtype;
	if (outfmt == DTrgbe) {
		if (mres->ncomp > 3)
			outfmt = DTspec;
		else if (mres->dtype == DTxyze)
			outfmt = DTxyze;
	}
	if (outfmt != DTascii)		/* write result to stdout */
		SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);

	return(rmx_write(mres, outfmt, stdout) ? 0 : 1);
userr:
	fprintf(stderr,
	"Usage: %s [-v][-f[adfc][-t][-s sf .. | -c ce ..][-rf|-rb] m1 [.+*/] .. > mres\n",
			argv[0]);
	return(1);
}
