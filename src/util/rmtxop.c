#ifndef lint
static const char RCSid[] = "$Id: rmtxop.c,v 2.15 2019/08/12 16:55:24 greg Exp $";
#endif
/*
 * General component matrix operations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "rtio.h"
#include "resolu.h"
#include "rmatrix.h"
#include "platform.h"

#define MAXCOMP		16		/* #components we support */

static const char	stdin_name[] = "<stdin>";

/* unary matrix operation(s) */
typedef struct {
	double		sca[MAXCOMP];		/* scalar coefficients */
	double		cmat[MAXCOMP*MAXCOMP];	/* component transformation */
	short		nsf;			/* number of scalars */
	short		clen;			/* number of coefficients */
	short		transpose;		/* do transpose? */
} RUNARYOP;

/* matrix input source and requested operation(s) */
typedef struct {
	const char	*inspec;		/* input specification */
	RUNARYOP	preop;			/* unary operation(s) */
	RMATRIX		*mtx;			/* original matrix if loaded */
	int		binop;			/* binary op with next (or 0) */
} ROPMAT;

int	verbose = 0;			/* verbose reporting? */

/* Load matrix */
static int
loadmatrix(ROPMAT *rop)
{
	if (rop->mtx != NULL)
		return(0);

	rop->mtx = rmx_load(rop->inspec == stdin_name ?
				(const char *)NULL : rop->inspec);
	if (rop->mtx == NULL) {
		fputs(rop->inspec, stderr);
		fputs(": cannot load matrix\n", stderr);
		return(-1);
	}
	return(1);
}

/* Get matrix and perform unary operations */
static RMATRIX *
loadop(ROPMAT *rop)
{
	RMATRIX	*mres;
	int	i;

	if (loadmatrix(rop) < 0)		/* make sure we're loaded */
		return(NULL);

	if (rop->preop.nsf > 0) {		/* apply scalar(s) */
		if (rop->preop.clen > 0) {
			fputs("Options -s and -c are exclusive\n", stderr);
			goto failure;
		}
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
	if (rop->preop.clen > 0) {	/* apply transform */
		if (rop->preop.clen % rop->mtx->ncomp) {
			fprintf(stderr, "%s: -c must have N x %d coefficients\n",
					rop->inspec, rop->mtx->ncomp);
			goto failure;
		}
		mres = rmx_transform(rop->mtx, rop->preop.clen/rop->mtx->ncomp,
					rop->preop.cmat);
		if (mres == NULL) {
			fprintf(stderr, "%s: matrix transform failed\n",
						rop->inspec);
			goto failure;
		}
		if (verbose)
			fprintf(stderr, "%s: applied %d x %d transform\n",
					rop->inspec, mres->ncomp,
					rop->mtx->ncomp);
		rmx_free(rop->mtx);
		rop->mtx = mres;
	}
	if (rop->preop.transpose) {	/* transpose matrix? */
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
		mres = rmx_multiply(mleft, mright);
		rmx_free(mleft);
		rmx_free(mright);
		if (mres == NULL) {
			fputs(inspec, stderr);
			if (mleft->ncols != mright->nrows)
				fputs(": mismatched dimensions for multiply\n",
						stderr);
			else
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
		mright = binaryop(mop[rpos].inspec,
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
				break;
			case 'c':
				if (n > MAXCOMP*MAXCOMP) n = MAXCOMP*MAXCOMP;
				i += mop[nmats].preop.clen =
					get_factors(mop[nmats].preop.cmat,
							n, argv+i+1);
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
					/* write result to stdout */
	if (outfmt == DTfromHeader)
		outfmt = mres->dtype;
	if (outfmt != DTascii)
		SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	if (!rmx_write(mres, outfmt, stdout)) {
		fprintf(stderr, "%s: error writing result matrix\n", argv[0]);
		return(1);
	}
	/* rmx_free(mres); free(mop); */
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-v][-f[adfc][-t][-s sf .. | -c ce ..] m1 [.+*/] .. > mres\n",
			argv[0]);
	return(1);
}
