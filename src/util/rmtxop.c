#ifndef lint
static const char RCSid[] = "$Id: rmtxop.c,v 2.4 2014/08/05 21:45:05 greg Exp $";
#endif
/*
 * General component matrix operations.
 */

#include <stdio.h>
#include <stdlib.h>
#include "rtio.h"
#include "resolu.h"
#include "rmatrix.h"

#define MAXCOMP		50		/* #components we support */

typedef struct {
	double		sca[MAXCOMP];		/* scalar coefficients */
	int		nsf;			/* number of scalars */
	double		cmat[MAXCOMP*MAXCOMP];	/* component transformation */
	int		clen;			/* number of coefficients */
	int		transpose;		/* do transpose? */
	int		op;			/* '*' or '+' */
} ROPERAT;				/* matrix operation */

int	outfmt = DTfromHeader;		/* output format */
int	verbose = 0;			/* verbose reporting? */

static void
op_default(ROPERAT *op)
{
	memset(op, 0, sizeof(ROPERAT));
	op->op = '*';
}

static RMATRIX *
operate(RMATRIX *mleft, ROPERAT *op, const char *fname)
{
	RMATRIX	*mright = rmx_load(fname);
	RMATRIX	*mtmp;
	int	i;

	if (fname == NULL)
		fname = "<stdin>";
	if (mright == NULL) {
		fputs(fname, stderr);
		fputs(": cannot load matrix\n", stderr);
		return(NULL);
	}
	if (op->transpose) {		/* transpose matrix? */
		mtmp = rmx_transpose(mright);
		if (mtmp == NULL) {
			fputs(fname, stderr);
			fputs(": transpose failed\n", stderr);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose) {
			fputs(fname, stderr);
			fputs(": transposed rows and columns\n", stderr);
		}
		rmx_free(mright);
		mright = mtmp;
	}
	if (op->nsf > 0) {		/* apply scalar(s) */
		if (op->clen > 0) {
			fputs("Options -s and -c are exclusive\n", stderr);
			rmx_free(mright);
			return(NULL);
		}
		if (op->nsf == 1) {
			for (i = mright->ncomp; --i; )
				op->sca[i] = op->sca[0];
		} else if (op->nsf != mright->ncomp) {
			fprintf(stderr, "%s: -s must have one or %d factors\n",
					fname, mright->ncomp);
			rmx_free(mright);
			return(NULL);
		}
		if ((mleft == NULL) | (op->op != '+') &&
				!rmx_scale(mright, op->sca)) {
			fputs(fname, stderr);
			fputs(": scalar operation failed\n", stderr);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose) {
			fputs(fname, stderr);
			fputs(": applied scalar (", stderr);
			for (i = 0; i < op->nsf; i++)
				fprintf(stderr, " %f", op->sca[i]);
			fputs(" )\n", stderr);
		}
	}
	if (op->clen > 0) {		/* apply transform */
		if (op->clen % mright->ncomp) {
			fprintf(stderr, "%s: -c must have N x %d coefficients\n",
					fname, mright->ncomp);
			rmx_free(mright);
			return(NULL);
		}
		mtmp = rmx_transform(mright, op->clen/mright->ncomp, op->cmat);
		if (mtmp == NULL) {
			fprintf(stderr, "%s: matrix transform failed\n", fname);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose)
			fprintf(stderr, "%s: applied %d x %d transform\n",
					fname, mtmp->ncomp, mright->ncomp);
		rmx_free(mright);
		mright = mtmp;
	}
	if (mleft == NULL)		/* just one matrix */
		return(mright);
	if (op->op == '*') {		/* concatenate */
		RMATRIX	*mres = rmx_multiply(mleft, mright);
		if (mres == NULL) {
			fputs(fname, stderr);
			if (mleft->ncols != mright->nrows)
				fputs(": mismatched dimensions for multiply\n",
						stderr);
			else
				fputs(": concatenation failed\n", stderr);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose) {
			fputs(fname, stderr);
			fputs(": concatenated matrix\n", stderr);
		}
		rmx_free(mright);
		rmx_free(mleft);
		mleft = mres;
	} else if (op->op == '+') {
		if (!rmx_sum(mleft, mright, op->nsf ? op->sca : (double *)NULL)) {
			fputs(fname, stderr);
			fputs(": matrix sum failed\n", stderr);
			rmx_free(mright);
			return(NULL);
		}
		if (verbose) {
			fputs(fname, stderr);
			fputs(": added in matrix\n", stderr);
		}
		rmx_free(mright);
	} else {
		fprintf(stderr, "%s: unknown operation '%c'\n", fname, op->op);
		rmx_free(mright);
		return(NULL);
	}
	return(mleft);
}

static int
get_factors(double da[], int n, char *av[])
{
	int	ac;

	for (ac = 0; ac < n && isflt(av[ac]); ac++)
		da[ac] = atof(av[ac]);
	return(ac);
}

/* Load one or more matrices and operate on them, sending results to stdout */
int
main(int argc, char *argv[])
{
	RMATRIX	*mres = NULL;
	ROPERAT	op;
	long	nbw;
	int	i;
					/* initialize */
	op_default(&op);
					/* get options and arguments */
	for (i = 1; i < argc; i++)
		if (argv[i][0] == '+' && !argv[i][1]) {
			op.op = '+';
		} else if (argv[i][0] != '-' || !argv[i][1]) {
			char	*fname = NULL;	/* load matrix */
			if (argv[i][0] != '-')
				fname = argv[i];
			mres = operate(mres, &op, fname);
			if (mres == NULL) {
				fprintf(stderr, "%s: operation failed on '%s'\n",
						argv[0], argv[i]);
				return(0);
			}
			op_default(&op);	/* reset operator */
		} else {
			int	n = argc-1 - i;
			switch (argv[i][1]) {	/* get option */
			case 'v':
				verbose = !verbose;
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
				op.transpose = 1;
				break;
			case 's':
				if (n > MAXCOMP) n = MAXCOMP;
				op.nsf = get_factors(op.sca, n, argv+i+1);
				i += op.nsf;
				break;
			case 'c':
				if (n > MAXCOMP*MAXCOMP) n = MAXCOMP*MAXCOMP;
				op.clen = get_factors(op.cmat, n, argv+i+1);
				i += op.clen;
				break;
			default:
				fprintf(stderr, "%s: unknown operation '%s'\n",
						argv[0], argv[i]);
				goto userr;
			}
		}
	if (mres == NULL)		/* check that we got something */
		goto userr;
					/* write result to stdout */
#ifdef getc_unlocked
	flockfile(stdout);
#endif
#ifdef _WIN32
	if (outfmt != DTascii)
		_setmode(fileno(stdout), _O_BINARY);
#endif
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	nbw = rmx_write(mres, outfmt, stdout);
	/* rmx_free(mres); mres = NULL; */
	if (nbw <= 0) {
		fprintf(stderr, "%s: error writing result matrix\n", argv[0]);
		return(1);
	}
	if (verbose)
		fprintf(stderr, "%s: %ld bytes written\n", argv[0], nbw);
	return(0);
userr:
	fprintf(stderr,
	"Usage: %s [-v][-f[adfc][-t][-s sf .. | -c ce ..] m1 [+] .. > mres\n",
			argv[0]);
	return(1);
}
