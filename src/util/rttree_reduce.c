#ifndef lint
static const char RCSid[] = "$Id: rttree_reduce.c,v 2.2 2011/05/31 20:50:26 greg Exp $";
#endif
/*
 *  A utility called by genBSDF.pl to reduce tensor tree samples and output
 *  in a standard format as required by XML specification for variable
 *  resolution BSDF data.  We are not meant to be run by the user directly.
 */

#include "rtio.h"
#include "rterror.h"
#include "platform.h"
#include <stdlib.h>

float	*datarr;		/* our loaded BSDF data array */
int	ttrank = 4;		/* tensor tree rank */
int	log2g = 4;		/* log2 of grid resolution */
int	infmt = 'a';		/* input format ('a','f','d') */
double	pctcull = 99.;		/* target culling percentile */

#define HISTLEN		300	/* histogram resolution */
#define HISTMAX		10.	/* maximum recorded value in histogram */

int	histo[HISTLEN];		/* histogram freq. of max-min BSDF */

double	tthresh;		/* acceptance threshold (TBD) */

#define dval3(ix,ox,oy)		datarr[((((ix)<<log2g)+(ox))<<log2g)+(oy)]
#define dval4(ix,iy,ox,oy)	datarr[((((((ix)<<log2g)+(iy))<<log2g)+(ox))<<log2g)+(oy)]

#define above_threshold(tp)	((tp)->vmax - (tp)->vmin > tthresh)

/* Tensor tree node */
typedef struct ttree_s {
	float		vmin, vmax;	/* value extrema */
	float		vavg;		/* average */
	struct ttree_s	*kid;		/* 2^ttrank children */
} TNODE;

/* Allocate a new set of children for the given node (no checks) */
static void
new_kids(TNODE *pn)
{
	pn->kid = (TNODE *)calloc(1<<ttrank, sizeof(TNODE));
	if (pn->kid == NULL)
		error(SYSTEM, "out of memory in new_kids");
}

/* Free children for this node */
static void
free_kids(TNODE *pn)
{
	int	i;

	if (pn->kid == NULL)
		return;
	for (i = 1<<ttrank; i--; )
		free_kids(pn->kid+i);
	free(pn->kid);
	pn->kid = NULL;
}

/* Build a tensor tree starting from the given hypercube */
static void
build_tree(TNODE *tp, const int bmin[], int l2s)
{
	int	bkmin[4];
	int	i, j;

	tp->vmin = 1e20;
	tp->vmax = 0;
	tp->vavg = 0;
	if (l2s <= 1) {			/* reached upper leaves */
		for (i = 1<<ttrank; i--; ) {
			float	val;
			for (j = ttrank; j--; )
				bkmin[j] = bmin[j] + (i>>j & 1);
			val = (ttrank == 3) ? dval3(bkmin[0],bkmin[1],bkmin[2])
				: dval4(bkmin[0],bkmin[1],bkmin[2],bkmin[3]);
			if (val < tp->vmin)
				tp->vmin = val;
			if (val > tp->vmax)
				tp->vmax = val;
			tp->vavg += val;
		}
		tp->vavg /= (float)(1<<ttrank);
					/* record stats */
		i = (HISTLEN/HISTMAX) * (tp->vmax - tp->vmin);
		if (i >= HISTLEN) i = HISTLEN-1;
		++histo[i];
		return;
	}
	--l2s;				/* else still branching */
	new_kids(tp);			/* grow recursively */
	for (i = 1<<ttrank; i--; ) {
		for (j = ttrank; j--; )
			bkmin[j] = bmin[j] + ((i>>j & 1)<<l2s);
		build_tree(tp->kid+i, bkmin, l2s);
		if (tp->kid[i].vmin < tp->vmin)
			tp->vmin = tp->kid[i].vmin;
		if (tp->kid[i].vmax > tp->vmax)
			tp->vmax = tp->kid[i].vmax;
		tp->vavg += tp->kid[i].vavg;
	}
	tp->vavg /= (float)(1<<ttrank);
}

/* Set our trimming threshold */
static void
set_threshold()
{
	int	hsum = 0;
	int	i;

	for (i = HISTLEN; i--; )
		hsum += histo[i];
	hsum = pctcull*.01 * (double)hsum;
	for (i = 0; hsum > 0; i++)
		hsum -= histo[i];
	tthresh = (HISTMAX/HISTLEN) * i;
}

/* Trim our tree according to the current threshold */
static void
trim_tree(TNODE *tp)
{
	if (tp->kid == NULL)
		return;
	if (above_threshold(tp)) {	/* keeping branches? */
		int	i = 1<<ttrank;
		while (i--)
			trim_tree(tp->kid+i);
		return;
	}
	free_kids(tp);			/* else trim at this point */
}

/* Print a tensor tree from the given hypercube */
static void
print_tree(const TNODE *tp, const int bmin[], int l2s)
{
	int	bkmin[4];
	int	i, j;

	for (j = log2g-l2s; j--; )	/* indent based on branch level */
		fputs("    ", stdout);
	fputc('{', stdout);		/* special case for upper leaves */
	if (l2s <= 1 && above_threshold(tp)) {
		for (i = 0; i < 1<<ttrank; i++) {
			float	val;
			for (j = ttrank; j--; )
				bkmin[j] = bmin[j] + (i>>j & 1);
			val = (ttrank == 3) ? dval3(bkmin[0],bkmin[1],bkmin[2])
				: dval4(bkmin[0],bkmin[1],bkmin[2],bkmin[3]);
			printf(" %.4e", val);
		}
		fputs(" }\n", stdout);
		return;
	}
	if (tp->kid == NULL) {		/* trimmed limb */
		printf(" %.6e }\n", tp->vavg);
		return;
	}
	--l2s;				/* else still branching */
	fputc('\n', stdout);
	for (i = 0; i < 1<<ttrank; i++) {
		for (j = ttrank; j--; )
			bkmin[j] = bmin[j] + ((i>>j & 1)<<l2s);
		print_tree(tp->kid+i, bkmin, l2s);
	}
	++l2s;
	for (j = log2g-l2s; j--; )
		fputs("    ", stdout);
	fputs("}\n", stdout);
}

/* Read a row of data in ASCII */
static int
read_ascii(float *rowp, int n)
{
	int	n2go;

	if ((rowp == NULL) | (n <= 0))
		return(0);
	for (n2go = n; n2go; n2go--)
		if (scanf("%f", rowp++) != 1)
			break;
	if (n2go)
		error(USER, "unexpected EOD on ascii input");
	return(n-n2go);
}

/* Read a row of float data */
static int
read_float(float *rowp, int n)
{
	int	nread;

	if ((rowp == NULL) | (n <= 0))
		return(0);
	nread = fread(rowp, sizeof(float), n, stdin);
	if (nread != n)
		error(USER, "unexpected EOF on float input");
	return(nread);
}

/* Read a row of double data */
static int
read_double(float *rowp, int n)
{
	static double	*rowbuf = NULL;
	static int	rblen = 0;
	int		nread, i;

	if ((rowp == NULL) | (n <= 0)) {
		if (rblen) {
			free(rowbuf);
			rowbuf = NULL; rblen = 0;
		}
		return(0);
	}
	if (rblen < n) {
		rowbuf = (double *)realloc(rowbuf, sizeof(double)*(rblen=n));
		if (rowbuf == NULL)
			error(SYSTEM, "out of memory in read_double");
	}
	nread = fread(rowbuf, sizeof(double), n, stdin);
	if (nread != n)
		error(USER, "unexpected EOF on double input");
	for (i = 0; i < nread; i++)
		*rowp++ = rowbuf[i];
	return(nread);
}

/* Load data array, filling zeroes for rank 3 demi-tensor */
static void
load_data()
{
	int	(*readf)(float *, int) = NULL;
	
	switch (infmt) {
	case 'a':
		readf = &read_ascii;
		break;
	case 'f':
		readf = &read_float;
		break;
	case 'd':
		readf = &read_double;
		break;
	default:
		error(COMMAND, "unsupported input format");
		break;
	}
	datarr = (float *)calloc(1<<(log2g*ttrank), sizeof(float));
	if (datarr == NULL)
		error(SYSTEM, "out of memory in load_data");
	if (ttrank == 3) {
		int	ix, ox;
		for (ix = 0; ix < 1<<log2g; ix++)
			for (ox = 0; ox < 1<<log2g; ox++)
				(*readf)(datarr+((((ix)<<log2g)+(ox))<<log2g),
						1<<(log2g-1));
	} else /* ttrank == 4 */ {
		int	ix, iy, ox;
		for (ix = 0; ix < 1<<log2g; ix++)
		    for (iy = 0; iy < 1<<log2g; iy++)
			for (ox = 0; ox < 1<<log2g; ox++)
				(*readf)(datarr +
				((((((ix)<<log2g)+(iy))<<log2g)+(ox))<<log2g),
						1<<log2g);
	}
	(*readf)(NULL, 0);	/* releases any buffers */
	if (infmt == 'a') {
		int	c;
		while ((c = getc(stdin)) != EOF) {
			switch (c) {
			case ' ':
			case '\t':
			case '\r':
			case '\n':
				continue;
			}
			error(WARNING, "data past end of expected input");
			break;
		}
	} else if (getc(stdin) != EOF)
		error(WARNING, "binary data past end of expected input");
}

/* Load BSDF array, coalesce uniform regions and format as tensor tree */
int
main(int argc, char *argv[])
{
	int	doheader = 1;
	int	bmin[4];
	TNODE	gtree;
	int	i;
					/* get options and parameters */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'h':
			doheader = !doheader;
			break;
		case 'r':
			ttrank = atoi(argv[++i]);
			if (ttrank != 3 && ttrank != 4)
				goto userr;
			break;
		case 'g':
			log2g = atoi(argv[++i]);
			if (log2g <= 1)
				goto userr;
			break;
		case 't':
			pctcull = atof(argv[++i]);
			if ((pctcull < 0) | (pctcull >= 100.))
				goto userr;
			break;
		case 'f':
			infmt = argv[i][2];
			if (!infmt || strchr("afd", infmt) == NULL)
				goto userr;
			break;
		default:
			goto userr;
		}
	if (i < argc-1)
		goto userr;
					/* load input data */
	if (i == argc-1 && freopen(argv[i], "rb", stdin) == NULL) {
		sprintf(errmsg, "cannot open input file '%s'", argv[i]);
		error(SYSTEM, errmsg);
	}
	if (infmt != 'a')
		SET_FILE_BINARY(stdin);
	load_data();
	if (doheader) {
		for (i = 0; i < argc; i++) {
			fputs(argv[i], stdout);
			fputc(i < argc-1 ? ' ' : '\n', stdout);
		}
		fputc('\n', stdout);
	}
	gtree.kid = NULL;		/* create our tree */
	bmin[0] = bmin[1] = bmin[2] = bmin[3] = 0;
	build_tree(&gtree, bmin, log2g);
					/* compute threshold & trim tree */
	set_threshold();
	trim_tree(&gtree);
					/* format to stdout */
	print_tree(&gtree, bmin, log2g);
	/* Clean up isn't necessary for main()...
	free_kids(&gtree);
	free(datarr);
	*/
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-h][-f{a|f|d}][-r {3|4}][-g log2grid][-t trim%%] [input]\n",
			argv[0]);
	return(1);
}
