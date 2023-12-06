#ifndef lint
static const char RCSid[] = "$Id: rmatrix.c,v 2.73 2023/12/06 17:57:34 greg Exp $";
#endif
/*
 * General matrix operations.
 */

#include <stdlib.h>
#include <errno.h>
#include "rtio.h"
#include "platform.h"
#include "resolu.h"
#include "paths.h"
#include "rmatrix.h"
#if !defined(_WIN32) && !defined(_WIN64)
#include <sys/mman.h>
#endif

static const char	rmx_mismatch_warn[] = "WARNING: data type mismatch\n";

#define array_size(rm)	(sizeof(double)*(rm)->nrows*(rm)->ncols*(rm)->ncomp)
#define mapped_size(rm)	((char *)(rm)->mtx + array_size(rm) - (char *)(rm)->mapped)

/* Initialize a RMATRIX struct but don't allocate array space */
RMATRIX *
rmx_new(int nr, int nc, int n)
{
	RMATRIX	*dnew;

	if (n <= 0)
		return(NULL);

	dnew = (RMATRIX *)calloc(1, sizeof(RMATRIX));
	if (!dnew)
		return(NULL);

	dnew->dtype = DTdouble;
	dnew->nrows = nr;
	dnew->ncols = nc;
	dnew->ncomp = n;
	setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	memcpy(dnew->wlpart, WLPART, sizeof(dnew->wlpart));

	return(dnew);
}

/* Prepare a RMATRIX for writing (allocate array if needed) */
int
rmx_prepare(RMATRIX *rm)
{
	if (!rm) return(0);
	if (rm->mtx)
		return(1);
	if ((rm->nrows <= 0) | (rm->ncols <= 0) | (rm->ncomp <= 0))
		return(0);
	rm->mtx = (double *)malloc(array_size(rm));
	return(rm->mtx != NULL);
}

/* Call rmx_new() and rmx_prepare() */
RMATRIX	*
rmx_alloc(int nr, int nc, int n)
{
	RMATRIX	*dnew = rmx_new(nr, nc, n);

	if (dnew && !rmx_prepare(dnew)) {
		rmx_free(dnew);
		dnew = NULL;
	}
	return(dnew);
}

/* Clear state by freeing info and matrix data */
void
rmx_reset(RMATRIX *rm)
{
	if (!rm) return;
	if (rm->info) {
		free(rm->info);
		rm->info = NULL;
	}
	if (rm->mtx) {
#ifdef MAP_FILE
		if (rm->mapped) {
			munmap(rm->mapped, mapped_size(rm));
			rm->mapped = NULL;
		} else
#endif
			free(rm->mtx);
		rm->mtx = NULL;
	}
}

/* Free an RMATRIX struct and data */
void
rmx_free(RMATRIX *rm)
{
	if (!rm) return;
	rmx_reset(rm);
	free(rm);
}

/* Resolve data type based on two input types (returns 0 for mismatch) */
int
rmx_newtype(int dtyp1, int dtyp2)
{
	if ((dtyp1==DTxyze) | (dtyp1==DTrgbe) | (dtyp1==DTspec) |
			(dtyp2==DTxyze) | (dtyp2==DTrgbe) | (dtyp2==DTspec)
			&& dtyp1 != dtyp2)
		return(0);
	if (dtyp1 < dtyp2)
		return(dtyp1);
	return(dtyp2);
}

/* Append header information associated with matrix data */
int
rmx_addinfo(RMATRIX *rm, const char *info)
{
	int	oldlen = 0;

	if (!rm || !info || !*info)
		return(0);
	if (!rm->info) {
		rm->info = (char *)malloc(strlen(info)+1);
		if (rm->info) rm->info[0] = '\0';
	} else {
		oldlen = strlen(rm->info);
		rm->info = (char *)realloc(rm->info,
				oldlen+strlen(info)+1);
	}
	if (!rm->info)
		return(0);
	strcpy(rm->info+oldlen, info);
	return(1);
}

static int
get_dminfo(char *s, void *p)
{
	RMATRIX	*ip = (RMATRIX *)p;
	char	fmt[MAXFMTLEN];
	int	i;

	if (headidval(NULL, s))
		return(0);
	if (isncomp(s)) {
		ip->ncomp = ncompval(s);
		return(0);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		ip->nrows = atoi(s+6);
		return(0);
	}
	if (!strncmp(s, "NCOLS=", 6)) {
		ip->ncols = atoi(s+6);
		return(0);
	}
	if ((i = isbigendian(s)) >= 0) {
		ip->swapin = (nativebigendian() != i);
		return(0);
	}
	if (isexpos(s)) {
		float	f = exposval(s);
		scalecolor(ip->cexp, f);
		return(0);
	}
	if (iscolcor(s)) {
		COLOR	ctmp;
		colcorval(ctmp, s);
		multcolor(ip->cexp, ctmp);
		return(0);
	}
	if (iswlsplit(s)) {
		wlsplitval(ip->wlpart, s);
		return(0);
	}
	if (!formatval(fmt, s)) {
		rmx_addinfo(ip, s);
		return(0);
	}			/* else check format */
	for (i = 1; i < DTend; i++)
		if (!strcmp(fmt, cm_fmt_id[i])) {
			ip->dtype = i;
			return(0);
		}
	return(-1);
}

static int
rmx_load_ascii(double *drp, const RMATRIX *rm, FILE *fp)
{
	int	j, k;

	for (j = 0; j < rm->ncols; j++)
	        for (k = rm->ncomp; k-- > 0; )
			if (fscanf(fp, "%lf", drp++) != 1)
				return(0);
	return(1);
}

static int
rmx_load_float(double *drp, const RMATRIX *rm, FILE *fp)
{
	int	j, k;
	float	val[100];

	if (rm->ncomp > 100) {
		fputs("Unsupported # components in rmx_load_float()\n", stderr);
		exit(1);
	}
	for (j = 0; j < rm->ncols; j++) {
		if (getbinary(val, sizeof(val[0]), rm->ncomp, fp) != rm->ncomp)
			return(0);
		if (rm->swapin)
			swap32((char *)val, rm->ncomp);
	        for (k = 0; k < rm->ncomp; k++)
			*drp++ = val[k];
	}
	return(1);
}

static int
rmx_load_double(double *drp, const RMATRIX *rm, FILE *fp)
{
	if (getbinary(drp, sizeof(*drp)*rm->ncomp, rm->ncols, fp) != rm->ncols)
		return(0);
	if (rm->swapin)
		swap64((char *)drp, rm->ncols*rm->ncomp);
	return(1);
}

static int
rmx_load_rgbe(double *drp, const RMATRIX *rm, FILE *fp)
{
	COLR	*scan;
	COLOR	col;
	int	j;

	if (rm->ncomp != 3)
		return(0);
	scan = (COLR *)tempbuffer(sizeof(COLR)*rm->ncols);
	if (!scan)
		return(0);
	if (freadcolrs(scan, rm->ncols, fp) < 0)
		return(0);
	for (j = 0; j < rm->ncols; j++) {
		colr_color(col, scan[j]);
		*drp++ = colval(col,RED);
		*drp++ = colval(col,GRN);
		*drp++ = colval(col,BLU);
	}
	return(1);
}

static int
rmx_load_spec(double *drp, const RMATRIX *rm, FILE *fp)
{
	uby8	*scan;
	SCOLOR	scol;
	int	j, k;

	if ((rm->ncomp < 3) | (rm->ncomp > MAXCSAMP))
		return(0);
	scan = (uby8 *)tempbuffer((rm->ncomp+1)*rm->ncols);
	if (!scan)
		return(0);
	if (freadscolrs(scan, rm->ncomp, rm->ncols, fp) < 0)
		return(0);
	for (j = 0; j < rm->ncols; j++) {
		scolr2scolor(scol, scan+j*(rm->ncomp+1), rm->ncomp);
		for (k = 0; k < rm->ncomp; k++)
			*drp++ = scol[k];
	}
	return(1);
}

/* Read matrix header from input stream (cannot be XML) */
int
rmx_load_header(RMATRIX *rm, FILE *fp)
{
	if (!rm | !fp)
		return(0);
	rmx_reset(rm);				/* clear state */
	if (rm->nrows | rm->ncols | !rm->dtype) {
		rm->nrows = rm->ncols = 0;
		rm->ncomp = 3;
		setcolor(rm->cexp, 1.f, 1.f, 1.f);
		memcpy(rm->wlpart, WLPART, sizeof(rm->wlpart));
		rm->swapin = 0;
	}
	rm->dtype = DTascii;			/* assumed w/o FORMAT */
	if (getheader(fp, get_dminfo, rm) < 0) {
		fputs("Unrecognized matrix format\n", stderr);
		return(0);
	}
	if ((rm->dtype == DTrgbe) | (rm->dtype == DTxyze) &&
			rm->ncomp != 3)
		return(0);
	if (rm->ncols <= 0 &&			/* resolution string? */
			!fscnresolu(&rm->ncols, &rm->nrows, fp))
		return(0);
	if (rm->dtype == DTascii)		/* set file type (WINDOWS) */
		SET_FILE_TEXT(fp);
	else
		SET_FILE_BINARY(fp);
	return(1);
}

/* Load next row as double (cannot be XML) */
int
rmx_load_row(double *drp, const RMATRIX *rm, FILE *fp)
{
	switch (rm->dtype) {
	case DTascii:
		return(rmx_load_ascii(drp, rm, fp));
	case DTfloat:
		return(rmx_load_float(drp, rm, fp));
	case DTdouble:
		return(rmx_load_double(drp, rm, fp));
	case DTrgbe:
	case DTxyze:
		return(rmx_load_rgbe(drp, rm, fp));
	case DTspec:
		return(rmx_load_spec(drp, rm, fp));
	default:
		fputs("Unsupported data type in rmx_load_row()\n", stderr);
	}
	return(0);
}

/* Allocate & load post-header data from stream given type set in rm->dtype */
int
rmx_load_data(RMATRIX *rm, FILE *fp)
{
	int	i;
#ifdef MAP_FILE
	long	pos;		/* map memory for file > 1MB if possible */
	if ((rm->dtype == DTdouble) & !rm->swapin && array_size(rm) >= 1L<<20 &&
			(pos = ftell(fp)) >= 0 && !(pos % sizeof(double))) {
		rm->mapped = mmap(NULL, array_size(rm)+pos, PROT_READ|PROT_WRITE,
					MAP_PRIVATE, fileno(fp), 0);
		if (rm->mapped != MAP_FAILED) {
			rm->mtx = (double *)rm->mapped + pos/sizeof(double);
			return(1);
		}		/* else fall back on reading into memory */
		rm->mapped = NULL;
	}
#endif
	if (!rmx_prepare(rm)) {	/* need in-core matrix array */
		fprintf(stderr, "Cannot allocate %g MByte matrix array\n",
				(1./(1L<<20))*(double)array_size(rm));
		return(0);
	}
	for (i = 0; i < rm->nrows; i++)
		if (!rmx_load_row(rmx_lval(rm,i,0), rm, fp))
			return(0);
	return(1);
}

/* Load matrix from supported file type */
RMATRIX *
rmx_load(const char *inspec, RMPref rmp)
{
	FILE		*fp;
	RMATRIX		*dnew;
	int		ok;

	if (!inspec)
		inspec = stdin_name;
	else if (!*inspec)
		return(NULL);
	if (inspec == stdin_name)		/* reading from stdin? */
		fp = stdin;
	else if (inspec[0] == '!')
		fp = popen(inspec+1, "r");
	else {
		const char	*sp = inspec;	/* check suffix */
		while (*sp)
			++sp;
		while (sp > inspec && sp[-1] != '.')
			--sp;
		if (!strcasecmp(sp, "XML")) {	/* assume it's a BSDF */
			CMATRIX	*cm = rmp==RMPnone ? (CMATRIX *)NULL :
					rmp==RMPtrans ? cm_loadBTDF(inspec) :
					cm_loadBRDF(inspec, rmp==RMPreflB) ;
			if (!cm)
				return(NULL);
			dnew = rmx_from_cmatrix(cm);
			cm_free(cm);
			dnew->dtype = DTascii;
			return(dnew);		/* return here */
		}				/* else open it ourselves */
		fp = fopen(inspec, "r");
	}
	if (!fp)
		return(NULL);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	SET_FILE_BINARY(fp);			/* load header info */
	if (!rmx_load_header(dnew = rmx_new(0,0,3), fp)) {
		fprintf(stderr, "Bad header in: %s\n", inspec);
		if (inspec[0] == '!') pclose(fp);
		else fclose(fp);
		rmx_free(dnew);
		return(NULL);
	}
	ok = rmx_load_data(dnew, fp);		/* allocate & load data */

	if (fp != stdin) {			/* close input stream */
		if (inspec[0] == '!')
			pclose(fp);
		else
			fclose(fp);
	}
#ifdef getc_unlocked
	else
		funlockfile(fp);
#endif
	if (!ok) {				/* load failure? */
		fprintf(stderr, "Error loading data from: %s\n", inspec);
		rmx_free(dnew);
		return(NULL);
	}
						/* undo exposure? */
	if ((dnew->cexp[0] != 1.f) |
			(dnew->cexp[1] != 1.f) | (dnew->cexp[2] != 1.f)) {
		double	cmlt[MAXCSAMP];
		int	i;
		cmlt[0] = 1./dnew->cexp[0];
		cmlt[1] = 1./dnew->cexp[1];
		cmlt[2] = 1./dnew->cexp[2];
		if (dnew->ncomp > MAXCSAMP) {
			fprintf(stderr, "Excess spectral components in: %s\n",
					inspec);
			rmx_free(dnew);
			return(NULL);
		}
		for (i = dnew->ncomp; i-- > 3; )
			cmlt[i] = cmlt[1];
		rmx_scale(dnew, cmlt);
		setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	}
	return(dnew);
}

static int
rmx_write_ascii(const double *dp, int nc, int len, FILE *fp)
{
	while (len-- > 0) {
		int	k = nc;
		while (k-- > 0)
			fprintf(fp, " %.7e", *dp++);
		fputc('\t', fp);
	}
	return(fputc('\n', fp) != EOF);
}

static int
rmx_write_float(const double *dp, int len, FILE *fp)
{
	float	val;

	while (len--) {
		val = *dp++;
		if (putbinary(&val, sizeof(float), 1, fp) != 1)
			return(0);
	}
	return(1);
}

static int
rmx_write_rgbe(const double *dp, int nc, int len, FILE *fp)
{
	COLR	*scan;
	int	j;

	if ((nc != 1) & (nc != 3)) return(0);
	scan = (COLR *)tempbuffer(sizeof(COLR)*len);
	if (!scan) return(0);

	for (j = 0; j < len; j++, dp += nc)
	    	if (nc == 1)
	    		setcolr(scan[j], dp[0], dp[0], dp[0]);
	    	else
	        	setcolr(scan[j], dp[0], dp[1], dp[2]);

	return(fwritecolrs(scan, len, fp) >= 0);
}

static int
rmx_write_spec(const double *dp, int nc, int len, FILE *fp)
{
	uby8	*scan;
	SCOLOR	scol;
	int	j, k;

	if (nc < 3) return(0);
	scan = (uby8 *)tempbuffer((nc+1)*len);
	if (!scan) return(0);
	for (j = len; j--; dp += nc) {
	    	for (k = nc; k--; )
	    		scol[k] = dp[k];
		scolor2scolr(scan+j*(nc+1), scol, nc);
	}
	return(fwritescolrs(scan, nc, len, fp) >= 0);
}

/* Check if CIE XYZ primaries were specified */
static int
findCIEprims(const char *info)
{
	RGBPRIMS	prims;

	if (!info)
		return(0);
	info = strstr(info, PRIMARYSTR);
	if (!info || !primsval(prims, info))
		return(0);

	return((prims[RED][CIEX] > .99) & (prims[RED][CIEY] < .01) &&
			(prims[GRN][CIEX] < .01) & (prims[GRN][CIEY] > .99) &&
			(prims[BLU][CIEX] < .01) & (prims[BLU][CIEY] < .01));
}

/* Finish writing header data with resolution and format, returning type used */
int
rmx_write_header(const RMATRIX *rm, int dtype, FILE *fp)
{
	if (!rm | !fp || rm->ncols <= 0)
		return(0);
	if (rm->info)
		fputs(rm->info, fp);
	if (dtype == DTfromHeader)
		dtype = rm->dtype;
	else if (dtype == DTrgbe && (rm->dtype == DTxyze ||
					findCIEprims(rm->info)))
		dtype = DTxyze;
	else if ((dtype == DTxyze) & (rm->dtype == DTrgbe))
		dtype = DTrgbe;
	if ((dtype == DTspec) & (rm->ncomp < 3))
		return(0);

	if (dtype == DTascii)			/* set file type (WINDOWS) */
		SET_FILE_TEXT(fp);
	else
		SET_FILE_BINARY(fp);
						/* write exposure? */
	if (rm->ncomp == 3 && (rm->cexp[RED] != rm->cexp[GRN]) |
			(rm->cexp[GRN] != rm->cexp[BLU]))
		fputcolcor(rm->cexp, fp);
	else if (rm->cexp[GRN] != 1.f)
		fputexpos(rm->cexp[GRN], fp);
						/* matrix size? */
	if ((dtype > DTspec) | (rm->nrows <= 0)) {
		if (rm->nrows > 0)
			fprintf(fp, "NROWS=%d\n", rm->nrows);
		fprintf(fp, "NCOLS=%d\n", rm->ncols);
	}
	if (dtype >= DTspec) {			/* # components & split? */
		fputncomp(rm->ncomp, fp);
		if (dtype == DTspec || (rm->ncomp > 3 &&
				memcmp(rm->wlpart, WLPART, sizeof(WLPART))))
			fputwlsplit(rm->wlpart, fp);
	} else if ((rm->ncomp != 3) & (rm->ncomp != 1))
		return(0);			/* wrong # components */
	if ((dtype == DTfloat) | (dtype == DTdouble))
		fputendian(fp);			/* important to record */
	fputformat(cm_fmt_id[dtype], fp);
	fputc('\n', fp);			/* end of header */
	if ((dtype <= DTspec) & (rm->nrows > 0))
		fprtresolu(rm->ncols, rm->nrows, fp);
	return(dtype);
}

/* Write out matrix data (usually by row) */
int
rmx_write_data(const double *dp, int nc, int len, int dtype, FILE *fp)
{
	switch (dtype) {
	case DTascii:
		return(rmx_write_ascii(dp, nc, len, fp));
	case DTfloat:
		return(rmx_write_float(dp, nc*len, fp));
	case DTdouble:
		return(putbinary(dp, sizeof(*dp)*nc, len, fp) == len);
	case DTrgbe:
	case DTxyze:
		return(rmx_write_rgbe(dp, nc, len, fp));
	case DTspec:
		return(rmx_write_spec(dp, nc, len, fp));
	}
	return(0);
}

/* Write matrix using file format indicated by dtype */
int
rmx_write(const RMATRIX *rm, int dtype, FILE *fp)
{
	int	ok = 0;
	int	i;
						/* complete header */
	dtype = rmx_write_header(rm, dtype, fp);
	if (dtype <= 0)
		return(0);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (dtype == DTdouble)			/* write all at once? */
		ok = rmx_write_data(rm->mtx, rm->ncomp,
				rm->nrows*rm->ncols, dtype, fp);
	else					/* else row by row */
		for (i = 0; i < rm->nrows; i++) {
			ok = rmx_write_data(rmx_val(rm,i,0), rm->ncomp,
					rm->ncols, dtype, fp);
			if (!ok) break;
		}

	if (ok) ok = (fflush(fp) == 0);
#ifdef getc_unlocked
	funlockfile(fp);
#endif
	if (!ok) fputs("Error writing matrix\n", stderr);
	return(ok);
}

/* Allocate and assign square identity matrix with n components */
RMATRIX *
rmx_identity(const int dim, const int n)
{
	RMATRIX	*rid = rmx_alloc(dim, dim, n);
	int	i, k;

	if (!rid)
		return(NULL);
	memset(rid->mtx, 0, array_size(rid));
	for (i = dim; i--; ) {
	    double	*dp = rmx_lval(rid,i,i);
	    for (k = n; k--; )
		dp[k] = 1.;
	}
	return(rid);
}

/* Duplicate the given matrix (may be unallocated) */
RMATRIX *
rmx_copy(const RMATRIX *rm)
{
	RMATRIX	*dnew;

	if (!rm)
		return(NULL);
	dnew = rmx_new(rm->nrows, rm->ncols, rm->ncomp);
	if (!dnew)
		return(NULL);
	if (rm->mtx) {
		if (!rmx_prepare(dnew)) {
			rmx_free(dnew);
			return(NULL);
		}
		memcpy(dnew->mtx, rm->mtx, array_size(dnew));
	}
	rmx_addinfo(dnew, rm->info);
	dnew->dtype = rm->dtype;
	copycolor(dnew->cexp, rm->cexp);
	memcpy(dnew->wlpart, rm->wlpart, sizeof(dnew->wlpart));
	return(dnew);
}

/* Replace data in first matrix with data from second */
int
rmx_transfer_data(RMATRIX *rdst, RMATRIX *rsrc, int dometa)
{
	if (!rdst | !rsrc || (rdst->nrows != rsrc->nrows) |
			(rdst->ncols != rsrc->ncols) |
			(rdst->ncomp != rsrc->ncomp))
		return(0);

	if (dometa) {		/* transfer everything? */
		rmx_reset(rdst);
		*rdst = *rsrc;
		rsrc->info = NULL; rsrc->mapped = NULL; rsrc->mtx = NULL;
		return(1);
	}
	if (rdst->mapped)
		return(0);	/* XXX can't handle this case */
				/* just matrix data -- leave metadata */
	if (rdst->mtx) free(rdst->mtx);
	rdst->mtx = rsrc->mtx;
	rsrc->mtx = NULL;
	return(1);
}

/* Allocate and assign transposed matrix */
RMATRIX *
rmx_transpose(const RMATRIX *rm)
{
	RMATRIX	*dnew;
	int	i, j;

	if (!rm || !rm->mtx)
		return(0);
	if ((rm->nrows == 1) | (rm->ncols == 1)) {
		dnew = rmx_copy(rm);
		if (!dnew)
			return(NULL);
		dnew->nrows = rm->ncols;
		dnew->ncols = rm->nrows;
		return(dnew);
	}
	dnew = rmx_alloc(rm->ncols, rm->nrows, rm->ncomp);
	if (!dnew)
		return(NULL);
	if (rm->info) {
		rmx_addinfo(dnew, rm->info);
		rmx_addinfo(dnew, "Transposed rows and columns\n");
	}
	dnew->dtype = rm->dtype;
	copycolor(dnew->cexp, rm->cexp);
	memcpy(dnew->wlpart, rm->wlpart, sizeof(dnew->wlpart));
	for (j = dnew->ncols; j--; )
	    for (i = dnew->nrows; i--; )
	    	memcpy(rmx_lval(dnew,i,j), rmx_val(rm,j,i),
	    			sizeof(double)*dnew->ncomp);
	return(dnew);
}

/* Multiply (concatenate) two matrices and allocate the result */
RMATRIX *
rmx_multiply(const RMATRIX *m1, const RMATRIX *m2)
{
	RMATRIX	*mres;
	int	i, j, k, h;

	if (!m1 | !m2 || !m1->mtx | !m2->mtx |
			(m1->ncomp != m2->ncomp) | (m1->ncols != m2->nrows))
		return(NULL);
	mres = rmx_alloc(m1->nrows, m2->ncols, m1->ncomp);
	if (!mres)
		return(NULL);
	i = rmx_newtype(m1->dtype, m2->dtype);
	if (i)
		mres->dtype = i;
	else
		rmx_addinfo(mres, rmx_mismatch_warn);
	for (i = mres->nrows; i--; )
	    for (j = mres->ncols; j--; )
	        for (k = mres->ncomp; k--; ) {
		    double	d = 0;
		    for (h = m1->ncols; h--; )
			d += rmx_val(m1,i,h)[k] * rmx_val(m2,h,j)[k];
		    rmx_lval(mres,i,j)[k] = d;
		}
	return(mres);
}

/* Element-wise multiplication (or division) of m2 into m1 */
int
rmx_elemult(RMATRIX *m1, const RMATRIX *m2, int divide)
{
	int	zeroDivides = 0;
	int	i, j, k;

	if (!m1 | !m2 || !m1->mtx | !m2->mtx |
			 (m1->ncols != m2->ncols) | (m1->nrows != m2->nrows))
		return(0);
	if ((m2->ncomp > 1) & (m2->ncomp != m1->ncomp))
		return(0);
	i = rmx_newtype(m1->dtype, m2->dtype);
	if (i)
		m1->dtype = i;
	else
		rmx_addinfo(m1, rmx_mismatch_warn);
	for (i = m1->nrows; i--; )
	    for (j = m1->ncols; j--; )
		if (divide) {
		    double	d;
		    if (m2->ncomp == 1) {
			d = rmx_val(m2,i,j)[0];
			if (d == 0) {
			    ++zeroDivides;
			    for (k = m1->ncomp; k--; )
				rmx_lval(m1,i,j)[k] = 0;
			} else {
			    d = 1./d;
			    for (k = m1->ncomp; k--; )
				rmx_lval(m1,i,j)[k] *= d;
			}
		    } else
		        for (k = m1->ncomp; k--; ) {
			    d = rmx_val(m2,i,j)[k];
			    if (d == 0) {
			        ++zeroDivides;
				rmx_lval(m1,i,j)[k] = 0;
			    } else
				rmx_lval(m1,i,j)[k] /= d;
			}
		} else {
		    if (m2->ncomp == 1) {
			const double	d = rmx_val(m2,i,j)[0];
		        for (k = m1->ncomp; k--; )
			    rmx_lval(m1,i,j)[k] *= d;
		    } else
		        for (k = m1->ncomp; k--; )
			    rmx_lval(m1,i,j)[k] *= rmx_val(m2,i,j)[k];
		}
	if (zeroDivides) {
		rmx_addinfo(m1, "WARNING: zero divide(s) corrupted results\n");
		errno = ERANGE;
	}
	return(1);
}

/* Sum second matrix into first, applying scale factor beforehand */
int
rmx_sum(RMATRIX *msum, const RMATRIX *madd, const double sf[])
{
	double	*mysf = NULL;
	int	i, j, k;

	if (!msum | !madd || !msum->mtx | !madd->mtx |
			(msum->nrows != madd->nrows) |
			(msum->ncols != madd->ncols) |
			(msum->ncomp != madd->ncomp))
		return(0);
	if (!sf) {
		mysf = (double *)malloc(sizeof(double)*msum->ncomp);
		if (!mysf)
			return(0);
		for (k = msum->ncomp; k--; )
			mysf[k] = 1;
		sf = mysf;
	}
	i = rmx_newtype(msum->dtype, madd->dtype);
	if (i)
		msum->dtype = i;
	else
		rmx_addinfo(msum, rmx_mismatch_warn);
	for (i = msum->nrows; i--; )
	    for (j = msum->ncols; j--; ) {
	    	const double	*da = rmx_val(madd,i,j);
	    	double		*ds = rmx_lval(msum,i,j);
		for (k = msum->ncomp; k--; )
		     ds[k] += sf[k] * da[k];
	    }
	if (mysf)
		free(mysf);
	return(1);
}

/* Scale the given matrix by the indicated scalar component vector */
int
rmx_scale(RMATRIX *rm, const double sf[])
{
	int	i, j, k;

	if (!rm | !sf || !rm->mtx)
		return(0);
	for (i = rm->nrows; i--; )
	    for (j = rm->ncols; j--; ) {
	    	double	*dp = rmx_lval(rm,i,j);
		for (k = rm->ncomp; k--; )
		    dp[k] *= sf[k];
	    }
	if (rm->info)
		rmx_addinfo(rm, "Applied scalar\n");
	/* XXX: should record as exposure for COLR and SCOLR types? */
	return(1);
}

/* Allocate new matrix and apply component transformation */
RMATRIX *
rmx_transform(const RMATRIX *msrc, int n, const double cmat[])
{
	int	i, j, ks, kd;
	RMATRIX	*dnew;

	if (!msrc | (n <= 0) | !cmat || !msrc->mtx)
		return(NULL);
	dnew = rmx_alloc(msrc->nrows, msrc->ncols, n);
	if (!dnew)
		return(NULL);
	if (msrc->info) {
		char	buf[128];
		sprintf(buf, "Applied %dx%d component transform\n",
				dnew->ncomp, msrc->ncomp);
		rmx_addinfo(dnew, msrc->info);
		rmx_addinfo(dnew, buf);
	}
	dnew->dtype = msrc->dtype;
	for (i = dnew->nrows; i--; )
	    for (j = dnew->ncols; j--; ) {
		const double	*ds = rmx_val(msrc,i,j);
	        for (kd = dnew->ncomp; kd--; ) {
		    double	d = 0;
		    for (ks = msrc->ncomp; ks--; )
		        d += cmat[kd*msrc->ncomp + ks] * ds[ks];
		    rmx_lval(dnew,i,j)[kd] = d;
		}
	    }
	return(dnew);
}

/* Convert a color matrix to newly allocated RMATRIX buffer */
RMATRIX *
rmx_from_cmatrix(const CMATRIX *cm)
{
	int	i, j;
	RMATRIX	*dnew;

	if (!cm)
		return(NULL);
	dnew = rmx_alloc(cm->nrows, cm->ncols, 3);
	if (!dnew)
		return(NULL);
	dnew->dtype = DTfloat;
	for (i = dnew->nrows; i--; )
	    for (j = dnew->ncols; j--; ) {
		const COLORV	*cv = cm_lval(cm,i,j);
		double		*dp = rmx_lval(dnew,i,j);
		dp[0] = cv[0];
		dp[1] = cv[1];
		dp[2] = cv[2];
	    }
	return(dnew);
}

/* Convert general matrix to newly allocated CMATRIX buffer */
CMATRIX *
cm_from_rmatrix(const RMATRIX *rm)
{
	int	i, j;
	CMATRIX	*cnew;

	if (!rm || !rm->mtx | (rm->ncomp == 2))
		return(NULL);
	cnew = cm_alloc(rm->nrows, rm->ncols);
	if (!cnew)
		return(NULL);
	for (i = cnew->nrows; i--; )
	    for (j = cnew->ncols; j--; ) {
		const double	*dp = rmx_val(rm,i,j);
		COLORV		*cv = cm_lval(cnew,i,j);
		switch (rm->ncomp) {
		case 3:
	    	    setcolor(cv, dp[0], dp[1], dp[2]);
	    	    break;
		case 1:
		    setcolor(cv, dp[0], dp[0], dp[0]);
		    break;
		default: {
			SCOLOR	scol;
			int	k;
			for (k = rm->ncomp; k--; )
				scol[k] = dp[k];
			scolor2color(cv, scol, rm->ncomp, rm->wlpart);
		    } break;
		}
	    }
	return(cnew);
}
