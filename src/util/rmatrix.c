#ifndef lint
static const char RCSid[] = "$Id$";
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

/* Initialize a RMATRIX struct but don't allocate array space */
RMATRIX *
rmx_new(int nr, int nc, int ncomp)
{
	RMATRIX	*dnew;

	if (ncomp <= 0)
		return(NULL);

	dnew = (RMATRIX *)calloc(1, sizeof(RMATRIX));
	if (!dnew)
		return(NULL);

	dnew->dtype = DTrmx_native;
	dnew->nrows = nr;
	dnew->ncols = nc;
	dnew->ncomp = ncomp;
	setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	memcpy(dnew->wlpart, WLPART, sizeof(dnew->wlpart));

	return(dnew);
}

/* Prepare a RMATRIX for writing (allocate array if needed) */
int
rmx_prepare(RMATRIX *rm)
{
	if (!rm) return(0);
	if (rm->mtx)			/* assume it's right size */
		return(1);
	if ((rm->nrows <= 0) | (rm->ncols <= 0) | (rm->ncomp <= 0))
		return(0);
	rm->mtx = (rmx_dtype *)malloc(rmx_array_size(rm));
	rm->pflags |= RMF_FREEMEM;
	return(rm->mtx != NULL);
}

/* Call rmx_new() and rmx_prepare() */
RMATRIX	*
rmx_alloc(int nr, int nc, int ncomp)
{
	RMATRIX	*dnew = rmx_new(nr, nc, ncomp);

	if (!rmx_prepare(dnew)) {
		rmx_free(dnew);
		return(NULL);
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
#ifdef MAP_FILE
	if (rm->mapped) {
		munmap(rm->mapped, rmx_mapped_size(rm));
		rm->mapped = NULL;
	} else
#endif
	if (rm->pflags & RMF_FREEMEM) {
		free(rm->mtx);
		rm->pflags &= ~RMF_FREEMEM;
	}
	rm->mtx = NULL;
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
	size_t	oldlen = 0;

	if (!rm || !info || !*info)
		return(0);
	if (!rm->info) {
		rm->info = (char *)malloc(strlen(info)+1);
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

	if (isheadid(s))
		return(0);
	if (isncomp(s)) {
		ip->ncomp = ncompval(s);
		return(ip->ncomp - 1);
	}
	if (!strncmp(s, "NROWS=", 6)) {
		ip->nrows = atoi(s+6);
		return(ip->nrows - 1);
	}
	if (!strncmp(s, "NCOLS=", 6)) {
		ip->ncols = atoi(s+6);
		return(ip->ncols - 1);
	}
	if ((i = isbigendian(s)) >= 0) {
		if (nativebigendian() != i)
			ip->pflags |= RMF_SWAPIN;
		else
			ip->pflags &= ~RMF_SWAPIN;
		return(0);
	}
	if (isexpos(s)) {
		float	f = exposval(s);
		scalecolor(ip->cexp, f);
		return(f > .0 ? 0 : -1);
	}
	if (iscolcor(s)) {
		COLOR	ctmp;
		if (!colcorval(ctmp, s)) return(-1);
		multcolor(ip->cexp, ctmp);
		return(0);
	}
	if (iswlsplit(s))
		return(wlsplitval(ip->wlpart, s) - 1);

	if (!formatval(fmt, s)) {
		rmx_addinfo(ip, s);
		return(0);
	}			/* else check format */
	for (i = 1; i < DTend; i++)
		if (!strcmp(fmt, cm_fmt_id[i])) {
			ip->dtype = i;
			return(0);
		}
	return(-1);		/* bad format */
}

static int
rmx_load_ascii(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
{
	int	j, k;

	for (j = 0; j < rm->ncols; j++)
	        for (k = rm->ncomp; k-- > 0; )
			if (fscanf(fp, rmx_scanfmt, drp++) != 1)
				return(0);
	return(1);
}

static int
rmx_load_float(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
{
#if DTrmx_native==DTfloat
	if (getbinary(drp, sizeof(*drp)*rm->ncomp, rm->ncols, fp) != rm->ncols)
		return(0);
	if (rm->pflags & RMF_SWAPIN)
		swap32((char *)drp, rm->ncols*rm->ncomp);
#else
	int	j, k;
	float	val[MAXCOMP];

	if (rm->ncomp > MAXCOMP) {
		fputs("Unsupported # components in rmx_load_float()\n", stderr);
		exit(1);
	}
	for (j = 0; j < rm->ncols; j++) {
		if (getbinary(val, sizeof(val[0]), rm->ncomp, fp) != rm->ncomp)
			return(0);
		if (rm->pflags & RMF_SWAPIN)
			swap32((char *)val, rm->ncomp);
	        for (k = 0; k < rm->ncomp; k++)
			*drp++ = val[k];
	}
#endif
	return(1);
}

static int
rmx_load_double(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
{
#if DTrmx_native==DTdouble
	if (getbinary(drp, sizeof(*drp)*rm->ncomp, rm->ncols, fp) != rm->ncols)
		return(0);
	if (rm->pflags & RMF_SWAPIN)
		swap64((char *)drp, rm->ncols*rm->ncomp);
#else
	int	j, k;
	double	val[MAXCOMP];

	if (rm->ncomp > MAXCOMP) {
		fputs("Unsupported # components in rmx_load_double()\n", stderr);
		exit(1);
	}
	for (j = 0; j < rm->ncols; j++) {
		if (getbinary(val, sizeof(val[0]), rm->ncomp, fp) != rm->ncomp)
			return(0);
		if (rm->pflags & RMF_SWAPIN)
			swap64((char *)val, rm->ncomp);
	        for (k = 0; k < rm->ncomp; k++)
			*drp++ = (float)val[k];
	}
#endif
	return(1);
}

static int
rmx_load_rgbe(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
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

#if DTrmx_native==DTfloat
static int
rmx_load_spec(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
{
	COLRV	*scan;
	int	j;

	if ((rm->ncomp < 3) | (rm->ncomp > MAXCOMP))
		return(0);
	scan = (COLRV *)tempbuffer((rm->ncomp+1)*rm->ncols);
	if (!scan)
		return(0);
	if (freadscolrs(scan, rm->ncomp, rm->ncols, fp) < 0)
		return(0);
	for (j = 0; j < rm->ncols; j++, drp += rm->ncomp)
		scolr2scolor(drp, scan+j*(rm->ncomp+1), rm->ncomp);
	return(1);
}
#else
static int
rmx_load_spec(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
{
	COLRV	*scan;
	COLORV	scol[MAXCOMP];
	int	j, k;

	if ((rm->ncomp < 3) | (rm->ncomp > MAXCOMP))
		return(0);
	scan = (COLRV *)tempbuffer((rm->ncomp+1)*rm->ncols);
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
#endif

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
		rm->pflags = 0;
	}
	rm->dtype = DTascii;			/* assumed w/o FORMAT */
	if (getheader(fp, get_dminfo, rm) < 0) {
		fputs("Bad matrix header\n", stderr);
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

/* Load next row as rmx_dtype (cannot be XML) */
int
rmx_load_row(rmx_dtype *drp, const RMATRIX *rm, FILE *fp)
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
	if ((rm->dtype == DTrmx_native) & !(rm->pflags & RMF_SWAPIN) &
			(rmx_array_size(rm) >= 1L<<20) &&
			(pos = ftell(fp)) >= 0 && !(pos % sizeof(rmx_dtype))) {
		rm->mapped = mmap(NULL, rmx_array_size(rm)+pos, PROT_READ|PROT_WRITE,
					MAP_PRIVATE, fileno(fp), 0);
		if (rm->mapped != MAP_FAILED) {
			if (rm->pflags & RMF_FREEMEM)
				free(rm->mtx);
			rm->mtx = (rmx_dtype *)rm->mapped + pos/sizeof(rmx_dtype);
			rm->pflags &= ~RMF_FREEMEM;
			return(1);
		}		/* else fall back on reading into memory */
		rm->mapped = NULL;
	}
#endif
	if (!rmx_prepare(rm)) {	/* need in-core matrix array */
		fprintf(stderr, "Cannot allocate %g MByte matrix array\n",
				(1./(1L<<20))*(double)rmx_array_size(rm));
		return(0);
	}
	for (i = 0; i < rm->nrows; i++)
		if (!rmx_load_row(rmx_lval(rm,i,0), rm, fp))
			return(0);
	return(1);
}

/* Load matrix from supported file type */
RMATRIX *
rmx_load(const char *inspec)
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
	else
		fp = fopen(inspec, "r");
	if (!fp) {
		fprintf(stderr, "Cannot open for reading: %s\n", inspec);
		return(NULL);
	}
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
			ok &= pclose(fp)==0;
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
		double	cmlt[MAXCOMP];
		int	i;
		if (dnew->ncomp > MAXCOMP) {
			fprintf(stderr, "Excess spectral components in: %s\n",
					inspec);
			rmx_free(dnew);
			return(NULL);
		}
		cmlt[0] = 1./dnew->cexp[0];
		cmlt[1] = 1./dnew->cexp[1];
		cmlt[2] = 1./dnew->cexp[2];
		for (i = dnew->ncomp; i-- > 3; )
			cmlt[i] = cmlt[1];	/* XXX hack! */
		rmx_scale(dnew, cmlt);
		setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	}
	return(dnew);
}

#if DTrmx_native==DTdouble
static int
rmx_write_float(const rmx_dtype *dp, int len, FILE *fp)
{
	float	val;

	while (len--) {
		val = (float)*dp++;
		if (putbinary(&val, sizeof(val), 1, fp) != 1)
			return(0);
	}
	return(1);
}
#else
static int
rmx_write_double(const rmx_dtype *dp, int len, FILE *fp)
{
	double	val;

	while (len--) {
		val = *dp++;
		if (putbinary(&val, sizeof(val), 1, fp) != 1)
			return(0);
	}
	return(1);
}
#endif

static int
rmx_write_ascii(const rmx_dtype *dp, int ncomp, int len, FILE *fp)
{
	while (len-- > 0) {
		int	k = ncomp;
		while (k-- > 0)
			fprintf(fp, " %.7e", *dp++);
		fputc('\t', fp);
	}
	return(fputc('\n', fp) != EOF);
}

static int
rmx_write_rgbe(const rmx_dtype *dp, int ncomp, int len, FILE *fp)
{
	COLR	*scan;
	int	j;

	if ((ncomp != 1) & (ncomp != 3)) return(0);
	scan = (COLR *)tempbuffer(sizeof(COLR)*len);
	if (!scan) return(0);

	for (j = 0; j < len; j++, dp += ncomp)
	    	if (ncomp == 1)
	    		setcolr(scan[j], dp[0], dp[0], dp[0]);
	    	else
	        	setcolr(scan[j], dp[0], dp[1], dp[2]);

	return(fwritecolrs(scan, len, fp) >= 0);
}

#if DTrmx_native==DTfloat
static int
rmx_write_spec(const rmx_dtype *dp, int ncomp, int len, FILE *fp)
{
	COLRV	*scan;
	int	j;

	if ((ncomp < 3) | (ncomp > MAXCOMP)) return(0);
	scan = (COLRV *)tempbuffer((ncomp+1)*len);
	if (!scan) return(0);
	for (j = 0; j < len; j++, dp += ncomp)
		scolor2scolr(scan+j*(ncomp+1), dp, ncomp);

	return(fwritescolrs(scan, ncomp, len, fp) >= 0);
}
#else
static int
rmx_write_spec(const rmx_dtype *dp, int ncomp, int len, FILE *fp)
{
	COLRV	*scan;
	COLORV	scol[MAXCOMP];
	int	j, k;

	if ((ncomp < 3) | (ncomp > MAXCOMP)) return(0);
	scan = (COLRV *)tempbuffer((ncomp+1)*len);
	if (!scan) return(0);
	for (j = 0; j < len; j++, dp += ncomp) {
	    	for (k = ncomp; k--; )
	    		scol[k] = dp[k];
		scolor2scolr(scan+j*(ncomp+1), scol, ncomp);
	}
	return(fwritescolrs(scan, ncomp, len, fp) >= 0);
}
#endif

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
	if (dtype == DTfromHeader) {
		dtype = rm->dtype;
#if DTrmx_native==DTfloat
		if (dtype == DTdouble)		/* but stored as float? */
			dtype = DTfloat;
#endif
	} else if (dtype == DTrgbe && (rm->dtype == DTxyze ||
					findCIEprims(rm->info)))
		dtype = DTxyze;
	else if ((dtype == DTxyze) & (rm->dtype == DTrgbe))
		dtype = DTrgbe;
	if ((dtype < DTspec) & (rm->ncomp > 3))
		dtype = DTspec;
	else if ((dtype == DTspec) & (rm->ncomp <= 3))
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
		if (rm->ncomp > 3 &&
				memcmp(rm->wlpart, WLPART, sizeof(WLPART)))
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
rmx_write_data(const rmx_dtype *dp, int ncomp, int len, int dtype, FILE *fp)
{
	switch (dtype) {
#if DTrmx_native==DTdouble
	case DTfloat:
		return(rmx_write_float(dp, ncomp*len, fp));
#else
	case DTdouble:
		return(rmx_write_double(dp, ncomp*len, fp));
#endif
	case DTrmx_native:
		return(putbinary(dp, sizeof(*dp)*ncomp, len, fp) == len);
	case DTascii:
		return(rmx_write_ascii(dp, ncomp, len, fp));
	case DTrgbe:
	case DTxyze:
		return(rmx_write_rgbe(dp, ncomp, len, fp));
	case DTspec:
		return(rmx_write_spec(dp, ncomp, len, fp));
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
	if (dtype == DTrmx_native)		/* write all at once? */
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
	memset(rid->mtx, 0, rmx_array_size(rid));
	for (i = dim; i--; ) {
	    rmx_dtype	*dp = rmx_lval(rid,i,i);
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
		memcpy(dnew->mtx, rm->mtx, rmx_array_size(dnew));
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
	if (!rdst | !rsrc)
		return(0);
	if (dometa) {		/* transfer everything? */
		rmx_reset(rdst);
		*rdst = *rsrc;
		rsrc->info = NULL; rsrc->mapped = NULL; rsrc->mtx = NULL;
		return(1);
	}
				/* just matrix data -- leave metadata */
	if ((rdst->nrows != rsrc->nrows) |
			(rdst->ncols != rsrc->ncols) |
			(rdst->ncomp != rsrc->ncomp))
		return(0);
#ifdef MAP_FILE
	if (rdst->mapped)
		munmap(rdst->mapped, rmx_mapped_size(rdst));
	else
#endif
	if (rdst->pflags & RMF_FREEMEM) {
		free(rdst->mtx);
		rdst->pflags &= ~RMF_FREEMEM;
	}
	rdst->mapped = rsrc->mapped;
	rdst->mtx = rsrc->mtx;
	rdst->pflags |= rsrc->pflags & RMF_FREEMEM;
	rsrc->mapped = NULL; rsrc->mtx = NULL;
	return(1);
}

/* Transpose the given matrix */
int
rmx_transpose(RMATRIX *rm)
{
	uby8		*bmap;
	rmx_dtype	val[MAXCOMP];
	RMATRIX		dold;
	int		i, j;

	if (!rm || !rm->mtx | (rm->ncomp > MAXCOMP))
		return(0);
	if (rm->info)
		rmx_addinfo(rm, "Transposed rows and columns\n");
	if ((rm->nrows == 1) | (rm->ncols == 1)) { /* vector? */
		j = rm->ncols;
		rm->ncols = rm->nrows;
		rm->nrows = j;
		return(1);
	}
	if (rm->nrows == rm->ncols) {	/* square matrix case */
	     for (i = rm->nrows; --i > 0; )
	     	for (j = i; j-- > 0; ) {
		    memcpy(val, rmx_val(rm,i,j),
				sizeof(rmx_dtype)*rm->ncomp);
		    memcpy(rmx_lval(rm,i,j), rmx_val(rm,j,i),
				sizeof(rmx_dtype)*rm->ncomp);
		    memcpy(rmx_lval(rm,j,i), val,
				sizeof(rmx_dtype)*rm->ncomp);
		}
	    return(1);
	}
#define	bmbyte(r,c)	bmap[((r)*rm->ncols+(c))>>3]
#define	bmbit(r,c)	(1 << ((r)*rm->ncols+(c) & 7))
#define	bmop(r,c, op)	(bmbyte(r,c) op bmbit(r,c))
#define	bmtest(r,c)	bmop(r,c,&)
#define	bmset(r,c)	bmop(r,c,|=)
					/* loop completion bitmap */
	bmap = (uby8 *)calloc(((size_t)rm->nrows*rm->ncols+7)>>3, 1);
	if (!bmap)
		return(0);
	dold = *rm;
	rm->ncols = dold.nrows; rm->nrows = dold.ncols;
	for (i = rm->nrows; i--; )	/* try every starting point */
	    for (j = rm->ncols; j--; ) {
	    	int	i0, j0;
	    	int	i1 = i;
	    	size_t	j1 = j;
		if (bmtest(i, j))
			continue;	/* traversed loop earlier */
		memcpy(val, rmx_val(rm,i,j),
			sizeof(rmx_dtype)*rm->ncomp);
		for ( ; ; ) {		/* new transpose loop */
		    const rmx_dtype	*ds;
		    i0 = i1; j0 = j1;
		    ds = rmx_val(&dold, j0, i0);
		    j1 = (ds - dold.mtx)/dold.ncomp;
		    i1 = j1 / rm->ncols;
		    j1 -= (size_t)i1*rm->ncols;
		    bmset(i1, j1);	/* mark as done */
		    if ((i1 == i) & (j1 == j))
		    	break;		/* back at start */
		    memcpy(rmx_lval(rm,i0,j0), ds,
				sizeof(rmx_dtype)*rm->ncomp);
		}			/* complete the loop */
		memcpy(rmx_lval(rm,i0,j0), val,
			sizeof(rmx_dtype)*rm->ncomp);
	    }
	free(bmap);			/* all done! */
	return(1);
#undef	bmbyte
#undef	bmbit
#undef	bmop
#undef	bmtest
#undef	bmset
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
			d += (double)rmx_val(m1,i,h)[k] *
					rmx_val(m2,h,j)[k];
		    rmx_lval(mres,i,j)[k] = (rmx_dtype)d;
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
		    rmx_dtype	d;
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
			const rmx_dtype	d = rmx_val(m2,i,j)[0];
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
	    	const rmx_dtype	*da = rmx_val(madd,i,j);
	    	rmx_dtype	*ds = rmx_lval(msum,i,j);
		for (k = msum->ncomp; k--; )
		     ds[k] += (rmx_dtype)sf[k] * da[k];
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
	    	rmx_dtype	*dp = rmx_lval(rm,i,j);
		for (k = rm->ncomp; k--; )
		    dp[k] *= (rmx_dtype)sf[k];
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
		const rmx_dtype	*ds = rmx_val(msrc,i,j);
	        for (kd = dnew->ncomp; kd--; ) {
		    double	d = 0;
		    for (ks = msrc->ncomp; ks--; )
		        d += cmat[kd*msrc->ncomp + ks] * ds[ks];
		    rmx_lval(dnew,i,j)[kd] = (rmx_dtype)d;
		}
	    }
	return(dnew);
}

