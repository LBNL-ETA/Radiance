#ifndef lint
static const char RCSid[] = "$Id: rmatrix.c,v 2.59 2023/07/01 15:25:26 greg Exp $";
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

static char	rmx_mismatch_warn[] = "WARNING: data type mismatch\n";

#define array_size(rm)	(sizeof(double)*(rm)->nrows*(rm)->ncols*(rm)->ncomp)
#define mapped_size(rm)	((char *)(rm)->mtx + array_size(rm) - (char *)(rm)->mapped)

/* Initialize a RMATRIX struct but don't allocate array space */
RMATRIX *
rmx_new(int nr, int nc, int n)
{
	RMATRIX	*dnew = (RMATRIX *)calloc(1, sizeof(RMATRIX));

	if (dnew) {
		dnew->dtype = DTdouble;
		dnew->nrows = nr;
		dnew->ncols = nc;
		dnew->ncomp = n;
		setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	}
	return(dnew);
}

/* Prepare a RMATRIX for writing (allocate array if needed) */
int
rmx_prepare(RMATRIX *rm)
{
	if (!rm) return(0);
	if (rm->mtx)
		return(1);
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

/* Free a RMATRIX array */
void
rmx_free(RMATRIX *rm)
{
	if (!rm) return;
	if (rm->info)
		free(rm->info);
#ifdef MAP_FILE
	if (rm->mapped)
		munmap(rm->mapped, mapped_size(rm));
	else
#endif
		free(rm->mtx);
	free(rm);
}

/* Resolve data type based on two input types (returns 0 for mismatch) */
int
rmx_newtype(int dtyp1, int dtyp2)
{
	if ((dtyp1==DTxyze) | (dtyp1==DTrgbe) |
			(dtyp2==DTxyze) | (dtyp2==DTrgbe)
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
	if (!strncmp(s, "NCOMP=", 6)) {
		ip->ncomp = atoi(s+6);
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
rmx_load_ascii(RMATRIX *rm, FILE *fp)
{
	int	i, j, k;

	if (!rmx_prepare(rm))
		return(0);
	for (i = 0; i < rm->nrows; i++)
	    for (j = 0; j < rm->ncols; j++) {
	    	double	*dp = rmx_lval(rm,i,j);
	        for (k = 0; k < rm->ncomp; k++)
		    if (fscanf(fp, "%lf", &dp[k]) != 1)
			return(0);
	    }
	return(1);
}

static int
rmx_load_float(RMATRIX *rm, FILE *fp)
{
	int	i, j, k;
	float	val[100];

	if (rm->ncomp > 100) {
		fputs("Unsupported # components in rmx_load_float()\n", stderr);
		exit(1);
	}
	if (!rmx_prepare(rm))
		return(0);
	for (i = 0; i < rm->nrows; i++)
	    for (j = 0; j < rm->ncols; j++) {
	    	double	*dp = rmx_lval(rm,i,j);
		if (getbinary(val, sizeof(val[0]), rm->ncomp, fp) != rm->ncomp)
		    return(0);
		if (rm->swapin)
		    swap32((char *)val, rm->ncomp);
	        for (k = rm->ncomp; k--; )
		     dp[k] = val[k];
	    }
	return(1);
}

static int
rmx_load_double(RMATRIX *rm, FILE *fp)
{
	int	i;
#ifdef MAP_FILE
	long	pos;		/* map memory for file > 1MB if possible */
	if (!rm->swapin && array_size(rm) >= 1L<<20 &&
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
	if (!rmx_prepare(rm))
		return(0);
	for (i = 0; i < rm->nrows; i++) {
		if (getbinary(rmx_lval(rm,i,0), sizeof(double)*rm->ncomp,
					rm->ncols, fp) != rm->ncols)
			return(0);
		if (rm->swapin)
			swap64((char *)rmx_lval(rm,i,0), rm->ncols*rm->ncomp);
	}
	return(1);
}

static int
rmx_load_rgbe(RMATRIX *rm, FILE *fp)
{
	COLOR	*scan = (COLOR *)malloc(sizeof(COLOR)*rm->ncols);
	int	i, j;

	if (!scan)
		return(0);
	if (!rmx_prepare(rm))
		return(0);
	for (i = 0; i < rm->nrows; i++) {
	    double	*dp = rmx_lval(rm,i,0);
	    if (freadscan(scan, rm->ncols, fp) < 0) {
		free(scan);
		return(0);
	    }
	    for (j = 0; j < rm->ncols; j++, dp += 3) {
	        dp[0] = colval(scan[j],RED);
	        dp[1] = colval(scan[j],GRN);
	        dp[2] = colval(scan[j],BLU);
	    }
	}
	free(scan);
	return(1);
}

/* Load matrix from supported file type */
RMATRIX *
rmx_load(const char *inspec, RMPref rmp)
{
	FILE		*fp;
	RMATRIX		*dnew;

	if (!inspec)
		inspec = stdin_name;
	else if (!*inspec)
		return(NULL);
	if (inspec == stdin_name) {		/* reading from stdin? */
		fp = stdin;
	} else if (inspec[0] == '!') {
		if (!(fp = popen(inspec+1, "r")))
			return(NULL);
	} else {
		const char	*sp = inspec;	/* check suffix */
		while (*sp)
			++sp;
		while (sp > inspec && sp[-1] != '.')
			--sp;
		if (!strcasecmp(sp, "XML")) {	/* assume it's a BSDF */
			CMATRIX	*cm = rmp==RMPtrans ? cm_loadBTDF(inspec) :
					cm_loadBRDF(inspec, rmp==RMPreflB) ;
			if (!cm)
				return(NULL);
			dnew = rmx_from_cmatrix(cm);
			cm_free(cm);
			dnew->dtype = DTascii;
			return(dnew);
		}
						/* else open it ourselves */
		if (!(fp = fopen(inspec, "r")))
			return(NULL);
	}
	SET_FILE_BINARY(fp);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (!(dnew = rmx_new(0,0,3))) {
		fclose(fp);
		return(NULL);
	}
	dnew->dtype = DTascii;			/* assumed w/o FORMAT */
	if (getheader(fp, get_dminfo, dnew) < 0) {
		fclose(fp);
		return(NULL);
	}
	if ((dnew->nrows <= 0) | (dnew->ncols <= 0)) {
		if (!fscnresolu(&dnew->ncols, &dnew->nrows, fp)) {
			fclose(fp);
			return(NULL);
		}
		if ((dnew->dtype == DTrgbe) | (dnew->dtype == DTxyze) &&
				dnew->ncomp != 3) {
			fclose(fp);
			return(NULL);
		}
	}
	switch (dnew->dtype) {
	case DTascii:
		SET_FILE_TEXT(fp);
		if (!rmx_load_ascii(dnew, fp))
			goto loaderr;
		dnew->dtype = DTascii;		/* should leave double? */
		break;
	case DTfloat:
		if (!rmx_load_float(dnew, fp))
			goto loaderr;
		dnew->dtype = DTfloat;
		break;
	case DTdouble:
		if (!rmx_load_double(dnew, fp))
			goto loaderr;
		dnew->dtype = DTdouble;
		break;
	case DTrgbe:
	case DTxyze:
		if (!rmx_load_rgbe(dnew, fp))
			goto loaderr;
		break;
	default:
		goto loaderr;
	}
	if (fp != stdin) {
		if (inspec[0] == '!')
			pclose(fp);
		else
			fclose(fp);
	}
#ifdef getc_unlocked
	else
		funlockfile(fp);
#endif
						/* undo exposure? */
	if (dnew->ncomp == 3 && (dnew->cexp[0] != 1.f) |
			(dnew->cexp[1] != 1.f) | (dnew->cexp[2] != 1.f)) {
		double	cmlt[3];
		cmlt[0] = 1./dnew->cexp[0];
		cmlt[1] = 1./dnew->cexp[1];
		cmlt[2] = 1./dnew->cexp[2];
		rmx_scale(dnew, cmlt);
		setcolor(dnew->cexp, 1.f, 1.f, 1.f);
	}
	return(dnew);
loaderr:					/* should report error? */
	if (inspec[0] == '!')
		pclose(fp);
	else
		fclose(fp);
	rmx_free(dnew);
	return(NULL);
}

static int
rmx_write_ascii(const RMATRIX *rm, FILE *fp)
{
	const char	*fmt = (rm->dtype == DTfloat) ? " %.7e" :
			(rm->dtype == DTrgbe) | (rm->dtype == DTxyze) ? " %.3e" :
				" %.15e" ;
	int	i, j, k;

	for (i = 0; i < rm->nrows; i++) {
	    for (j = 0; j < rm->ncols; j++) {
		const double	*dp = rmx_lval(rm,i,j);
	        for (k = 0; k < rm->ncomp; k++)
		    fprintf(fp, fmt, dp[k]);
		fputc('\t', fp);
	    }
	    fputc('\n', fp);
	}
	return(1);
}

static int
rmx_write_float(const RMATRIX *rm, FILE *fp)
{
	int	i, j, k;
	float	val[100];

	if (rm->ncomp > 100) {
		fputs("Unsupported # components in rmx_write_float()\n", stderr);
		exit(1);
	}
	for (i = 0; i < rm->nrows; i++)
	    for (j = 0; j < rm->ncols; j++) {
	    	const double	*dp = rmx_lval(rm,i,j);
	        for (k = rm->ncomp; k--; )
		    val[k] = (float)dp[k];
		if (putbinary(val, sizeof(float), rm->ncomp, fp) != rm->ncomp)
			return(0);
	    }
	return(1);
}

static int
rmx_write_double(const RMATRIX *rm, FILE *fp)
{
	int	i;

	for (i = 0; i < rm->nrows; i++)
		if (putbinary(rmx_lval(rm,i,0), sizeof(double)*rm->ncomp,
					rm->ncols, fp) != rm->ncols)
			return(0);
	return(1);
}

static int
rmx_write_rgbe(const RMATRIX *rm, FILE *fp)
{
	COLR	*scan = (COLR *)malloc(sizeof(COLR)*rm->ncols);
	int	i, j;

	if (!scan)
		return(0);
	for (i = 0; i < rm->nrows; i++) {
	    for (j = rm->ncols; j--; ) {
	    	const double	*dp = rmx_lval(rm,i,j);
	    	if (rm->ncomp == 1)
	    		setcolr(scan[j], dp[0], dp[0], dp[0]);
	    	else
	        	setcolr(scan[j], dp[0], dp[1], dp[2]);
	    }
	    if (fwritecolrs(scan, rm->ncols, fp) < 0) {
		free(scan);
		return(0);
	    }
	}
	free(scan);
	return(1);
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

/* Write matrix to file type indicated by dtype */
int
rmx_write(const RMATRIX *rm, int dtype, FILE *fp)
{
	int	ok = 1;

	if (!rm | !fp || !rm->mtx)
		return(0);
#ifdef getc_unlocked
	flockfile(fp);
#endif
						/* complete header */
	if (rm->info)
		fputs(rm->info, fp);
	if (dtype == DTfromHeader)
		dtype = rm->dtype;
	else if (dtype == DTrgbe && (rm->dtype == DTxyze ||
					findCIEprims(rm->info)))
		dtype = DTxyze;
	else if ((dtype == DTxyze) & (rm->dtype == DTrgbe))
		dtype = DTrgbe;
	if (rm->ncomp == 3) {			/* write exposure? */
		if ((rm->cexp[RED] != rm->cexp[GRN]) |
				(rm->cexp[GRN] != rm->cexp[BLU]))
			fputcolcor(rm->cexp, fp);
		else if (rm->cexp[GRN] != 1.f)
			fputexpos(rm->cexp[GRN], fp);
	}
	if ((dtype != DTrgbe) & (dtype != DTxyze)) {
		fprintf(fp, "NROWS=%d\n", rm->nrows);
		fprintf(fp, "NCOLS=%d\n", rm->ncols);
		fprintf(fp, "NCOMP=%d\n", rm->ncomp);
	} else if ((rm->ncomp != 3) & (rm->ncomp != 1))
		return(0);			/* wrong # components */
	if ((dtype == DTfloat) | (dtype == DTdouble))
		fputendian(fp);			/* important to record */
	fputformat(cm_fmt_id[dtype], fp);
	fputc('\n', fp);
	switch (dtype) {			/* write data */
	case DTascii:
		ok = rmx_write_ascii(rm, fp);
		break;
	case DTfloat:
		ok = rmx_write_float(rm, fp);
		break;
	case DTdouble:
		ok = rmx_write_double(rm, fp);
		break;
	case DTrgbe:
	case DTxyze:
		fprtresolu(rm->ncols, rm->nrows, fp);
		ok = rmx_write_rgbe(rm, fp);
		break;
	default:
		return(0);
	}
	ok &= (fflush(fp) == 0);
#ifdef getc_unlocked
	funlockfile(fp);
#endif
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

/* Duplicate the given matrix */
RMATRIX *
rmx_copy(const RMATRIX *rm)
{
	RMATRIX	*dnew;

	if (!rm)
		return(NULL);
	dnew = rmx_alloc(rm->nrows, rm->ncols, rm->ncomp);
	if (!dnew)
		return(NULL);
	rmx_addinfo(dnew, rm->info);
	dnew->dtype = rm->dtype;
	memcpy(dnew->mtx, rm->mtx, array_size(dnew));
	return(dnew);
}

/* Allocate and assign transposed matrix */
RMATRIX *
rmx_transpose(const RMATRIX *rm)
{
	RMATRIX	*dnew;
	int	i, j;

	if (!rm)
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
	for (j = dnew->ncols; j--; )
	    for (i = dnew->nrows; i--; )
	    	memcpy(rmx_lval(dnew,i,j), rmx_lval(rm,j,i),
	    			sizeof(double)*dnew->ncomp);
	return(dnew);
}

/* Multiply (concatenate) two matrices and allocate the result */
RMATRIX *
rmx_multiply(const RMATRIX *m1, const RMATRIX *m2)
{
	RMATRIX	*mres;
	int	i, j, k, h;

	if (!m1 | !m2 || (m1->ncomp != m2->ncomp) | (m1->ncols != m2->nrows))
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
			d += rmx_lval(m1,i,h)[k] * rmx_lval(m2,h,j)[k];
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

	if (!m1 | !m2 || (m1->ncols != m2->ncols) | (m1->nrows != m2->nrows))
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
			d = rmx_lval(m2,i,j)[0];
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
			    d = rmx_lval(m2,i,j)[k];
			    if (d == 0) {
			        ++zeroDivides;
				rmx_lval(m1,i,j)[k] = 0;
			    } else
				rmx_lval(m1,i,j)[k] /= d;
			}
		} else {
		    if (m2->ncomp == 1) {
			const double	d = rmx_lval(m2,i,j)[0];
		        for (k = m1->ncomp; k--; )
			    rmx_lval(m1,i,j)[k] *= d;
		    } else
		        for (k = m1->ncomp; k--; )
			    rmx_lval(m1,i,j)[k] *= rmx_lval(m2,i,j)[k];
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

	if (!msum | !madd ||
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
	    	const double	*da = rmx_lval(madd,i,j);
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

	if (!rm | !sf)
		return(0);
	for (i = rm->nrows; i--; )
	    for (j = rm->ncols; j--; ) {
	    	double	*dp = rmx_lval(rm,i,j);
		for (k = rm->ncomp; k--; )
		    dp[k] *= sf[k];
	    }
	if (rm->info)
		rmx_addinfo(rm, "Applied scalar\n");
	return(1);
}

/* Allocate new matrix and apply component transformation */
RMATRIX *
rmx_transform(const RMATRIX *msrc, int n, const double cmat[])
{
	int	i, j, ks, kd;
	RMATRIX	*dnew;

	if (!msrc | (n <= 0) | !cmat)
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
		const double	*ds = rmx_lval(msrc,i,j);
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

	if (!rm || !rm->mtx | ((rm->ncomp != 3) & (rm->ncomp != 1)))
		return(NULL);
	cnew = cm_alloc(rm->nrows, rm->ncols);
	if (!cnew)
		return(NULL);
	for (i = cnew->nrows; i--; )
	    for (j = cnew->ncols; j--; ) {
		const double	*dp = rmx_lval(rm,i,j);
		COLORV		*cv = cm_lval(cnew,i,j);
		if (rm->ncomp == 1)
		    setcolor(cv, dp[0], dp[0], dp[0]);
		else
	    	    setcolor(cv, dp[0], dp[1], dp[2]);
	    }
	return(cnew);
}
