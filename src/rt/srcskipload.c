#ifndef lint
static const char RCSid[] = "$Id: srcskipload.c,v 2.1 2024/11/15 20:47:42 greg Exp $";
#endif
/*
 * Load source exclusion maps (BMP and MTX)
 *
 *  External symbols declared in source.h
 */

#ifdef SSKIPOPT		/* currently an optional compile */

#include "copyright.h"

#include  "platform.h"
#include  "paths.h"
#include  "ray.h"
#include  "source.h"
#include  "bmpfile.h"
#include  "resolu.h"

int		sskip_dim[2] = {0,0};		/* source skip image size */

/* correction file types */
enum {CFnone=0, CFfloatY, CFbmpY8};

/* struct for bitmap file loading */
typedef struct skipbmp {
	struct skipbmp	*next;			/* next in skip list */
	BMPReader	*bmp;			/* BMP file reader */
	uby8		*sfl;			/* corresponding sources */
	char		bname[1];		/* BMP file name */
} SKIPBMP;

/* holder for source exclusion images */
struct {
	SKIPBMP	*sbmplist;		/* list of BMP inputs */
	int	cftype;			/* correction file type */
	union {
		struct {
			int	fd;
			int	dstart;
		} Y;				/* float Y input */
		BMPReader	*bmp;		/* BMP input pointer */
	}	cf;			/* correction file input */
	int	*ndxmap;		/* allocated index map */
	float	*corrmap;		/* allocated correction map */
	char	cfname[MAXSTR];		/* correction file name */
} skipin;		/* just need the one */

/* call-back for matrix header line */
static int
check_mtx(char *hl, void *p)
{
	int	*dim = (int *)p;
	char	fmt[MAXFMTLEN];
	int	rv;

	if ((rv = isbigendian(hl)) >= 0) {
		if (rv != nativebigendian()) {
			eputs("cannot handle byte-swapped data\n");
			return(-1);
		}
		return(1);
	}
	if (!strncmp(hl, "NCOLS=", 6)) {
		dim[0] = atoi(hl+6);
		if (dim[0] <= 0)
			return(-1);
		return(1);
	}
	if (!strncmp(hl, "NROWS=", 6)) {
		dim[1] = atoi(hl+6);
		if (dim[1] <= 0)
			return(-1);
		return(1);
	}
	if (isncomp(hl)) {
		if (ncompval(hl) != 1) {
			eputs("require single component\n");
			return(-1);
		}
		return(1);
	}
	if (formatval(fmt, hl)) {
		if (strcmp(fmt, "float")) {
			eputs("require binary float format\n");
			return(-1);
		}
		return(1);
	}
	return(0);
}

/* private call to open floating-point matrix and check dimensions */
static void
open_float_mtx()
{
	int	mydim[2];
	FILE	*fp;

	skipin.cf.Y.fd = open(skipin.cfname, O_RDONLY);
	if (skipin.cf.Y.fd < 0) {
		sprintf(errmsg, "cannot open matrix '%s'", skipin.cfname);
		error(SYSTEM, errmsg);
	}
	SET_FD_BINARY(skipin.cf.Y.fd);	/* get temporary FILE pointer */
	fp = fdopen(dup(skipin.cf.Y.fd), "rb");
	if (fp == NULL)
		error(SYSTEM, "out of memory in open_float_mtx()");
	mydim[0] = mydim[1] = 0;	/* check header format, etc. */
	if (getheader(fp, check_mtx, mydim) < 0)
		goto badheader;
	if (!mydim[0] | !mydim[1] &&
			!fscnresolu(&mydim[0], &mydim[1], fp))
		goto badheader;
	if ((mydim[0] == sskip_dim[0]) & (mydim[1] == sskip_dim[1])) {
		skipin.cf.Y.dstart = ftell(fp);
		fclose(fp);
		return;			/* lookin' good! */
	}
badheader:
	sprintf(errmsg, "incompatible header for matrix '%s'", skipin.cfname);
	error(USER, errmsg);
	fclose(fp);
}

/* Open a set of bitmaps corresponding to loaded sources w/ corrections */
int
srcskip_open(char *bmpspec, char *scorrimg)
{
	int	bmcnt = 0;
	char	fname[MAXSTR];
	char	*sfx;
	int	sn;

	srcskip_free_maps();		/* clear previous */
	srcskip_close();
	sskip_dim[0] = sskip_dim[1] = 0;
	if (!nsources) {
		sskip_rsi(NULL);
		return(0);
	}
	if (bmpspec == NULL)
		return(0);
	if (strstr(bmpspec, "%s") == NULL)
		error(USER, "missing '%s' in source skip BMP file spec");
	if (ssf_select == NULL)		/* starting fresh? */
		ssf_select = sskip_new();
	if (scorrimg == NULL)		/* no correction map? */
		skipin.cfname[0] = '\0';
	else
		strcpy(skipin.cfname, scorrimg);
	skipin.cftype = CFnone;		/* open input BMPs */
	for (sn = 0; sn < nsources; sn++) {
		OBJECT		mod;
		const char	*modname;
		SKIPBMP		*sbmp;
		int		sn1;
		if (source[sn].sflags & (SSKIP|SVIRTUAL) ||
				sskip_chk(ssf_select, sn))
			continue;
		mod = source[sn].so->omod;
		modname = objptr(mod)->oname;
		sprintf(fname, bmpspec, modname);
		if (access(fname, R_OK) < 0)
			continue;	/* none such input */
					/* else check it out */
		sbmp = (SKIPBMP *)emalloc(sizeof(SKIPBMP)+strlen(fname));
		strcpy(sbmp->bname, fname);
		sbmp->bmp = BMPopenInputFile(fname);
		if (sbmp->bmp == NULL) {
			sprintf(errmsg, "cannot open bitmap '%s'", fname);
			error(SYSTEM, errmsg);
		}
		if (sbmp->bmp->hdr->nColors > 2) {
			sprintf(errmsg, "expected bilevel bitmap in '%s'", fname);
			error(USER, errmsg);
		}
		if (!sbmp->bmp->hdr->yIsDown) {
			sprintf(errmsg, "bad orientation for '%s'", fname);
			error(INTERNAL, errmsg);
		}
		if (skipin.sbmplist == NULL) {
			sskip_dim[0] = sbmp->bmp->hdr->width;
			sskip_dim[1] = sbmp->bmp->hdr->height;
		} else if ((sbmp->bmp->hdr->width != sskip_dim[0]) |
				(sbmp->bmp->hdr->height != sskip_dim[1])) {
			sprintf(errmsg, "dimensions do not match for '%s'",
					fname);
			error(USER, errmsg);
		}
					/* flag this light source */
		sbmp->sfl = sskip_new();
		sskip_set(sbmp->sfl, sn);
					/* flag others w/ same modifier */
		for (sn1 = sn; ++sn1 < nsources; ) {
			OBJECT	mod1;
			if (source[sn1].sflags & (SSKIP|SVIRTUAL) ||
					sskip_chk(ssf_select, sn1))
				continue;
			mod1 = source[sn1].so->omod;
			if (mod1 == mod || !strcmp(objptr(mod1)->oname, modname))
				sskip_set(sbmp->sfl, sn1);
		}
		sskip_addflags(ssf_select, sbmp->sfl);
		sbmp->next = skipin.sbmplist;
		skipin.sbmplist = sbmp;
		bmcnt++;
	}
	if (!bmcnt) {
		sprintf(errmsg, "no matching BMP input files for '%s'", bmpspec);
		error(WARNING, errmsg);
		return(0);
	}
	if (!skipin.cfname[0])		/* no correction image? */
		return(bmcnt);
					/* else open correction image */
	sfx = skipin.cfname;		/* find file type */
	while (*sfx) sfx++;
	while (--sfx > skipin.cfname && !ISDIRSEP(*sfx))
		if (*sfx == '.') {
			if (!strcasecmp(sfx, ".mtx"))
				skipin.cftype = CFfloatY;
			else if (!strcasecmp(sfx, ".bmp"))
				skipin.cftype = CFbmpY8;
			break;
		}
	switch (skipin.cftype) {
	case CFfloatY:
		open_float_mtx();
		break;
	case CFbmpY8:
		skipin.cf.bmp = BMPopenInputFile(skipin.cfname);
		if (skipin.cf.bmp == NULL) {
			sprintf(errmsg, "cannot open image '%s'", skipin.cfname);
			error(SYSTEM, errmsg);
		}
		if (!skipin.cf.bmp->hdr->yIsDown |
				(skipin.cf.bmp->hdr->bpp != 8) |
				(skipin.cf.bmp->hdr->width != sskip_dim[0]) |
				(skipin.cf.bmp->hdr->height != sskip_dim[1])) {
			sprintf(errmsg, "bad type/size/orientation for '%s'",
					skipin.cfname);
			error(USER, errmsg);
		}
		break;
	case CFnone:
		sprintf(errmsg, "unsupported image type for '%s'", skipin.cfname);
		error(USER, errmsg);
	}
	return(bmcnt);			/* ready to roll */
}

/* private function to convert 8-bit to float correction multiplier */
static void
srcskip_cvtY8(float *scorr, const uby8 *y8, int n)
{
	static float	gamtab[256];

	if (gamtab[0] < 0.5f) {		/* initialize lookup */
		int	i = 256;
		while (i--)
			gamtab[i] = 1./(1. - pow((i+.5)*(1./255.), 2.2));
	}
	while (n-- > 0)
		*scorr++ = gamtab[*y8++];
}

/* Read and convert the specified source skip scanline */
int
srcskip_getrow(int row, int *sndx, float *scorr)
{
	int	err;

	if ((0 > row) | (row >= sskip_dim[1]))
		goto badargum;
	if (sndx != NULL) {		/* read bitmap flags & convert? */
		uby8	*scanflags;
		SKIPBMP	*sbmp;
		int	x;
		if (skipin.sbmplist == NULL)
			goto inpclosed;
		errno = 0;
		for (sbmp = skipin.sbmplist; sbmp != NULL; sbmp = sbmp->next)
			if ((err = BMPseekScanline(row, sbmp->bmp)) != BIR_OK) {
				sprintf(errmsg, "%s '%s'",
						BMPerrorMessage(err), sbmp->bname);
				error(SYSTEM, errmsg);
			}
					/* per-column source skip flags */
		scanflags = (uby8 *)ecalloc(sskip_dim[0], SSKIPFLSIZ);
		for (sbmp = skipin.sbmplist; sbmp != NULL; sbmp = sbmp->next) {
			const uby8	*bscn = sbmp->bmp->scanline;
			for (x = 0; x < sskip_dim[0]; bscn += !(++x & 7))
				if (!(*bscn & 0x80>>(x&7)))
					sskip_addflags(scanflags + x*SSKIPFLSIZ,
							sbmp->sfl);
		}
					/* convert to lookup indices */
		for (x = sskip_dim[0]; x-- > 0; )
			sndx[x] = sskip_rsi(scanflags + x*SSKIPFLSIZ);
		efree(scanflags);
	}
	if (scorr == NULL)		/* all done? */
		return(row);
	switch (skipin.cftype) {	/* else read correction row */
	case CFfloatY:
		if (skipin.cf.Y.fd < 0)
			goto inpclosed;
		if (pread(skipin.cf.Y.fd, scorr, sizeof(float)*sskip_dim[0],
				skipin.cf.Y.dstart + sizeof(float)*sskip_dim[0]*row)
					!= sizeof(float)*sskip_dim[0]) {
			sprintf(errmsg, "read error from '%s'", skipin.cfname);
			error(SYSTEM, errmsg);
		}
		return(row);
	case CFbmpY8:
		if (skipin.cf.bmp == NULL)
			goto inpclosed;
		err = BMPseekScanline(row, skipin.cf.bmp);
		if (err != BIR_OK) {
			sprintf(errmsg, "%s '%s'",
					BMPerrorMessage(err), skipin.cfname);
			error(SYSTEM, errmsg);
		}
		srcskip_cvtY8(scorr, skipin.cf.bmp->scanline, sskip_dim[0]);
		return(row);
	case CFnone:			/* caller asking for missing input */
		break;
	}
inpclosed:
	error(CONSISTENCY, "call to srcskip_getrow() on closed input");
badargum:
	error(CONSISTENCY, "bad argument in srcskip_readrow()");
	return(EOF);	/* pro forma */
}

/* Close input images and free memory (leaving any maps) */
void
srcskip_close()
{
	while (skipin.sbmplist != NULL) {
		SKIPBMP	*sbmp = skipin.sbmplist;
		skipin.sbmplist = sbmp->next;
		BMPcloseInput(sbmp->bmp);
		sskip_free(sbmp->sfl);
		efree(sbmp);
	}
	switch (skipin.cftype) {
	case CFfloatY:
		if (skipin.cf.Y.fd >= 0) {
			close(skipin.cf.Y.fd);
			skipin.cf.Y.fd = -1;
		}
		break;
	case CFbmpY8:
		if (skipin.cf.bmp != NULL) {
			BMPcloseInput(skipin.cf.bmp);
			skipin.cf.bmp = NULL;
		}
		break;
	}
	skipin.cfname[0] = '\0';
	if (ssf_select == NULL) {	/* freed lookup table already? */
		srcskip_free_maps();
		skipin.cftype = CFnone;
		sskip_dim[0] = sskip_dim[1] = 0;
	}
}

#if defined(_WIN32) || defined(_WIN64)

/* Allocate and load entire source skip index array */
int *
srcskip_ndxmap()
{
	int	y;

	if (ssf_select == NULL) {	/* rug pulled from under us? */
		srcskip_free_maps();
		return(NULL);
	}
	if (skipin.ndxmap != NULL)
		return(skipin.ndxmap);

	skipin.ndxmap = (int *)emalloc(sizeof(int) *
					sskip_dim[0]*sskip_dim[1]);

	for (y = 0; y < sskip_dim[1]; y++)
		srcskip_getrow(y, skipin.ndxmap + sskip_dim[0]*y, NULL);

	return(skipin.ndxmap);
}

/* Allocate and load entire source correction array */
float *
srcskip_corrmap()
{
	int	y;
	ssize_t	nbytes;

	if (ssf_select == NULL) {	/* rug pulled from under us? */
		srcskip_free_maps();
		return(NULL);
	}
	if (skipin.corrmap != NULL)
		return(skipin.corrmap);

	if (skipin.cftype == CFnone)
		return(NULL);

	nbytes = sizeof(float)*sskip_dim[0]*sskip_dim[1];

	skipin.corrmap = (float *)emalloc(nbytes);

	switch (skipin.cftype) {
	case CFfloatY:
		if (pread(skipin.cf.Y.fd, skipin.corrmap, nbytes,
				skipin.cf.Y.dstart) != nbytes) {
			sprintf(errmsg, "read error from '%s'", skipin.cfname);
			error(SYSTEM, errmsg);
		}
		break;
	case CFbmpY8:
		for (y = 0; y < sskip_dim[1]; y++)
			srcskip_getrow(y, NULL, skipin.corrmap + sskip_dim[0]*y);
		break;
	}
	skipin.cref = 1;
	return(skipin.corrmap);
}

/* Free allocated memory for source skip index and correction maps */
void
srcskip_free_maps()
{
	efree(skipin.ndxmap); skipin.ndxmap = NULL;
	efree(skipin.corrmap); skipin.corrmap = NULL;
}

#else	/* ! Windows */

#include  <sys/mman.h>			/* mmap() support */

/* Create memory map with entire source skip index array */
int *
srcskip_ndxmap()
{
	int	y;

	if (ssf_select == NULL) {	/* rug pulled from under us? */
		srcskip_free_maps();
		return(NULL);
	}
	if (skipin.ndxmap != NULL)
		return(skipin.ndxmap);

	skipin.ndxmap = (int *)mmap(NULL, sizeof(int)*sskip_dim[0]*sskip_dim[1],
				PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE, -1, 0);

	if ((void *)skipin.ndxmap == MAP_FAILED)
		error(SYSTEM, "out of memory in srcskip_ndxmap()");

	for (y = 0; y < sskip_dim[1]; y++)
		srcskip_getrow(y, skipin.ndxmap + sskip_dim[0]*y, NULL);

	return(skipin.ndxmap);
}

/* Create memory map with entire source skip correction array */
float *
srcskip_corrmap()
{
	int	y;
	ssize_t	nbytes;

	if (ssf_select == NULL) {	/* rug pulled from under us? */
		srcskip_free_maps();
		return(NULL);
	}
	if (skipin.corrmap != NULL)
		return(skipin.corrmap);

	if (skipin.cftype == CFnone)
		return(NULL);

	nbytes = sizeof(float)*sskip_dim[0]*sskip_dim[1];

	if (skipin.cftype != CFfloatY || skipin.cf.Y.dstart % sizeof(float)) {
		skipin.corrmap = (float *)mmap(NULL, nbytes,
					PROT_READ|PROT_WRITE,
					MAP_ANON|MAP_PRIVATE, -1, 0);
		if ((void *)skipin.corrmap == MAP_FAILED)
			error(SYSTEM, "out of memory in srcskip_corrmap()");
	}
	switch (skipin.cftype) {
	case CFfloatY:
		if (skipin.corrmap != NULL) {
			if (pread(skipin.cf.Y.fd, skipin.corrmap, nbytes,
					skipin.cf.Y.dstart) != nbytes) {
				sprintf(errmsg, "read error from '%s'",
						skipin.cfname);
				error(SYSTEM, errmsg);
			}
			break;		/* best we could do */
		}			/* else map directly to file data */
		skipin.corrmap = (float *)mmap(NULL, skipin.cf.Y.dstart + nbytes,
					PROT_READ|PROT_WRITE,
					MAP_FILE|MAP_PRIVATE, skipin.cf.Y.fd, 0);
		if ((void *)skipin.corrmap == MAP_FAILED)
			error(SYSTEM, "cannot map file in srcskip_corrmap()");
		skipin.corrmap += skipin.cf.Y.dstart/sizeof(float);
		break;
	case CFbmpY8:
		for (y = 0; y < sskip_dim[1]; y++)
			srcskip_getrow(y, NULL, skipin.corrmap + sskip_dim[0]*y);
		break;
	}
	return(skipin.corrmap);
}

/* Unmap memory allocated for source skip index and correction arrays */
void
srcskip_free_maps()
{
	if (skipin.ndxmap != NULL) {
		munmap(skipin.ndxmap, sizeof(int)*sskip_dim[0]*sskip_dim[1]);
		skipin.ndxmap = NULL;
	}
	if (skipin.corrmap != NULL) {
		off_t	headlen = 0;
		if (skipin.cftype == CFfloatY &&
				!(skipin.cf.Y.dstart % sizeof(float)))
			headlen = skipin.cf.Y.dstart;

		munmap((char *)skipin.corrmap - headlen,
				headlen + sizeof(float)*sskip_dim[0]*sskip_dim[1]);
		skipin.corrmap = NULL;
	}
}

#endif 	/* ! Windows */
#endif	/* SSKIPOPT */
