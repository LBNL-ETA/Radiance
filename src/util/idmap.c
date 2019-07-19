#ifndef lint
static const char RCSid[] = "$Id: idmap.c,v 2.1 2019/07/19 22:25:03 greg Exp $";
#endif
/*
 * Routines for loading identifier maps
 *
 */

#include <stdlib.h>
#include "platform.h"
#include "rtio.h"
#include "idmap.h"


/* Open ID map file for reading, copying info to stdout based on hflags */
IDMAP *
idmap_ropen(const char *fname, int hflags)
{
	IDMAP	idinit, *idmp = NULL;
	char	infmt[MAXFMTLEN];
	int	gotfmt;
	int	tablength;
	int	i;

	if (!fname || !*fname)
		return NULL;

	memset(&idinit, 0, sizeof(IDMAP));

	idinit.finp = fopen(fname, "rb");

	if (!idinit.finp) {
		fputs(fname, stderr);
		fputs(": cannot open for reading\n", stderr);
		return NULL;
	}
	SET_FILE_BINARY(idinit.finp);
#ifdef getc_unlocked			/* avoid stupid semaphores */
	flockfile(idinit.finp);
#endif
	strcpy(infmt, IDMAPFMT);	/* read/copy header */
	gotfmt = checkheader(idinit.finp, infmt,
			(hflags & HF_HEADOUT) ? stdout : (FILE *)NULL);
	if (gotfmt <= 0) {
		fputs(fname, stderr);
		fputs(": bad or missing header format\n", stderr);
		fclose(idinit.finp);
		return NULL;
	}
	if (hflags & HF_HEADOUT) {
		fputformat("ascii", stdout);
		fputc('\n', stdout);
	}
	switch (atoi(infmt)) {		/* get index size */
	case 8:
		idinit.bytespi = 1;
		break;
	case 16:
		idinit.bytespi = 2;
		break;
	case 24:
		idinit.bytespi = 3;
		break;
	default:
		fputs(fname, stderr);
		fputs(": unsupported index size: ", stderr);
		fputs(infmt, stderr);
		fputc('\n', stderr);
		fclose(idinit.finp);
		return NULL;
	}
					/* get/copy resolution string */
	if (!fgetsresolu(&idinit.res, idinit.finp)) {
		fputs(fname, stderr);
		fputs(": missing map resolution\n", stderr);
		fclose(idinit.finp);
		return NULL;
	}
	if (hflags & HF_RESOUT)
		fputsresolu(&idinit.res, stdout);
					/* set up identifier table */
	idinit.dstart = ftell(idinit.finp);
	if (fseek(idinit.finp, 0, SEEK_END) < 0)
		goto seekerr;
	tablength = ftell(idinit.finp) - idinit.dstart -
			(long)idinit.res.xr*idinit.res.yr*idinit.bytespi;
	if (tablength < 2) {
		fputs(fname, stderr);
		fputs(": truncated file\n", stderr);
		fclose(idinit.finp);
		return NULL;
	}
					/* allocate and load table */
	idmp = (IDMAP *)malloc(sizeof(IDMAP)+tablength);
	if (idmp == NULL) {
		fputs("Out of memory allocating ID table!\n", stderr);
		fclose(idinit.finp);
		return NULL;
	}
	*idmp = idinit;			/* initialize allocated struct */
	idmp->idtable[tablength] = '\0';
	if (fseek(idmp->finp, -tablength, SEEK_END) < 0)
		goto seekerr;
	if (fread(idmp->idtable, 1, tablength, idmp->finp) != tablength) {
		fputs(fname, stderr);
		fputs(": error reading identifier table\n", stderr);
		fclose(idmp->finp);
		free(idmp);
		return NULL;
	}
	if (fseek(idmp->finp, idmp->curpos = idmp->dstart, SEEK_SET) < 0)
		goto seekerr;
	while (tablength > 0)		/* create index */
		idmp->nids += !idmp->idtable[--tablength];
	idmp->idoffset = (int *)malloc(sizeof(int)*idmp->nids);
	if (!idmp->idoffset) {
		fputs("Out of memory in idmap_ropen!\n", stderr);
		fclose(idmp->finp);
		free(idmp);
		return NULL;
	}
	for (i = 0; i < idmp->nids; i++) {
		idmp->idoffset[i] = tablength;
		while (idmp->idtable[tablength++])
			;
	}
	return idmp;
seekerr:
	fputs(fname, stderr);
	fputs(": cannot seek on file\n", stderr);
	fclose(idinit.finp);
	if (idmp) free(idmp);
	return NULL;
}


/* Read the next ID from input */
const char *
idmap_next(IDMAP *idmp)
{
	int	ndx = getint(idmp->bytespi, idmp->finp);

	if (ndx == EOF && feof(idmp->finp))
		return NULL;

	idmp->curpos += idmp->bytespi;

	ndx &= (1<<(idmp->bytespi<<3)) - 1;	/* undo sign extension */

	return mapID(idmp, ndx);
}


/* Seek to a specific pixel position in ID map */
int
idmap_seek(IDMAP *idmp, int x, int y)
{
	long	seekto;

	if ((x < 0) | (y < 0) | (x >= idmp->res.xr) | (y >= idmp->res.yr))
		return 0;
	
	seekto = idmp->dstart + ((long)y*idmp->res.xr + x)*idmp->bytespi;

	if (seekto != idmp->curpos && fseek(idmp->finp, seekto, SEEK_SET) < 0)
		return -1;

	idmp->curpos = seekto;
	
	return 1;
}


/* Read ID at a pixel position */
const char *
idmap_pix(IDMAP *idmp, int x, int y)
{
	if (idmap_seek(idmp, x, y) <= 0)
		return NULL;

	return idmap_next(idmp);
}


/* Close ID map and free resources */
void
idmap_close(IDMAP *idmp)
{
	fclose(idmp->finp);
	free(idmp->idoffset);
	free(idmp);
}
