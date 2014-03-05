#ifndef lint
static const char RCSid[] = "$Id: eplus_idf.c,v 2.12 2014/03/05 19:49:39 greg Exp $";
#endif
/*
 *  eplus_idf.c
 *
 *  EnergyPlus Input Data File i/o routines
 *
 *  Created by Greg Ward on 1/31/14.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "eplus_idf.h"

#ifdef getc_unlocked            /* avoid horrendous overhead of flockfile */
#undef getc
#define getc    getc_unlocked
#endif

/* Create a new object with empty field list (comment optional) */
IDF_OBJECT *
idf_newobject(IDF_LOADED *idf, const char *pname, const char *comm,
			IDF_OBJECT *prev)
{
	LUENT		*pent;
	IDF_OBJECT	*pnew;

	if ((idf == NULL) | (pname == NULL))
		return(NULL);
	if (comm == NULL) comm = "";
	pent = lu_find(&idf->ptab, pname);
	if (pent == NULL)
		return(NULL);
	if (pent->key == NULL) {	/* new object name/type? */
		pent->key = (char *)malloc(strlen(pname)+1);
		if (pent->key == NULL)
			return(NULL);
		strcpy(pent->key, pname);
	}
	pnew = (IDF_OBJECT *)malloc(sizeof(IDF_OBJECT)+strlen(comm));
	if (pnew == NULL)
		return(NULL);
	strcpy(pnew->rem, comm);
	pnew->nfield = 0;
	pnew->flist = NULL;
	pnew->pname = pent->key;	/* add to table */
	pnew->pnext = (IDF_OBJECT *)pent->data;
	pent->data = (char *)pnew;
	pnew->dnext = NULL;		/* add to file list */
	if (prev != NULL || (prev = idf->plast) != NULL) {
		pnew->dnext = prev->dnext;
		if (prev == idf->plast)
			idf->plast = pnew;
	}
	if (idf->pfirst == NULL)
		idf->pfirst = idf->plast = pnew;
	else
		prev->dnext = pnew;
	return(pnew);
}

/* Add a field to the given object and follow with the given text */
int
idf_addfield(IDF_OBJECT *param, const char *fval, const char *comm)
{
	int		fnum = 1;	/* returned argument number */
	IDF_FIELD	*fnew, *flast;
	char		*cp;

	if ((param == NULL) | (fval == NULL))
		return(0);
	if (comm == NULL) comm = "";
	fnew = (IDF_FIELD *)malloc(sizeof(IDF_FIELD)+strlen(fval)+strlen(comm));
	if (fnew == NULL)
		return(0);
	fnew->next = NULL;
	cp = fnew->val;			/* copy value and comments */
	while ((*cp++ = *fval++))
		;
	fnew->rem = cp;
	while ((*cp++ = *comm++))
		;
					/* add to object's field list */
	if ((flast = param->flist) != NULL) {
		++fnum;
		while (flast->next != NULL) {
			flast = flast->next;
			++fnum;
		}
	}
	if (flast == NULL)
		param->flist = fnew;
	else
		flast->next = fnew;
	param->nfield++;
	return(fnum);
}

/* Retrieve the indexed field from object (first field index is 1) */
IDF_FIELD *
idf_getfield(IDF_OBJECT *param, int fn)
{
	IDF_FIELD	*fld;

	if ((param == NULL) | (fn <= 0))
		return(NULL);
	fld = param->flist;
	while ((--fn > 0) & (fld != NULL))
		fld = fld->next;
	return(fld);
}

/* Delete the specified object from loaded IDF */
int
idf_delobject(IDF_LOADED *idf, IDF_OBJECT *param)
{
	LUENT		*pent;
	IDF_OBJECT	*pptr, *plast;

	if ((idf == NULL) | (param == NULL))
		return(0);
					/* remove from object table */
	pent = lu_find(&idf->ptab, param->pname);
	for (plast = NULL, pptr = (IDF_OBJECT *)pent->data;
				pptr != NULL; plast = pptr, pptr = pptr->pnext)
		if (pptr == param)
			break;
	if (pptr == NULL)
		return(0);
	if (plast == NULL)
		pent->data = (char *)param->pnext;
	else
		plast->pnext = param->pnext;
					/* remove from global list */
	for (plast = NULL, pptr = idf->pfirst;
				pptr != param; plast = pptr, pptr = pptr->dnext)
		if (pptr == NULL)
			return(0);
	if (plast == NULL)
		idf->pfirst = param->dnext;
	else
		plast->dnext = param->dnext;
	if (idf->plast == param)
		idf->plast = plast;
					/* free field list */
	while (param->flist != NULL) {
		IDF_FIELD	*fdel = param->flist;
		param->flist = fdel->next;
		free(fdel);
	}
	free(param);			/* free object struct */
	return(1);
}

/* Move the specified object to the given position in the IDF */
int
idf_movobject(IDF_LOADED *idf, IDF_OBJECT *param, IDF_OBJECT *prev)
{
	IDF_OBJECT	*pptr, *plast;

	if ((idf == NULL) | (param == NULL))
		return(0);
					/* quick check if already there */
	if (param == (prev==NULL ? idf->pfirst : prev->dnext))
		return(1);
					/* find in IDF list, first*/
	for (plast = NULL, pptr = idf->pfirst;
				pptr != param; plast = pptr, pptr = pptr->dnext)
		if (pptr == NULL)
			return(0);
	if (plast == NULL)
		idf->pfirst = param->dnext;
	else
		plast->dnext = param->dnext;
	if (idf->plast == param)
		idf->plast = plast;
	if (prev == NULL) {		/* means they want it at beginning */
		param->dnext = idf->pfirst;
		idf->pfirst = param;
	} else {
		param->dnext = prev->dnext;
		prev->dnext = param;
	}
	return(1);
}

/* Get a named object list */
IDF_OBJECT *
idf_getobject(IDF_LOADED *idf, const char *pname)
{
	if ((idf == NULL) | (pname == NULL))
		return(NULL);

	return((IDF_OBJECT *)lu_find(&idf->ptab,pname)->data);
}

/* Read an argument including terminating ',' or ';' -- return which */
static int
idf_read_argument(char *buf, FILE *fp, int trim)
{
	int	skipwhite = trim;
	char	*cp = buf;
	int	c;

	while ((c = getc(fp)) != EOF && (c != ',') & (c != ';')) {
		if (skipwhite && isspace(c))
			continue;
		skipwhite = 0;
		if (cp-buf < IDF_MAXARGL-1)
			*cp++ = c;
	}
	if (trim)
		while (cp > buf && isspace(cp[-1]))
			--cp;
	*cp = '\0';
	return(c);
}

/* Read a comment, including all white space up to next alpha character */
static void
idf_read_comment(char *buf, int len, FILE *fp)
{
	int	incomm = 0;
	char	*cp = buf;
	char	dummys[2];
	int	c;

	if ((buf == NULL) | (len <= 0)) {
		buf = dummys;
		len = sizeof(dummys);
	}
	while ((c = getc(fp)) != EOF &&
			(isspace(c) || (incomm += (c == '!')))) {
		if (c == '\n')
			incomm = 0;
		if (cp-buf < len-2)
			*cp++ = c;
		else if (cp-buf == len-2)
			*cp++ = '\n';
	}
	*cp = '\0';
	if (c != EOF)
		ungetc(c, fp);
}

/* Read a object and fields from an open file and add to end of list */
IDF_OBJECT *
idf_readobject(IDF_LOADED *idf, FILE *fp)
{
	char		abuf[IDF_MAXARGL], cbuf[100*IDF_MAXLINE];
	int		delim;
	IDF_OBJECT	*pnew;
	
	if ((delim = idf_read_argument(abuf, fp, 1)) == EOF)
		return(NULL);
	idf_read_comment(cbuf, sizeof(cbuf), fp);
	pnew = idf_newobject(idf, abuf, cbuf, NULL);
	while (delim == ',')
		if ((delim = idf_read_argument(abuf, fp, 1)) != EOF) {
			idf_read_comment(cbuf, sizeof(cbuf), fp);
			idf_addfield(pnew, abuf, cbuf);
		}
	if (delim != ';')
		fputs("Expected ';' at end of object list!\n", stderr);
	return(pnew);
}

/* Upper-case string hashing function */
static unsigned long
strcasehash(const char *s)
{
	char		strup[IDF_MAXARGL];
	char		*cdst = strup;

	while ((*cdst++ = toupper(*s++)))
		if (cdst >= strup+(sizeof(strup)-1)) {
			*cdst = '\0';
			break;
		}
	return(lu_shash(strup));
}

/* Initialize an IDF struct */
IDF_LOADED *
idf_create(const char *hdrcomm)
{
	IDF_LOADED	*idf = (IDF_LOADED *)calloc(1, sizeof(IDF_LOADED));

	if (idf == NULL)
		return(NULL);
	idf->ptab.hashf = &strcasehash;
	idf->ptab.keycmp = &strcasecmp;
	idf->ptab.freek = &free;
	lu_init(&idf->ptab, 200);
	if (hdrcomm != NULL && *hdrcomm) {
		idf->hrem = (char *)malloc(strlen(hdrcomm)+1);
		if (idf->hrem != NULL)
			strcpy(idf->hrem, hdrcomm);
	}
	return(idf);
}

/* Add comment(s) to header */
int
idf_add2hdr(IDF_LOADED *idf, const char *hdrcomm)
{
	int	olen, len;

	if ((idf == NULL) | (hdrcomm == NULL))
		return(0);
	len = strlen(hdrcomm);
	if (!len)
		return(0);
	if (idf->hrem == NULL)
		olen = 0;
	else
		olen = strlen(idf->hrem);
	if (olen)
		idf->hrem = (char *)realloc(idf->hrem, olen+len+1);
	else
		idf->hrem = (char *)malloc(len+1);
	if (idf->hrem == NULL)
		return(0);
	strcpy(idf->hrem+olen, hdrcomm);
	return(1);
}

/* Load an Input Data File */
IDF_LOADED *
idf_load(const char *fname)
{
	char		hdrcomm[256*IDF_MAXLINE];
	FILE		*fp;
	IDF_LOADED	*idf;

	if (fname == NULL)
		fp = stdin;		/* open file if not stdin */
	else if ((fp = fopen(fname, "r")) == NULL)
		return(NULL);
					/* read header comments */
	idf_read_comment(hdrcomm, sizeof(hdrcomm), fp);
	idf = idf_create(hdrcomm);	/* create IDF struct */
	if (idf == NULL)
		return(NULL);
					/* read each object */
	while (idf_readobject(idf, fp) != NULL)
		;
	if (fp != stdin)		/* close file if not stdin */
		fclose(fp);
	return(idf);			/* success! */
}

/* Check string for end-of-line */
static int
idf_hasEOL(const char *s)
{
	while (*s)
		if (*s++ == '\n')
			return(1);
	return(0);
}

/* Write a object and fields to an open file */
int
idf_writeparam(IDF_OBJECT *param, FILE *fp, int incl_comm)
{
	IDF_FIELD	*fptr;

	if ((param == NULL) | (fp == NULL))
		return(0);
	fputs(param->pname, fp);
	fputc(',', fp);
	if (incl_comm)
		fputs(param->rem, fp);
	else
		fputc('\n', fp);
	for (fptr = param->flist; fptr != NULL; fptr = fptr->next) {
		if (!incl_comm)
			fputs("    ", fp);
		fputs(fptr->val, fp);
		if (fptr->next == NULL) {
			fputc(';', fp);
			if (incl_comm && !idf_hasEOL(fptr->rem))
				fputc('\n', fp);
		} else
			fputc(',', fp);
		if (incl_comm)
			fputs(fptr->rem, fp);
		else
			fputc('\n', fp);
	}
	if (!incl_comm)
		fputc('\n', fp);
	return(!ferror(fp));
}

/* Write out an Input Data File */
int
idf_write(IDF_LOADED *idf, const char *fname, int incl_comm)
{
	FILE		*fp;
	IDF_OBJECT	*pptr;

	if (idf == NULL)
		return(0);
	if (fname == NULL)
		fp = stdout;		/* open file if not stdout */
	else if ((fp = fopen(fname, "w")) == NULL)
		return(0);
	if (incl_comm)
		fputs(idf->hrem, fp);	/* write header then parameters */
	for (pptr = idf->pfirst; pptr != NULL; pptr = pptr->dnext)
		if (!idf_writeparam(pptr, fp, incl_comm>0)) {
			fclose(fp);
			return(0);
		}
	if (fp == stdout)		/* flush/close file & check status */
		return(fflush(fp) == 0);
	return(fclose(fp) == 0);
}

/* Free a loaded IDF */
void
idf_free(IDF_LOADED *idf)
{
	if (idf == NULL)
		return;
	while (idf->pfirst != NULL) {
		IDF_OBJECT	*pdel = idf->pfirst;
		idf->pfirst = pdel->dnext;
		while (pdel->flist != NULL) {
			IDF_FIELD	*fdel = pdel->flist;
			pdel->flist = fdel->next;
			free(fdel);
		}
		free(pdel);
	}
	lu_done(&idf->ptab);
	if (idf->hrem != NULL)
		free(idf->hrem);
	free(idf);
}
