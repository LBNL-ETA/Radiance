#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Read a CSV file into general struct
 *
 *	2/9/2026	G. Ward
 */

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "readCSV.h"

#ifdef getc_unlocked            /* avoid horrendous overhead of flockfile */
#undef getc
#define getc    getc_unlocked
#endif

/* Allocate a CVS record with the specified number of fields */
CSVREC	*
alloc_csvrec(int len)
{
	CSVREC	*rp;

	if (len <= 0) return NULL;
	rp = (CSVREC *)calloc(1,
			sizeof(CSVREC) + (len-1)*sizeof(char *));
	if (!rp) return NULL;
	rp->nf = len;
	return rp;
}

/* Reallocate a CVS record to the specified number of fields */
CSVREC	*
realloc_csvrec(CSVREC *rp, int newlen)
{
	if (!rp)
		return alloc_csvrec(newlen);
	if (newlen == rp->nf)
		return rp;
	if (newlen <= 0) {
		if (rp->next)
			fputs("Warning - deleting list with first record\n",
					stderr);
		free_csv(rp);
		return NULL;
	}
	while (rp->nf > newlen)
		if (rp->f[--rp->nf] && *rp->f[rp->nf])
			free(rp->f[rp->nf]);
	rp = (CSVREC *)realloc(rp,
			sizeof(CSVREC) + (newlen-1)*sizeof(char *));
	if (!rp) return NULL;
	memset(rp->f+rp->nf, 0, sizeof(char *)*(newlen-rp->nf));
	rp->nf = newlen;
	return rp;
}

/* Allocate and copy a string field (or clear if NULL) */
int
set_csvfield(CSVREC *rp, int i, const char *fstr)
{
	int	nbytes;

	if (!rp | (i < 0) || i >= rp->nf)
		return 0;
	if (rp->f[i] && *rp->f[i])
		free(rp->f[i]);
	if (!fstr) {
		rp->f[i] = NULL;
		return -1;
	}
	nbytes = strlen(fstr);
	if (nbytes++) {
		rp->f[i] = (char *)malloc(nbytes);
		if (!rp->f[i])
			return 0;
		memcpy(rp->f[i], fstr, nbytes);
	} else
		rp->f[i] = "";

	return nbytes;
}

/* Free a CVS record or list of records */
void
free_csv(CSVREC *rp)
{
	while (rp) {
		CSVREC	*rpn = rp->next;
		while (rp->nf > 0)
			if (rp->f[--rp->nf] && *rp->f[rp->nf])
				free(rp->f[rp->nf]);
		free(rp);
		rp = rpn;
	}
}

/* Read a CVS record from the given stream */
CSVREC	*
read_csvrec(FILE *fp, CSVREC *toappend)
{
	static int	nguess = 64;
	int		nfread = 0;
	int		flen=0, fterm=0;
	char		field[1024];
	int		c;
	CSVREC		*rp;

	if (!fp || (c = getc(fp)) == EOF)
		return NULL;
	if (toappend) {			/* go to end of list */
		while (toappend->next)
			toappend = toappend->next;
		nguess = toappend->nf;
	}
	rp = alloc_csvrec(nguess);	/* start with probable length */
	if (!rp) return NULL;
	for ( ; ; c = getc(fp)) {	/* loop to end of line */
		switch (c) {
		case '\t':
		case ' ':		/* skip leading white space */
			if (flen) {	/* trailing white leaves fterm */
				field[flen++] = c;
				if (flen >= sizeof(field))
					goto field2big;
			}
			continue;
		case '"':		/* quoted string */
			while ((c = getc(fp)) != EOF && c != '"') {
				if (c == '\\' && (c = getc(fp)) == EOF)
					break;
				field[flen++] = c;
				if (flen >= sizeof(field))
					goto field2big;
			}
			if (c != '"') {
				fputs("Missing end-quote\n", stderr);
				free_csv(rp);
				errno = EINVAL;
				return NULL;
			}
			fterm = flen;	/* don't trim white space */
			continue;
		case '\r':		/* non-Unix EOL */
			if ((c = getc(fp)) != EOF && c != '\n') {
				ungetc(c, fp);
				c = '\n';	/* still counts as EOL */
			}
			/* fall through */
		case EOF:
		case '\n':
		case ',':		/* end of this field for sure */
			if (nfread >= rp->nf) {
				rp = realloc_csvrec(rp, 2*rp->nf);
				if (!rp) return NULL;
			}
			field[fterm] = '\0';
			if (!set_csvfield(rp, nfread++, field)) {
				free_csv(rp);
				return NULL;
			}
			if (c != ',')	/* no more fields? */
				break;
			flen = fterm = 0;
			continue;
		case '\\':		/* literal next */
			if ((c = getc(fp)) == EOF)
				break;
			/* fall through */
		default:		/* ordinary or escaped char */
			field[flen++] = c;
			if (flen >= sizeof(field))
				goto field2big;
			fterm = flen;	/* doesn't count as white space */
			continue;
		}
		break;
	}				/* fix record length */
	rp = realloc_csvrec(rp, nfread);
	if (!rp) return NULL;
	if (toappend)			/* append list if given */
		toappend->next = rp;
	else if (nguess != nfread)
		nguess = nfread;
	return rp;			/* return the new record */
field2big:
	fprintf(stderr, "Field exceeds %ld bytes\n", sizeof(field));
	free_csv(rp);
	errno = EOVERFLOW;
	return NULL;
}

/* Read the named CVS file into an allocated record list */
CSVREC	*
read_csvfile(char *fname)
{
	FILE	*fp = stdin;
	int	lineno = 0;
	CSVREC	*csv = NULL;
	CSVREC	*rp;
	int	ok;

	if (fname && *fname) {
		errno = 0;
		fp = fopen(fname, "r");
		if (!fp) goto loaderr;
	} else
		fname = "<stdin>";
#ifdef getc_unlocked
	flockfile(fp);			/* avoid stupid semaphores */
#endif
	lineno = 1;
	errno = 0;
	rp = csv = read_csvrec(fp, NULL);
	while (rp) {
		++lineno;
		errno = 0;
		rp = read_csvrec(fp, rp);
	}
	ok = csv && feof(fp);
	if (fp != stdin)
		fclose(fp);
#ifdef getc_unlocked
	else if (ok)
		funlockfile(stdin);
#endif
	if (ok) return csv;
loaderr:
	if (lineno) fprintf(stderr, "At line %d in ", lineno);
	perror(fname);
	free_csv(csv);
	return NULL;
}
