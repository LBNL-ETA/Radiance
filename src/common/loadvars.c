#ifndef lint
static const char	RCSid[] = "$Id: loadvars.c,v 2.9 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  Routines for loading and checking variables from file.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "vars.h"

#define NOCHAR	127		/* constant for character to delete */

extern char  *fgetline();


void
loadvars(rfname)		/* load variables into vv from file */
char	*rfname;
{
	FILE	*fp;
	char	buf[512];
	register char	*cp;

	if (rfname == NULL)
		fp = stdin;
	else if ((fp = fopen(rfname, "r")) == NULL) {
		perror(rfname);
		quit(1);
	}
	while (fgetline(buf, sizeof(buf), fp) != NULL) {
		for (cp = buf; *cp; cp++) {
			switch (*cp) {
			case '\\':
				*cp++ = NOCHAR;
				continue;
			case '#':
				*cp = '\0';
				break;
			default:
				continue;
			}
			break;
		}
		if (setvariable(buf, matchvar) < 0) {
			fprintf(stderr, "%s: unknown variable: %s\n",
					rfname, buf);
			quit(1);
		}
	}
	fclose(fp);
}


int
setvariable(ass, mv)		/* assign variable according to string */
register char	*ass;
VARIABLE	*(*mv)();
{
	char	varname[32];
	int	n;
	register char	*cp;
	register VARIABLE	*vp;
	register int	i;

	while (isspace(*ass))		/* skip leading space */
		ass++;
	cp = varname;			/* extract name */
	while (cp < varname+sizeof(varname)-1
			&& *ass && !isspace(*ass) && *ass != '=')
		*cp++ = *ass++;
	*cp = '\0';
	if (!varname[0])
		return(0);	/* no variable name! */
					/* trim value */
	while (isspace(*ass) || *ass == '=')
		ass++;
	for (n = strlen(ass); n > 0; n--)
		if (!isspace(ass[n-1]))
			break;
	if (!n)
		return(0);	/* no assignment */
					/* match variable from list */
	vp = (*mv)(varname);
	if (vp == NULL)
		return(-1);
					/* assign new value */
	if (i = vp->nass) {
		cp = vp->value;
		while (i--)
			while (*cp++)
				;
		i = cp - vp->value;
		vp->value = (char *)realloc(vp->value, i+n+1);
	} else
		vp->value = (char *)malloc(n+1);
	if (vp->value == NULL) {
		perror(progname);
		quit(1);
	}
	cp = vp->value+i;		/* copy value, squeezing spaces */
	*cp = *ass;
	for (i = 1; i <= n; i++) {
		if (ass[i] == NOCHAR)
			continue;
		if (isspace(*cp))
			while (isspace(ass[i]))
				i++;
		*++cp = ass[i];
	}
	if (isspace(*cp))		/* remove trailing space */
		*cp = '\0';
	return(++vp->nass);
}


VARIABLE *
matchvar(nam)			/* match a variable by its name */
char	*nam;
{
	int	n = strlen(nam);
	register int	i;

	for (i = 0; i < NVARS; i++)
		if (n >= vv[i].nick && !strncmp(nam, vv[i].name, n))
			return(vv+i);
	return(NULL);
}


char *
nvalue(vn, n)			/* return nth variable value */
register int	vn;
register int	n;
{
	register char	*cp;

	if (vval(vn) == NULL | n < 0 | n >= vdef(vn))
		return(NULL);
	cp = vval(vn);
	while (n--)
		while (*cp++)
			;
	return(cp);
}


void
checkvalues()			/* check assignments */
{
	register int	i;

	for (i = 0; i < NVARS; i++)
		if (vv[i].fixval != NULL)
			(*vv[i].fixval)(vv+i);
}


void
onevalue(vp)			/* only one assignment for this variable */
register VARIABLE	*vp;
{
	if (vp->nass < 2)
		return;
	if (!nowarn)
		fprintf(stderr,
		"%s: warning - multiple assignment of variable '%s'\n",
			progname, vp->name);
	do
		vp->value += strlen(vp->value)+1;
	while (--vp->nass > 1);
}


void
catvalues(vp)			/* concatenate variable values */
register VARIABLE	*vp;
{
	register char	*cp;

	if (vp->nass < 2)
		return;
	for (cp = vp->value; vp->nass > 1; vp->nass--) {
		while (*cp)
			cp++;
		*cp++ = ' ';
	}
}


int
badmatch(tv, cv)		/* case insensitive truncated comparison */
register char	*tv, *cv;
{
	if (!*tv) return(1);		/* null string cannot match */
	do
		if (UPPER(*tv) != *cv++)
			return(1);
	while (*++tv);
	return(0);		/* OK */
}


void
boolvalue(vp)			/* check boolean for legal values */
register VARIABLE	*vp;
{
	if (!vp->nass) return;
	onevalue(vp);
	switch (UPPER(vp->value[0])) {
	case 'T':
		if (badmatch(vp->value, "TRUE")) break;
		return;
	case 'F':
		if (badmatch(vp->value, "FALSE")) break;
		return;
	}
	fprintf(stderr, "%s: illegal value for boolean variable '%s'\n",
			progname, vp->name);
	quit(1);
}


void
qualvalue(vp)			/* check qualitative var. for legal values */
register VARIABLE	*vp;
{
	if (!vp->nass) return;
	onevalue(vp);
	switch (UPPER(vp->value[0])) {
	case 'L':
		if (badmatch(vp->value, "LOW")) break;
		return;
	case 'M':
		if (badmatch(vp->value, "MEDIUM")) break;
		return;
	case 'H':
		if (badmatch(vp->value, "HIGH")) break;
		return;
	}
	fprintf(stderr, "%s: illegal value for qualitative variable '%s'\n",
			progname, vp->name);
	quit(1);
}


void
intvalue(vp)			/* check integer variable for legal values */
register VARIABLE	*vp;
{
	if (!vp->nass) return;
	onevalue(vp);
	if (isint(vp->value)) return;
	fprintf(stderr, "%s: illegal value for integer variable '%s'\n",
			progname, vp->name);
	quit(1);
}


void
fltvalue(vp)			/* check float variable for legal values */
register VARIABLE	*vp;
{
	if (!vp->nass) return;
	onevalue(vp);
	if (isflt(vp->value)) return;
	fprintf(stderr, "%s: illegal value for real variable '%s'\n",
			progname, vp->name);
	quit(1);
}


void
printvars(fp)			/* print variable values */
register FILE	*fp;
{
	int	i, j, k, clipline;
	register char	*cp;

	for (i = 0; i < NVARS; i++)		/* print each variable */
	    for (j = 0; j < vdef(i); j++) {	/* print each assignment */
		fputs(vnam(i), fp);
		fputs("= ", fp);
		k = clipline = ( vv[i].fixval == catvalues ? 64 : 236 )
				- strlen(vnam(i)) ;
		cp = nvalue(i, j);
		while (*cp) {
		    putc(*cp++, fp);
		    if (--k <= 0) {		/* line too long */
			while (*cp && !isspace(*cp))
			    fputc(*cp++, fp);	/* finish this word */
			if (*cp) {		/* start new line */
			    if (vv[i].fixval == catvalues) {
				fputc('\n', fp);
				fputs(vnam(i), fp);
				fputc('=', fp);
			    } else
				fputs(" \\\n", fp);
			    k = clipline;
			}
		    }
		}
	        fputc('\n', fp);
	    }
	fflush(fp);
}
