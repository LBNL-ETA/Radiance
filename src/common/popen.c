#ifndef lint
static const char	RCSid[] = "$Id: popen.c,v 2.5 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * popen() and pclose() calls for systems without pipe facilities
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
#include <string.h>
#include <ctype.h>
#include "paths.h"

#ifndef NFILE
#define NFILE		32
#endif

static struct {
	char	*f;		/* file name, NULL if inactive */
	char	*c;		/* command if writing */
} pips[NFILE];


FILE *
popen(cmd, mode)		/* open command for reading or writing */
register char	*cmd;
char	*mode;
{
	FILE	*fp;
	char	*newcmd, *fname;
	register char	*cp, *cp2 = NULL;
	int	quote = '\0', paren = 0;
	
	while (isspace(*cmd))
		cmd++;
	if (!*cmd)
		return(NULL);
	fname = (char *)malloc(TEMPLEN+1);
	newcmd = (char *)malloc(TEMPLEN+6+strlen(cmd));
	if (fname == NULL | newcmd == NULL)
		return(NULL);
	mktemp(strcpy(fname,TEMPLATE));
				/* build our command */
	for (cp = newcmd; ; cmd++) {
		switch (*cmd) {
		case '\'':
		case '"':
			if (!quote)
				quote = *cmd;
			else if (quote == *cmd)
				quote = '\0';
#ifdef MSDOS
			else
				break;
			*cp++ = '"';		/* double quotes only */
			continue;
		case '\\':
			if (!quote && cmd[1] == '\n') {
				*cp++ = ' ';
				cmd++;
				continue;
			}
			cmd++;
			break;
#else
			break;
		case '\\':
			*cp++ = *cmd++;
			break;
		case '(':
			if (!quote)
				paren++;
			break;
		case ')':
			if (!quote && paren)
				paren--;
			break;
		case ';':
#endif
		case '|':
			if (!quote && !paren && cp2 == NULL)
				cp2 = fname;
			break;
		case ' ':
		case '\t':
			if (!quote)
				while (isspace(cmd[1]))
					cmd++;
			break;
		case '\n':
		case '\0':
			if (cp2 == NULL)
				cp2 = fname;
			break;
		}
		if (cp2 == fname && *mode == 'w') {	/* add input file */
			*cp++ = ' '; *cp++ = '<';
			while (*cp2)
				*cp++ = *cp2++;
			*cp++ = ' ';
		}
		if (!(*cp = *cmd))
			break;
		cp++;
	}
	if (*mode == 'r') {			/* add output file */
		while (isspace(cp[-1]) || cp[-1] == ';')
			cp--;
		*cp++ = ' '; *cp++ = '>';
		while (*cp2)
			*cp++ = *cp2++;
		*cp = '\0';
		system(newcmd);			/* execute command first */
	}
	if ((fp = fopen(fname, mode)) == NULL) {	/* open file */
		free(fname);
		free(newcmd);
		return(NULL);
	}
	if (fileno(fp) >= NFILE) {
		eputs("popen: too many open files\n");
		quit(1);
	}
	pips[fileno(fp)].f = fname;
	if (*mode == 'r') {
		free(newcmd);
		pips[fileno(fp)].c = NULL;
	} else
		pips[fileno(fp)].c = newcmd;
	return(fp);
}


pclose(fp)			/* close pipe and cleanup */
FILE	*fp;
{
	register int	pn = fileno(fp);

	if (pn >= NFILE || pips[pn].f == NULL || fclose(fp) == -1)
		return(-1);
	if (pips[pn].c != NULL) {
		system(pips[pn].c);		/* execute command last */
		free(pips[pn].c);
	}
	unlink(pips[pn].f);		/* remove temp file */
	free(pips[pn].f);
	pips[pn].f = NULL;
	return(0);
}
