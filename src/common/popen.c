#ifndef lint
static const char	RCSid[] = "$Id: popen.c,v 2.6 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * popen() and pclose() calls for systems without pipe facilities
 */
 
#include "copyright.h"

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
