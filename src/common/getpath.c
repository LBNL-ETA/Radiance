#ifndef lint
static const char	RCSid[] = "$Id: getpath.c,v 2.9 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  getpath.c - function to search for file in a list of directories
 *
 *  External symbols declared in standard.h
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

#include  "standard.h"

#include  "paths.h"

#ifndef	 NIX
#include  <pwd.h>
extern struct passwd  *getpwnam();
#endif


char *
getpath(fname, searchpath, mode)	/* expand fname, return full path */
register char  *fname;
register char  *searchpath;
int  mode;
{
#ifndef	 NIX
	struct passwd  *pwent;
#endif
	static char  pname[MAXPATH];
	register char  *cp;

	if (fname == NULL)
		return(NULL);

	pname[0] = '\0';		/* check for full specification */
	switch (*fname) {
	CASEDIRSEP:				/* relative to root */
	case '.':				/* relative to cwd */
		strcpy(pname, fname);
		break;
#ifndef NIX
	case '~':				/* relative to home directory */
		fname++;
		if (*fname == '\0' || ISDIRSEP(*fname)) {	/* ours */
			if ((cp = getenv("HOME")) == NULL)
				return(NULL);
			strcpy(pname, cp);
			strcat(pname, fname);
			break;
		}
		cp = pname;					/* user */
		do
			*cp++ = *fname++;
		while (*fname && !ISDIRSEP(*fname));
		*cp = '\0';
		if ((pwent = getpwnam(pname)) == NULL)
			return(NULL);
		strcpy(pname, pwent->pw_dir);
		strcat(pname, fname);
		break;
#endif
	}
	if (pname[0])		/* got it, check access if search requested */
		return(searchpath==NULL||access(pname,mode)==0 ? pname : NULL);

	if (searchpath == NULL) {			/* don't search */
		strcpy(pname, fname);
		return(pname);
	}
							/* check search path */
	do {
		cp = pname;
		while (*searchpath && (*cp = *searchpath++) != PATHSEP)
			cp++;
		if (cp > pname && !ISDIRSEP(cp[-1]))
			*cp++ = DIRSEP;
		strcpy(cp, fname);
		if (access(pname, mode) == 0)		/* file accessable? */
			return(pname);
	} while (*searchpath);
							/* not found */
	return(NULL);
}
