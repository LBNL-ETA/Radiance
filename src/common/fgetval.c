#ifndef lint
static const char	RCSid[] = "$Id: fgetval.c,v 2.3 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Read white space separated values from stream
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

#include  <stdio.h>

#include  <math.h>

#include  <stdlib.h>

#include  <string.h>

#include  <ctype.h>



int
fgetval(fp, ty, vp)			/* get specified data word */
register FILE	*fp;
int	ty;
char	*vp;
{
	char	wrd[64];
	register char	*cp;
	register int	c;
					/* elide comments (# to EOL) */
	do {
		while ((c = getc(fp)) != EOF && isspace(c))
			;
		if (c == '#')
			while ((c = getc(fp)) != EOF && c != '\n')
				;
	} while (c == '\n');
	if (c == EOF)
		return(EOF);
					/* get word */
	cp = wrd;
	do {
		*cp++ = c;
		if (cp - wrd >= sizeof(wrd))
			return(0);
	} while ((c = getc(fp)) != EOF && !isspace(c) && c != '#');
	if (c != EOF)
		ungetc(c, fp);
	*cp = '\0';
	switch (ty) {			/* check and convert it */
	case 'h':			/* short */
		if (!isint(wrd))
			return(0);
		*(short *)vp = c = atoi(wrd);
		if (*(short *)vp != c)
			return(0);
		return(1);
	case 'i':			/* integer */
		if (!isint(wrd))
			return(0);
		*(int *)vp = atoi(wrd);
		return(1);
	case 'l':			/* long */
		if (!isint(wrd))
			return(0);
		*(long *)vp = atol(wrd);
		return(1);
	case 'f':			/* float */
		if (!isflt(wrd))
			return(0);
		*(float *)vp = atof(wrd);
		return(1);
	case 'd':			/* double */
		if (!isflt(wrd))
			return(0);
		*(double *)vp = atof(wrd);
		return(1);
	case 's':			/* string */
		strcpy(vp, wrd);
		return(1);
	default:			/* unsupported type */
		return(0);
	}
}
