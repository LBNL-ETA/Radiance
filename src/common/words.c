#ifndef lint
static const char	RCSid[] = "$Id: words.c,v 2.5 2003/02/22 02:07:23 greg Exp $";
#endif
/*
 * Routines for recognizing and moving about words in strings.
 *
 * External symbols declared in standard.h
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

#include  <ctype.h>
#include  <string.h>

#ifdef  BSD
#define  strchr		index
#endif


char *
atos(rs, nb, s)			/* get word from string, returning rs */
char  *rs;
register int  nb;
register char  *s;
{
	register char  *cp = rs;

	while (isspace(*s))
		s++;
	while (--nb > 0 && *s && !isspace(*s))
		*cp++ = *s++;
	*cp = '\0';
	return(rs);
}


char *
nextword(cp, nb, s)		/* get (quoted) word, returning new s */
register char  *cp;
register int  nb;
register char  *s;
{
	int	quote = 0;

	if (s == NULL) return(NULL);
	while (isspace(*s))
		s++;
	switch (*s) {
	case '\0':
		return(NULL);
	case '"':
	case '\'':
		quote = *s++;
	}
	while (--nb > 0 && *s && (quote ? *s!=quote : !isspace(*s)))
		*cp++ = *s++;
	*cp = '\0';
	if (quote && *s==quote)
		s++;
	return(s);
}


char *
sskip(s)			/* skip word in string, leaving on space */
register char  *s;
{
	while (isspace(*s))
		s++;
	while (*s && !isspace(*s))
		s++;
	return(s);
}


char *
sskip2(s, n)			/* skip word(s) in string, leaving on word */
register char  *s;
register int	n;
{
	while (isspace(*s))
		s++;
	while (n-- > 0) {
		while (*s && !isspace(*s))
			s++;
		while (isspace(*s))
			s++;
	}
	return(s);
}


char *
iskip(s)			/* skip integer in string */
register char  *s;
{
	while (isspace(*s))
		s++;
	if (*s == '-' || *s == '+')
		s++;
	if (!isdigit(*s))
		return(NULL);
	do
		s++;
	while (isdigit(*s));
	return(s);
}


char *
fskip(s)			/* skip float in string */
register char  *s;
{
	register char  *cp;

	while (isspace(*s))
		s++;
	if (*s == '-' || *s == '+')
		s++;
	cp = s;
	while (isdigit(*cp))
		cp++;
	if (*cp == '.') {
		cp++; s++;
		while (isdigit(*cp))
			cp++;
	}
	if (cp == s)
		return(NULL);
	if (*cp == 'e' || *cp == 'E')
		return(iskip(cp+1));
	return(cp);
}


int
isint(s)			/* check integer format */
char  *s;
{
	register char  *cp;

	cp = iskip(s);
	return(cp != NULL && *cp == '\0');
}


int
isintd(s, ds)			/* check integer format with delimiter set */
char  *s, *ds;
{
	register char  *cp;

	cp = iskip(s);
	return(cp != NULL && strchr(ds, *cp) != NULL);
}


int
isflt(s)			/* check float format */
char  *s;
{
	register char  *cp;

	cp = fskip(s);
	return(cp != NULL && *cp == '\0');
}


int
isfltd(s, ds)			/* check integer format with delimiter set */
char  *s, *ds;
{
	register char  *cp;

	cp = fskip(s);
	return(cp != NULL && strchr(ds, *cp) != NULL);
}
