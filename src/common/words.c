#ifndef lint
static const char RCSid[] = "$Id: words.c,v 2.13 2021/01/01 02:10:32 greg Exp $";
#endif
/*
 * Routines for recognizing and moving about words in strings.
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include  <ctype.h>

#include  "rtio.h"

char *
atos(char *rs, int nb, char *s)		/* get word from string, returning rs */
{
	char  *cp = rs;

	while (isspace(*s))
		s++;
	while (--nb > 0 && *s && !isspace(*s))
		*cp++ = *s++;
	*cp = '\0';
	return(rs);
}


char *
nextword(char *cp, int nb, char *s)	/* get (quoted) word, returning new s */
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
sskip(char *s)			/* skip word in string, leaving on space */
{
	while (isspace(*s))
		s++;
	while (*s && !isspace(*s))
		s++;
	return(s);
}


char *
sskip2(char *s, int n)		/* skip word(s) in string, leaving on word */
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
iskip(char *s)			/* skip integer in string */
{
	while (isspace(*s))
		s++;
	s += (*s == '-') | (*s == '+');
	if (!isdigit(*s))
		return(NULL);
	do
		s++;
	while (isdigit(*s));
	return(s);
}


char *
fskip(char *s)			/* skip float in string */
{
	char  *cp;

	while (isspace(*s))
		s++;
	s += (*s == '-') | (*s == '+');
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
	if ((*cp == 'e') | (*cp == 'E'))
		return(isspace(*++cp) ? NULL : iskip(cp));
	return(cp);
}


int
isint(char *s)			/* check integer format */
{
	char  *cp;

	cp = iskip(s);
	return(cp != NULL && *cp == '\0');
}


int
isintd(char *s, char *ds)	/* check integer format with delimiter set */
{
	char  *cp;

	cp = iskip(s);
	return(cp != NULL && strchr(ds, *cp) != NULL);
}


int
isflt(char *s)			/* check float format */
{
	char  *cp;

	cp = fskip(s);
	return(cp != NULL && *cp == '\0');
}


int
isfltd(char *s, char *ds)	/* check integer format with delimiter set */
{
	char  *cp;

	cp = fskip(s);
	return(cp != NULL && strchr(ds, *cp) != NULL);
}


int
isname(char *s)			/* check for legal identifier name */
{
	while (*s == '_')			/* skip leading underscores */
		s++;
	if (!isascii(*s) | !isalpha(*s))	/* start with a letter */
		return(0);
	while (isascii(*++s) && isgraph(*s))	/* all visible characters */
		;
	return(*s == '\0');			/* ending in nul */
}
