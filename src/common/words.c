#ifndef lint
static const char	RCSid[] = "$Id: words.c,v 2.6 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Routines for recognizing and moving about words in strings.
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

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
