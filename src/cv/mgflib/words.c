#ifndef lint
static const char	RCSid[] = "$Id: words.c,v 1.3 2003/02/28 20:11:30 greg Exp $";
#endif
/*
 * Routines for recognizing and moving about words in strings.
 */

#include  <string.h>
#include  <ctype.h>

#ifdef  BSD
#define  strchr		index
#endif

#define  NULL		0


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


int
isname(s)			/* check for legal identifier name */
register char  *s;
{
	while (*s == '_')			/* skip leading underscores */
		s++;
	if (!isascii(*s) || !isalpha(*s))	/* start with a letter */
		return(0);
	while (isascii(*++s) && isgraph(*s))	/* all visible characters */
		;
	return(*s == '\0');			/* ending in nul */
}
