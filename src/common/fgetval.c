#ifndef lint
static const char	RCSid[] = "$Id: fgetval.c,v 2.6 2003/06/27 06:53:21 greg Exp $";
#endif
/*
 * Read white space separated values from stream
 *
 *  External symbols declared in rtio.h
 */

#include  "rtio.h"

#include  <math.h>
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
