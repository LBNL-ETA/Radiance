#ifndef lint
static const char	RCSid[] = "$Id: wordfile.c,v 2.16 2016/03/21 19:06:08 greg Exp $";
#endif
/*
 * Load whitespace separated words from a file into an array.
 * Assume the passed pointer array is big enough to hold them all.
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "platform.h"
#include "standard.h"


#ifndef MAXFLEN
#define MAXFLEN		204800	/* file must be smaller than this */
#endif


int
wordfile(			/* get words from fname, put in words */
	char	**words,
	int	nargs,
	char	*fname
)
{
	int	fd;
	char	buf[MAXFLEN];
	int	n;
					/* load file into buffer */
	if (fname == NULL)
		return(-1);			/* no filename */
	if ((fd = open(fname, 0)) < 0)
		return(-1);			/* open error */
	n = read(fd, buf, MAXFLEN);
	close(fd);
	if (n < 0)				/* read error */
		return(-1);
	if (n == MAXFLEN)		/* file too big, take what we can */
		while (!isspace(buf[--n]))
			if (n <= 0)		/* one long word! */
				return(-1);
	buf[n] = '\0';				/* terminate */
	return(wordstring(words, nargs, buf));	/* wordstring does the rest */
}


int
wordstring(				/* allocate and load argument list */
	char	**avl,
	int	nargs,
	char	*str
)
{
	char	*cp, **ap;
	
	if (str == NULL)
		return(-1);
	cp = bmalloc(strlen(str)+1);
	if (cp == NULL)			/* ENOMEM */
		return(-1);
	strcpy(cp, str);
					/* parse into words */
	for (ap = avl; ap-avl < nargs-1; ap++) {
		while (isspace(*cp))	/* nullify spaces */
			*cp++ = '\0';
		if (!*cp)		/* all done? */
			break;
		*ap = cp;		/* add argument to list */
		while (*++cp && !isspace(*cp))
			;
	}
	*ap = NULL;
	return(ap - avl);
}
