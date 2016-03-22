#ifndef lint
static const char	RCSid[] = "$Id: wordfile.c,v 2.17 2016/03/22 15:14:39 greg Exp $";
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


#define MAXWLEN		4096	/* words must be shorter than this */


int
wordfile(			/* get words from fname, put in words */
	char	**words,
	int	nargs,
	char	*fname
)
{
	int	wrdcnt = 0;
	int	n = 0;
	int	fd;
	char	buf[MAXWLEN];
					/* load file into buffer */
	if (fname == NULL || !*fname)
		return(-1);			/* no filename */
	if (nargs <= 1)
		return(-1);
	if ((fd = open(fname, 0)) < 0)
		return(-1);			/* open error */
	words[0] = NULL;
	while (nargs > 1 && (n += read(fd, buf+n, MAXWLEN-n)) > 0) {
		int	crem = 0;
		if (n >= MAXWLEN)		/* still something left? */
			while (!isspace(buf[--n])) {
				if (n <= 0)	/* one long word! */
					goto done;
				++crem;
			}
		buf[n] = '\0';			/* terminate & parse */
		n = wordstring(words, nargs, buf);
		if (n < 0) {
			wrdcnt = -1;		/* memory error */
			break;
		}
		words += n;
		nargs -= n;
		wrdcnt += n;
		if ((n = crem) > 0)		/* move remainder */
			strncpy(buf, buf+MAXWLEN-crem, crem);
	}
done:
	close(fd);
	return(wrdcnt);
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
	*cp = '\0';			/* terminates overflow */
	*ap = NULL;
	return(ap - avl);
}
