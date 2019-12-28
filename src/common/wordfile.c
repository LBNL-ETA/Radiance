#ifndef lint
static const char	RCSid[] = "$Id: wordfile.c,v 2.24 2019/12/28 18:05:14 greg Exp $";
#endif
/*
 * Load whitespace separated words from a file into an array.
 *
 * External symbols declared in rtio.h
 */

#include "copyright.h"

#include <ctype.h>

#include "platform.h"
#include "rtio.h"
#include "rtmisc.h"


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
		int	dlen = n;
		int	crem = 0;
		if (n > MAXWLEN/2)		/* check for mid-word end */
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
		words += n; nargs -= n;
		wrdcnt += n;
		if ((n = crem) > 0)		/* move remainder */
			memmove(buf, buf+dlen-crem, crem);
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
