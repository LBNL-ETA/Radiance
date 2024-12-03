#ifndef lint
static const char RCSid[] = "$Id: strnstr.c,v 3.1 2024/12/03 01:11:43 greg Exp $";
#endif
/*
 * String in string (limited search length)
 */

#include "copyright.h"
#include "rtio.h"

char *
strnstr(const char *haystack, const char *needle, size_t len)
{
	const size_t	nlen = strlen(needle);
	const char	*spos = haystack;
	int		i;

	if (!nlen) {		/* special case */
		while (len--) {
			if (!*spos)
				return((char *)spos);
			++spos;
		}
		return(NULL);	/* haystack is longer than len */
	}
	if (len < nlen)		/* not long enough for match! */
		return(NULL);
	len -= nlen;
	do {
		for (i = 0; i < nlen; i++)
			if (spos[i] != needle[i])
				break;

		if (i == nlen)	/* found match? */
			return((char *)spos);

	} while (((spos++)[i] != '\0') & (len-- != 0));

	return(NULL);		/* no match found */
}
