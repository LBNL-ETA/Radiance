#ifndef lint
static const char	RCSid[] = "$Id: strlcpy.c,v 2.2 2018/05/04 21:22:44 greg Exp $";
#endif
/*
 * String copy routines similar to strncpy
 */

#include "copyright.h"
#include "rtio.h"

size_t
strlcpy(char *dst, const char *src, size_t siz)
{
	size_t	n = siz;

	while (--n > 0)
		if (!(*dst++ = *src++))
			return(siz-1-n);
	*dst = '\0';
	return(siz-1);
}

size_t
strlcat(char *dst, const char *src, size_t siz)
{
	size_t	n = siz;

	while (*dst && --n > 0)
		++dst;
	if (n <= 0)
		return(siz);
	return(siz-n + strlcpy(dst, src, n));
}
