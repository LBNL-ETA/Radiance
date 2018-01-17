#ifndef lint
static const char	RCSid[] = "$Id: strlcpy.c,v 2.1 2018/01/17 22:36:44 greg Exp $";
#endif
/*
 * String copy routines similar to strncpy
 */

#include "copyright.h"

int
strlcpy(char *dst, const char *src, int siz)
{
	int	n = siz;

	while (--n > 0)
		if (!(*dst++ = *src++))
			return(siz-1-n);
	*dst = '\0';
	return(siz-1);
}

int
strlcat(char *dst, const char *src, int siz)
{
	int	n = siz;

	while (*dst && --n > 0)
		++dst;
	if (n <= 0)
		return(siz);
	return(siz-n + strlcpy(dst, src, n));
}
