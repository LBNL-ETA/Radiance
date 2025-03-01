#ifndef lint
static const char RCSid[] = "$Id: strncasecmp.c,v 3.1 2025/03/01 00:51:26 greg Exp $";
#endif
/*
 * Replacement for strncasecmp() library call
 */

#include "rtio.h"
#include <ctype.h>

int
strncasecmp(const char *s1, const char *s2, size_t n)
{
	while (n-- > 0) {
		int	d = tolower(*s1++) - tolower(*s2);

		if (d | !*s2++) return(d);
	}
	return(0);
}
