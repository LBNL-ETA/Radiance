#ifndef lint
static const char	RCSid[] = "$Id: myhostname.c,v 2.4 2003/07/14 20:02:29 schorsch Exp $";
#endif
/*
 * Query system for host name
 */

#include "copyright.h"

#include <unistd.h>

char *
myhostname()
{
	static char	hostname[65];

	if (!hostname[0])
		gethostname(hostname, sizeof(hostname));
	return(hostname);
}

