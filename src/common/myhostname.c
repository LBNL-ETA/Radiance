#ifndef lint
static const char	RCSid[] = "$Id: myhostname.c,v 2.5 2003/07/17 09:21:29 schorsch Exp $";
#endif
/*
 * Query system for host name
 */

#include "copyright.h"

#include <unistd.h>

#include "rtmisc.h"

char *
myhostname()
{
	static char	hostname[65];

	if (!hostname[0])
		gethostname(hostname, sizeof(hostname));
	return(hostname);
}

