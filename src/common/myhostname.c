#ifndef lint
static const char	RCSid[] = "$Id: myhostname.c,v 2.6 2003/10/21 02:02:31 schorsch Exp $";
#endif
/*
 * Query system for host name
 */

#include "copyright.h"

#ifdef _WIN32
  #include <winsock.h>
#else
  #include <unistd.h>
#endif

#include "rtmisc.h"

char *
myhostname()
{
	static char	hostname[65];

	if (!hostname[0])
		gethostname(hostname, sizeof(hostname));
	return(hostname);
}

