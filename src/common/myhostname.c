#ifndef lint
static const char	RCSid[] = "$Id: myhostname.c,v 2.7 2003/10/27 10:19:31 schorsch Exp $";
#endif
/*
 * Query system for host name
 */

#include "copyright.h"

#ifdef _WIN32
  #include <winsock2.h>
#else
  #include <unistd.h>
#endif

#include "rtmisc.h"

extern char *
myhostname()
{
	static char	hostname[65];

	if (!hostname[0])
		gethostname(hostname, sizeof(hostname));
	return(hostname);
}

