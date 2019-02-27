#ifndef lint
static const char	RCSid[] = "$Id: fdate.c,v 2.10 2019/02/27 21:30:01 greg Exp $";
#endif
/*
 * Return file date (UNIX seconds as returned by time(2) call)
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  "rtio.h"
#include  <sys/stat.h>
#if defined(_WIN32) || defined(_WIN64)
#include  <sys/utime.h>
#else
#include  <utime.h>
#endif


time_t
fdate(				/* get file date */
char  *fname
)
{
	struct stat  sbuf;

	if (stat(fname, &sbuf) == -1)
		return(0);

	return(sbuf.st_mtime);
}


int
setfdate(			/* set file date */
char  *fname,
long  ftim
)
{
	struct utimbuf utb;

	utb.actime = utb.modtime = ftim;

	return(utime(fname, &utb));
}
