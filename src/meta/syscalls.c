#ifndef lint
static const char	RCSid[] = "$Id: syscalls.c,v 1.5 2023/06/09 15:25:49 greg Exp $";
#endif
/*
 *  System calls for meta-file routines
 */

#include "rtprocess.h" /* getpid() */
#include "rterror.h"
#include  "meta.h"


FILE *
efopen(		/* open a file, report errors */
const char  *fname,
const char  *mode
)
{
 FILE  *fp;
 FILE  *fopen();

 if ((fp = fopen(fname, mode)) == NULL)  {
    sprintf(errmsg, "cannot open file \"%s\", mode \"%s\"", fname, mode);
    error(USER, errmsg);
    }

 return(fp);
 }



FILE *
mfopen(		/* open a program metafile */
const char  *fname,
const char  *mode
)
{
    char  *mdir, stemp[MAXFNAME];
    char  *getenv();

    if ((mdir = getenv("MDIR")) == NULL)
	mdir = MDIR;
    sprintf(stemp, "%s%s", mdir, fname);

    return(efopen(stemp, mode));
}


