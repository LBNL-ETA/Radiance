#ifndef lint
static const char	RCSid[] = "$Id: syscalls.c,v 1.2 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
 *  System calls for meta-file routines
 */

#ifdef _WIN32
  #include <process.h> /* getpid() */
#endif

#include  "meta.h"





FILE *
efopen(fname, mode)		/* open a file, report errors */

char  *fname, *mode;

{
 register FILE  *fp;
 FILE  *fopen();

 if ((fp = fopen(fname, mode)) == NULL)  {
    sprintf(errmsg, "cannot open file \"%s\", mode \"%s\"", fname, mode);
    error(USER, errmsg);
    }

 return(fp);
 }



FILE *
mfopen(fname, mode)		/* open a program metafile */

char  *fname;
char  *mode;

{
    char  *mdir, stemp[MAXFNAME];
#if  UNIX || MAC || _WIN32
    char  *getenv();

    if ((mdir = getenv("MDIR")) == NULL)
#endif
	mdir = MDIR;
    sprintf(stemp, "%s%s", mdir, fname);

    return(efopen(stemp, mode));
}



#ifdef  CPM
getpid()			/* for CP/M, get user number */

{

 return(getusr());
 }
#endif



#ifdef  MAC
getpid()			/* dummy value for MacIntosh */

{

 return(0);
 }
#endif
