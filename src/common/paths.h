/* RCSid $Id: paths.h,v 2.14 2003/06/26 00:58:09 schorsch Exp $ */
/*
 * Definitions for paths on different machines
 */
#ifndef _RAD_PATHS_H_
#define _RAD_PATHS_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
  #include <io.h>
#else /* _WIN32 */
  #include <unistd.h>
  #include <sys/param.h>
#endif /* _WIN32 */

#ifdef _WIN32

  #define access _access
  #define PATH_MAX _MAX_PATH
  #define DIRSEP		'/'
  #define ISDIRSEP(c)	((c)=='/' || (c)=='\\')
  #define ISABS(s)	((s)!=NULL \
		  && (s[0])!='\0' \
          && (   ISDIRSEP(s[0]) \
              || (   (s[1])!='\0' \
                  && (isupper(s[0])||islower(s[0])) \
                  && (s[1])==':')))
  #define CASEDIRSEP	case '/': case '\\'
  #define PATHSEP		';'
  #define CURDIR		'.'
  #define MAXPATH		512 /* XXX obsoleted by posix PATH_MAX */
  #define DEFAULT_TEMPDIRS {"C:/TMP", "C:/TEMP", ".", NULL}
  #define TEMPLATE	"rtXXXXXX"
  #define TEMPLEN		8
  #define ULIBVAR		"RAYPATH"
  #ifndef DEFPATH
    #define DEFPATH		";c:/ray/lib"
  #endif
  /* <io.h> only does half the work for access() */
  #ifndef F_OK
    #define  F_OK 00
    #define  X_OK 01
    #define  W_OK 02
    #define  R_OK 04
  #endif
extern char  *fixargv0();

#else
  #ifndef PATH_MAX /* doesn't implement posix? */
    #define PATH_MAX 512 /* quite conservative (I think...) */
  #endif

  #ifdef AMIGA

    #define DIRSEP		'/'
    #define ISABS(s) ((s)!=NULL && (ISDIRSEP(s[0])))
    #define PATHSEP		';'
	#define CURDIR		'.'
    #define MAXPATH		128
    #define DEFAULT_TEMPDIRS {"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE	"/tmp/rtXXXXXX"
    #define TEMPLEN		13
    #define ULIBVAR		"RAYPATH"
    #ifndef DEFPATH
      #define DEFPATH		";/ray/lib"
    #endif
    #define	 fixargv0(a0)	(a0)

  #else

    #define DIRSEP		'/'
    #define ISABS(s) ((s)!=NULL && (ISDIRSEP(s[0])))
    #define PATHSEP		':'
	#define CURDIR		'.'
    #define MAXPATH		256
    #define DEFAULT_TEMPDIRS {"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE	"/usr/tmp/rtXXXXXX"
    #define TEMPLEN		17
    #define ULIBVAR		"RAYPATH"
    #ifndef DEFPATH
      #define DEFPATH		":/usr/local/lib/ray"
    #endif
    #define	 fixargv0(a0)	(a0)

  #endif
#endif

#ifndef ISDIRSEP
  #define ISDIRSEP(c)	((c)==DIRSEP)
#endif
#ifndef CASEDIRSEP
  #define CASEDIRSEP	case DIRSEP
#endif

extern char  *mktemp(), *getenv();

#ifdef BSD
extern char  *getwd();
#else
extern char  *getcwd();
#define  getwd(p)	getcwd(p, sizeof(p))
#endif


/* Find a writeable directory for temporary files */
/* If s is NULL, we return a static string */
extern char *temp_directory(char *s, size_t len);

/* Compose a *currently* unique name within a temporary directory */
/* If s is NULL, we return a static string */
/* If templ is NULL, we take our default template */
/* WARNING: On Windows, there's a maximum of 27 unique names within
            one process for the same template. */
extern char *temp_filename(char *s, size_t len, char *templ);

/* Same as above, but also open the file and return the descriptor */
/* This one is supposed to protect against race conditions on unix */
/* WARNING: On Windows, there's no protection against race conditions */
extern int temp_file(char *s, size_t len, char *templ);

/* Concatenate two strings, leaving exactly one DIRSEP in between */
extern char *append_filepath(char *s1, char *s2, size_t len);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PATHS_H_ */

