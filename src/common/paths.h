/* RCSid $Id: paths.h,v 2.22 2005/09/19 11:30:10 schorsch Exp $ */
/*
 * Definitions for paths on different machines
 */
#ifndef _RAD_PATHS_H_
#define _RAD_PATHS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
  #include <io.h>
  #include <direct.h> /* getcwd(), chdir(), _mkdir(), etc. */

  #define access _access
  #define mkdir(dirname,perms) _mkdir(dirname)
  #define PATH_MAX _MAX_PATH
  #define NULL_DEVICE	"NUL"
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
  /* to make the permissions user specific we'd need to use CreateFile() */
  #ifndef S_IRUSR
    #define S_IRUSR _S_IREAD
    #define S_IWUSR _S_IWRITE
  #endif

  #ifdef __cplusplus
    extern "C" {
  #endif
  extern char  *fixargv0();
  #ifdef __cplusplus
    }
  #endif

#else /* everything but Windows */
  #include <unistd.h>
  #include <sys/param.h>

  #define RMAX_PATH_MAX 4096 /* our own maximum */
  #ifndef PATH_MAX
    #define PATH_MAX 512
  #elif PATH_MAX > RMAX_PATH_MAX /* the OS is exaggerating */
    #undef PATH_MAX
    #define PATH_MAX RMAX_PATH_MAX
  #endif

  #ifdef AMIGA

	#define NULL_DEVICE	"NIL:"
    #define DIRSEP		'/'
    #define ISABS(s) ((s)!=NULL && (ISDIRSEP(s[0])))
    #define PATHSEP		';'
	#define CURDIR		'.'
    #define DEFAULT_TEMPDIRS {"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE	"/tmp/rtXXXXXX"
    #define TEMPLEN		13
    #define ULIBVAR		"RAYPATH"
    #ifndef DEFPATH
      #define DEFPATH		";/ray/lib"
    #endif
    #define	 fixargv0(a0)	(a0)

  #else

    /* posix */

	/* this is defined as _PATH_DEVNULL in /usr/include/paths.h on Linux */
	#define NULL_DEVICE	"/dev/null"
    #define DIRSEP		'/'
    #define ISABS(s) ((s)!=NULL && (ISDIRSEP(s[0])))
    #define PATHSEP		':'
	#define CURDIR		'.'
    #define DEFAULT_TEMPDIRS {"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE	"/tmp/rtXXXXXX"
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

#ifdef __cplusplus
extern "C" {
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
extern int temp_fd(char *s, size_t len, char *templ);

/* Same as above, but return a file pointer instead of a descriptor */
extern FILE *temp_fp(char *s, size_t len, char *templ);

/* Concatenate two strings, leaving exactly one DIRSEP in between */
extern char *append_filepath(char *s1, char *s2, size_t len);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PATHS_H_ */

