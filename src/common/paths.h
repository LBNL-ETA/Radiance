/* RCSid $Id: paths.h,v 2.30 2017/07/10 20:15:18 greg Exp $ */
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

#if defined(_WIN32) || defined(_WIN64)
  #include <io.h>
  #include <direct.h> /* getcwd(), chdir(), _mkdir(), etc. */
  #define getcwd _getcwd
  #define chdir _chdir

  #define access		_access
  #define mkdir(dirname,perms)	_mkdir(dirname)
  /* The windows _popen with the native shell breaks '\\\n' escapes.
   * RT_WINPROC (used by SCons) enables our replacement functions to fix that.
   * XXX This should really not depend on the compiler used! */
  #if defined(_MSC_VER) && !defined(RT_WINPROC)
    #define popen(cmd,mode)	_popen(cmd,mode)
    #define pclose(p)		_pclose(p)
  #else
    #define popen(cmd,mode)	win_popen(cmd,mode)
    #define pclose(p)		win_pclose(p)
  #endif
  #define kill(pid,sig)		win_kill(pid,sig)
  #define nice(inc)		win_nice(inc)
  #define PATH_MAX		_MAX_PATH
  #define NULL_DEVICE		"NUL"
  #define DIRSEP		'/'
  #define ISDIRSEP(c)		(((c)=='/') | ((c)=='\\'))
  #define ISABS(s)		( ISDIRSEP((s)[0]) \
		|| ( isupper((s)[0]) | islower((s)[0]) \
			&& (s)[1]==':' && ISDIRSEP((s)[2]) ) )
  #define CASEDIRSEP		case '/': case '\\'
  #define PATHSEP		';'
  #define CURDIR		'.'
  #define DEFAULT_TEMPDIRS	{"C:/TMP", "C:/TEMP", ".", NULL}
  #define TEMPLATE		"rtXXXXXX"
  #define TEMPLEN		8
  #define ULIBVAR		"RAYPATH"
  #ifndef DEFPATH
    #define DEFPATH		";c:/ray/lib"
  #endif
  #define SPECIALS		" \t\"$*?|"
  #define QUOTCHAR		'"'
  /* <io.h> only does half the work for access() */
  #ifndef F_OK
    #define  F_OK 00
    #define  W_OK 02
  #endif
  #ifndef R_OK
    #define X_OK 01
    #define R_OK 04
  #endif
  /* to make the permissions user specific we'd need to use CreateFile() */
  #ifndef S_IRUSR
    #define S_IRUSR _S_IREAD
    #define S_IWUSR _S_IWRITE
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
    #define ISABS(s)		ISDIRSEP((s)[0])
    #define PATHSEP		';'
	#define CURDIR		'.'
    #define DEFAULT_TEMPDIRS	{"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE		"/tmp/rtXXXXXX"
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
    #define ISABS(s)		ISDIRSEP((s)[0])
    #define PATHSEP		':'
	#define CURDIR		'.'
    #define DEFAULT_TEMPDIRS	{"/var/tmp", "/usr/tmp", "/tmp", ".", NULL}
    #define TEMPLATE		"/tmp/rtXXXXXX"
    #define TEMPLEN		13
    #define ULIBVAR		"RAYPATH"
    #ifndef DEFPATH
      #define DEFPATH		":/usr/local/lib/ray"
    #endif
    #define SPECIALS		" \t\n'\"()${}*?[];|&"
    #define QUOTCHAR		'\''
    #define ALTQUOT		'"'
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

#if defined(_WIN32) || defined(_WIN64)
  extern FILE *win_popen(char *command, char *type);
  extern int win_pclose(FILE *p);
  extern char  *fixargv0(char *arg0);
#endif

/* Check if any of the characters in str2 are found in str1 */
extern int matchany(const char *str1, const char *str2);

/* Convert a set of arguments into a command line for pipe() or system() */
extern char *convert_commandline(char *cmd, const int len, char *av[]);

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

