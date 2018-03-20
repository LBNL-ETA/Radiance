/* RCSid $Id: platform.h,v 3.16 2018/03/20 22:45:29 greg Exp $ */
/*
 *  platform.h - header file for platform specific definitions
 */
#ifndef _RAD_PLATFORM_H_
#define _RAD_PLATFORM_H_

#if defined(_WIN32) || defined(_WIN64)

  #include <io.h>     /* _setmode() and stuff from unistd.h */
  #include <stdio.h>
  typedef long off_t;
  #undef fdopen
  #define fdopen _fdopen
  #undef read
  #define read _read
  #undef open
  #define open _open
  #undef close
  #define close _close
  #undef write
  #define write _write
  #undef ftruncate
  #define ftruncate _chsize_s
  #undef unlink
  #define unlink _unlink
  #undef fileno
  #define fileno _fileno
  #undef snprintf
  #define snprintf _snprintf
  #undef vsnprintf
  #define vsnprintf _vsnprintf
  /* XXX should we check first if size_t is 32 bit? */
  #undef fseeko
  #define fseeko _fseeki64
  #undef lseek
  #define lseek _lseek
  #undef access
  #define access _access
  #undef mktemp
  #define mktemp _mktemp

  #include <string.h>
  #undef strcasecmp
  #define strcasecmp _stricmp
  #undef strncasecmp
  #define strncasecmp _strnicmp
  #undef strdup
  #define strdup _strdup

  #include <windows.h>
  /* really weird defines by Microsoft in <resource.h>
	 generating lots of name collisions in Radiance. */
  #if defined(rad1)
    #undef rad1
    #undef rad2
    #undef rad3
    #undef rad4
    #undef rad5
    #undef rad6
    #undef rad7
    #undef rad8
    #undef rad9
  #endif
  #define sleep(s) Sleep((DWORD)((s)*1000))

  #define NON_POSIX

  #include <sys/types.h>
  #include <sys/stat.h>
  #define RHAS_STAT
  #define S_IFREG _S_IFREG
  #define W_IFDIR _S_IFDIR

  #include <fcntl.h>  /* _O_BINARY, _O_TEXT */
  #include <stdlib.h> /* _fmode */
  #define SET_DEFAULT_BINARY() (_fmode = _O_BINARY)
  #define SET_DEFAULT_TEXT() (_fmode = _O_TEXT)
  #define SET_FILE_BINARY(fp) _setmode(_fileno(fp),_O_BINARY)
  #define SET_FILE_TEXT(fp) _setmode(_fileno(fp),_O_TEXT)
  #define SET_FD_BINARY(fd) _setmode(fd,_O_BINARY)
  #define SET_FD_TEXT(fd) _setmode(fd,_O_TEXT)
  #define putenv _putenv

#else /* _WIN32 || _WIN64 */

  #ifdef AMIGA
    #define NON_POSIX
  #else
    /* assumedly posix systems */
	#include <unistd.h>
    #define RHAS_STAT
    #define RHAS_FORK_EXEC
  #endif

  /* everybody except Windows */

  /* NOPs */
  #define SET_DEFAULT_BINARY()
  #define SET_FILE_BINARY(fp)
  #define SET_FD_BINARY(fd)
  #define SET_DEFAULT_TEXT()
  #define SET_FILE_TEXT(fp)
  #define SET_FD_TEXT(fd)

#endif /* _WIN32 || _WIN64 */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(_WIN64)

extern	int	usleep(__int64 usec);

#endif /* _WIN32 || _WIN64 */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PLATFORM_H_ */

