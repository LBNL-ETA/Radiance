/* RCSid $Id: platform.h,v 3.4 2003/06/27 06:53:21 greg Exp $ */
/*
 *  platform.h - header file for platform specific definitions
 */
#ifndef _RAD_PLATFORM_H_
#define _RAD_PLATFORM_H_
#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32

#include <stdio.h>
#define popen _popen
#define pclose _pclose
#include <fcntl.h>  /* _O_BINARY, _O_TEXT */
#include <io.h>     /* _setmode() */
#include <stdlib.h> /* _fmode */

#define SET_DEFAULT_BINARY() _fmode = _O_BINARY
#define SET_FILE_BINARY(fp) _setmode(fileno(fp),_O_BINARY)
#define SET_FD_BINARY(fd) _setmode(fd,_O_BINARY)







#else /* _WIN32 */

/* NOPs on unix */
#define SET_DEFAULT_BINARY()
#define SET_FILE_BINARY(fp)
#define SET_FD_BINARY(fd)







#endif /* _WIN32 */


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PLATFORM_H_ */

