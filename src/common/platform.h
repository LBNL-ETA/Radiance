/* RCSid $Id: platform.h,v 3.1 2003/06/05 19:29:34 schorsch Exp $ */
/*
 *  platform.h - header file for platform specific definitions
 */

#include "copyright.h"


#ifndef _RAD_PLATFORM_H_
#define _RAD_PLATFORM_H_


#ifdef _WIN32

#include <stdio.h>  /* fileno() */
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

#endif /* _RAD_PLATFORM_H_ */

