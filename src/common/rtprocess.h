/* RCSid $Id: rtprocess.h,v 3.1 2003/06/26 00:58:09 schorsch Exp $ */
/*
 *   rtprocess.h 
 *   Routines to communicate with separate process via dual pipes
 *
 *   WARNING: On Windows, there's a system header named <process.h>.
 */
#ifndef _RAD_PROCESS_H_
#define _RAD_PROCESS_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "copyright.h"

#include  <sys/types.h>
#ifdef _WIN32
  #include <windows.h>
#else
  #include <sys/param.h>
  #include <unistd.h>
#endif
#ifndef BSD
#include  <errno.h>
#endif

#include "paths.h"


/* On Windows, a process ID is a DWORD. That might actually be the
   same thing as an int, but it's better not to assume anything.

   This means that we shouldn't rely on PIDs and file descriptors
   being the same type, so we have to describe processes with a struct,
   instead of the original int[3]. To keep things simple, we typedef
   the posix pid_t on those systems that don't have it already.

   Some older Windows systems use negative PIDs. Open_process() and
   close_process() will convert those to positive values during
   runtime, so that client modules can still use -1 as invalid PID.
*/

#ifdef _WIN32
  typedef DWORD pid_t;
#endif

#ifndef PIPE_BUF
  #ifdef PIPSIZ
    #define PIPE_BUF	PIPSIZ
  #else
    #ifdef PIPE_MAX
      #define PIPE_BUF	PIPE_MAX
    #else
      #define PIPE_BUF	512		/* hyperconservative */
    #endif
  #endif
#endif

typedef struct {
	int r; /* read handle */
	int w; /* write handle */
	int running; /* doing something */
	pid_t pid; /* process ID */
} SUBPROC;

#define SP_INACTIVE {-1,-1,0,0} /* for static initializations */

extern int open_process(SUBPROC *pd, char *av[]);
extern int close_process(SUBPROC *pd);
extern int process(SUBPROC *pd, char *recvbuf, char *sendbuf, int nbr, int nbs);
extern int readbuf(int fd, char *bpos, int siz);
extern int writebuf(int fd, char *bpos, int siz);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PROCESS_H_ */

