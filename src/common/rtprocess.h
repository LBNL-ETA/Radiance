/* RCSid $Id: rtprocess.h,v 3.5 2003/10/21 19:19:28 schorsch Exp $ */
/*
 *   rtprocess.h 
 *   Routines to communicate with separate process via dual pipes
 *
 *   WARNING: On Windows, there's a system header named <process.h>.
 */
#ifndef _RAD_PROCESS_H_
#define _RAD_PROCESS_H_

#include  <errno.h>
#ifdef _WIN32
  #include <windows.h> /* DWORD etc. */
  typedef DWORD pid_t;
  #include <process.h> /* getpid() and others */
  #define nice(inc) win_nice(inc)
#else
  #include <sys/param.h>
#endif

#include "paths.h"

#ifdef __cplusplus
extern "C" {
#endif

/* On Windows, a process ID is a DWORD. That might actually be the
   same thing as an int, but it's better not to assume anything.

   This means that we shouldn't rely on PIDs and file descriptors
   being the same type, so we have to describe processes with a struct,
   instead of the original int[3]. To keep things simple, we typedef
   the posix pid_t on those systems that don't have it already.
*/


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

#ifdef _WIN32
/* any non-negative increment will send the process to IDLE_PRIORITY_CLASS. */
extern int win_nice(int inc);
#endif


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PROCESS_H_ */

