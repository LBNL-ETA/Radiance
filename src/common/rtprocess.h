/* RCSid $Id: rtprocess.h,v 3.18 2020/02/28 05:18:49 greg Exp $ */
/*
 *   rtprocess.h 
 *   Routines to communicate with separate process via dual pipes
 *
 *   WARNING: On Windows, there's a system header named <process.h>.
 */
#ifndef _RAD_PROCESS_H_
#define _RAD_PROCESS_H_

#include  <errno.h>
#include <stdio.h>
#if defined(_WIN32) || defined(_WIN64)
  #include <windows.h> /* DWORD etc. */
  typedef DWORD RT_PID;
  #include <process.h> /* getpid() and others */
  #define getpid _getpid
  #define execv _execv
  #define execvp _execvp
#else
  #include <sys/param.h>
  #include <sys/types.h>
  typedef pid_t RT_PID;
#endif

#include "paths.h"

#ifdef __cplusplus
extern "C" {
#endif

/* On Windows, a process ID is a DWORD. That might actually be the
   same thing as an int, but it's better not to assume anything.

   This means that we shouldn't rely on PIDs and file descriptors
   being the same type, so we have to describe processes with a struct,
   instead of the original int[3]. For that purpose, we typedef a
   platform independent RT_PID.
*/

/* On Unix, we can set flags and assign descriptors before opening a
   process, coupling an existing input or output to the new process rather
   than opening both pipes.  If PF_FILT_INP is passed in the flags member of
   SUBPROC, then the given r stream will be attached to the standard input
   of the child process, and subsequent reads from that descriptor in the
   parent get data from the standard output of the child, instead.  The
   returned w descriptor is set to -1, since there is no longer any way
   to write to the input of the child.  The default r descriptor of 0 will
   compel the child to act as a filter on the standard input of the parent.
   Whatever r handle you specify, the child will filter its read operations.
   Note that this should be called before anything has been buffered using r.
   If PF_FILT_OUT is set in flags, then the given w stream will be
   attached to the standard output of the child, and subsequent writes
   to that descriptor in the parent send data to the standard input
   of the child. The returned r descriptor is set to -1, since
   there is no output to read from any longer in the child.  The
   default w descriptor of 1 will cause the child to act as a filter
   on the output of the parent.  Make sure to call fflush(stdout) first
   if any data was buffered.  It is illegal to set both PF_FILT_INP and
   PF_FILT_OUT, as a circular process is guaranteed to hang.
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
				/* process flags */
#define PF_RUNNING	1		/* process is running */
#define PF_FILT_INP	2		/* use assigned read descriptor */
#define PF_FILT_OUT	4		/* use assigned write descriptor */

typedef struct {
	int	flags;		/* what is being done */
	int	r;		/* read handle */
	int	w;		/* write handle */
	RT_PID	pid;		/* process ID */
} SUBPROC;

#define SP_INACTIVE {0,0,1,-1}	/* for static initializations */

#define close_process(pd)	close_processes(pd,1)

extern int open_process(SUBPROC *pd, char *av[]);
extern int close_processes(SUBPROC pd[], int nproc);
extern int process(SUBPROC *pd, char *recvbuf, char *sendbuf, int nbr, int nbs);
extern int readbuf(int fd, char *bpos, int siz);
extern int writebuf(int fd, char *bpos, int siz);

#if defined(_WIN32) || defined(_WIN64)
/* any non-negative increment will send the process to IDLE_PRIORITY_CLASS. */
extern int win_kill(RT_PID pid, int sig /* ignored */);
extern int win_nice(int inc);
#endif

extern SUBPROC	sp_inactive;

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PROCESS_H_ */

