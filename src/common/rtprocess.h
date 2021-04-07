/* RCSid $Id: rtprocess.h,v 3.23 2021/04/07 01:15:53 greg Exp $ */
/*
 *   rtprocess.h 
 *   Routines to communicate with separate process via dual pipes
 *
 *   WARNING: On Windows, there's a system header named <process.h>.
 */
#ifndef _RAD_PROCESS_H_
#define _RAD_PROCESS_H_

#include <errno.h>
#include "paths.h"
#if defined(_WIN32) || defined(_WIN64)
  #include <windows.h> /* DWORD etc. */
  typedef DWORD RT_PID;
  #include <process.h> /* getpid() and others */
  #define getpid _getpid
  #define execv _execv
  #define execvp _execvp
  #ifdef _MSC_VER
    #include <BaseTsd.h>
    typedef SSIZE_T ssize_t;
  #endif
#else
  typedef pid_t RT_PID;
#endif


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
   there is no longer any output to read from the child.  The
   default w descriptor of 1 will cause the child to act as a filter
   on the output of the parent.  Make sure to call fflush(stdout) first
   if any data was buffered.  It is illegal to set both PF_FILT_INP and
   PF_FILT_OUT, as a circular process is guaranteed to hang.
   
   If you want behavior similar to popen(cmd, "w") (again Unix-only),
   keeping stdout open in parent, use a duplicate descriptor like so:
   {
	SUBPROC	rtp = sp_inactive;
	FILE	*fout;
	fflush(stdout);
	rtp.w = dup(fileno(stdout));
	rtp.flags |= PF_FILT_OUT;
	if (open_process(&rtp, cmd_argv) <= 0) {
		perror(cmd_argv[0]); exit(1);
	}
	fout = fdopen(rtp.w, "w");
	...write data to filter using fout until finished...
	fclose(fout);
	if (close_process(&rtp)) {
		perror(cmd_argv[0]); exit(1);
	}
	...can continue sending data directly to stdout...
    }
    We could also have called open_process() after fdopen() above, or after
    using fopen() on a file if we wanted to insert our filter before it.
    A similar sequence may be used to filter from stdin without closing
    it, though process termination becomes more difficult with two readers.
    Filtering input from a file works better, since the file is then read by
    the child only, as in:
    {
	SUBPROC rtp = sp_inactive;
	FILE	*fin = fopen(fname, "r");
	if (fin == NULL) {
		open_error(fname); exit(1);
	}
	rtp.r = fileno(fin);
	rtp.flags |= PF_FILT_INP;
	if (open_process(&rtp, cmd_argv) <= 0) {
		perror(cmd_argv[0]); fclose(fin); exit(1);
	}
	...read filtered file data from fin until EOF...
	fclose(fin);
	if (close_process(&rtp)) {
		perror(cmd_argv[0]); exit(1);
	}
    }
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
extern ssize_t readbuf(int fd, char *bpos, ssize_t siz);
extern ssize_t writebuf(int fd, char *bpos, ssize_t siz);

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

