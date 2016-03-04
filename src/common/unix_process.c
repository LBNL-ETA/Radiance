#ifndef lint
static const char	RCSid[] = "$Id: unix_process.c,v 3.12 2016/03/04 19:16:49 greg Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes
 * Unix version
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <stdlib.h>

#include "rtprocess.h"
#include "rtio.h"


int
open_process(		/* open communication to separate process */
SUBPROC *pd,
char	*av[]
)
{
	char	*compath;
	int	p0[2], p1[2];

	pd->pid = -1;
	pd->running = 0;		/* not going yet */
	
	if (av == NULL)			/* cloning operation? */
		compath = NULL;
	else if ((compath = getpath(av[0], getenv("PATH"), X_OK)) == NULL)
		return(0);
	if (pipe(p0) < 0 || pipe(p1) < 0)
		return(-1);
	if ((pd->pid = fork()) == 0) {	/* if child... */
		close(p0[1]);
		close(p1[0]);
		if (p0[0] != 0) {	/* connect p0 to stdin */
			dup2(p0[0], 0);
			close(p0[0]);
		}
		if (p1[1] != 1) {	/* connect p1 to stdout */
			dup2(p1[1], 1);
			close(p1[1]);
		}
		if (compath == NULL)	/* just cloning? */
			return(0);
		execv(compath, av);	/* else exec command */
		perror(compath);
		_exit(127);
	}
	if (pd->pid == -1)
		return(-1);
	close(p0[0]);
	close(p1[1]);
	pd->r = p1[0];
	pd->w = p0[1];
	/*
	 * Close write stream on exec to avoid multiprocessing deadlock.
	 * No use in read stream without it, so set flag there as well.
	 * GW: This bug took me two days to figure out!!
	 */
	fcntl(pd->r, F_SETFD, FD_CLOEXEC);
	fcntl(pd->w, F_SETFD, FD_CLOEXEC);
	pd->running = 1;
	return(PIPE_BUF);
}


int
close_processes(	/* close pipes and wait for processes to finish */
SUBPROC pd[],
int nproc
)
{
	int	togo = nproc;
	int	status, rtn_status = 0;
	RT_PID	pid;
	int	i;

	for (i = 0; i < nproc; i++)		/* close pipes, first */
		if (pd[i].running) {
			close(pd[i].w);
			close(pd[i].r);
			pd[i].running = 0;
		}
	if (nproc == 1) {			/* await specific process? */
		if (waitpid(pd->pid, &status, 0) != pd->pid)
			return(-1);
		pd->pid = -1;
		return(status>>8 & 0xff);
	}
						/* else unordered wait */
	while (togo > 0 && (pid = wait(&status)) >= 0) {
		for (i = nproc; i-- > 0; )
			if (pd[i].pid == pid) {
				pd[i].pid = -1;
				--togo;
				break;
			}
		if (i < 0)
			continue;		/* child we don't know? */
		status = status>>8 & 0xff;
		if (status)			/* record non-zero status */
			rtn_status = status;
	}
	if (togo)				/* child went missing? */
		return(-1);
	return(rtn_status);
}
