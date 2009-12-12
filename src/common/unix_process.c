#ifndef lint
static const char	RCSid[] = "$Id: unix_process.c,v 3.9 2009/12/12 23:08:13 greg Exp $";
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

	pd->running = 0; /* not yet */
					/* find executable */
	compath = getpath(av[0], getenv("PATH"), 1);
	if (compath == 0)
		return(0);
	if (pipe(p0) < 0 || pipe(p1) < 0)
		return(-1);
	if ((pd->pid = fork()) == 0) {		/* if child */
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
		execv(compath, av);	/* exec command */
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
close_process(		/* close pipes and wait for process */
SUBPROC *pd
)
{
	int	status;

	if (!pd->running)
		return(0);
	close(pd->w);
	pd->running = 0;
	if (waitpid(pd->pid, &status, 0) == pd->pid) {
		close(pd->r);
		return(status>>8 & 0xff);
	}
	return(-1);		/* ? unknown status */
}


