#ifndef lint
static const char	RCSid[] = "$Id: unix_process.c,v 3.2 2003/07/17 09:21:29 schorsch Exp $";
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

#include "rtprocess.h"
#include  "vfork.h"


int
open_process(		/* open communication to separate process */
SUBPROC *pd,
char	*av[]
)
{
	extern char	*getpath(), *getenv();
	char	*compath;
	int	p0[2], p1[2];

	pd->running = 0; /* not yet */
					/* find executable */
	compath = getpath(av[0], getenv("PATH"), 1);
	if (compath == 0)
		return(0);
	if (pipe(p0) < 0 || pipe(p1) < 0)
		return(-1);
	if ((pd->pid = vfork()) == 0) {		/* if child */
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
	pd->running = 1;
	return(PIPE_BUF);
}



int
close_process(		/* close pipes and wait for process */
SUBPROC *pd
)
{
	int	pid, status;

	close(pd->r);
	close(pd->w);
	pd->running = 0;
	while ((pid = wait(&status)) != -1)
		if (pid == pd->pid)
			return(status>>8 & 0xff);
	return(-1);		/* ? unknown status */
}


