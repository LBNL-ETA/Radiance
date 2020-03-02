#ifndef lint
static const char	RCSid[] = "$Id: unix_process.c,v 3.17 2020/03/02 19:06:48 greg Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes
 * Unix version
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <sys/wait.h>
#include <stdlib.h>

#include "rtprocess.h"
#include "rtio.h"


SUBPROC		sp_inactive = SP_INACTIVE;


int
open_process(		/* open communication to separate process */
	SUBPROC *pd,
	char	*av[]
)
{
	char	*compath;
	int	p0[2], p1[2];

	if (pd->pid > 0)
		return(-1);		/* need to close, first */
	
	if ((pd->flags&(PF_FILT_INP|PF_FILT_OUT)) == (PF_FILT_INP|PF_FILT_OUT))
		return(-1);		/* circular process not supported */

	pd->flags &= ~PF_RUNNING;	/* not running for sure */
	pd->pid = -1;

	if (av == NULL)			/* cloning operation? */
		compath = NULL;
	else if ((compath = getpath(av[0], getenv("PATH"), X_OK)) == NULL)
		return(0);

	if (pd->flags & PF_FILT_INP) {	/* filterning input stream? */
		if ((pd->r < 0) | (pd->r == 1))
			return(-1);
		p0[0] = pd->r;
		p0[1] = -1;
	} else if (pipe(p0) < 0)
		return(-1);

	if (pd->flags & PF_FILT_OUT) {	/* filtering output stream? */
		if (pd->w < 1)
			return(-1);
		p1[0] = -1;
		p1[1] = pd->w;
	} else if (pipe(p1) < 0)
		return(-1);
#ifdef BSD
	if (compath != NULL)
		pd->pid = vfork();	/* more efficient with exec() */
	else
#endif
	pd->pid = fork();
	if (pd->pid == 0) {		/* if child... */
		close(p0[1]);
		close(p1[0]);
		if (p0[0] != 0) {	/* connect p0 to stdin */
			if (dup2(p0[0], 0) < 0)
				return(-1);
			close(p0[0]);
		}
		if (p1[1] != 1) {	/* connect p1 to stdout */
			if (dup2(p1[1], 1) < 0)
				return(-1);
			close(p1[1]);
		}
		if (compath == NULL)	/* just cloning? */
			return(0);
					/* else clear streams' FD_CLOEXEC */
		if (p0[0] == 0)
			fcntl(0, F_SETFD, 0);
		if (p1[1] == 1)
			fcntl(1, F_SETFD, 0);
		execv(compath, av);	/* exec command */
		perror(compath);
		_exit(127);
	}
	if (pd->pid == -1)
		return(-1);
					/* connect parent's streams */
	if (!(pd->flags & PF_FILT_INP)) {
		close(p0[0]);
		pd->r = p1[0];
	} else if (p1[0] != pd->r) {
		if (dup2(p1[0], pd->r) < 0)
			return(-1);
		close(p1[0]);
	}
	if (!(pd->flags & PF_FILT_OUT)) {
		close(p1[1]);
		pd->w = p0[1];
	} else if (p0[1] != pd->w) {
		if (dup2(p0[1], pd->w) < 0)
			return(-1);
		close(p0[1]);
	}
	/*
	 * Close write stream on exec to avoid multiprocessing deadlock.
	 * No use in read stream without it, so set flag there as well.
	 * GW: This bug took me two days to figure out!!
	 */
	if (pd->r > 0)
		fcntl(pd->r, F_SETFD, FD_CLOEXEC);
	if (pd->w > 1)
		fcntl(pd->w, F_SETFD, FD_CLOEXEC);
	pd->flags |= PF_RUNNING;
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
		if (pd[i].flags & PF_RUNNING) {
			close(pd[i].w);
			close(pd[i].r);
			pd[i].flags &= ~PF_RUNNING;
		} else 
			togo -= (pd[i].pid < 0);
	if (nproc == 1) {			/* await specific process? */
		status = 0;
		if (waitpid(pd->pid, &status, 0) != pd->pid)
			return(-1);
		*pd = sp_inactive;
		return(status>>8 & 0xff);
	}
						/* else unordered wait */
	while (togo > 0 && (pid = wait(&status)) >= 0) {
		for (i = nproc; i-- > 0; )
			if (pd[i].pid == pid) {
				pd[i] = sp_inactive;
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
