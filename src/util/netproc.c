#ifndef lint
static const char	RCSid[] = "$Id: netproc.c,v 2.13 2003/11/10 16:52:24 greg Exp $";
#endif
/*
 * Parallel network process handling routines
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <fcntl.h>
#include <unistd.h>

#include "selcall.h"
#include "netproc.h"
#include "paths.h"

PSERVER	*pslist = NULL;		/* global process server list */

static NETPROC	*pindex[FD_SETSIZE];	/* process index table */

static char	ourhost[64];	/* this host name */
static char	ourdir[PATH_MAX];	/* our working directory */
static char	ouruser[32];	/* our user name */
static char	*ourshell;	/* our user's shell */

static fd_set	errdesc;	/* error file descriptors */
static int	maxfd;		/* maximum assigned descriptor */

extern char	*remsh;		/* externally defined remote shell program */


PSERVER *
addpserver(host, dir, usr, np)		/* add a new process server */
char	*host, *dir, *usr;
int	np;
{
	register PSERVER	*ps;
					/* allocate the struct */
	if (np < 1)
		return(NULL);
	ps = (PSERVER *)malloc(sizeof(PSERVER)+(np-1)*sizeof(NETPROC));
	if (ps == NULL)
		return(NULL);
	if (!ourhost[0]) {		/* initialize */
		char	dirtmp[PATH_MAX];
		register char	*cp;
		register int	len;

		strcpy(ourhost, myhostname());
		getcwd(dirtmp, sizeof(dirtmp));
		if ((cp = getenv("HOME")) != NULL) {
			if (!strcmp(cp, dirtmp))
				ourdir[0] = '\0';
			else if (!strncmp(cp, dirtmp, len=strlen(cp)) &&
					dirtmp[len] == '/')
				strcpy(ourdir, dirtmp+len+1);
			else
				strcpy(ourdir, dirtmp);
		} else
			strcpy(ourdir, dirtmp);
		if ((cp = getenv("USER")) != NULL)
			strcpy(ouruser, cp);
		if ((ourshell = getenv("SHELL")) == NULL)
			ourshell = "/bin/sh";
		FD_ZERO(&errdesc);
		maxfd = -1;
	}
					/* assign host, directory, user */
	if (host == NULL || !strcmp(host, ourhost) ||
			!strcmp(host, LHOSTNAME))
		ps->hostname[0] = '\0';
	else
		strcpy(ps->hostname, host);
	if (dir == NULL)
		strcpy(ps->directory, ourdir);
	else
		strcpy(ps->directory, dir);
	if (usr == NULL || !strcmp(usr, ouruser))
		ps->username[0] = '\0';
	else
		strcpy(ps->username, usr);
					/* clear process slots */
	ps->nprocs = np;
	while (np--) {
		ps->proc[np].com = NULL;
		ps->proc[np].pid = -1;
		ps->proc[np].efd = -1;
		ps->proc[np].errs = NULL;
		ps->proc[np].elen = 0;
		ps->proc[np].cf = NULL;
	}
					/* insert in our list */
	ps->next = pslist;
	pslist = ps;
					/* check for signs of life */
	if (!pserverOK(ps)) {
		delpserver(ps);			/* failure -- abort */
		return(NULL);
	}
	return(ps);
}


delpserver(ps)				/* delete a process server */
PSERVER	*ps;
{
	PSERVER	pstart;
	register PSERVER	*psp;
	register int	i;
					/* find server in our list */
	pstart.next = pslist;
	for (psp = &pstart; ps != psp->next; psp = psp->next)
		if (psp->next == NULL)
			return;			/* not in our list! */
					/* kill any running jobs */
	for (i = 0; i < ps->nprocs; i++)
		if (ps->proc[i].com != NULL) {
			kill(SIGTERM, ps->proc[i].pid);
			wait4job(ps, ps->proc[i].pid);
		}
					/* remove server from list */
	psp->next = ps->next;
	pslist = pstart.next;
	free((void *)ps);		/* free associated memory */
}


PSERVER *
findjob(pnp)			/* find out where process is running */
register int	*pnp;	/* modified */
{
	register PSERVER	*ps;
	register int	i;

	for (ps = pslist; ps != NULL; ps = ps->next)
		for (i = 0; i < ps->nprocs; i++)
			if (ps->proc[i].pid == *pnp) {
				*pnp = i;
				return(ps);
			}
	return(NULL);		/* not found */
}


int
startjob(ps, command, compf)	/* start a job on a process server */
register PSERVER	*ps;
char	*command;
int	(*compf)();
{
	char	udirt[PATH_MAX];
	char	*av[16];
	int	pfd[2], pid;
	register int	i;

	if (ps == NULL) {		/* find a server */
		for (ps = pslist; ps != NULL; ps = ps->next)
			if ((i = startjob(ps, command, compf)) != -1)
				return(i);		/* got one */
		return(-1);	/* no slots anywhere */
	}
	for (i = 0; i < ps->nprocs; i++)
		if (ps->proc[i].com == NULL)
			break;
	if (i >= ps->nprocs)
		return(-1);		/* out of process slots */
					/* open pipe */
	if (pipe(pfd) < 0) {
		perror("cannot open pipe");
		exit(1);
	}
					/* start child process */
	if ((pid = vfork()) == 0) {
		close(pfd[0]);			/* connect stderr to pipe */
		if (pfd[1] != 2) {
			dup2(pfd[1], 2);
			close(pfd[1]);
		}
		if (ps->hostname[0]) {		/* rsh command */
			av[i=0] = remsh;
			av[++i] = ps->hostname;
			av[++i] = "-n";			/* no stdin */
			if (ps->username[0]) {		/* different user */
				av[++i] = "-l";
				av[++i] = ps->username;
				av[++i] = "cd";
				udirt[0] = '~';
				strcpy(udirt+1, ouruser);
				av[++i] = udirt;
				av[++i] = ";";
			}
			if (ps->directory[0]) {		/* change directory */
				av[++i] = "cd";
				av[++i] = ps->directory;
				av[++i] = ";";
			}
			av[++i] = command;
			av[++i] = NULL;
		} else {			/* shell command */
			av[0] = ourshell;
			av[1] = "-c";
			av[2] = command;
			av[3] = NULL;
		}
		execv(av[0], av);
		_exit(1);
	}
	if (pid == -1) {
		perror("fork failed");
		exit(1);
	}
	ps->proc[i].com = command;	/* assign process slot */
	ps->proc[i].cf = compf;
	ps->proc[i].pid = pid;
	close(pfd[1]);			/* get piped stderr file descriptor */
	ps->proc[i].efd = pfd[0];
	fcntl(pfd[0], F_SETFD, 1);	/* set close on exec flag */
	pindex[pfd[0]] = ps->proc + i;	/* assign error fd index */
	FD_SET(pfd[0], &errdesc);	/* add to select call parameter */
	if (pfd[0] > maxfd)
		maxfd = pfd[0];
	return(pid);			/* return to parent process */
}


static int
readerrs(fd)			/* read error output from fd */
int	fd;
{
	char	errbuf[BUFSIZ];
	int	nr;
	register NETPROC	*pp;
				/* look up associated process */
	if ((pp = pindex[fd]) == NULL)
		abort();		/* serious consistency error */
	nr = read(fd, errbuf, BUFSIZ-1);
	if (nr < 0) {
		perror("read error");
		exit(1);
	}
	if (nr == 0)		/* stream closed (process finished) */
		return(0);
	errbuf[nr] = '\0';	/* add to error buffer */
	if (pp->elen == 0)
		pp->errs = (char *)malloc(nr+1);
	else
		pp->errs = (char *)realloc((void *)pp->errs, pp->elen+nr+1);
	if (pp->errs == NULL) {
		perror("malloc failed");
		exit(1);
	}
	strcpy(pp->errs+pp->elen, errbuf);
	pp->elen += nr;
	return(nr);
}


static
wait4end()			/* read error streams until someone is done */
{
	fd_set	readfds, excepfds;
	register int	i;
				/* find end of descriptor set */
	for ( ; maxfd >= 0; maxfd--)
		if (FD_ISSET(maxfd, &errdesc))
			break;
	if (maxfd < 0)
		return;		/* nothing to read */
	readfds = excepfds = errdesc;
	while (select(maxfd+1, &readfds, NULL, &excepfds, NULL) > 0)
		for (i = 0; i <= maxfd; i++)		/* get pending i/o */
			if (FD_ISSET(i, &readfds) || FD_ISSET(i, &excepfds))
				if (readerrs(i) == 0)
					return;		/* finished process */
	perror("select call failed");
	exit(1);
}


static int
finishjob(ps, pn, status)	/* clean up finished process */
PSERVER	*ps;
int	pn;
int	status;
{
	register NETPROC	*pp;

	pp = ps->proc + pn;
	if (pp->cf != NULL)			/* client cleanup */
		status = (*pp->cf)(ps, pn, status);
	close(pp->efd);				/* close error stream */
	pindex[pp->efd] = NULL;
	FD_CLR(pp->efd, &errdesc);
	free((void *)pp->errs);
	pp->com = NULL;				/* clear settings */
	pp->pid = -1;
	pp->efd = -1;
	pp->errs = NULL;
	pp->elen = 0;
	pp->cf = NULL;
	return(status);
}


int
wait4job(ps, pid)		/* wait for process to finish */
PSERVER	*ps;
int	pid;
{
	int	status, psn, psn2;
	PSERVER	*ps2;

	if (pid == -1) {			/* wait for first job */
		if (ps != NULL) {
			for (psn = ps->nprocs; psn--; )
				if (ps->proc[psn].com != NULL)
					break;
			if (psn < 0)
				return(-1);	/* no processes this server */
		}
		do {
			wait4end();		/* wait for something to end */
			if ((psn2 = wait(&status)) == -1)
				return(-1);	/* none left */
			ps2 = findjob(&psn2);
			if (ps2 != NULL)	/* clean up job if ours */
				status = finishjob(ps2, psn2, status);
		} while (ps2 == NULL || (ps != NULL && ps2 != ps));
		return(status);			/* return job status */
	}
	psn = pid;				/* else find specific job */
	ps2 = findjob(&psn);			/* find process slot */
	if (ps2 == NULL || (ps != NULL && ps2 != ps))
		return(-1);		/* inconsistent target */
	ps = ps2;
	do {
		wait4end();			/* wait for something to end */
		if ((psn2 = wait(&status)) == -1)
			return(-1);		/* none left */
		ps2 = findjob(&psn2);
		if (ps2 != NULL)		/* clean up job if ours */
			status = finishjob(ps2, psn2, status);
	} while (ps2 != ps || psn2 != psn);
	return(status);				/* return job status */
}
