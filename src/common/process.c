#ifndef lint
static const char	RCSid[] = "$Id: process.c,v 2.6 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

				/* find pipe buffer limit */
#include  <sys/param.h>

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

#include  "vfork.h"

#ifndef BSD
#include  <errno.h>
#endif


int
open_process(pd, av)		/* open communication to separate process */
int	pd[3];
char	*av[];
{
	extern char	*getpath(), *getenv();
	char	*compath;
	int	p0[2], p1[2];
					/* find executable */
	compath = getpath(av[0], getenv("PATH"), 1);
	if (compath == 0)
		return(0);
	if (pipe(p0) < 0 || pipe(p1) < 0)
		return(-1);
	if ((pd[2] = vfork()) == 0) {		/* if child */
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
	if (pd[2] == -1)
		return(-1);
	close(p0[0]);
	close(p1[1]);
	pd[0] = p1[0];
	pd[1] = p0[1];
	return(PIPE_BUF);
}


int
process(pd, recvbuf, sendbuf, nbr, nbs)		/* process data through pd */
int	pd[3];
char	*recvbuf, *sendbuf;
int	nbr, nbs;
{
	if (nbs > PIPE_BUF)
		return(-1);
	if (writebuf(pd[1], sendbuf, nbs) < nbs)
		return(-1);
	return(readbuf(pd[0], recvbuf, nbr));
}


int
close_process(pd)		/* close pipes and wait for process */
int	pd[3];
{
	int	pid, status;

	close(pd[1]);
	close(pd[0]);
	while ((pid = wait(&status)) != -1)
		if (pid == pd[2])
			return(status>>8 & 0xff);
	return(-1);		/* ? unknown status */
}


int
readbuf(fd, bpos, siz)		/* read all of requested buffer */
int	fd;
char	*bpos;
int	siz;
{
	register int	cc = 0, nrem = siz;
retry:
	while (nrem > 0 && (cc = read(fd, bpos, nrem)) > 0) {
		bpos += cc;
		nrem -= cc;
	}
	if (cc < 0) {
#ifndef BSD
		if (errno == EINTR)	/* we were interrupted! */
			goto retry;
#endif
		return(cc);
	}
	return(siz-nrem);
}


int
writebuf(fd, bpos, siz)		/* write all of requested buffer */
int	fd;
char	*bpos;
int	siz;
{
	register int	cc = 0, nrem = siz;
retry:
	while (nrem > 0 && (cc = write(fd, bpos, nrem)) > 0) {
		bpos += cc;
		nrem -= cc;
	}
	if (cc < 0) {
#ifndef BSD
		if (errno == EINTR)	/* we were interrupted! */
			goto retry;
#endif
		return(cc);
	}
	return(siz-nrem);
}
