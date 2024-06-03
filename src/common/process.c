#ifndef lint
static const char	RCSid[] = "$Id: process.c,v 2.12 2024/06/03 18:55:51 greg Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes
 *
 * External symbols declared in rtprocess.h
 */

#include "copyright.h"

#include "rtprocess.h"

/*

The functions open_process() and close_process() exist in
(currently) two versions, which are found in the files:

	win_process.c
	unix_process.c

*/

int
process(		/* process data through pd */
	SUBPROC *pd,
	void	*recvbuf, void *sendbuf,
	int	nbr, int nbs
)
{
	if (!(pd->flags & PF_RUNNING))
		return(-1);
	if (writebuf(pd->w, sendbuf, nbs) < nbs)
		return(-1);
	return(readbuf(pd->r, recvbuf, nbr));
}



ssize_t
readbuf(		/* read all of requested buffer */
	int	fd,
	void	*buf,
	ssize_t	siz
)
{
	char	*bpos = (char *)buf;
	ssize_t	cc = 0, nrem = siz;
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


ssize_t
writebuf(		/* write all of requested buffer */
int	fd,
void	*buf,
ssize_t	siz
)
{
	char	*bpos = (char *)buf;
	ssize_t	cc = 0, nrem = siz;
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
