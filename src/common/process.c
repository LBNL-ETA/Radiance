#ifndef lint
static const char	RCSid[] = "$Id: process.c,v 2.8 2004/09/17 21:43:49 greg Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes
 *
 * External symbols declared in standard.h
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
char	*recvbuf, char *sendbuf,
int	nbr, int nbs
)
{
	if (nbs > PIPE_BUF)
		return(-1);
	if (!pd->running)
		return(-1);
	if (writebuf(pd->w, sendbuf, nbs) < nbs)
		return(-1);
	return(readbuf(pd->r, recvbuf, nbr));
}



int
readbuf(		/* read all of requested buffer */
int	fd,
char	*bpos,
int	siz
)
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
writebuf(		/* write all of requested buffer */
int	fd,
char	*bpos,
int	siz
)
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
