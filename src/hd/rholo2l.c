#ifndef lint
static const char	RCSid[] = "$Id: rholo2l.c,v 3.13 2003/06/26 00:58:10 schorsch Exp $";
#endif
/*
 * Routines for local rtrace execution
 */

#include <signal.h>
#include <sys/time.h>

#include "rholo.h"
#include "random.h"
#include "paths.h"
#include "selcall.h"
#include "rtprocess.h"

#ifndef MAXPROC
#define MAXPROC		64
#endif

int	nprocs = 0;				/* running process count */

static char	pfile[] = TEMPLATE;		/* persist file name */

static SUBPROC	rtpd[MAXPROC];		/* process descriptors */
static float	*rtbuf = NULL;			/* allocated i/o buffer */
static int	maxqlen = 0;			/* maximum packets per queue */

static PACKET	*pqueue[MAXPROC];		/* packet queues */
static int	pqlen[MAXPROC];			/* packet queue lengths */


int
start_rtrace()			/* start rtrace process */
{
	static char	buf1[8];
	int	rmaxpack = 0;
	int	psiz, n;
					/* get number of processes */
	if (ncprocs <= 0 || nprocs > 0)
		return(0);
	if (ncprocs > MAXPROC) {
		sprintf(errmsg,
			"number of rtrace processes reduced from %d to %d",
				ncprocs, MAXPROC);
		error(WARNING, errmsg);
		ncprocs = MAXPROC;
	}
	if (rtargv[rtargc-1] != vval(OCTREE)) {
						/* add compulsory options */
		rtargv[rtargc++] = "-i-";
		rtargv[rtargc++] = "-I-";
		rtargv[rtargc++] = "-h-";
		rtargv[rtargc++] = "-ld-";
		sprintf(buf1, "%d", RPACKSIZ);
		rtargv[rtargc++] = "-x"; rtargv[rtargc++] = buf1;
		rtargv[rtargc++] = "-y"; rtargv[rtargc++] = "0";
		rtargv[rtargc++] = "-fff";
		rtargv[rtargc++] = vbool(VDIST) ? "-ovl" : "-ovL";
		if (nowarn)
			rtargv[rtargc++] = "-w-";
		if (ncprocs > 1) {
			mktemp(pfile);
			rtargv[rtargc++] = "-PP"; rtargv[rtargc++] = pfile;
		}
		rtargv[rtargc++] = vval(OCTREE);
		rtargv[rtargc] = NULL;
	}
	maxqlen = 0;
	for (nprocs = 0; nprocs < ncprocs; nprocs++) {	/* spawn children */
		psiz = open_process(&rtpd[nprocs], rtargv);
		if (psiz <= 0)
			error(SYSTEM, "cannot start rtrace process");
		n = psiz/(RPACKSIZ*6*sizeof(float));
		if (maxqlen == 0) {
			if (!(maxqlen = n))
				error(INTERNAL,
					"bad pipe buffer size assumption");
			sleep(2);
		} else if (n != maxqlen)
			error(INTERNAL, "varying pipe buffer size!");
		rmaxpack += n;
	}
	rtbuf = (float *)malloc(RPACKSIZ*6*sizeof(float)*maxqlen);
	if (rtbuf == NULL)
		error(SYSTEM, "malloc failure in start_rtrace");
	return(rmaxpack);
}


static int
bestout()			/* get best process to process packet */
{
	int	cnt;
	register int	pn, i;

	pn = 0;			/* find shortest queue */
	for (i = 1; i < nprocs; i++)
		if (pqlen[i] < pqlen[pn])
			pn = i;
				/* sanity check */
	if (pqlen[pn] == maxqlen)
		return(-1);
	cnt = 0;		/* count number of ties */
	for (i = pn; i < nprocs; i++)
		if (pqlen[i] == pqlen[pn])
			cnt++;
				/* break ties fairly */
	if ((cnt = random() % cnt))
		for (i = pn; i < nprocs; i++)
			if (pqlen[i] == pqlen[pn] && !cnt--)
				return(i);
	return(pn);
}


int
slots_avail()			/* count packet slots available */
{
	register int	nslots = 0;
	register int	i;

	for (i = nprocs; i--; )
		nslots += maxqlen - pqlen[i];
	return(nslots);
}


queue_packet(p)			/* queue up a beam packet */
register PACKET	*p;
{
	int	pn, n;
				/* determine process to write to */
	if ((pn = bestout()) < 0)
		error(INTERNAL, "rtrace input queues are full!");
				/* write out the packet */
	packrays(rtbuf, p);
	if ((n = p->nr) < RPACKSIZ)	/* add flush block? */
		bzero((char *)(rtbuf+6*n++), 6*sizeof(float));
	if (writebuf(rtpd[pn].w, (char *)rtbuf, 6*sizeof(float)*n) < 0)
		error(SYSTEM, "write error in queue_packet");
	p->next = NULL;
	if (!pqlen[pn]++)	/* add it to the end of the queue */
		pqueue[pn] = p;
	else {
		register PACKET	*rpl = pqueue[pn];
		while (rpl->next != NULL)
			rpl = rpl->next;
		rpl->next = p;
	}
}


PACKET *
get_packets(poll)		/* read packets from rtrace processes */
int	poll;
{
	static struct timeval	tpoll;	/* zero timeval struct */
	fd_set	readset, errset;
	PACKET	*pldone = NULL, *plend;
	register PACKET	*p;
	int	n, nr;
	register int	pn;
	float	*bp;
					/* prepare select call */
	FD_ZERO(&readset); FD_ZERO(&errset); n = 0;
	for (pn = nprocs; pn--; ) {
		if (pqlen[pn])
			FD_SET(rtpd[pn].r, &readset);
		FD_SET(rtpd[pn].r, &errset);
		if (rtpd[pn].r >= n)
			n = rtpd[pn].r + 1;
	}
					/* make the call */
	n = select(n, &readset, (fd_set *)NULL, &errset,
			poll ? &tpoll : (struct timeval *)NULL);
	if (n < 0) {
		if (errno == EINTR)	/* interrupted select call */
			return(NULL);
		error(SYSTEM, "select call failure in get_packets");
	}
	if (n == 0)			/* is nothing ready? */
		return(NULL);
					/* make read call(s) */
	for (pn = 0; pn < nprocs; pn++) {
		if (!FD_ISSET(rtpd[pn].r, &readset) &&
				!FD_ISSET(rtpd[pn].r, &errset))
			continue;
	reread:
		n = read(rtpd[pn].r, (char *)rtbuf,
				4*sizeof(float)*RPACKSIZ*pqlen[pn]);
		if (n < 0) {
			if (errno == EINTR | errno == EAGAIN)
				goto reread;
			error(SYSTEM, "read error in get_packets");
		}
		if (n == 0)
			goto eoferr;
		bp = rtbuf;				/* finish processing */
		for (p = pqueue[pn]; n && p != NULL; p = p->next) {
			if ((nr = p->nr) < RPACKSIZ)
				nr++;			/* add flush block */
			n -= 4*sizeof(float)*nr;
			if (n < 0) {			/* get remainder */
				n += readbuf(rtpd[pn].r,
						(char *)(bp+4*nr)+n, -n);
				if (n)
					goto eoferr;
			}
			donerays(p, bp);
			bp += 4*nr;
			pqlen[pn]--;
		}
		if (n)					/* read past end? */
			error(INTERNAL, "packet sync error in get_packets");
							/* take from queue */
		if (pldone == NULL)
			pldone = plend = pqueue[pn];
		else
			plend->next = pqueue[pn];
		while (plend->next != p)
			plend = plend->next;
		plend->next = NULL;
		pqueue[pn] = p;
	}
	return(pldone);				/* return finished packets */
eoferr:
	error(USER, "rtrace process died");
}


PACKET *
do_packets(pl)			/* queue a packet list, return finished */
register PACKET	*pl;
{
	register PACKET	*p;
					/* consistency check */
	if (nprocs < 1)
		error(CONSISTENCY, "do_packets called with no active process");
					/* queue each new packet */
	while (pl != NULL) {
		p = pl; pl = p->next; p->next = NULL;
		queue_packet(p);
	}
	return(get_packets(slots_avail()));	/* return processed packets */
}


PACKET *
flush_queue()			/* empty all rtrace queues */
{
	PACKET	*rpdone = NULL;
	register PACKET	*rpl;
	float	*bp;
	register PACKET	*p;
	int	i, n, nr;

	for (i = 0; i < nprocs; i++)
		if (pqlen[i]) {
			if (rpdone == NULL) {		/* tack on queue */
				rpdone = rpl = pqueue[i];
				if ((nr = rpl->nr) < RPACKSIZ) nr++;
			} else {
				rpl->next = pqueue[i];
				nr = 0;
			}
			while (rpl->next != NULL) {
				nr += (rpl = rpl->next)->nr;
				if (rpl->nr < RPACKSIZ)
					nr++;		/* add flush block */
			}
			n = readbuf(rtpd[i].r, (char *)rtbuf,
					4*sizeof(float)*nr);
			if (n < 0)
				error(SYSTEM, "read failure in flush_queue");
			bp = rtbuf;			/* process packets */
			for (p = pqueue[i]; p != NULL; p = p->next) {
				if ((nr = p->nr) < RPACKSIZ)
					nr++;		/* add flush block */
				n -= 4*sizeof(float)*nr;
				if (n >= 0)
					donerays(p, bp);
				else
					p->nr = 0;	/* short data error */
				bp += 4*nr;
			}
			pqueue[i] = NULL;		/* zero this queue */
			pqlen[i] = 0;
		}
	return(rpdone);		/* return all packets completed */
}


static
killpersist()			/* kill persistent process */
{
	FILE	*fp;
	int	pid;

	if ((fp = fopen(pfile, "r")) == NULL)
		return;
	if (fscanf(fp, "%*s %d", &pid) != 1 || kill(pid, SIGALRM) < 0)
		unlink(pfile);
	fclose(fp);
}


int
end_rtrace()			/* close rtrace process(es) */
{
	int	status = 0, rv;

	if (nprocs > 1)
		killpersist();
	while (nprocs > 0) {
		rv = close_process(&rtpd[--nprocs]);
		if (rv > 0)
			status = rv;
	}
	free((void *)rtbuf);
	rtbuf = NULL;
	maxqlen = 0;
	return(status);
}
