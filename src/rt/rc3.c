#ifndef lint
static const char RCSid[] = "$Id: rc3.c,v 2.1 2012/06/09 07:16:47 greg Exp $";
#endif
/*
 * Accumulate ray contributions for a set of materials
 * Controlling process for multiple children
 */

#include "rcontrib.h"
#include "platform.h"
#include "rtprocess.h"
#include "selcall.h"

/* Modifier contribution queue (results waiting to be output) */
typedef struct s_binq {
	int		ndx;		/* index for this entry */
	int		n2add;		/* number left to add */
	struct s_binq	*next;		/* next in queue */
	MODCONT		*mca[1];	/* contrib. array (extends struct) */
} BINQ;

static BINQ	*out_bq = NULL;		/* output bin queue */
static BINQ	*free_bq = NULL;	/* free queue entries */

static SUBPROC	kida[MAXPROCESS];	/* child processes */


/* Get new (empty) bin queue entry */
static BINQ *
new_binq()
{
	BINQ	*bp = free_bq;
	int	i;

	if (bp != NULL) {		/* something already available? */
		free_bq = bp->next;
		bp->next = NULL;
		bp->n2add = accumulate-1;
		return(bp);
	}
					/* else allocate fresh */
	bp = (BINQ *)malloc(sizeof(BINQ)+(nmods-1)*sizeof(MODCONT *));
	if (bp == NULL)
		goto memerr;
	for (i = nmods; i--; ) {
		MODCONT	*mp = (MODCONT *)lu_find(&modconttab,modname[i])->data;
		bp->mca[i] = (MODCONT *)malloc(sizeof(MODCONT) +
						sizeof(DCOLOR)*(mp->nbins-1));
		if (bp->mca[i] == NULL)
			goto memerr;
		memcpy(bp->mca[i], mp, sizeof(MODCONT)-sizeof(DCOLOR));
		/* memset(bp->mca[i]->cbin, 0, sizeof(DCOLOR)*mp->nbins); */
	}
	bp->ndx = 0;
	bp->n2add = accumulate-1;
	bp->next = NULL;
	return(bp);
memerr:
	error(SYSTEM, "out of memory in new_binq()");
	return(NULL);
}


/* Free a bin queue entry */
static void
free_binq(BINQ *bp)
{
	int	i;

	if (bp == NULL) {		/* signal to release our free list */
		while ((bp = free_bq) != NULL) {
			free_bq = bp->next;
			for (i = nmods; i--; )
				free(bp->mca[i]);
				/* Note: we don't own bp->mca[i]->binv */
			free(bp);
		}
		return;
	}
					/* clear sums for next use */
/*	for (i = nmods; i--; )
		memset(bp->mca[i]->cbin, 0, sizeof(DCOLOR)*bp->mca[i]->nbins);
*/
	bp->ndx = 0;
	bp->next = free_bq;		/* push onto free list */
	free_bq = bp;
}


/* Add modifier values to accumulation record in queue and clear */
void
queue_modifiers()
{
	MODCONT	*mpin, *mpout;
	int	i, j;

	if ((accumulate > 0) | (out_bq == NULL))
		error(CONSISTENCY, "bad call to queue_modifiers()");

	for (i = nmods; i--; ) {
		mpin = (MODCONT *)lu_find(&modconttab,modname[i])->data;
		mpout = out_bq->mca[i];
		for (j = mpout->nbins; j--; )
			addcolor(mpout->cbin[j], mpin->cbin[j]);
		memset(mpin->cbin, 0, sizeof(DCOLOR)*mpin->nbins);
	}
}


/* Sum one modifier record into another (doesn't update n2add) */
static void
add_modbin(BINQ *dst, BINQ *src)
{
	int	i, j;
	
	for (i = nmods; i--; ) {
		MODCONT	*mpin = src->mca[i];
		MODCONT	*mpout = dst->mca[i];
		for (j = mpout->nbins; j--; )
			addcolor(mpout->cbin[j], mpin->cbin[j]);
	}
}


/* Queue output, catching up with and freeing FIFO entries when possible */
static int
queue_output(BINQ *bp)
{
	int	nout = 0;
	BINQ	*b_last, *b_cur;
	int	i;

	if (accumulate <= 0) {		/* just accumulating? */
		if (out_bq == NULL) {
			bp->next = NULL;
			out_bq = bp;
		} else {
			add_modbin(out_bq, bp);
			free_binq(bp);
		}
		return(0);
	}
	b_last = NULL;			/* else insert in output queue */
	for (b_cur = out_bq; b_cur != NULL && b_cur->ndx < bp->ndx;
				b_cur = b_cur->next)
		b_last = b_cur;
	if (b_last != NULL) {
		bp->next = b_cur;
		b_last->next = bp;
	} else {
		bp->next = out_bq;
		out_bq = bp;
	}
	if (accumulate > 1) {		/* merge accumulation entries */
		b_cur = out_bq;
		while (b_cur->next != NULL) {
			if (b_cur->n2add <= 0 ||
					(b_cur->ndx-1)/accumulate !=
					(b_cur->next->ndx-1)/accumulate) {
				b_cur = b_cur->next;
				continue;
			}
			add_modbin(b_cur, b_cur->next);
			b_cur->n2add--;
			b_last = b_cur->next;
			b_cur->next = b_last->next;
			b_last->next = NULL;
			free_binq(b_last);
		}
	}
					/* output ready results */
	while (out_bq != NULL && (out_bq->ndx == lastdone+1) & !out_bq->n2add) {
		b_cur = out_bq;			/* pop off first entry */
		out_bq = b_cur->next;
		b_cur->next = NULL;
		for (i = 0; i < nmods; i++)	/* output record */
			mod_output(b_cur->mca[i]);
		end_record();
		free_binq(b_cur);		/* free this entry */
		lastdone += accumulate;
		++nout;
	}
	return(nout);
}


/* Put a zero record in results queue */
void
zero_record(int ndx)
{
	BINQ	*bp = new_binq();
	int	i;

	for (i = nmods; i--; )
		memset(bp->mca[i]->cbin, 0, sizeof(DCOLOR)*bp->mca[i]->nbins);
	bp->ndx = ndx;
	queue_output(bp);
}


/* callback to set output spec to NULL (stdout) */
static int
set_stdout(const LUENT *le, void *p)
{
	(*(MODCONT *)le->data).outspec = NULL;
	return(0);
}


/* Start child processes if we can */
int
in_rchild()
{
#ifdef _WIN32
	error(WARNING, "multiprocessing unsupported -- running solo");
	nproc = 1;
	return(1);
#else
					/* try to fork ourselves */
	while (nchild < nproc) {
		int	p0[2], p1[2];
		int	pid;
					/* prepare i/o pipes */
		if (pipe(p0) < 0 || pipe(p1) < 0)
			error(SYSTEM, "pipe() call failed!");
		pid = fork();		/* fork parent process */
		if (pid == 0) {		/* if in child, set up & return true */
			close(p0[1]); close(p1[0]);
			dup2(p0[0], 0); close(p0[0]);
			dup2(p1[1], 1); close(p1[1]);
			inpfmt = (sizeof(RREAL)==sizeof(double)) ? 'd' : 'f';
			outfmt = 'd';
			header = 0;
			waitflush = xres = 1;
			yres = 0;
			raysleft = 0;
			account = accumulate = 1;
			lu_doall(&modconttab, set_stdout, NULL);
			nchild = -1;
			return(1);	/* child return value */
		}
		if (pid < 0)
			error(SYSTEM, "fork() call failed!");
					/* connect our pipes */
		close(p0[0]); close(p1[1]);
		kida[nchild].r = p1[0];
		kida[nchild].w = p0[1];
		kida[nchild].pid = pid;
		kida[nchild].running = -1;
		++nchild;
	}
	return(0);			/* parent return value */
#endif
}


/* Close child processes */
void
end_children()
{
	int	status;
	
	while (nchild-- > 0)
		if ((status = close_process(&kida[nchild])) > 0) {
			sprintf(errmsg,
				"rendering process returned bad status (%d)",
					status);
			error(WARNING, errmsg);
		}
}


/* Wait for the next available child */
static int
next_child_nq(int force_wait)
{
	fd_set	readset, errset;
	int	i, j, n, nr;

	if (!force_wait)		/* see if there's one free */
		for (i = nchild; i--; )
			if (kida[i].running < 0)
				return(i);
					/* prepare select() call */
	FD_ZERO(&readset); FD_ZERO(&errset);
	n = nr = 0;
	for (i = nchild; i--; ) {
		if (kida[i].running > 0) {
			FD_SET(kida[i].r, &readset);
			++nr;
		}
		FD_SET(kida[i].r, &errset);
		if (kida[i].r >= n)
			n = kida[i].r + 1;
	}
	if (!nr)			/* nothing going on */
		return(-1);
	if (nr > 1) {			/* call select if multiple busy */
		errno = 0;
		if (select(n, &readset, NULL, &errset, NULL) < 0)
			error(SYSTEM, "select call error in next_child_nq()");
	} else
		FD_ZERO(&errset);
	n = -1;				/* read results from child(ren) */
	for (i = nchild; i--; ) {
		BINQ	*bq;
		if (FD_ISSET(kida[i].r, &errset))
			error(USER, "rendering process died");
		if (!FD_ISSET(kida[i].r, &readset))
			continue;
		bq = new_binq();	/* get results holder */
		bq->ndx = kida[i].running;
					/* read from child */
		for (j = 0; j < nmods; j++) {
			n = sizeof(DCOLOR)*bq->mca[j]->nbins;
			nr = readbuf(kida[i].r, (char *)bq->mca[j]->cbin, n);
			if (nr != n)
				error(SYSTEM, "read error from render process");
		}
		queue_output(bq);	/* put results in output queue */
		kida[i].running = -1;	/* mark child as available */
		n = i;
	}
	return(n);			/* last available child */
}


/* Run parental oversight loop */
void
parental_loop()
{
	static int	ignore_warning_given = 0;
	FVECT		orgdir[2];
	double		d;
					/* load rays from stdin & process */
#ifdef getc_unlocked
	flockfile(stdin);		/* avoid lock/unlock overhead */
#endif
	while (getvec(orgdir[0]) == 0 && getvec(orgdir[1]) == 0) {
		d = normalize(orgdir[1]);
							/* asking for flush? */
		if ((d == 0.0) & (accumulate != 1)) {
			if (!ignore_warning_given++)
				error(WARNING,
				"dummy ray(s) ignored during accumulation\n");
			continue;
		}
		if ((d == 0.0) | (lastray+1 < lastray)) {
			while (next_child_nq(1) >= 0)
				;
			lastdone = lastray = 0;
		}
		if (d == 0.0) {
			if ((yres <= 0) | (xres <= 0))
				waitflush = 1;		/* flush right after */
			zero_record(++lastray);
		} else {				/* else assign */
			int	avail = next_child_nq(0);
			if (writebuf(kida[avail].w, (char *)orgdir,
					sizeof(FVECT)*2) != sizeof(FVECT)*2)
				error(SYSTEM, "pipe write error");
			kida[avail].running = ++lastray;
		}
		if (raysleft && !--raysleft)
			break;				/* preemptive EOI */
	}
	while (next_child_nq(1) >= 0)		/* empty results queue */
		;
						/* output accumulated record */
	if (accumulate <= 0 || account < accumulate) {
		int	i;
		if (account < accumulate) {
			error(WARNING, "partial accumulation in final record");
			accumulate -= account;
		}
		for (i = 0; i < nmods; i++)
			mod_output(out_bq->mca[i]);
		end_record();
	}
	if (raysleft)
		error(USER, "unexpected EOF on input");
	free_binq(NULL);			/* clean up */
	lu_done(&ofiletab);
}
