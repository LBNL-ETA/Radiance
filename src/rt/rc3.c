#ifndef lint
static const char RCSid[] = "$Id: rc3.c,v 2.4 2012/06/11 05:07:55 greg Exp $";
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
	int		nadded;		/* accumulated so far */
	struct s_binq	*next;		/* next in queue */
	MODCONT		*mca[1];	/* contrib. array (extends struct) */
} BINQ;

static BINQ	*out_bq = NULL;		/* output bin queue */
static BINQ	*free_bq = NULL;	/* free queue entries */

static SUBPROC	kida[MAXPROCESS];	/* child processes */
static FILE	*inq_fp[MAXPROCESS];	/* input streams */


/* Get new bin queue entry */
static BINQ *
new_binq()
{
	BINQ	*bp;
	int	i;

	if (free_bq != NULL) {		/* something already available? */
		bp = free_bq;
		free_bq = bp->next;
		bp->next = NULL;
		bp->nadded = 1;
		return(bp);
	}
					/* else allocate fresh */
	bp = (BINQ *)malloc(sizeof(BINQ) + sizeof(MODCONT *)*(nmods-1));
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
	bp->nadded = 1;
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


/* Sum one modifier record into another (updates nadded) */
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
	dst->nadded += src->nadded;
}


/* Queue values for later output */
static void
queue_output(BINQ *bp)
{
	BINQ	*b_last, *b_cur;

	if (accumulate <= 0) {		/* just accumulating? */
		if (out_bq == NULL) {
			bp->next = NULL;
			out_bq = bp;
		} else {
			add_modbin(out_bq, bp);
			free_binq(bp);
		}
		return;
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
	if (accumulate == 1)		/* no accumulation? */
		return;
	b_cur = out_bq;			/* else merge accumulation entries */
	while (b_cur->next != NULL) {
		if (b_cur->nadded >= accumulate ||
				(b_cur->ndx-1)/accumulate !=
				(b_cur->next->ndx-1)/accumulate) {
			b_cur = b_cur->next;
			continue;
		}
		add_modbin(b_cur, b_cur->next);
		b_last = b_cur->next;
		b_cur->next = b_last->next;
		b_last->next = NULL;
		free_binq(b_last);
	}
}


/* Count number of records ready for output */
static int
queue_ready()
{
	int	nready = 0;
	BINQ	*bp;

	if (accumulate <= 0)		/* just accumulating? */
		return(0);

	for (bp = out_bq; bp != NULL && bp->nadded >= accumulate &&
				bp->ndx == lastdone+nready*accumulate+1;
				bp = bp->next)
		++nready;

	return(nready);
}


/* Catch up with output queue by producing ready results */
static int
output_catchup(int nmax)
{
	int	nout = 0;
	BINQ	*bp;
	int	i;

	if (accumulate <= 0)		/* just accumulating? */
		return(0);
					/* else output ready results */
	while (out_bq != NULL && out_bq->nadded >= accumulate
				&& out_bq->ndx == lastdone+1) {
		if ((nmax > 0) & (nout >= nmax))
			break;
		bp = out_bq;			/* pop off first entry */
		out_bq = bp->next;
		bp->next = NULL;
		for (i = 0; i < nmods; i++)	/* output record */
			mod_output(bp->mca[i]);
		end_record();
		free_binq(bp);			/* free this entry */
		lastdone += accumulate;
		++nout;
	}
	return(nout);
}


/* Put a zero record in results queue & output */
void
put_zero_record(int ndx)
{
	BINQ	*bp = new_binq();
	int	i;

	for (i = nmods; i--; )
		memset(bp->mca[i]->cbin, 0, sizeof(DCOLOR)*bp->mca[i]->nbins);
	bp->ndx = ndx;
	queue_output(bp);
	output_catchup(0);
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
		errno = 0;
		if (pipe(p0) < 0 || pipe(p1) < 0)
			error(SYSTEM, "pipe() call failed!");
		pid = fork();		/* fork parent process */
		if (pid == 0) {		/* if in child, set up & return true */
			close(p0[1]); close(p1[0]);
			lu_doall(&modconttab, set_stdout, NULL);
			lu_done(&ofiletab);
			while (nchild--) {
				close(kida[nchild].w);
				fclose(inq_fp[nchild]);
			}
			dup2(p0[0], 0); close(p0[0]);
			dup2(p1[1], 1); close(p1[1]);
			inpfmt = (sizeof(RREAL)==sizeof(double)) ? 'd' : 'f';
			outfmt = 'd';
			header = 0;
			waitflush = xres = 1;
			yres = 0;
			raysleft = 0;
			account = accumulate = 1;
			return(1);	/* child return value */
		}
		if (pid < 0)
			error(SYSTEM, "fork() call failed!");
					/* connect parent's pipes */
		close(p0[0]); close(p1[1]);
		kida[nchild].r = p1[0];
		kida[nchild].w = p0[1];
		kida[nchild].pid = pid;
		kida[nchild].running = -1;
		inq_fp[nchild] = fdopen(p1[0], "rb");
		if (inq_fp[nchild] == NULL)
			error(SYSTEM, "out of memory in in_rchild()");
#ifdef getc_unlocked
		flockfile(inq_fp[nchild]);	/* avoid mutex overhead */
#endif
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
	
	while (nchild-- > 0) {
		kida[nchild].r = -1;	/* close(-1) error is ignored */
		if ((status = close_process(&kida[nchild])) > 0) {
			sprintf(errmsg,
				"rendering process returned bad status (%d)",
					status);
			error(WARNING, errmsg);
		}
		fclose(inq_fp[nchild]);	/* performs actual close() */
	}
}


/* Wait for the next available child, managing output queue as well */
static int
next_child_nq(int force_wait)
{
	static struct timeval	polling;
	struct timeval		*pmode;
	fd_set			readset, errset;
	int			i, j, n, nr, nqr;

	if (!force_wait)		/* see if there's one free */
		for (i = nchild; i--; )
			if (kida[i].running < 0)
				return(i);

	nqr = queue_ready();		/* wait mode or polling? */
	if (!nqr | force_wait | (accumulate <= 0))
		pmode = NULL;
	else
		pmode = &polling;
tryagain:
	n = 0;				/* catch up with output? */
	if ((pmode == &polling) & (nqr > nchild))
		n = nqr - nchild;
	if ((pmode == NULL) & (nqr > 0) | (n > 0))
		nqr -= output_catchup(n);
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
	if (!nr)			/* nothing to wait for? */
		return(-1);
	if ((nr > 1) | (pmode == &polling)) {
		errno = 0;
		i = select(n, &readset, NULL, &errset, pmode);
		if (!i) {
			pmode = NULL;	/* try again, blocking this time */
			goto tryagain;
		}
		if (i < 0)
			error(SYSTEM, "select() error in next_child_nq()");
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
			nr = bq->mca[j]->nbins;
			if (fread(bq->mca[j]->cbin, sizeof(DCOLOR), nr,
							inq_fp[i]) != nr)
				error(SYSTEM, "read error from render process");
		}
		queue_output(bq);	/* add results to output queue */
		kida[i].running = -1;	/* mark child as available */
		n = i;
	}
	return(n);			/* first available child */
}


/* Run parental oversight loop */
void
parental_loop()
{
	static int	ignore_warning_given = 0;
	FVECT		orgdir[2];
	double		d;
	int		i;
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
				;			/* clear the queue */
			lastdone = lastray = 0;
		}
		if (d == 0.0) {
			if ((yres <= 0) | (xres <= 0))
				waitflush = 1;		/* flush right after */
			put_zero_record(++lastray);
		} else {				/* else assign ray */
			i = next_child_nq(0);
			if (writebuf(kida[i].w, (char *)orgdir,
					sizeof(orgdir)) != sizeof(orgdir))
				error(SYSTEM, "pipe write error");
			kida[i].running = ++lastray;
		}
		if (raysleft && !--raysleft)
			break;				/* preemptive EOI */
	}
	while (next_child_nq(1) >= 0)		/* empty results queue */
		;
						/* output accumulated record */
	if (accumulate <= 0 || account < accumulate) {
		if (account < accumulate) {
			error(WARNING, "partial accumulation in final record");
			accumulate -= account;
		}
		for (i = 0; i < nmods; i++)
			mod_output(out_bq->mca[i]);
		end_record();
		free_binq(out_bq);
		out_bq = NULL;
	}
	if (raysleft)
		error(USER, "unexpected EOF on input");
	free_binq(NULL);			/* clean up */
	lu_done(&ofiletab);
}
