#ifndef lint
static const char	RCSid[] = "$Id: raypcalls.c,v 2.9 2004/09/20 16:26:58 greg Exp $";
#endif
/*
 *  raypcalls.c - interface for parallel rendering using Radiance
 *
 *  External symbols declared in ray.h
 */

#include "copyright.h"

/*
 *  These calls are designed similarly to the ones in raycalls.c,
 *  but allow for multiple rendering processes on the same host
 *  machine.  There is no sense in specifying more child processes
 *  than you have processors, but one child may help by allowing
 *  asynchronous ray computation in an interactive program, and
 *  will protect the caller from fatal rendering errors.
 *
 *  You should first read and undrstand the header in raycalls.c,
 *  as some things are explained there that are not repated here.
 *
 *  The first step is opening one or more rendering processes
 *  with a call to ray_pinit(oct, nproc).  Before calling fork(),
 *  ray_pinit() loads the octree and data structures into the
 *  caller's memory.  This permits all sorts of queries that
 *  wouldn't be possible otherwise, without causing any real
 *  memory overhead, since all the static data are shared
 *  between processes.  Rays are then traced using a simple
 *  queuing mechanism, explained below.
 *
 *  The ray queue holds as many rays as there are rendering
 *  processes.  Rays are queued and returned by a single
 *  ray_pqueue() call.  A ray_pqueue() return
 *  value of 0 indicates that no rays are ready
 *  and the queue is not yet full.  A return value of 1
 *  indicates that a ray was returned, though it is probably
 *  not the one you just requested.  Rays may be identified by
 *  the rno member of the RAY struct, which is incremented
 *  by the rayorigin() call, or may be set explicitly by
 *  the caller.  Below is an example call sequence:
 *
 *	myRay.rorg = ( ray origin point )
 *	myRay.rdir = ( normalized ray direction )
 *	myRay.rmax = ( maximum length, or zero for no limit )
 *	rayorigin(&myRay, NULL, PRIMARY, 1.0);
 *	myRay.rno = ( my personal ray identifier )
 *	if (ray_pqueue(&myRay) == 1)
 *		{ do something with results }
 *
 *  Note the differences between this and the simpler ray_trace()
 *  call.  In particular, the call may or may not return a value
 *  in the passed ray structure.  Also, you need to call rayorigin()
 *  yourself, which is normally called for you by ray_trace().  The
 *  benefit is that ray_pqueue() will trace rays faster in
 *  proportion to the number of CPUs you have available on your
 *  system.  If the ray queue is full before the call, ray_pqueue()
 *  will block until a result is ready so it can queue this one.
 *  The global int ray_pnidle indicates the number of currently idle
 *  children.  If you want to check for completed rays without blocking,
 *  or get the results from rays that have been queued without
 *  queuing any new ones, the ray_presult() call is for you:
 *
 *	if (ray_presult(&myRay, 1) == 1)
 *		{ do something with results }
 *
 *  If the second argument is 1, the call won't block when
 *  results aren't ready, but will immediately return 0.
 *  If the second argument is 0, the call will block
 *  until a value is available, returning 0 only if the
 *  queue is completely empty.  A negative return value
 *  indicates that a rendering process died.  If this
 *  happens, ray_close(0) is automatically called to close
 *  all child processes, and ray_pnprocs is set to zero.
 *
 *  If you just want to fill the ray queue without checking for
 *  results, check ray_pnidle and call ray_psend():
 *
 *	while (ray_pnidle) {
 *		( set up ray )
 *		ray_psend(&myRay);
 *	}
 *
 *  Note that it is a fatal error to call ra_psend() when
 *  ray_pnidle is zero.  The ray_presult() and/or ray_pqueue()
 *  functions may be called subsequently to read back the results.
 *
 *  When you are done, you may call ray_pdone(1) to close
 *  all child processes and clean up memory used by Radiance.
 *  Any queued ray calculations will be awaited and discarded.
 *  As with ray_done(), ray_pdone(0) hangs onto data files
 *  and fonts that are likely to be used in subsequent renderings.
 *  Whether you want to bother cleaning up memory or not, you
 *  should at least call ray_pclose(0) to clean the child processes.
 *
 *  Warning:  You cannot affect any of the rendering processes
 *  by changing global parameter values onece ray_pinit() has
 *  been called.  Changing global parameters will have no effect
 *  until the next call to ray_pinit(), which restarts everything.
 *  If you just want to reap children so that you can alter the
 *  rendering parameters without reloading the scene, use the
 *  ray_pclose(0) and ray_popen(nproc) calls to close
 *  then restart the child processes after the changes are made.
 *
 *  Note:  These routines are written to coordinate with the
 *  definitions in raycalls.c, and in fact depend on them.
 *  If you want to trace a ray and get a result synchronously,
 *  use the ray_trace() call to compute it in the parent process
 *  This will not interfere with any subprocess calculations,
 *  but beware that a fatal error may end with a call to quit().
 *
 *  Note:  One of the advantages of using separate processes
 *  is that it gives the calling program some immunity from
 *  fatal rendering errors.  As discussed in raycalls.c,
 *  Radiance tends to throw up its hands and exit at the
 *  first sign of trouble, calling quit() to return control
 *  to the top level.  Although you can avoid exit() with
 *  your own longjmp() in quit(), the cleanup afterwards
 *  is always suspect.  Through the use of subprocesses,
 *  we avoid this pitfall by closing the processes and
 *  returning a negative value from ray_pqueue() or
 *  ray_presult().  If you get a negative value from either
 *  of these calls, you can assume that the processes have
 *  been cleaned up with a call to ray_close(), though you
 *  will have to call ray_pdone() yourself if you want to
 *  free memory.  Obviously, you cannot continue rendering
 *  without risking further errors, but otherwise your
 *  process should not be compromised.
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h> /* XXX platform */

#include  "rtprocess.h"
#include  "ray.h"
#include  "ambient.h"
#include  "selcall.h"

#ifndef RAYQLEN
#define RAYQLEN		16		/* # rays to send at once */
#endif

#ifndef MAX_RPROCS
#if (FD_SETSIZE/2-4 < 64) 
#define MAX_NPROCS	(FD_SETSIZE/2-4)
#else
#define MAX_NPROCS	64		/* max. # rendering processes */
#endif
#endif

extern char	*shm_boundary;		/* boundary of shared memory */

int		ray_pnprocs = 0;	/* number of child processes */
int		ray_pnidle = 0;		/* number of idle children */

static struct child_proc {
	int	pid;				/* child process id */
	int	fd_send;			/* write to child here */
	int	fd_recv;			/* read from child here */
	int	npending;			/* # rays in process */
	unsigned long  rno[RAYQLEN];		/* working on these rays */
} r_proc[MAX_NPROCS];			/* our child processes */

static RAY	r_queue[2*RAYQLEN];	/* ray i/o buffer */
static int	r_send_next;		/* next send ray placement */
static int	r_recv_first;		/* position of first unreported ray */
static int	r_recv_next;		/* next receive ray placement */

#define sendq_full()	(r_send_next >= RAYQLEN)

static int ray_pflush(void);
static void ray_pchild(int	fd_in, int	fd_out);


extern void
ray_pinit(		/* initialize ray-tracing processes */
	char	*otnm,
	int	nproc
)
{
	if (nobjects > 0)		/* close old calculation */
		ray_pdone(0);

	ray_init(otnm);			/* load the shared scene */

	preload_objs();			/* preload auxiliary data */

					/* set shared memory boundary */
	shm_boundary = (char *)malloc(16);
	strcpy(shm_boundary, "SHM_BOUNDARY");

	r_send_next = 0;		/* set up queue */
	r_recv_first = r_recv_next = RAYQLEN;

	ray_popen(nproc);		/* fork children */
}


static int
ray_pflush(void)			/* send queued rays to idle children */
{
	int	nc, n, nw, i, sfirst;

	if ((ray_pnidle <= 0) | (r_send_next <= 0))
		return(0);		/* nothing we can send */
	
	sfirst = 0;			/* divvy up labor */
	nc = ray_pnidle;
	for (i = ray_pnprocs; nc && i--; ) {
		if (r_proc[i].npending > 0)
			continue;	/* child looks busy */
		n = (r_send_next - sfirst)/nc--;
		if (!n)
			continue;
					/* smuggle set size in crtype */
		r_queue[sfirst].crtype = n;
		nw = writebuf(r_proc[i].fd_send, (char *)&r_queue[sfirst],
				sizeof(RAY)*n);
		if (nw != sizeof(RAY)*n)
			return(-1);	/* write error */
		r_proc[i].npending = n;
		while (n--)		/* record ray IDs */
			r_proc[i].rno[n] = r_queue[sfirst+n].rno;
		sfirst += r_proc[i].npending;
		ray_pnidle--;		/* now she's busy */
	}
	if (sfirst != r_send_next)
		error(CONSISTENCY, "code screwup in ray_pflush");
	r_send_next = 0;
	return(sfirst);			/* return total # sent */
}


extern void
ray_psend(			/* add a ray to our send queue */
	RAY	*r
)
{
	if (r == NULL)
		return;
					/* flush output if necessary */
	if (sendq_full() && ray_pflush() <= 0)
		error(INTERNAL, "ray_pflush failed in ray_psend");

	r_queue[r_send_next] = *r;
	r_send_next++;
}


extern int
ray_pqueue(			/* queue a ray for computation */
	RAY	*r
)
{
	if (r == NULL)
		return(0);
					/* check for full send queue */
	if (sendq_full()) {
		RAY	mySend;
		int	rval;
		mySend = *r;
					/* wait for a result */
		rval = ray_presult(r, 0);
					/* put new ray in queue */
		r_queue[r_send_next] = mySend;
		r_send_next++;
		return(rval);		/* done */
	}
					/* add ray to send queue */
	r_queue[r_send_next] = *r;
	r_send_next++;
					/* check for returned ray... */
	if (r_recv_first >= r_recv_next)
		return(0);
					/* ...one is sitting in queue */
	*r = r_queue[r_recv_first];
	r_recv_first++;
	return(1);
}


extern int
ray_presult(		/* check for a completed ray */
	RAY	*r,
	int	poll
)
{
	static struct timeval	tpoll;	/* zero timeval struct */
	static fd_set	readset, errset;
	int	n, ok;
	register int	pn;

	if (r == NULL)
		return(0);
					/* check queued results first */
	if (r_recv_first < r_recv_next) {
		*r = r_queue[r_recv_first];
		r_recv_first++;
		return(1);
	}
	n = ray_pnprocs - ray_pnidle;	/* pending before flush? */

	if (ray_pflush() < 0)		/* send new rays to process */
		return(-1);
					/* reset receive queue */
	r_recv_first = r_recv_next = RAYQLEN;

	if (!poll)			/* count newly sent unless polling */
		n = ray_pnprocs - ray_pnidle;
	if (n <= 0)			/* return if nothing to await */
		return(0);
getready:				/* any children waiting for us? */
	for (pn = ray_pnprocs; pn--; )
		if (FD_ISSET(r_proc[pn].fd_recv, &readset) ||
				FD_ISSET(r_proc[pn].fd_recv, &errset))
			break;
					/* call select if we must */
	if (pn < 0) {
		FD_ZERO(&readset); FD_ZERO(&errset); n = 0;
		for (pn = ray_pnprocs; pn--; ) {
			if (r_proc[pn].npending > 0)
				FD_SET(r_proc[pn].fd_recv, &readset);
			FD_SET(r_proc[pn].fd_recv, &errset);
			if (r_proc[pn].fd_recv >= n)
				n = r_proc[pn].fd_recv + 1;
		}
					/* find out who is ready */
		while ((n = select(n, &readset, (fd_set *)NULL, &errset,
				poll ? &tpoll : (struct timeval *)NULL)) < 0)
			if (errno != EINTR) {
				error(WARNING,
					"select call failed in ray_presult");
				ray_pclose(0);
				return(-1);
			}
		if (n > 0)		/* go back and get it */
			goto getready;
		return(0);		/* else poll came up empty */
	}
	if (r_recv_next + r_proc[pn].npending > sizeof(r_queue)/sizeof(RAY))
		error(CONSISTENCY, "buffer shortage in ray_presult()");

					/* read rendered ray data */
	n = readbuf(r_proc[pn].fd_recv, (char *)&r_queue[r_recv_next],
			sizeof(RAY)*r_proc[pn].npending);
	if (n > 0) {
		r_recv_next += n/sizeof(RAY);
		ok = (n == sizeof(RAY)*r_proc[pn].npending);
	} else
		ok = 0;
					/* reset child's status */
	FD_CLR(r_proc[pn].fd_recv, &readset);
	if (n <= 0)
		FD_CLR(r_proc[pn].fd_recv, &errset);
	r_proc[pn].npending = 0;
	ray_pnidle++;
					/* check for rendering errors */
	if (!ok) {
		ray_pclose(0);		/* process died -- clean up */
		return(-1);
	}
					/* preen returned rays */
	for (n = r_recv_next - r_recv_first; n--; ) {
		register RAY	*rp = &r_queue[r_recv_first + n];
		rp->rno = r_proc[pn].rno[n];
		rp->parent = NULL;
		rp->newcset = rp->clipset = NULL;
		rp->rox = NULL;
		rp->slights = NULL;
	}
					/* return first ray received */
	*r = r_queue[r_recv_first];
	r_recv_first++;
	return(1);
}


extern void
ray_pdone(		/* reap children and free data */
	int	freall
)
{
	ray_pclose(0);			/* close child processes */

	if (shm_boundary != NULL) {	/* clear shared memory boundary */
		free((void *)shm_boundary);
		shm_boundary = NULL;
	}
	ray_done(freall);		/* free rendering data */
}


static void
ray_pchild(	/* process rays (never returns) */
	int	fd_in,
	int	fd_out
)
{
	int	n;
	register int	i;
					/* read each ray request set */
	while ((n = read(fd_in, (char *)r_queue, sizeof(r_queue))) > 0) {
		int	n2;
		if (n % sizeof(RAY))
			break;
		n /= sizeof(RAY);
					/* get smuggled set length */
		n2 = r_queue[0].crtype - n;
		if (n2 < 0)
			error(INTERNAL, "buffer over-read in ray_pchild");
		if (n2 > 0) {		/* read the rest of the set */
			i = readbuf(fd_in, (char *)(r_queue+n),
					sizeof(RAY)*n2);
			if (i != sizeof(RAY)*n2)
				break;
			n += n2;
		}
					/* evaluate rays */
		for (i = 0; i < n; i++) {
			r_queue[i].crtype = r_queue[i].rtype;
			r_queue[i].parent = NULL;
			r_queue[i].clipset = NULL;
			r_queue[i].slights = NULL;
			r_queue[i].revf = raytrace;
			samplendx++;
			rayclear(&r_queue[i]);
			rayvalue(&r_queue[i]);
		}
					/* write back our results */
		i = writebuf(fd_out, (char *)r_queue, sizeof(RAY)*n);
		if (i != sizeof(RAY)*n)
			error(SYSTEM, "write error in ray_pchild");
	}
	if (n)
		error(SYSTEM, "read error in ray_pchild");
	ambsync();
	quit(0);			/* normal exit */
}


extern void
ray_popen(			/* open the specified # processes */
	int	nadd
)
{
					/* check if our table has room */
	if (ray_pnprocs + nadd > MAX_NPROCS)
		nadd = MAX_NPROCS - ray_pnprocs;
	if (nadd <= 0)
		return;
	fflush(stderr);			/* clear pending output */
	fflush(stdout);
	while (nadd--) {		/* fork each new process */
		int	p0[2], p1[2];
		if (pipe(p0) < 0 || pipe(p1) < 0)
			error(SYSTEM, "cannot create pipe");
		if ((r_proc[ray_pnprocs].pid = fork()) == 0) {
			int	pn;	/* close others' descriptors */
			for (pn = ray_pnprocs; pn--; ) {
				close(r_proc[pn].fd_send);
				close(r_proc[pn].fd_recv);
			}
			close(p0[0]); close(p1[1]);
					/* following call never returns */
			ray_pchild(p1[0], p0[1]);
		}
		if (r_proc[ray_pnprocs].pid < 0)
			error(SYSTEM, "cannot fork child process");
		close(p1[0]); close(p0[1]);
		/*
		 * Close write stream on exec to avoid multiprocessing deadlock.
		 * No use in read stream without it, so set flag there as well.
		 */
		fcntl(p1[1], F_SETFD, FD_CLOEXEC);
		fcntl(p0[0], F_SETFD, FD_CLOEXEC);
		r_proc[ray_pnprocs].fd_send = p1[1];
		r_proc[ray_pnprocs].fd_recv = p0[0];
		r_proc[ray_pnprocs].npending = 0;
		ray_pnprocs++;
		ray_pnidle++;
	}
}


extern void
ray_pclose(		/* close one or more child processes */
	int	nsub
)
{
	static int	inclose = 0;
	RAY	res;
					/* check recursion */
	if (inclose)
		return;
	inclose++;
					/* check argument */
	if ((nsub <= 0) | (nsub > ray_pnprocs))
		nsub = ray_pnprocs;
					/* clear our ray queue */
	while (ray_presult(&res,0) > 0)
		;
					/* clean up children */
	while (nsub--) {
		int	status;
		ray_pnprocs--;
		close(r_proc[ray_pnprocs].fd_recv);
		close(r_proc[ray_pnprocs].fd_send);
		if (waitpid(r_proc[ray_pnprocs].pid, &status, 0) < 0)
			status = 127<<8;
		if (status) {
			sprintf(errmsg,
				"rendering process %d exited with code %d",
					r_proc[ray_pnprocs].pid, status>>8);
			error(WARNING, errmsg);
		}
		ray_pnidle--;
	}
	inclose--;
}


void
quit(ec)			/* make sure exit is called */
int	ec;
{
	exit(ec);
}
