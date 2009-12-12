#ifndef lint
static const char RCSid[] = "$Id: rayfifo.c,v 2.3 2009/12/12 23:08:13 greg Exp $";
#endif
/*
 *  rayfifo.c - parallelize ray queue that respects order
 *
 *  External symbols declared in ray.h
 */

#include "copyright.h"

/*
 *  These routines are essentially an adjunct to raypcalls.c, providing
 *  a convenient means to get first-in/first-out behavior from multiple
 *  processor cores.  The interface is quite simple, with two functions
 *  and a callback, which must be defined by the calling program.  The
 *  hand-off for finished rays is assigned to ray_fifo_out, which takes
 *  a single pointer to the finished ray and returns a non-negative
 *  integer.  If there is an exceptional condition where termination
 *  is desired, a negative value may be returned.
 *
 *  The ray_fifo_in() call takes a ray that has been initialized in
 *  the same manner as for the ray_trace() call, i.e., the origin,
 *  direction, and maximum length have been assigned. The ray number
 *  will be set according to the number of rays traced since the
 *  last call to ray_fifo_flush().  This final call completes all
 *  pending ray calculations and frees the FIFO buffer.  If any of
 *  the automatic calls to the ray_fifo_out callback return a
 *  negative value, processing stops and -1 is returned.
 *
 *  Note:  The ray passed to ray_fifo_in() may be overwritten
 *  arbitrarily, since it is passed to ray_pqueue().
 */

#include  "ray.h"
#include  <string.h>

int		(*ray_fifo_out)(RAY *r) = NULL;	/* ray output callback */

static RAY	*r_fifo_buf = NULL;		/* circular FIFO out buffer */
static int	r_fifo_len = 0;			/* allocated FIFO length */
static RNUMBER	r_fifo_start = 1;		/* first awaited ray */
static RNUMBER	r_fifo_end = 1;			/* one past FIFO last */
static RNUMBER	r_fifo_next = 1;		/* next ray assignment */

#define r_fifo(rn)	(&r_fifo_buf[(rn)&(r_fifo_len-1)])


static void
ray_fifo_growbuf(void)	/* double buffer size (or set to minimum if NULL) */
{
	RAY	*old_buf = r_fifo_buf;
	int	old_len = r_fifo_len;
	int	i;

	if (r_fifo_buf == NULL)
		r_fifo_len = 1<<4;
	else
		r_fifo_len <<= 1;
						/* allocate new */
	r_fifo_buf = (RAY *)calloc(r_fifo_len, sizeof(RAY));
	if (r_fifo_buf == NULL)
		error(SYSTEM, "out of memory in ray_fifo_growbuf");
	if (old_buf == NULL)
		return;
						/* copy old & free */
	for (i = r_fifo_start; i < r_fifo_end; i++)
		*r_fifo(i) = old_buf[i&(old_len-1)];

	free(old_buf);
}


static int
ray_fifo_push(		/* send finished ray to output (or queue it) */
	RAY *r
)
{
	int	rv, nsent = 0;

	if (ray_fifo_out == NULL)
		error(INTERNAL, "ray_fifo_out is NULL");
	if ((r->rno < r_fifo_start) | (r->rno >= r_fifo_next))
		error(INTERNAL, "unexpected ray number in ray_fifo_push()");

	if (r->rno > r_fifo_start) {		/* insert into output queue */
		if (r->rno - r_fifo_start >= r_fifo_len)
			ray_fifo_growbuf();	/* need more space */
		*r_fifo(r->rno) = *r;
		if (r->rno >= r_fifo_end)
			r_fifo_end = r->rno + 1;
		return(0);
	}
			/* r->rno == r_fifo_start, so transfer ray(s) */
	do {
		rv = (*ray_fifo_out)(r);
		r->rno = 0;			/* flag this entry complete */
		if (rv < 0)
			return(-1);
		nsent += rv;
		if (++r_fifo_start < r_fifo_end)
			r = r_fifo(r_fifo_start);
		else if (r_fifo_start > r_fifo_end)
			r_fifo_end = r_fifo_start;
	} while (r->rno == r_fifo_start);

	return(nsent);
}


int
ray_fifo_in(		/* add ray to FIFO */
	RAY	*r
)
{
	static int	incall = 0;		/* prevent recursion */
	int		rv, rval = 0;

	if (incall++)
		error(INTERNAL, "recursive call to ray_fifo_in()");

	if (r_fifo_start >= 1L<<30) {		/* reset our counters */
		if ((rv = ray_fifo_flush()) < 0)
			{--incall; return(-1);}
		rval += rv;
	}
						/* queue ray */
	rayorigin(r, PRIMARY, NULL, NULL);
	r->rno = r_fifo_next++;
	if ((rv = ray_pqueue(r)) < 0)
		{--incall; return(-1);}

	if (!rv)				/* no result this time? */
		{--incall; return(rval);}
	
	do {					/* else send/queue result */
		if ((rv = ray_fifo_push(r)) < 0)
			{--incall; return(-1);}
		rval += rv;

	} while (ray_presult(r, -1) > 0);	/* empty in-core queue */

	--incall; return(rval);
}


int
ray_fifo_flush(void)	/* flush everything and release buffer */
{
	RAY	myRay;
	int	rv, rval = 0;
						/* clear parallel queue */
	while ((rv = ray_presult(&myRay, 0)) > 0 &&
			(rv = ray_fifo_push(&myRay)) >= 0)
		rval += rv;

	if (rv < 0)				/* check for exception */
		return(-1);

	if (r_fifo_start != r_fifo_end)
		error(INTERNAL, "could not empty queue in ray_fifo_flush()");

	if (r_fifo_buf != NULL) {
		free(r_fifo_buf);
		r_fifo_buf = NULL; r_fifo_len = 0;
	}
	r_fifo_next = r_fifo_end = r_fifo_start = 1;

	return(rval);
}
