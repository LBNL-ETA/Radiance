#ifndef lint
static const char RCSid[] = "$Id: raypwin.c,v 2.8 2010/10/26 03:45:35 greg Exp $";
#endif
/*
 *  raypwin.c - interface for parallel rendering using Radiance (Windows ver)
 *
 *  External symbols declared in ray.h
 */

#include "copyright.h"

/*
 * See raypcalls.c for an explanation of these routines.
 */

/***** XXX CURRENTLY, THIS IS JUST A COLLECTION OF IMPOTENT STUBS XXX *****/

#include  "ray.h"

int		ray_pnprocs = 0;	/* number of child processes */
int		ray_pnidle = 0;		/* number of idle children */

static RAY	queued_ray;

void
ray_pinit(		/* initialize ray-tracing processes */
	char	*otnm,
	int	nproc
)
{
	ray_pdone(0);
	ray_init(otnm);
	ray_popen(nproc);
}


int
ray_psend(			/* add a ray to our send queue */
	RAY	*r
)
{
	if (r == NULL)
		return(0);
	if (ray_pnidle <= 0)
		return(0);
	queued_ray = *r;
	ray_pnidle = 0;
	return(1);
}


int
ray_pqueue(			/* queue a ray for computation */
	RAY	*r
)
{
	RNUMBER	rno;

	if (r == NULL)
		return(0);
	if (ray_pnidle <= 0) {
		RAY	new_ray = *r;
		*r = queued_ray;
		queued_ray = new_ray;
	}
	rno = r->rno;
	r->rno = raynum++;
	samplendx++;
	rayvalue(r);
	r->rno = rno;
	return(1);
}


int
ray_presult(		/* check for a completed ray */
	RAY	*r,
	int	poll
)
{
	if (r == NULL)
		return(0);
	if (ray_pnidle <= 0) {
		*r = queued_ray;
		r->rno = raynum++;
		samplendx++;
		rayvalue(r);
		r->rno = queued_ray.rno;
		ray_pnidle = 1;
		return(1);
	}
	return(0);
}


void
ray_pdone(		/* reap children and free data */
	int	freall
)
{
	ray_done(freall);
	ray_pnprocs = ray_pnidle = 0;
}


void
ray_popen(			/* open the specified # processes */
	int	nadd
)
{
	if (ray_pnprocs + nadd > 1) {
		error(WARNING, "only single process supported");
		nadd = 1 - ray_pnprocs;
	}
	ray_pnprocs += nadd;
	ray_pnidle += nadd;
}


void
ray_pclose(		/* close one or more child processes */
	int	nsub
)
{
	if (nsub > ray_pnprocs)
		nsub = ray_pnprocs;
	ray_pnprocs -= nsub;
	if ((ray_pnidle -= nsub) < 0)
		ray_pnidle = 0;
}
