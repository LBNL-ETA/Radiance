#ifndef lint
static const char RCSid[] = "$Id: raypwin.c,v 2.4 2008/02/13 01:06:10 greg Exp $";
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

extern void
ray_pinit(		/* initialize ray-tracing processes */
	char	*otnm,
	int	nproc
)
{
	ray_pdone(0);
	ray_init(otnm);
	ray_popen(nproc);
}


extern void
ray_psend(			/* add a ray to our send queue */
	RAY	*r
)
{
	if (r == NULL)
		return;
	if (ray_pnidle <= 0)
		error(INTERNAL, "illegal call to ray_psend");
	queued_ray = *r;
	ray_pnidle = 0;
}


extern int
ray_pqueue(			/* queue a ray for computation */
	RAY	*r
)
{
	if (r == NULL)
		return(0);
	samplendx++;
	rayvalue(r);
	return(1);
}


extern int
ray_presult(		/* check for a completed ray */
	RAY	*r,
	int	poll
)
{
	if (r == NULL)
		return(0);
	if (ray_pnidle <= 0) {
		samplendx++;
		rayvalue(r);
		*r = queued_ray;
		ray_pnidle = 1;
		return(1);
	}
	return(0);
}


extern void
ray_pdone(		/* reap children and free data */
	int	freall
)
{
	ray_done(freall);
	ray_pnprocs = ray_pnidle = 0;
}


extern void
ray_popen(			/* open the specified # processes */
	int	nadd
)
{
	if (ray_pnprocs + nadd > 1) {
		error(WARNING, "Only single process supported");
		nadd = 1 - ray_pnprocs;
	}
	ray_pnprocs += nadd;
	ray_pnidle += nadd;
}


extern void
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


void
quit(ec)			/* make sure exit is called */
int	ec;
{
	exit(ec);
}
