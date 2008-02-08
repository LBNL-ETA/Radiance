#ifndef lint
static const char RCSid[] = "$Id: raypwin.c,v 2.1 2008/02/08 18:27:31 greg Exp $";
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

/***** XXX CURRENTLY, THIS IS JUST A COLLECTION OF FATAL STUBS XXX *****/

#include  "ray.h"


extern void
ray_pinit(		/* initialize ray-tracing processes */
	char	*otnm,
	int	nproc
)
{
	error(CONSISTENCY, "parallel ray processing unimplemented");
}


extern void
ray_psend(			/* add a ray to our send queue */
	RAY	*r
)
{
	error(CONSISTENCY, "parallel ray processing unimplemented");
}


extern int
ray_pqueue(			/* queue a ray for computation */
	RAY	*r
)
{
	return(0);
}


extern int
ray_presult(		/* check for a completed ray */
	RAY	*r,
	int	poll
)
{
	return(-1);
}


extern void
ray_pdone(		/* reap children and free data */
	int	freall
)
{
}


extern void
ray_popen(			/* open the specified # processes */
	int	nadd
)
{
	error(CONSISTENCY, "parallel ray processing unimplemented");
}


extern void
ray_pclose(		/* close one or more child processes */
	int	nsub
)
{
}


void
quit(ec)			/* make sure exit is called */
int	ec;
{
	exit(ec);
}
