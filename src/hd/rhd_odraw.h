/* RCSid $Id: rhd_odraw.h,v 3.10 2003/06/30 14:59:11 schorsch Exp $ */
/*
 * Header for OpenGL cone drawing routines with depth buffer checks.
 *
 * Include after "standard.h"
 */

#include "color.h"
#include "tonemap.h"
#include "rhdriver.h"

extern struct ODview {
	int	sfirst, snext;	/* first sample and first in next view */
	short	hhi, vhi;	/* screen image resolution */
	short	hlow, vlow;	/* block resolution */
	int	n2redraw;	/* approx. number of samples needing redraw */
	struct ODblock {
		short	nsamp;		/* number of samples in block */
		short	nused;		/* number actually allocated */
		int	first;		/* first sample in this block */
		int	free;		/* index for block free list */
		float	pthresh;	/* proximity threshold */
	}	*bmap;		/* low resolution image map */
	int32	*emap;		/* low resolution edge presence map */
	int32	*pmap;		/* high resolution presence map */
	GLfloat	*dmap;		/* high resolution depth map */
} *odView;		/* our view list */

extern int	odNViews;	/* number of views in our list */

extern struct ODsamp {
	union ODfunion {
		float	prox;			/* viewpoint proximity */
		int32	next;			/* next in free list */
		
	} *f;				/* free list next or proximity */
	short		(*ip)[2];	/* image position array */
	TMbright	*brt;		/* encoded brightness array */
	BYTE		(*chr)[3];	/* encoded chrominance array */
	BYTE		(*rgb)[3];	/* tone-mapped color array */
	int32		*redraw;	/* redraw flags */
	int		nsamp;		/* total number of samples */
	char		*base;		/* base of allocated memory */
} odS;			/* sample values */

#ifndef FL4OP
#define FL4OP(f,i,op)	((f)[(i)>>5] op (1L<<((i)&0x1f)))
#define CHK4(f,i)	FL4OP(f,i,&)
#define SET4(f,i)	FL4OP(f,i,|=)
#define CLR4(f,i)	FL4OP(f,i,&=~)
#define TGL4(f,i)	FL4OP(f,i,^=)
#define FL4NELS(n)	(((n)+0x1f)>>5)
#define CLR4ALL(f,n)	memset((char *)(f),'\0',FL4NELS(n)*sizeof(int32))
#endif

#define OMAXDEPTH	32000			/* maximum depth value */

#define nextfree(i)	f[i].next		/* free pointers */
#define closeness(i)	f[i].prox		/* viewpoint proximity */
#define ENDFREE		(-1)			/* free list terminator */

#define odClean()	odInit(odS.nsamp)	/* clear samples */
#define odDone()	odInit(0)		/* free samples */


/*****************************************************************************
 * Interface routines:


int
odInit(nsamps)			: allocate and initialize memory
int	nsamps;			: number of samples to make available

If nsamps is zero, then this becomes a deallocation routine.  If nsamps
is the same as last time, then this merely clears all data.  The dev_auxview()
function may be called to determine the current view(s).  The odAlloc()
function returns the number of samples actually allocated.


void
odSample(c, d, p)		: register new sample value
COLR	c;			: pixel color (RGBE)
FVECT	d;			: ray direction vector
FVECT	p;			: world intersection point

If p is NULL, then the point is at infinity.


void
odDepthMap(vn, dm)		: set depth map for the given view
int	vn;			: view number
GLfloat	*dm;			: depth map

Assign the depth map associated with view number vn.  The map has
been preallocated, and won't be freed until it is no longer needed
(after an odAlloc() call).  All depth values are the projected
distance along the view direction from the view origin.  If dm
is NULL, then there is no depth map associated with this view.


void
odRedraw(vn, x0, y0, x1, y1)	: region needs to be redrawn
int	vn;			: view number
int	x0, y0, x1, y1;		: rectangle to redraw

This call indicates that the given rectangular region in view vn
needs to be redrawn in the next call to odUpdate().


void
odRedrawAll()			: everything needs to be redrawn

Redraw everything in all views on the next call to odUpdate().
Unless odRemap() is called, no new tone mapping will be done.


void
odUpdate(vn)			: update the current view
int	vn;			: view number

Draw all new and undrawn sample values since last call for this view.


void
odRemap(newhist)		: recompute tone mapping
int	newhist;		: flag whether to clear history

Recompute the tone mapping for all the samples in all the views
and redraw them on the next call(s) to odUpdate().  If newhist
is non-zero, then clear the previous sample history.

 **********************************************************************/
