#ifndef lint
static const char RCSid[] = "$Id: o_instance.c,v 2.6 2004/03/30 16:13:00 schorsch Exp $";
#endif
/*
 *  o_instance.c - routines for creating octrees for other octrees
 */

#include  "standard.h"

#include  "object.h"

#include  "octree.h"

#include  "instance.h"

#include  "mesh.h"

#include  "plocate.h"

/*
 *	To determine if two cubes intersect:
 *
 *	1) Check to see if any vertices of first cube are inside the
 *	   second (intersection).
 *
 *	2) Check to see if all vertices of first are to one side of
 *	   second (no intersection).
 *
 *	3) Perform 1 and 2 with roles reversed.
 *
 *	4) Check to see if any portion of any edge of second is inside
 *	   first (intersection).
 *
 *	5) If test 4 fails, we have no intersection.
 *
 *	Note that if we were testing two boxes, we would need
 *  to check that neither had any edges inside the other to be sure.
 *	Since an octree is a volume rather than a surface, we will
 *  return a value of 2 if the cube is entirely within the octree.
 */


static int o_cube(CUBE  *cu1, FULLXF  *fxf, CUBE  *cu);


static int
o_cube(			/* determine if cubes intersect */
	CUBE  *cu1,
	FULLXF  *fxf,
	CUBE  *cu
)
{
	static int  vstart[4] = {0, 3, 5, 6};
	FVECT  cumin, cumax;
	FVECT  vert[8];
	FVECT  v1, v2;
	int  vloc, vout;
	register int  i, j;
					/* check if cube vertex in octree */
	for (j = 0; j < 3; j++)
		cumax[j] = (cumin[j] = cu1->cuorg[j]) + cu1->cusize;
	vloc = ABOVE | BELOW;
	vout = 0;
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 3; j++) {
			v1[j] = cu->cuorg[j];
			if (i & 1<<j)
				v1[j] += cu->cusize;
		}
		multp3(v2, v1, fxf->b.xfm);
		if ( (j = plocate(v2, cumin, cumax)) )
			vout++;
		vloc &= j;
	}
	if (vout == 0)			/* all inside */
		return(O_IN);
	if (vout < 8)			/* some inside */
		return(O_HIT);
	if (vloc)			/* all to one side */
		return(O_MISS);
					/* octree vertices in cube? */
	for (j = 0; j < 3; j++)
		cumax[j] = (cumin[j] = cu->cuorg[j]) + cu->cusize;
	vloc = ABOVE | BELOW;
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 3; j++) {
			v1[j] = cu1->cuorg[j];
			if (i & 1<<j)
				v1[j] += cu1->cusize;
		}
		multp3(vert[i], v1, fxf->f.xfm);
		if ( (j = plocate(vert[i], cumin, cumax)) )
			vloc &= j;
		else
			return(O_HIT);	/* vertex inside */
	}
	if (vloc)			/* all to one side */
		return(O_MISS);
					/* check edges */
	for (i = 0; i < 4; i++)
		for (j = 0; j < 3; j++) {
						/* clip modifies vertices! */
			VCOPY(v1, vert[vstart[i]]);
			VCOPY(v2, vert[vstart[i] ^ 1<<j]);
			if (clip(v1, v2, cumin, cumax))
				return(O_HIT);		/* edge inside */
		}

	return(O_MISS);			/* no intersection */
}


int
o_instance(			/* determine if instance intersects */
	OBJREC  *o,
	CUBE  *cu
)
{
	INSTANCE  *ins;
					/* get octree bounds */
	ins = getinstance(o, IO_BOUNDS);
					/* call o_cube to do the work */
	return(o_cube(&ins->obj->scube, &ins->x, cu));
}


int
o_mesh(				/* determine if mesh intersects */
	OBJREC  *o,
	CUBE  *cu
)
{
	MESHINST	*mip;
					/* get mesh bounds */
	mip = getmeshinst(o, IO_BOUNDS);
					/* call o_cube to do the work */
	return(o_cube(&mip->msh->mcube, &mip->x, cu));
}
