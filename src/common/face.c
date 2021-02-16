#ifndef lint
static const char RCSid[] = "$Id: face.c,v 2.16 2021/02/16 16:48:11 greg Exp $";
#endif
/*
 *  face.c - routines dealing with polygonal faces.
 */

#include "copyright.h"

#include  "standard.h"

#include  "object.h"

#include  "face.h"

/*
 *	A face is given as a list of 3D vertices.  The normal
 *  direction and therefore the surface orientation is determined
 *  by the ordering of the vertices.  Looking in the direction opposite
 *  the normal (at the front of the face), the vertices will be
 *  listed in counter-clockwise order.
 *	There is no checking done to insure that the edges do not cross
 *  one another.  This was considered too expensive and should be unnecessary.
 *  The last vertex is automatically connected to the first.
 */

#ifdef  SMLFLT
#define  VERTEPS	1e-3		/* allowed vertex error */
#else
#define  VERTEPS	1e-5		/* allowed vertex error */
#endif


FACE *
getface(				/* get arguments for a face */
	OBJREC  *o
)
{
	double  d1;
	int  smalloff, badvert;
	FVECT  v1, v2, v3;
	FACE  *f;
	int  i;

	if ((f = (FACE *)o->os) != NULL)
		return(f);			/* already done */

	f = (FACE *)malloc(sizeof(FACE));
	if (f == NULL)
		error(SYSTEM, "out of memory in makeface");

	if (o->oargs.nfargs < 9 || o->oargs.nfargs % 3)
		objerror(o, USER, "bad # arguments");

	o->os = (char *)f;			/* save face */

	f->va = o->oargs.farg;
	f->nv = o->oargs.nfargs / 3;
						/* check for last==first */
	if (f->nv > 3 && dist2(VERTEX(f,0),VERTEX(f,f->nv-1)) <= FTINY*FTINY)
		f->nv--;
						/* compute area and normal */
	f->norm[0] = f->norm[1] = f->norm[2] = 0.0;
	v1[0] = VERTEX(f,1)[0] - VERTEX(f,0)[0];
	v1[1] = VERTEX(f,1)[1] - VERTEX(f,0)[1];
	v1[2] = VERTEX(f,1)[2] - VERTEX(f,0)[2];
	for (i = 2; i < f->nv; i++) {
		v2[0] = VERTEX(f,i)[0] - VERTEX(f,0)[0];
		v2[1] = VERTEX(f,i)[1] - VERTEX(f,0)[1];
		v2[2] = VERTEX(f,i)[2] - VERTEX(f,0)[2];
		fcross(v3, v1, v2);
		f->norm[0] += v3[0];
		f->norm[1] += v3[1];
		f->norm[2] += v3[2];
		VCOPY(v1, v2);
	}
	f->area = normalize(f->norm);
	if (f->area == 0.0) {
		objerror(o, WARNING, "zero area");	/* used to be fatal */
		f->offset = 0.0;
		f->ax = 0;
		return(f);
	}
	f->area *= 0.5;
						/* compute offset */
	badvert = 0;
	f->offset = DOT(f->norm, VERTEX(f,0));
	smalloff = fabs(f->offset) <= VERTEPS;
	for (i = 1; i < f->nv; i++) {
		d1 = DOT(f->norm, VERTEX(f,i));
		if (smalloff)
			badvert += fabs(d1 - f->offset/i) > VERTEPS;
		else
			badvert += fabs(1.0 - d1*i/f->offset) > VERTEPS;
		f->offset += d1;
	}
	f->offset /= (double)f->nv;
	if (f->nv > 3 && badvert)
		objerror(o, WARNING, "non-planar vertex");
						/* find axis */
	f->ax = (fabs(f->norm[1]) > fabs(f->norm[0]));
	if (fabs(f->norm[2]) > fabs(f->norm[f->ax]))
		f->ax = 2;

	return(f);
}


void
freeface(			/* free memory associated with face */
	OBJREC  *o
)
{
	if (o->os == NULL)
		return;
	free(o->os);
	o->os = NULL;
}


int
inface(				/* determine if point is in face */
	FVECT  p,
	FACE  *f
)
{
	int  ncross, n;
	double  x, y;
	int  tst;
	int  xi, yi;
	RREAL  *p0, *p1;

	if ((xi = f->ax + 1) >= 3) xi -= 3;
	if ((yi = xi + 1) >= 3) yi -= 3;
	x = p[xi];
	y = p[yi];
	n = f->nv;
	p0 = f->va + 3*(n-1);		/* connect last to first */
	p1 = f->va;
	ncross = 0;
					/* positive x axis cross test */
	while (n--) {
		if ((p0[yi] > y) ^ (p1[yi] > y)) {
			tst = (p0[xi] > x) + (p1[xi] > x);
			if (tst == 2)
				ncross++;
			else if (tst) {
				double	prodA = (p0[yi]-y)*(p1[xi]-x);
				double	prodB = (p0[xi]-x)*(p1[yi]-y);
				if (FRELEQ(prodA, prodB))
					return(1);	/* edge case #1 */
				ncross += (p1[yi] > p0[yi]) ^ (prodA > prodB);
			} else if (FRELEQ(p0[xi], x) && FRELEQ(p1[xi], x))
				return(1);		/* edge case #2 */
		} else if (FRELEQ(p0[yi], y) && FRELEQ(p1[yi], y) &&
				(p0[xi] > x) ^ (p1[xi] > x))
			return(1);			/* edge case #3 */
		p0 = p1;
		p1 += 3;
	}
	return(ncross & 01);
}
