#ifndef lint
static const char	RCSid[] = "$Id: face.c,v 2.6 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  face.c - routines dealing with polygonal faces.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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
#define  VERTEPS	1e-2		/* allowed vertex error */
#else
#define  VERTEPS	1e-4		/* allowed vertex error */
#endif


FACE *
getface(o)			/* get arguments for a face */
OBJREC  *o;
{
	double  d1;
	int  badvert;
	FVECT  v1, v2, v3;
	register FACE  *f;
	register int  i;

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
	if (dist2(VERTEX(f,0),VERTEX(f,f->nv-1)) <= FTINY*FTINY)
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
	for (i = 1; i < f->nv; i++) {
		d1 = DOT(f->norm, VERTEX(f,i));
		badvert += fabs(d1 - f->offset/i) > VERTEPS;
		f->offset += d1;
	}
	f->offset /= (double)f->nv;
	if (badvert)
		objerror(o, WARNING, "non-planar vertex");
						/* find axis */
	f->ax = fabs(f->norm[0]) > fabs(f->norm[1]) ? 0 : 1;
	if (fabs(f->norm[2]) > fabs(f->norm[f->ax]))
		f->ax = 2;

	return(f);
}


void
freeface(o)			/* free memory associated with face */
OBJREC  *o;
{
	if (o->os == NULL)
		return;
	free(o->os);
	o->os = NULL;
}


int
inface(p, f)			/* determine if point is in face */
FVECT  p;
FACE  *f;
{
	int  ncross, n;
	double  x, y;
	register int  xi, yi;
	register FLOAT  *p0, *p1;

	xi = (f->ax+1)%3;
	yi = (f->ax+2)%3;
	x = p[xi];
	y = p[yi];
	n = f->nv;
	p0 = f->va + 3*(n-1);		/* connect last to first */
	p1 = f->va;
	ncross = 0;
					/* positive x axis cross test */
	while (n--) {
		if ((p0[yi] > y) ^ (p1[yi] > y))
			if (p0[xi] > x && p1[xi] > x)
				ncross++;
			else if (p0[xi] > x || p1[xi] > x)
				ncross += (p1[yi] > p0[yi]) ^
						((p0[yi]-y)*(p1[xi]-x) >
						(p0[xi]-x)*(p1[yi]-y));
		p0 = p1;
		p1 += 3;
	}
	return(ncross & 01);
}
