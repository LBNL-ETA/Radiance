#ifndef lint
static const char	RCSid[] = "$Id: bbox.c,v 2.2 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  bbox.c - routines for bounding box computation.
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

#include  "octree.h"

#include  "otypes.h"

#include  "face.h"

#include  "cone.h"

#include  "instance.h"


add2bbox(o, bbmin, bbmax)		/* expand bounding box to fit object */
register OBJREC  *o;
FVECT  bbmin, bbmax;
{
	CONE  *co;
	FACE  *fo;
	INSTANCE  *io;
	FVECT  v;
	register int  i, j;

	switch (o->otype) {
	case OBJ_SPHERE:
	case OBJ_BUBBLE:
		if (o->oargs.nfargs != 4)
			objerror(o, USER, "bad arguments");
		for (i = 0; i < 3; i++) {
			VCOPY(v, o->oargs.farg);
			v[i] -= o->oargs.farg[3];
			point2bbox(v, bbmin, bbmax);
			v[i] += 2.0 * o->oargs.farg[3];
			point2bbox(v, bbmin, bbmax);
		}
		break;
	case OBJ_FACE:
		fo = getface(o);
		j = fo->nv;
		while (j--)
			point2bbox(VERTEX(fo,j), bbmin, bbmax);
		break;
	case OBJ_CONE:
	case OBJ_CUP:
	case OBJ_CYLINDER:
	case OBJ_TUBE:
	case OBJ_RING:
		co = getcone(o, 0);
		if (o->otype != OBJ_RING)
			circle2bbox(CO_P0(co), co->ad, CO_R0(co), bbmin, bbmax);
		circle2bbox(CO_P1(co), co->ad, CO_R1(co), bbmin, bbmax);
		break;
	case OBJ_INSTANCE:
		io = getinstance(o, IO_BOUNDS);
		for (j = 0; j < 8; j++) {
			for (i = 0; i < 3; i++) {
				v[i] = io->obj->scube.cuorg[i];
				if (j & 1<<i)
					v[i] += io->obj->scube.cusize;
			}
			multp3(v, v, io->x.f.xfm);
			point2bbox(v, bbmin, bbmax);
		}
		break;
	}
}


point2bbox(p, bbmin, bbmax)		/* expand bounding box to fit point */
register FVECT  p, bbmin, bbmax;
{
	register int  i;

	for (i = 0; i < 3; i++) {
		if (p[i] < bbmin[i])
			bbmin[i] = p[i];
		if (p[i] > bbmax[i])
			bbmax[i] = p[i];
	}
}


circle2bbox(cent, norm, rad, bbmin, bbmax)	/* expand bbox to fit circle */
FVECT  cent, norm;
double  rad;
FVECT  bbmin, bbmax;
{
	FVECT  v1, v2;
	register int  i, j;

	for (i = 0; i < 3; i++) {
		v1[0] = v1[1] = v1[2] = 0;
		v1[i] = 1.0;
		fcross(v2, norm, v1);
		if (normalize(v2) == 0.0)
			continue;
		for (j = 0; j < 3; j++)
			v1[j] = cent[j] + rad*v2[j];
		point2bbox(v1, bbmin, bbmax);
		for (j = 0; j < 3; j++)
			v1[j] = cent[j] - rad*v2[j];
		point2bbox(v1, bbmin, bbmax);
	}
}
