#ifndef lint
static const char	RCSid[] = "$Id: cone.c,v 2.5 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  cone.c - routines for making cones
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

#include  "otypes.h"

#include  "cone.h"

/*
 *     In general, a cone may be any one of a cone, a cylinder, a ring,
 *  a cup (inverted cone), or a tube (inverted cylinder).
 *     Most cones are specified with a starting point and radius and
 *  an ending point and radius.  In the cases of a cylinder or tube,
 *  only one radius is needed.  In the case of a ring, a normal direction
 *  is specified instead of a second endpoint.
 *
 *	mtype (cone|cup) name
 *	0
 *	0
 *	8 P0x P0y P0z P1x P1y P1z R0 R1
 *
 *	mtype (cylinder|tube) name
 *	0
 *	0
 *	7 P0x P0y P0z P1x P1y P1z R
 *
 *	mtype ring name
 *	0
 *	0
 *	8 Px Py Pz Nx Ny Nz R0 R1
 */


CONE *
getcone(o, getxf)			/* get cone structure */
register OBJREC  *o;
int  getxf;
{
	int  sgn0, sgn1;
	register CONE  *co;

	if ((co = (CONE *)o->os) == NULL) {

		co = (CONE *)malloc(sizeof(CONE));
		if (co == NULL)
			error(SYSTEM, "out of memory in makecone");

		co->ca = o->oargs.farg;
						/* get radii */
		if (o->otype == OBJ_CYLINDER | o->otype == OBJ_TUBE) {
			if (o->oargs.nfargs != 7)
				goto argcerr;
			if (co->ca[6] < -FTINY) {
				objerror(o, WARNING, "negative radius");
				o->otype = o->otype == OBJ_CYLINDER ?
						OBJ_TUBE : OBJ_CYLINDER;
				co->ca[6] = -co->ca[6];
			} else if (co->ca[6] <= FTINY)
				goto raderr;
			co->p0 = 0; co->p1 = 3;
			co->r0 = co->r1 = 6;
		} else {
			if (o->oargs.nfargs != 8)
				goto argcerr;
			if (co->ca[6] < -FTINY) sgn0 = -1;
			else if (co->ca[6] > FTINY) sgn0 = 1;
			else sgn0 = 0;
			if (co->ca[7] < -FTINY) sgn1 = -1;
			else if (co->ca[7] > FTINY) sgn1 = 1;
			else sgn1 = 0;
			if (sgn0+sgn1 == 0)
				goto raderr;
			if (sgn0 < 0 | sgn1 < 0) {
				objerror(o, o->otype==OBJ_RING?USER:WARNING,
					"negative radii");
				o->otype = o->otype == OBJ_CONE ?
						OBJ_CUP : OBJ_CONE;
			}
			co->ca[6] = co->ca[6]*sgn0;
			co->ca[7] = co->ca[7]*sgn1;
			if (co->ca[7] - co->ca[6] > FTINY) {
				if (o->otype == OBJ_RING)
					co->p0 = co->p1 = 0;
				else {
					co->p0 = 0; co->p1 = 3;
				}
				co->r0 = 6; co->r1 = 7;
			} else if (co->ca[6] - co->ca[7] > FTINY) {
				if (o->otype == OBJ_RING)
					co->p0 = co->p1 = 0;
				else {
					co->p0 = 3; co->p1 = 0;
				}
				co->r0 = 7; co->r1 = 6;
			} else {
				if (o->otype == OBJ_RING)
					goto raderr;
				o->otype = o->otype == OBJ_CONE ?
						OBJ_CYLINDER : OBJ_TUBE;
				o->oargs.nfargs = 7;
				co->p0 = 0; co->p1 = 3;
				co->r0 = co->r1 = 6;
			}
		}
						/* get axis orientation */
		if (o->otype == OBJ_RING)
			VCOPY(co->ad, o->oargs.farg+3);
		else {
			co->ad[0] = CO_P1(co)[0] - CO_P0(co)[0];
			co->ad[1] = CO_P1(co)[1] - CO_P0(co)[1];
			co->ad[2] = CO_P1(co)[2] - CO_P0(co)[2];
		}
		co->al = normalize(co->ad);
		if (co->al == 0.0)
			objerror(o, USER, "zero orientation");
					/* compute axis and side lengths */
		if (o->otype == OBJ_RING) {
			co->al = 0.0;
			co->sl = CO_R1(co) - CO_R0(co);
		} else if (o->otype == OBJ_CONE | o->otype == OBJ_CUP) {
			co->sl = co->ca[7] - co->ca[6];
			co->sl = sqrt(co->sl*co->sl + co->al*co->al);
		} else { /* OBJ_CYLINDER or OBJ_TUBE */
			co->sl = co->al;
		}
		co->tm = NULL;
		o->os = (char *)co;
	}
	if (getxf && co->tm == NULL)
		conexform(co);
	return(co);

argcerr:
	objerror(o, USER, "bad # arguments");
raderr:
	objerror(o, USER, "illegal radii");
}


void
freecone(o)			/* free memory associated with cone */
OBJREC  *o;
{
	register CONE  *co = (CONE *)o->os;

	if (co == NULL)
		return;
	if (co->tm != NULL)
		free((void *)co->tm);
	free((void *)co);
	o->os = NULL;
}


void
conexform(co)			/* get cone transformation matrix */
register CONE  *co;
{
	MAT4  m4;
	register double  d;
	register int  i;

	co->tm = (FLOAT (*)[4])malloc(sizeof(MAT4));
	if (co->tm == NULL)
		error(SYSTEM, "out of memory in conexform");

				/* translate to origin */
	setident4(co->tm);
	if (co->r0 == co->r1)
		d = 0.0;
	else
		d = CO_R0(co) / (CO_R1(co) - CO_R0(co));
	for (i = 0; i < 3; i++)
		co->tm[3][i] = d*(CO_P1(co)[i] - CO_P0(co)[i])
				- CO_P0(co)[i];
	
				/* rotate to positive z-axis */
	setident4(m4);
	d = co->ad[1]*co->ad[1] + co->ad[2]*co->ad[2];
	if (d <= FTINY*FTINY) {
		m4[0][0] = 0.0;
		m4[0][2] = co->ad[0];
		m4[2][0] = -co->ad[0];
		m4[2][2] = 0.0;
	} else {
		d = sqrt(d);
		m4[0][0] = d;
		m4[1][0] = -co->ad[0]*co->ad[1]/d;
		m4[2][0] = -co->ad[0]*co->ad[2]/d;
		m4[1][1] = co->ad[2]/d;
		m4[2][1] = -co->ad[1]/d;
		m4[0][2] = co->ad[0];
		m4[1][2] = co->ad[1];
		m4[2][2] = co->ad[2];
	}
	multmat4(co->tm, co->tm, m4);

				/* scale z-axis */
	if (co->p0 != co->p1 & co->r0 != co->r1) {
		setident4(m4);
		m4[2][2] = (CO_R1(co) - CO_R0(co)) / co->al;
		multmat4(co->tm, co->tm, m4);
	}
}
