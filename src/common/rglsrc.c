#ifndef lint
static const char	RCSid[] = "$Id: rglsrc.c,v 3.3 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Routines for handling OpenGL light sources
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

#include "radogl.h"

double	expval = 0.;			/* global exposure value */

COLOR	ambval = {0.2, 0.2, 0.2};	/* global ambient value */

int	dolights = MAXLIGHTS;		/* do how many more light sources? */

int	glightid[MAXLIGHTS] = {GL_LIGHT0, GL_LIGHT1, GL_LIGHT2,
			GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7};

static int	lightlist;		/* light list id */

				/* source types */
#define L_NONE		0
#define L_SOURCE	1
#define L_FLAT		2
#define L_SPHERE	3

static struct {
	int	type;		/* light type (0 if none) */
	MATREC	*m;		/* light material */
	FVECT	pos;		/* light position (or direction) */
	FVECT	norm;		/* flat source normal */
	double	area;		/* source area (or solid angle) */
} lightrec[MAXLIGHTS];	/* light source list */

static int	nlights;		/* number of defined lights */

static void	l_flatsrc(), l_sphsrc(), l_source();


void
lightinit()			/* initialize lighting */
{
	GLfloat	ambv[4];
	register int	i;

	if (!dolights)
		return;
	glPushAttrib(GL_LIGHTING_BIT);
	if (expval <= FTINY && bright(ambval) > FTINY)
		expval = 0.2/bright(ambval);
	ambv[0] = expval*colval(ambval,RED);
	ambv[1] = expval*colval(ambval,GRN);
	ambv[2] = expval*colval(ambval,BLU);
	ambv[3] = 1.;
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambv);
	glCallList(lightlist = newglist());
	rgl_checkerr("in lightinit");
	nlights = 0;
}


void
lightclean()			/* clean up light source commands */
{
	if ((dolights += nlights) <= 0)
		return;
	glPopAttrib();
}


void
lightdefs()			/* define light source list */
{
	register int	i;

	if (!nlights)
		return;
	glNewList(lightlist, GL_COMPILE);
	while (nlights--) {
		switch (lightrec[nlights].type) {
		case L_FLAT:
			l_flatsrc(nlights);
			break;
		case L_SPHERE:
			l_sphsrc(nlights);
			break;
		case L_SOURCE:
			l_source(nlights);
			break;
		default:
			error(CONSISTENCY, "botched light type in lightdefs");
		}
		freemtl(lightrec[nlights].m);
		lightrec[nlights].type = L_NONE;
	}
	glEndList();
	rgl_checkerr("defining lights");
}


int
o_source(o)			/* record a distant source */
register OBJREC	*o;
{
	if (!dolights || !issrcmat((MATREC *)o->os))
		return(0);
	if (o->oargs.nfargs != 4)
		objerror(o, USER, "bad # real arguments");
						/* record type & material */
	lightrec[nlights].type = L_SOURCE;
	(lightrec[nlights].m = (MATREC *)o->os)->nlinks++;
						/* assign direction */
	VCOPY(lightrec[nlights].pos, o->oargs.farg);
						/* compute solid angle */
	if (o->oargs.farg[3] <= FTINY)
		objerror(o, USER, "zero size");
	lightrec[nlights].area = 2.*PI*(1. - cos(PI/180./2.*o->oargs.farg[3]));
	nlights++; dolights--;
	return(1);
}


int
doflatsrc(m, pos, norm, area)	/* record a flat source */
MATREC	*m;
FVECT	pos, norm;
double	area;
{
	if (!dolights || !issrcmat(m) || area <= FTINY)
		return(0);
						/* record type & material */
	lightrec[nlights].type = L_FLAT;
	(lightrec[nlights].m = m)->nlinks++;
						/* assign geometry */
	VCOPY(lightrec[nlights].pos, pos);
	VCOPY(lightrec[nlights].norm, norm);
	lightrec[nlights].area = area;
	nlights++; dolights--;
	return(1);
}


int
dosphsrc(m, pos, area)		/* record a spherical source */
register MATREC	*m;
FVECT	pos;
double	area;
{
	if (!dolights || !issrcmat(m) || area <= FTINY)
		return(0);
						/* record type & material */
	lightrec[nlights].type = L_SPHERE;
	(lightrec[nlights].m = m)->nlinks++;
						/* assign geometry */
	VCOPY(lightrec[nlights].pos, pos);
	lightrec[nlights].area = area;
	nlights++; dolights--;
	return(1);
}


static void
l_source(n)			/* convert a distant source */
register int	n;
{
	register MATREC	*m = lightrec[n].m;
	int	thislight = glightid[n];
	GLfloat	vec[4];
						/* assign direction */
	VCOPY(vec, lightrec[n].pos);
	vec[3] = 0.;
	glLightfv(thislight, GL_POSITION, vec);
						/* assign color */
	vec[0] = expval*lightrec[n].area*colval(m->u.l.emission,RED);
	vec[1] = expval*lightrec[n].area*colval(m->u.l.emission,GRN);
	vec[2] = expval*lightrec[n].area*colval(m->u.l.emission,BLU);
	vec[3] = 1.;
	glLightfv(thislight, GL_SPECULAR, vec);
	glLightfv(thislight, GL_DIFFUSE, vec);
	vec[0] = vec[1] = vec[2] = 0.; vec[3] = 1.;
	glLightfv(thislight, GL_AMBIENT, vec);
	glEnable(thislight);
}


static void
l_flatsrc(n)			/* convert a flat source */
register int	n;
{
	GLfloat	vec[4];
	register MATREC	*m = lightrec[n].m;
	int	thislight = glightid[n];
						/* assign position */
	VCOPY(vec, lightrec[n].pos); vec[3] = 1.;
	glLightfv(thislight, GL_POSITION, vec);
						/* assign color */
	vec[0] = expval*lightrec[n].area*colval(m->u.l.emission,RED);
	vec[1] = expval*lightrec[n].area*colval(m->u.l.emission,GRN);
	vec[2] = expval*lightrec[n].area*colval(m->u.l.emission,BLU);
	vec[3] = 1.;
	glLightfv(thislight, GL_SPECULAR, vec);
	glLightfv(thislight, GL_DIFFUSE, vec);
	vec[0] = vec[1] = vec[2] = 0.; vec[3] = 1.;
	glLightfv(thislight, GL_AMBIENT, vec);
	glLightf(thislight, GL_SPOT_EXPONENT, 1.);
	glLightf(thislight, GL_CONSTANT_ATTENUATION, 0.);
	glLightf(thislight, GL_LINEAR_ATTENUATION, 0.);
	glLightf(thislight, GL_QUADRATIC_ATTENUATION, 1.);
	if (m->type == MAT_SPOT && m->u.l.spotang < 90.) {
		glLightf(thislight, GL_SPOT_CUTOFF, m->u.l.spotang);
		glLightfv(thislight, GL_SPOT_DIRECTION, m->u.l.spotdir);
	} else {
		glLightf(thislight, GL_SPOT_CUTOFF, 90.);
		VCOPY(vec, lightrec[n].norm);
		glLightfv(thislight, GL_SPOT_DIRECTION, vec);
	}
	glEnable(thislight);
}


static void
l_sphsrc(n)			/* convert a spherical source */
register int	n;
{
	GLfloat	vec[4];
	register MATREC	*m = lightrec[n].m;
	int	thislight = glightid[n];
						/* assign position */
	VCOPY(vec, lightrec[n].pos); vec[3] = 1.;
	glLightfv(thislight, GL_POSITION, vec);
						/* assign color */
	vec[0] = expval*lightrec[n].area*colval(m->u.l.emission,RED);
	vec[1] = expval*lightrec[n].area*colval(m->u.l.emission,GRN);
	vec[2] = expval*lightrec[n].area*colval(m->u.l.emission,BLU);
	vec[3] = 1.;
	glLightfv(thislight, GL_SPECULAR, vec);
	glLightfv(thislight, GL_DIFFUSE, vec);
	vec[0] = vec[1] = vec[2] = 0.; vec[3] = 1.;
	glLightfv(thislight, GL_AMBIENT, vec);
	glLightf(thislight, GL_SPOT_EXPONENT, 0.);
	glLightf(thislight, GL_CONSTANT_ATTENUATION, 0.);
	glLightf(thislight, GL_LINEAR_ATTENUATION, 0.);
	glLightf(thislight, GL_QUADRATIC_ATTENUATION, 1.);
	if (m->type == MAT_SPOT && m->u.l.spotang <= 90.) {
		glLightf(thislight, GL_SPOT_CUTOFF, m->u.l.spotang);
		glLightfv(thislight, GL_SPOT_DIRECTION, m->u.l.spotdir);
	} else
		glLightf(thislight, GL_SPOT_CUTOFF, 180.);
	glEnable(thislight);
}
