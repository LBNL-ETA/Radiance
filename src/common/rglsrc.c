#ifndef lint
static const char	RCSid[] = "$Id: rglsrc.c,v 3.6 2003/11/14 17:22:06 schorsch Exp $";
#endif
/*
 * Routines for handling OpenGL light sources
 */

#include "copyright.h"

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

static void	l_flatsrc(int), l_sphsrc(int), l_source(int);


void
lightinit()			/* initialize lighting */
{
	GLfloat	ambv[4];

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
