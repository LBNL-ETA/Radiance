#ifndef lint
static const char RCSid[] = "$Id: rglsurf.c,v 3.12 2004/03/30 20:40:03 greg Exp $";
#endif
/*
 * Convert Radiance -> OpenGL surfaces.
 */

#include "copyright.h"

#include "radogl.h"

#ifndef NSLICES
#define NSLICES		18		/* number of quadric slices */
#endif
#ifndef NSTACKS
#define NSTACKS		10		/* number of quadric stacks */
#endif

MATREC	*curmat = NULL;			/* current material */

static int	curpolysize = 0;	/* outputting triangles/quads */

static GLUquadricObj	*gluqo;		/* shared quadric object */
static GLUtesselator	*gluto;		/* shared tessallation object */

static char	*glu_rout = "unk";	/* active GLU routine */

#define NOPOLY()	if (curpolysize) {glEnd(); curpolysize = 0;} else


void
setmaterial(mp, cent, ispoly)	/* prepare for new material */
register MATREC	*mp;
FVECT	cent;
int	ispoly;
{
	if (mp != curmat && domats) {
		NOPOLY();
		domatobj(curmat = mp, cent);
	} else if (!ispoly) {
		NOPOLY();
	}
}


double
polyarea(cent, norm, n, v)	/* compute polygon area & normal */
FVECT	cent, norm;	/* returned center and normal */
int	n;		/* number of vertices */
register FVECT	v[];	/* vertex list */
{
	FVECT	v1, v2, v3;
	double	d;
	register int	i;

	norm[0] = norm[1] = norm[2] = 0.;
	v1[0] = v[1][0] - v[0][0];
	v1[1] = v[1][1] - v[0][1];
	v1[2] = v[1][2] - v[0][2];
	for (i = 2; i < n; i++) {
		v2[0] = v[i][0] - v[0][0];
		v2[1] = v[i][1] - v[0][1];
		v2[2] = v[i][2] - v[0][2];
		fcross(v3, v1, v2);
		norm[0] += v3[0];
		norm[1] += v3[1];
		norm[2] += v3[2];
		VCOPY(v1, v2);
	}
	if (cent != NULL) {			/* compute center also */
		cent[0] = cent[1] = cent[2] = 0.;
		for (i = n; i--; ) {
			cent[0] += v[i][0];
			cent[1] += v[i][1];
			cent[2] += v[i][2];
		}
		d = 1./n;
		cent[0] *= d; cent[1] *= d; cent[2] *= d;
	}
	return(normalize(norm)*.5);
}


static void
glu_error(en)			/* report an error as a warning */
GLenum	en;
{
	sprintf(errmsg, "GLU error %s: %s", glu_rout, gluErrorString(en));
	error(WARNING, errmsg);
}


static void
myCombine(coords, vertex_data, weight, dataOut)
register GLdouble	coords[3];
GLdouble	*vertex_data[4];
GLfloat	weight[4];
GLdouble	**dataOut;
{
	register GLdouble	*newvert;

	newvert = (GLdouble *)malloc(3*sizeof(GLdouble));
	if (newvert == NULL)
		error(SYSTEM, "out of memory in myCombine");
	VCOPY(newvert, coords);		/* no data, just coordinates */
	*dataOut = newvert;
}


static void
newtess()			/* allocate GLU tessellation object */
{
	if ((gluto = gluNewTess()) == NULL)
		error(INTERNAL, "gluNewTess failed");
	gluTessCallback(gluto, GLU_TESS_BEGIN, glBegin);
	gluTessCallback(gluto, GLU_TESS_VERTEX, glVertex3dv);
	gluTessCallback(gluto, GLU_TESS_END, glEnd);
	gluTessCallback(gluto, GLU_TESS_COMBINE, myCombine);
	gluTessCallback(gluto, GLU_TESS_ERROR, glu_error);
	gluTessProperty(gluto, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
}


static void
newquadric()			/* allocate GLU quadric structure */
{
	if ((gluqo = gluNewQuadric()) == NULL)
		error(INTERNAL, "gluNewQuadric failed");
	gluQuadricDrawStyle(gluqo, GLU_FILL);
	gluQuadricCallback(gluqo, GLU_ERROR, glu_error);
}


int
o_face(o)			/* convert a face */
register OBJREC	*o;
{
	double	area;
	FVECT	norm, cent;
	register int	i;

	if ((o->oargs.nfargs < 9) | (o->oargs.nfargs % 3))
		objerror(o, USER, "bad # real arguments");
	area = polyarea(cent, norm, o->oargs.nfargs/3, (FVECT *)o->oargs.farg);
	if (area <= FTINY)
		return(0);
	if (dolights)					/* check for source */
		doflatsrc((MATREC *)o->os, cent, norm, area);
	setmaterial((MATREC *)o->os, cent, 1);		/* set material */
	if (o->oargs.nfargs/3 != curpolysize) {
		if (curpolysize) glEnd();
		curpolysize = o->oargs.nfargs/3;
		if (curpolysize == 3)
			glBegin(GL_TRIANGLES);
		else if (curpolysize == 4)
			glBegin(GL_QUADS);
	}
	glNormal3d((GLdouble)norm[0], (GLdouble)norm[1], (GLdouble)norm[2]);
	if (curpolysize > 4) {
		if (gluto == NULL) newtess();
		glu_rout = "tessellating polygon";
		gluTessNormal(gluto, (GLdouble)norm[0],
				(GLdouble)norm[1], (GLdouble)norm[2]);
		gluTessBeginPolygon(gluto, NULL);
		gluTessBeginContour(gluto);
#ifdef SMLFLT
		error(INTERNAL, "bad code segment in o_face");
#endif
		for (i = 0; i < curpolysize; i++)
			gluTessVertex(gluto, (GLdouble *)(o->oargs.farg+3*i),
					(void *)(o->oargs.farg+3*i));
		gluTessEndContour(gluto);
		gluTessEndPolygon(gluto);
		curpolysize = 0;
	} else {
		for (i = 0; i < curpolysize; i++)
			glVertex3d((GLdouble)o->oargs.farg[3*i],
					(GLdouble)o->oargs.farg[3*i+1],
					(GLdouble)o->oargs.farg[3*i+2]);
	}
	return(0);
}


void
surfclean()			/* clean up surface routines */
{
	setmaterial(NULL, NULL, 0);
	if (gluqo != NULL) {
		gluDeleteQuadric(gluqo);
		gluqo = NULL;
	}
	if (gluto != NULL) {
		gluDeleteTess(gluto);
		gluto = NULL;
	}
	rgl_checkerr("in surfclean");
}


int
o_sphere(o)			/* convert a sphere */
register OBJREC	*o;
{
					/* check arguments */
	if (o->oargs.nfargs != 4)
		objerror(o, USER, "bad # real arguments");
	if (o->oargs.farg[3] < -FTINY) {
		o->otype = o->otype==OBJ_SPHERE ? OBJ_BUBBLE : OBJ_SPHERE;
		o->oargs.farg[3] = -o->oargs.farg[3];
	} else if (o->oargs.farg[3] <= FTINY)
		return(0);
	if (dolights)
		dosphsrc((MATREC *)o->os, o->oargs.farg,
				PI*o->oargs.farg[3]*o->oargs.farg[3]);
	setmaterial((MATREC *)o->os, o->oargs.farg, 0);
	if (gluqo == NULL) newquadric();
	glu_rout = "making sphere";
	gluQuadricOrientation(gluqo,
			o->otype==OBJ_BUBBLE ? GLU_INSIDE : GLU_OUTSIDE);
	gluQuadricNormals(gluqo, GLU_SMOOTH);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslated((GLdouble)o->oargs.farg[0], (GLdouble)o->oargs.farg[1],
			(GLdouble)o->oargs.farg[2]);
	gluSphere(gluqo, (GLdouble)o->oargs.farg[3], NSLICES, NSTACKS);
	glPopMatrix();
	return(0);
}


int
o_cone(o)			/* convert a cone or cylinder */
register OBJREC *o;
{
	double	x1, y1, h, d;
	FVECT	cent;
	register int	iscyl;

	iscyl = (o->otype==OBJ_CYLINDER) | (o->otype==OBJ_TUBE);
	if (o->oargs.nfargs != (iscyl ? 7 : 8))
		objerror(o, USER, "bad # real arguments");
	if (o->oargs.farg[6] < -FTINY) {
		o->oargs.farg[6] = -o->oargs.farg[6];
		if (iscyl)
			o->otype = o->otype==OBJ_CYLINDER ?
					OBJ_TUBE : OBJ_CYLINDER;
		else {
			if ((o->oargs.farg[7] = -o->oargs.farg[7]) < -FTINY)
				objerror(o, USER, "illegal radii");
			o->otype = o->otype==OBJ_CONE ? OBJ_CUP : OBJ_CONE;
		}
	} else if (!iscyl && o->oargs.farg[7] < -FTINY)
		objerror(o, USER, "illegal radii");
	if (o->oargs.farg[6] <= FTINY && (iscyl || o->oargs.farg[7] <= FTINY))
		return(0);
	if (!iscyl) {
		if (o->oargs.farg[6] < 0.)	/* complains for tiny neg's */
			o->oargs.farg[6] = 0.;
		if (o->oargs.farg[7] < 0.)
			o->oargs.farg[7] = 0.;
	}
	h = sqrt(dist2(o->oargs.farg,o->oargs.farg+3));
	if (h <= FTINY)
		return(0);
	cent[0] = .5*(o->oargs.farg[0] + o->oargs.farg[3]);
	cent[1] = .5*(o->oargs.farg[1] + o->oargs.farg[4]);
	cent[2] = .5*(o->oargs.farg[2] + o->oargs.farg[5]);
	setmaterial((MATREC *)o->os, cent, 0);
	if (gluqo == NULL) newquadric();
	glu_rout = "making cylinder";
	gluQuadricOrientation(gluqo, (o->otype==OBJ_CUP) | (o->otype==OBJ_TUBE) ?
			GLU_INSIDE : GLU_OUTSIDE);
	gluQuadricNormals(gluqo, GLU_SMOOTH);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
					/* do base translation */
	glTranslated((GLdouble)o->oargs.farg[0], (GLdouble)o->oargs.farg[1],
			(GLdouble)o->oargs.farg[2]);
					/* compute height & rotation angle */
	x1 = o->oargs.farg[1] - o->oargs.farg[4];
	y1 = o->oargs.farg[3] - o->oargs.farg[0];
	/* z1 = 0; */
	d = x1*x1 + y1*y1;
	if (d <= FTINY*FTINY)
		x1 = 1.;
	else
		d = 180./PI * asin(sqrt(d) / h);
	if (o->oargs.farg[5] < o->oargs.farg[2])
		d = 180. - d;
	if (d > FTINY)
		glRotated(d, (GLdouble)x1, (GLdouble)y1, 0.);
	gluCylinder(gluqo, o->oargs.farg[6], o->oargs.farg[iscyl ? 6 : 7],
			h, NSLICES, 1);
	glPopMatrix();
	return(0);
}


int
o_ring(o)			/* convert a ring */
register OBJREC	*o;
{
	double	x1, y1, d, h;

	if (o->oargs.nfargs != 8)
		objerror(o, USER, "bad # real arguments");
	if (o->oargs.farg[7] < o->oargs.farg[6]) {
		register double	d = o->oargs.farg[7];
		o->oargs.farg[7] = o->oargs.farg[6];
		o->oargs.farg[6] = d;
	}
	if (o->oargs.farg[6] < -FTINY)
		objerror(o, USER, "negative radius");
	if (o->oargs.farg[6] < 0.)		/* complains for tiny neg's */
		o->oargs.farg[6] = 0.;
	if (o->oargs.farg[7] - o->oargs.farg[6] <= FTINY)
		return(0);
	h = VLEN(o->oargs.farg+3);
	if (h <= FTINY)
		return(0);
	if (dolights)
		doflatsrc((MATREC *)o->os, o->oargs.farg, o->oargs.farg+3, 
				PI*(o->oargs.farg[7]*o->oargs.farg[7] -
					o->oargs.farg[6]*o->oargs.farg[6]));
	setmaterial((MATREC *)o->os, o->oargs.farg, 0);
	if (gluqo == NULL) newquadric();
	glu_rout = "making disk";
	gluQuadricOrientation(gluqo, GLU_OUTSIDE);
	gluQuadricNormals(gluqo, GLU_FLAT);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslated((GLdouble)o->oargs.farg[0], (GLdouble)o->oargs.farg[1],
			(GLdouble)o->oargs.farg[2]);
					/* compute rotation angle */
	x1 = -o->oargs.farg[4];
	y1 = o->oargs.farg[3];
	/* z1 = 0; */
	d = x1*x1 + y1*y1;
	if (d <= FTINY*FTINY)
		x1 = 1.;
	else
		d = 180./PI * asin(sqrt(d) / h);
	if (o->oargs.farg[5] < 0.)
		d = 180. - d;
	if (d > FTINY)
		glRotated(d, (GLdouble)x1, (GLdouble)y1, 0.);
	gluDisk(gluqo, o->oargs.farg[6], o->oargs.farg[7], NSLICES, 1);
	glPopMatrix();
	return(0);
}
