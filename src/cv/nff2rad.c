#ifndef lint
static const char	RCSid[] = "$Id: nff2rad.c,v 2.6 2003/07/27 22:12:01 schorsch Exp $";
#endif
/*
 * Convert Neutral File Format input to Radiance scene description.
 *
 *	12/9/90		Greg Ward
 *	02/7/92         Peter Averkamp added X11(MTV)color names & 
 *                      fixed some lf's for direct import of MTV 
 *                      source files
 */

/******************************************************************

Since Eric Haines wrote such excellent documentation of his
Neutral File Format, I am just going to reprint it here with
my added comments in braces {}.

Neutral File Format (NFF), by Eric Haines

Draft document #1, 10/3/88

The NFF (Neutral File Format) is designed as a minimal scene description
language.  The language was designed in order to test various rendering
algorithms and efficiency schemes.  It is meant to describe the geometry and
basic surface characteristics of objects, the placement of lights, and the
viewing frustum for the eye.  Some additional information is provided for
esthetic reasons (such as the color of the objects, which is not strictly
necessary for testing rendering algorithms).

Future enhancements include:  circle and torus objects, spline surfaces
with trimming curves, directional lights, characteristics for positional
lights, CSG descriptions, and probably more by the time you read this.
Comments, suggestions, and criticisms are all welcome.

At present the NFF file format is used in conjunction with the SPD (Standard
Procedural Database) software, a package designed to create a variety of
databases for testing rendering schemes.  The SPD package is available
from Netlib and via ftp from drizzle.cs.uoregon.edu.  For more information
about SPD see "A Proposal for Standard Graphics Environments," IEEE Computer
Graphics and Applications, vol. 7, no. 11, November 1987, pp. 3-5.

By providing a minimal interface, NFF is meant to act as a simple format to
allow the programmer to quickly write filters to move from NFF to the
local file format.  Presently the following entities are supported:
     A simple perspective frustum
     A positional (vs. directional) light source description
     A background color description
     A surface properties description
     Polygon, polygonal patch, cylinder/cone, and sphere descriptions

Files are output as lines of text.  For each entity, the first line
defines its type.  The rest of the first line and possibly other lines
contain further information about the entity.  Entities include:

"v"  - viewing vectors and angles	{ optionally creates view file }
"l"  - positional light location	{ it's there, but bad to use }
"b"  - background color			{ ditto }
"f"  - object material properties	{ this is flakey }
"c"  - cone or cylinder primitive
"s"  - sphere primitive
"p"  - polygon primitive
"pp" - polygonal patch primitive	{ interpreted same as p for now }

These are explained in depth below:	{ see conversion routines }

***********************************************************************/

#include <stdio.h>

char	*viewfile = NULL;	/* view parameters file */

char	*progname;


main(argc, argv)		/* convert NFF file to Radiance */
int	argc;
char	*argv[];
{
	int	i;
	
	progname = argv[0];
	for (i = 1; i < argc; i++)
		if (argc-i > 1 && !strcmp(argv[i], "-vf"))
			viewfile = argv[++i];
		else if (!strncmp(argv[i], "-h",2))
			goto userr;
		else
			break;
	if (argc-i > 1)
		goto userr;
	if (argc-i == 1 && freopen(argv[i], "r", stdin) == NULL) {
		perror(argv[i]);
		exit(1);
	}
	init();
	nff2rad();
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-vf viewfile] [input]\n", progname);
	exit(1);
}


init()			/* spit out initial definitions */
{
	printf("# File created by %s\n", progname);
	printf("\nvoid light light\n");
	printf("0\n0\n3 1e6 1e6 1e6\n");
	printf("\nvoid plastic fill\n");
	printf("0\n0\n5 .5 .5 .5 0 0\n");
}


nff2rad()		/* convert NFF on stdin to Radiance on stdout */
{
	register int	c;
	
	while ((c = getchar()) != EOF)
		switch (c) {
		case ' ':			/* white space */
		case '\t':
		case '\n':
		case '\f':
		case '\r':
			continue;
		case '#':			/* comment */
			comment();
			break;
		case 'v':			/* view point */
			view();
			break;
		case 'l':			/* light source */
			light();
			break;
		case 'b':			/* background color */
			background();
			break;
		case 'f':			/* fill material */
			fill();
			break;
		case 'c':			/* cylinder or cone */
			cone();
			break;
		case 's':			/* sphere */
			sphere();
			break;
		case 'p':			/* polygon or patch */
			poly();
			break;
		default:			/* unknown */
			fprintf(stderr, "%c: unknown NFF primitive\n", c);
			exit(1);
		}
}


/*******************************************

Comment.  Description:
    "#" [ string ]

Format:
    # [ string ]

    As soon as a "#" character is detected, the rest of the line is considered
    a comment.
    
******************/

comment()
{
	register int	c;
	
	putchar('#');
	while ((c = getchar()) != EOF) {
		putchar(c);
		if (c == '\n')
			break;
	}
}


/***************************************************

Viewpoint location.  Description:
    "v"
    "from" Fx Fy Fz
    "at" Ax Ay Az
    "up" Ux Uy Uz
    "angle" angle
    "hither" hither
    "resolution" xres yres

Format:

    v
    from %g %g %g
    at %g %g %g
    up %g %g %g
    angle %g
    hither %g
    resolution %d %d

The parameters are:

    From:  the eye location in XYZ.
    At:    a position to be at the center of the image, in XYZ world
	   coordinates.  A.k.a. "lookat".
    Up:    a vector defining which direction is up, as an XYZ vector.
    Angle: in degrees, defined as from the center of top pixel row to
	   bottom pixel row and left column to right column.
    Resolution: in pixels, in x and in y.

  Note that no assumptions are made about normalizing the data (e.g. the
  from-at distance does not have to be 1).  Also, vectors are not
  required to be perpendicular to each other.

  For all databases some viewing parameters are always the same:
    Yon is "at infinity."
    Aspect ratio is 1.0.

  A view entity must be defined before any objects are defined (this
  requirement is so that NFF files can be used by hidden surface machines).

***************/

view()
{
	static FILE	*fp = NULL;
	float	from[3], at[3], up[3], angle;
	
	if (scanf(" from %f %f %f", &from[0], &from[1], &from[2]) != 3)
		goto fmterr;
	if (scanf(" at %f %f %f", &at[0], &at[1], &at[2]) != 3)
		goto fmterr;
	if (scanf(" up %f %f %f", &up[0], &up[1], &up[2]) != 3)
		goto fmterr;
	if (scanf(" angle %f", &angle) != 1)
		goto fmterr;
	scanf(" hither %*f");
	scanf(" resolution %*d %*d");
	if (viewfile != NULL) {
		if (fp == NULL && (fp = fopen(viewfile, "a")) == NULL) {
			perror(viewfile);
			exit(1);
		}
		fprintf(fp,
	"VIEW= -vp %g %g %g -vd %g %g %g -vu %g %g %g -vh %g -vv %g\n",
				from[0], from[1], from[2],
				at[0]-from[0], at[1]-from[1], at[2]-from[2],
				up[0], up[1], up[2],
				angle, angle);
	}
	return;
fmterr:
	fprintf(stderr, "%s: view syntax error\n", progname);
	exit(1);
}


/********************************

Positional light.  A light is defined by XYZ position.  Description:
    "l" X Y Z

Format:
    l %g %g %g

    All light entities must be defined before any objects are defined (this
    requirement is so that NFF files can be used by hidden surface machines).
    Lights have a non-zero intensity of no particular value [this definition
    may change soon, with the addition of an intensity and/or color].

**************************/

light()
{
	static int	nlights = 0;
	register int	c;
	float	x, y, z;

	if (scanf("%f %f %f",&x, &y, &z) != 3) {
	    fprintf(stderr, "%s: light source syntax error\n", progname);
	    exit(1);
	}
	while ((c = getchar()) != EOF && c != '\n')
		;
	printf("\nlight sphere l%d \n", ++nlights);
	printf("0\n0\n4 %g %g %g .01\n", x, y, z);
}


/**************************************************

Background color.  A color is simply RGB with values between 0 and 1:
    "b" R G B

Format:
    b %g %g %g

    If no background color is set, assume RGB = {0,0,0}.

********************/

background()
{
	float	r, g, b;
	char colname[50];
	double cvec[3];

	if (scanf("%s", colname) != 1) {
	    fprintf(stderr,"%s: background syntax error\n",progname);exit(1);
	}
        if(LookupColorByName(colname,cvec)==1){
	    r=cvec[0];g=cvec[1];b=cvec[2];
	}else{
	    if(sscanf(colname,"%f",&r)!=1 ||
	       scanf("%f %f", &g, &b) !=2) {
		fprintf(stderr, "%s: background syntax error\n", progname);
		exit(1);
	    }
	}
	printf("\nvoid glow backg_color\n");
	printf("0\n0\n4 %g %g %g 0\n", r, g, b);
	printf("\nbackg_color source background\n");
	printf("0\n0\n4 0 0 1 360\n");
}


/****************************************************

Fill color and shading parameters.  Description:
     "f" red green blue Kd Ks Shine T index_of_refraction

Format:
    f %g %g %g %g %g %g %g %g

    RGB is in terms of 0.0 to 1.0.

    Kd is the diffuse component, Ks the specular, Shine is the Phong cosine
    power for highlights, T is transmittance (fraction of light passed per
    unit).  Usually, 0 <= Kd <= 1 and 0 <= Ks <= 1, though it is not required
    that Kd + Ks == 1.  Note that transmitting objects ( T > 0 ) are considered
    to have two sides for algorithms that need these (normally objects have
    one side).
  
    The fill color is used to color the objects following it until a new color
    is assigned.

*********************/

fill()
{
	float	r, g, b, d, s, p, t, n;
	char colname[50];
	double cvec[3];

	if (scanf("%s", colname) != 1) {
	    fprintf(stderr,"%s: fill syntax error\n",progname);exit(1);
	}
        if(LookupColorByName(colname,cvec)==1){
	    r=cvec[0];g=cvec[1];b=cvec[2];
	}else{
	    if(sscanf(colname,"%f",&r)!=1 ||
	       scanf("%f %f", &g, &b) !=2) {
		fprintf(stderr, "%s: fill syntax error\n", progname);
		exit(1);
	    }
	}
	if (scanf("%f %f %f %f %f", &d, &s, &p, &t, &n) != 5) {
		fprintf(stderr, "%s: fill material syntax error\n", progname);
		exit(1);
	}
	if (p > 1.)
		p = 1./p;
	if (t > .001) {		/* has transmission */
		if (n > 1.1) {		/* has index of refraction */
			printf("\nvoid dielectric fill\n");
			printf("0\n0\n5 %g %g %g %g 0\n", r, g, b, n);
		} else {		/* transmits w/o refraction */
			printf("\nvoid trans fill\n");
			printf("0\n0\n7 %g %g %g %g 0 %g 1\n",
					r*d, g*d, b*d, s, t);
		}
	} else {		/* no transmission */
		printf("\nvoid plastic fill\n");
		printf("0\n0\n5 %g %g %g %g %g\n", r*d, g*d, b*d, s, p);
	}
}


/*****************************************************

Cylinder or cone.  A cylinder is defined as having a radius and an axis
    defined by two points, which also define the top and bottom edge of the
    cylinder.  A cone is defined similarly, the difference being that the apex
    and base radii are different.  The apex radius is defined as being smaller
    than the base radius.  Note that the surface exists without endcaps.  The
    cone or cylinder description:

    "c"
    base.x base.y base.z base_radius
    apex.x apex.y apex.z apex_radius

Format:
    c
    %g %g %g %g
    %g %g %g %g

    A negative value for both radii means that only the inside of the object is
    visible (objects are normally considered one sided, with the outside
    visible).  Note that the base and apex cannot be coincident for a cylinder
    or cone.

************************/

cone()
{
	static int	ncs = 0;
	int	invert;
	float	x0, y0, z0, x1, y1, z1, r0, r1;
	
	if (scanf("%f %f %f %f %f %f %f %f", &x0, &y0, &z0, &r0,
			&x1, &y1, &z1, &r1) != 8) {
		fprintf(stderr, "%s: cylinder or cone syntax error\n",
				progname);
		exit(1);
	}
	if ( (invert = r0 < 0.) ) {
		r0 = -r0;
		r1 = -r1;
	}
	if (r0-r1 < .001 && r1-r0 < .001) {	/* cylinder */
		printf("\nfill %s c%d \n", invert?"tube":"cylinder", ++ncs);
		printf("0\n0\n7\n");
		printf("\t%g\t%g\t%g\n", x0, y0, z0);
		printf("\t%g\t%g\t%g\n", x1, y1, z1);
		printf("\t%g\n", r0);
	} else {				/* cone */
		printf("\nfill %s c%d \n", invert?"cup":"cone", ++ncs);
		printf("0\n0\n8\n");
		printf("\t%g\t%g\t%g\n", x0, y0, z0);
		printf("\t%g\t%g\t%g\n", x1, y1, z1);
		printf("\t%g\t%g\n", r0, r1);
	}
}


/*****************************************

Sphere.  A sphere is defined by a radius and center position:
    "s" center.x center.y center.z radius

Format:
    s %g %g %g %g

    If the radius is negative, then only the sphere's inside is visible
    (objects are normally considered one sided, with the outside visible).

******************/

sphere()
{
	static int	nspheres = 0;
	float	x, y, z, r;
	
	if (scanf("%f %f %f %f", &x, &y, &z, &r) != 4) {
		fprintf(stderr, "%s: sphere syntax error\n", progname);
		exit(1);
	}
	if (r < 0.) {
		printf("\nfill bubble s%d \n", ++nspheres);
		printf("0\n0\n4 %g %g %g %g\n", x, y, z, -r);
	} else {
		printf("\nfill sphere s%d \n", ++nspheres);
		printf("0\n0\n4 %g %g %g %g\n", x, y, z, r);
	}
}


/*********************************************

Polygon.  A polygon is defined by a set of vertices.  With these databases,
    a polygon is defined to have all points coplanar.  A polygon has only
    one side, with the order of the vertices being counterclockwise as you
    face the polygon (right-handed coordinate system).  The first two edges
    must form a non-zero convex angle, so that the normal and side visibility
    can be determined.  Description:

    "p" total_vertices
    vert1.x vert1.y vert1.z
    [etc. for total_vertices vertices]

Format:
    p %d
    [ %g %g %g ] <-- for total_vertices vertices

--------

Polygonal patch.  A patch is defined by a set of vertices and their normals.
    With these databases, a patch is defined to have all points coplanar.
    A patch has only one side, with the order of the vertices being
    counterclockwise as you face the patch (right-handed coordinate system).
    The first two edges must form a non-zero convex angle, so that the normal
    and side visibility can be determined.  Description:

    "pp" total_vertices
    vert1.x vert1.y vert1.z norm1.x norm1.y norm1.z
    [etc. for total_vertices vertices]

Format:
    pp %d
    [ %g %g %g %g %g %g ] <-- for total_vertices vertices

*******************/

poly()
{
	static int	npolys = 0;
	int	ispatch;
	int	nverts;
	float	x, y, z;
	
	ispatch = getchar();
	if (ispatch != 'p') {
		ungetc(ispatch, stdin);
		ispatch = 0;
	}
	if (scanf("%d", &nverts) != 1)
		goto fmterr;
	printf("\nfill polygon p%d \n", ++npolys);
	printf("0\n0\n%d\n", 3*nverts);
	while (nverts-- > 0) {
		if (scanf("%f %f %f", &x, &y, &z) != 3)
			goto fmterr;
		if (ispatch)
			scanf("%*f %*f %*f");
		printf("\t%g\t%g\t%g\n", x, y, z);
	}
	return;
fmterr:
	fprintf(stderr, "%s: polygon or patch syntax error\n", progname);
	exit(1);
}
/***********************************************************************
 * $Author: schorsch $ (Mark VandeWettering, drizzle.cs.uoregon.edu)
 * $Revision: 2.6 $
 * $Date: 2003/07/27 22:12:01 $
 * $Log: nff2rad.c,v $
 * Revision 2.6  2003/07/27 22:12:01  schorsch
 * Added grouping parens to reduce ambiguity warnings.
 *
 * Revision 2.5  2003/02/22 02:07:23  greg
 * Changes and check-in for 3.5 release
 * Includes new source files and modifications not recorded for many years
 * See ray/doc/notes/ReleaseNotes for notes between 3.1 and 3.5 release
 *
 * Revision 1.2  88/09/12  12:53:47  markv
 * Fixed problem in LookupColorbyName, had return ; and return(0).
 * [ Thank you lint! ]
 * 
 * Revision 1.1  88/09/11  11:00:37  markv
 * Initial revision
 * 
 * Peter Averkamp 92/02/01
 * added complete X11R5 rgb.txt-table, hacked standalone version 
 * for nff2rad
 * 
 ***********************************************************************/

typedef double Flt ;
typedef Flt Vec[3] ;
typedef Vec Point ;
typedef Vec Color ;

#define VecCopy(a,b)     (b)[0]=(a)[0];(b)[1]=(a)[1];(b)[2]=(a)[2];
#define		NCOLORS		(738)

typedef struct t_color_entry {
	char *	ce_name ;
	Vec 	ce_color ;
} ColorEntry ;

#define LESS_THAN -1
#define GREATER_THAN 1
#define EQUAL_TO 0

/*
 * Note: These colors must be in sorted order, because we binary search
 * for them.
 *
 * They were swiped from the X-11 distribution.  Sorry....
 */

ColorEntry Colors[] = {
{"AliceBlue",  {0.941176 , 0.972549 , 1.000000 }},
{"AntiqueWhite",  {0.980392 , 0.921569 , 0.843137 }},
{"AntiqueWhite1",  {1.000000 , 0.937255 , 0.858824 }},
{"AntiqueWhite2",  {0.933333 , 0.874510 , 0.800000 }},
{"AntiqueWhite3",  {0.803922 , 0.752941 , 0.690196 }},
{"AntiqueWhite4",  {0.545098 , 0.513725 , 0.470588 }},
{"BlanchedAlmond",  {1.000000 , 0.921569 , 0.803922 }},
{"BlueViolet",  {0.541176 , 0.168627 , 0.886275 }},
{"CadetBlue",  {0.372549 , 0.619608 , 0.627451 }},
{"CadetBlue1",  {0.596078 , 0.960784 , 1.000000 }},
{"CadetBlue2",  {0.556863 , 0.898039 , 0.933333 }},
{"CadetBlue3",  {0.478431 , 0.772549 , 0.803922 }},
{"CadetBlue4",  {0.325490 , 0.525490 , 0.545098 }},
{"CornflowerBlue",  {0.392157 , 0.584314 , 0.929412 }},
{"DarkGoldenrod",  {0.721569 , 0.525490 , 0.043137 }},
{"DarkGoldenrod1",  {1.000000 , 0.725490 , 0.058824 }},
{"DarkGoldenrod2",  {0.933333 , 0.678431 , 0.054902 }},
{"DarkGoldenrod3",  {0.803922 , 0.584314 , 0.047059 }},
{"DarkGoldenrod4",  {0.545098 , 0.396078 , 0.031373 }},
{"DarkGreen",  {0.000000 , 0.392157 , 0.000000 }},
{"DarkKhaki",  {0.741176 , 0.717647 , 0.419608 }},
{"DarkOliveGreen",  {0.333333 , 0.419608 , 0.184314 }},
{"DarkOliveGreen1",  {0.792157 , 1.000000 , 0.439216 }},
{"DarkOliveGreen2",  {0.737255 , 0.933333 , 0.407843 }},
{"DarkOliveGreen3",  {0.635294 , 0.803922 , 0.352941 }},
{"DarkOliveGreen4",  {0.431373 , 0.545098 , 0.239216 }},
{"DarkOrange",  {1.000000 , 0.549020 , 0.000000 }},
{"DarkOrange1",  {1.000000 , 0.498039 , 0.000000 }},
{"DarkOrange2",  {0.933333 , 0.462745 , 0.000000 }},
{"DarkOrange3",  {0.803922 , 0.400000 , 0.000000 }},
{"DarkOrange4",  {0.545098 , 0.270588 , 0.000000 }},
{"DarkOrchid",  {0.600000 , 0.196078 , 0.800000 }},
{"DarkOrchid1",  {0.749020 , 0.243137 , 1.000000 }},
{"DarkOrchid2",  {0.698039 , 0.227451 , 0.933333 }},
{"DarkOrchid3",  {0.603922 , 0.196078 , 0.803922 }},
{"DarkOrchid4",  {0.407843 , 0.133333 , 0.545098 }},
{"DarkSalmon",  {0.913725 , 0.588235 , 0.478431 }},
{"DarkSeaGreen",  {0.560784 , 0.737255 , 0.560784 }},
{"DarkSeaGreen1",  {0.756863 , 1.000000 , 0.756863 }},
{"DarkSeaGreen2",  {0.705882 , 0.933333 , 0.705882 }},
{"DarkSeaGreen3",  {0.607843 , 0.803922 , 0.607843 }},
{"DarkSeaGreen4",  {0.411765 , 0.545098 , 0.411765 }},
{"DarkSlateBlue",  {0.282353 , 0.239216 , 0.545098 }},
{"DarkSlateGray",  {0.184314 , 0.309804 , 0.309804 }},
{"DarkSlateGray1",  {0.592157 , 1.000000 , 1.000000 }},
{"DarkSlateGray2",  {0.552941 , 0.933333 , 0.933333 }},
{"DarkSlateGray3",  {0.474510 , 0.803922 , 0.803922 }},
{"DarkSlateGray4",  {0.321569 , 0.545098 , 0.545098 }},
{"DarkSlateGrey",  {0.184314 , 0.309804 , 0.309804 }},
{"DarkTurquoise",  {0.000000 , 0.807843 , 0.819608 }},
{"DarkViolet",  {0.580392 , 0.000000 , 0.827451 }},
{"DeepPink",  {1.000000 , 0.078431 , 0.576471 }},
{"DeepPink1",  {1.000000 , 0.078431 , 0.576471 }},
{"DeepPink2",  {0.933333 , 0.070588 , 0.537255 }},
{"DeepPink3",  {0.803922 , 0.062745 , 0.462745 }},
{"DeepPink4",  {0.545098 , 0.039216 , 0.313725 }},
{"DeepSkyBlue",  {0.000000 , 0.749020 , 1.000000 }},
{"DeepSkyBlue1",  {0.000000 , 0.749020 , 1.000000 }},
{"DeepSkyBlue2",  {0.000000 , 0.698039 , 0.933333 }},
{"DeepSkyBlue3",  {0.000000 , 0.603922 , 0.803922 }},
{"DeepSkyBlue4",  {0.000000 , 0.407843 , 0.545098 }},
{"DimGray",  {0.411765 , 0.411765 , 0.411765 }},
{"DimGrey",  {0.411765 , 0.411765 , 0.411765 }},
{"DodgerBlue",  {0.117647 , 0.564706 , 1.000000 }},
{"DodgerBlue1",  {0.117647 , 0.564706 , 1.000000 }},
{"DodgerBlue2",  {0.109804 , 0.525490 , 0.933333 }},
{"DodgerBlue3",  {0.094118 , 0.454902 , 0.803922 }},
{"DodgerBlue4",  {0.062745 , 0.305882 , 0.545098 }},
{"FloralWhite",  {1.000000 , 0.980392 , 0.941176 }},
{"ForestGreen",  {0.133333 , 0.545098 , 0.133333 }},
{"GhostWhite",  {0.972549 , 0.972549 , 1.000000 }},
{"GreenYellow",  {0.678431 , 1.000000 , 0.184314 }},
{"HotPink",  {1.000000 , 0.411765 , 0.705882 }},
{"HotPink1",  {1.000000 , 0.431373 , 0.705882 }},
{"HotPink2",  {0.933333 , 0.415686 , 0.654902 }},
{"HotPink3",  {0.803922 , 0.376471 , 0.564706 }},
{"HotPink4",  {0.545098 , 0.227451 , 0.384314 }},
{"IndianRed",  {0.803922 , 0.360784 , 0.360784 }},
{"IndianRed1",  {1.000000 , 0.415686 , 0.415686 }},
{"IndianRed2",  {0.933333 , 0.388235 , 0.388235 }},
{"IndianRed3",  {0.803922 , 0.333333 , 0.333333 }},
{"IndianRed4",  {0.545098 , 0.227451 , 0.227451 }},
{"LavenderBlush",  {1.000000 , 0.941176 , 0.960784 }},
{"LavenderBlush1",  {1.000000 , 0.941176 , 0.960784 }},
{"LavenderBlush2",  {0.933333 , 0.878431 , 0.898039 }},
{"LavenderBlush3",  {0.803922 , 0.756863 , 0.772549 }},
{"LavenderBlush4",  {0.545098 , 0.513725 , 0.525490 }},
{"LawnGreen",  {0.486275 , 0.988235 , 0.000000 }},
{"LemonChiffon",  {1.000000 , 0.980392 , 0.803922 }},
{"LemonChiffon1",  {1.000000 , 0.980392 , 0.803922 }},
{"LemonChiffon2",  {0.933333 , 0.913725 , 0.749020 }},
{"LemonChiffon3",  {0.803922 , 0.788235 , 0.647059 }},
{"LemonChiffon4",  {0.545098 , 0.537255 , 0.439216 }},
{"LightBlue",  {0.678431 , 0.847059 , 0.901961 }},
{"LightBlue1",  {0.749020 , 0.937255 , 1.000000 }},
{"LightBlue2",  {0.698039 , 0.874510 , 0.933333 }},
{"LightBlue3",  {0.603922 , 0.752941 , 0.803922 }},
{"LightBlue4",  {0.407843 , 0.513725 , 0.545098 }},
{"LightCoral",  {0.941176 , 0.501961 , 0.501961 }},
{"LightCyan",  {0.878431 , 1.000000 , 1.000000 }},
{"LightCyan1",  {0.878431 , 1.000000 , 1.000000 }},
{"LightCyan2",  {0.819608 , 0.933333 , 0.933333 }},
{"LightCyan3",  {0.705882 , 0.803922 , 0.803922 }},
{"LightCyan4",  {0.478431 , 0.545098 , 0.545098 }},
{"LightGoldenrod",  {0.933333 , 0.866667 , 0.509804 }},
{"LightGoldenrod1",  {1.000000 , 0.925490 , 0.545098 }},
{"LightGoldenrod2",  {0.933333 , 0.862745 , 0.509804 }},
{"LightGoldenrod3",  {0.803922 , 0.745098 , 0.439216 }},
{"LightGoldenrod4",  {0.545098 , 0.505882 , 0.298039 }},
{"LightGoldenrodYellow",  {0.980392 , 0.980392 , 0.823529 }},
{"LightGray",  {0.827451 , 0.827451 , 0.827451 }},
{"LightGrey",  {0.827451 , 0.827451 , 0.827451 }},
{"LightPink",  {1.000000 , 0.713725 , 0.756863 }},
{"LightPink1",  {1.000000 , 0.682353 , 0.725490 }},
{"LightPink2",  {0.933333 , 0.635294 , 0.678431 }},
{"LightPink3",  {0.803922 , 0.549020 , 0.584314 }},
{"LightPink4",  {0.545098 , 0.372549 , 0.396078 }},
{"LightSalmon",  {1.000000 , 0.627451 , 0.478431 }},
{"LightSalmon1",  {1.000000 , 0.627451 , 0.478431 }},
{"LightSalmon2",  {0.933333 , 0.584314 , 0.447059 }},
{"LightSalmon3",  {0.803922 , 0.505882 , 0.384314 }},
{"LightSalmon4",  {0.545098 , 0.341176 , 0.258824 }},
{"LightSeaGreen",  {0.125490 , 0.698039 , 0.666667 }},
{"LightSkyBlue",  {0.529412 , 0.807843 , 0.980392 }},
{"LightSkyBlue1",  {0.690196 , 0.886275 , 1.000000 }},
{"LightSkyBlue2",  {0.643137 , 0.827451 , 0.933333 }},
{"LightSkyBlue3",  {0.552941 , 0.713725 , 0.803922 }},
{"LightSkyBlue4",  {0.376471 , 0.482353 , 0.545098 }},
{"LightSlateBlue",  {0.517647 , 0.439216 , 1.000000 }},
{"LightSlateGray",  {0.466667 , 0.533333 , 0.600000 }},
{"LightSlateGrey",  {0.466667 , 0.533333 , 0.600000 }},
{"LightSteelBlue",  {0.690196 , 0.768627 , 0.870588 }},
{"LightSteelBlue1",  {0.792157 , 0.882353 , 1.000000 }},
{"LightSteelBlue2",  {0.737255 , 0.823529 , 0.933333 }},
{"LightSteelBlue3",  {0.635294 , 0.709804 , 0.803922 }},
{"LightSteelBlue4",  {0.431373 , 0.482353 , 0.545098 }},
{"LightYellow",  {1.000000 , 1.000000 , 0.878431 }},
{"LightYellow1",  {1.000000 , 1.000000 , 0.878431 }},
{"LightYellow2",  {0.933333 , 0.933333 , 0.819608 }},
{"LightYellow3",  {0.803922 , 0.803922 , 0.705882 }},
{"LightYellow4",  {0.545098 , 0.545098 , 0.478431 }},
{"LimeGreen",  {0.196078 , 0.803922 , 0.196078 }},
{"MediumAquamarine",  {0.400000 , 0.803922 , 0.666667 }},
{"MediumBlue",  {0.000000 , 0.000000 , 0.803922 }},
{"MediumOrchid",  {0.729412 , 0.333333 , 0.827451 }},
{"MediumOrchid1",  {0.878431 , 0.400000 , 1.000000 }},
{"MediumOrchid2",  {0.819608 , 0.372549 , 0.933333 }},
{"MediumOrchid3",  {0.705882 , 0.321569 , 0.803922 }},
{"MediumOrchid4",  {0.478431 , 0.215686 , 0.545098 }},
{"MediumPurple",  {0.576471 , 0.439216 , 0.858824 }},
{"MediumPurple1",  {0.670588 , 0.509804 , 1.000000 }},
{"MediumPurple2",  {0.623529 , 0.474510 , 0.933333 }},
{"MediumPurple3",  {0.537255 , 0.407843 , 0.803922 }},
{"MediumPurple4",  {0.364706 , 0.278431 , 0.545098 }},
{"MediumSeaGreen",  {0.235294 , 0.701961 , 0.443137 }},
{"MediumSlateBlue",  {0.482353 , 0.407843 , 0.933333 }},
{"MediumSpringGreen",  {0.000000 , 0.980392 , 0.603922 }},
{"MediumTurquoise",  {0.282353 , 0.819608 , 0.800000 }},
{"MediumVioletRed",  {0.780392 , 0.082353 , 0.521569 }},
{"MidnightBlue",  {0.098039 , 0.098039 , 0.439216 }},
{"MintCream",  {0.960784 , 1.000000 , 0.980392 }},
{"MistyRose",  {1.000000 , 0.894118 , 0.882353 }},
{"MistyRose1",  {1.000000 , 0.894118 , 0.882353 }},
{"MistyRose2",  {0.933333 , 0.835294 , 0.823529 }},
{"MistyRose3",  {0.803922 , 0.717647 , 0.709804 }},
{"MistyRose4",  {0.545098 , 0.490196 , 0.482353 }},
{"NavajoWhite",  {1.000000 , 0.870588 , 0.678431 }},
{"NavajoWhite1",  {1.000000 , 0.870588 , 0.678431 }},
{"NavajoWhite2",  {0.933333 , 0.811765 , 0.631373 }},
{"NavajoWhite3",  {0.803922 , 0.701961 , 0.545098 }},
{"NavajoWhite4",  {0.545098 , 0.474510 , 0.368627 }},
{"NavyBlue",  {0.000000 , 0.000000 , 0.501961 }},
{"OldLace",  {0.992157 , 0.960784 , 0.901961 }},
{"OliveDrab",  {0.419608 , 0.556863 , 0.137255 }},
{"OliveDrab1",  {0.752941 , 1.000000 , 0.243137 }},
{"OliveDrab2",  {0.701961 , 0.933333 , 0.227451 }},
{"OliveDrab3",  {0.603922 , 0.803922 , 0.196078 }},
{"OliveDrab4",  {0.411765 , 0.545098 , 0.133333 }},
{"OrangeRed",  {1.000000 , 0.270588 , 0.000000 }},
{"OrangeRed1",  {1.000000 , 0.270588 , 0.000000 }},
{"OrangeRed2",  {0.933333 , 0.250980 , 0.000000 }},
{"OrangeRed3",  {0.803922 , 0.215686 , 0.000000 }},
{"OrangeRed4",  {0.545098 , 0.145098 , 0.000000 }},
{"PaleGoldenrod",  {0.933333 , 0.909804 , 0.666667 }},
{"PaleGreen",  {0.596078 , 0.984314 , 0.596078 }},
{"PaleGreen1",  {0.603922 , 1.000000 , 0.603922 }},
{"PaleGreen2",  {0.564706 , 0.933333 , 0.564706 }},
{"PaleGreen3",  {0.486275 , 0.803922 , 0.486275 }},
{"PaleGreen4",  {0.329412 , 0.545098 , 0.329412 }},
{"PaleTurquoise",  {0.686275 , 0.933333 , 0.933333 }},
{"PaleTurquoise1",  {0.733333 , 1.000000 , 1.000000 }},
{"PaleTurquoise2",  {0.682353 , 0.933333 , 0.933333 }},
{"PaleTurquoise3",  {0.588235 , 0.803922 , 0.803922 }},
{"PaleTurquoise4",  {0.400000 , 0.545098 , 0.545098 }},
{"PaleVioletRed",  {0.858824 , 0.439216 , 0.576471 }},
{"PaleVioletRed1",  {1.000000 , 0.509804 , 0.670588 }},
{"PaleVioletRed2",  {0.933333 , 0.474510 , 0.623529 }},
{"PaleVioletRed3",  {0.803922 , 0.407843 , 0.537255 }},
{"PaleVioletRed4",  {0.545098 , 0.278431 , 0.364706 }},
{"PapayaWhip",  {1.000000 , 0.937255 , 0.835294 }},
{"PeachPuff",  {1.000000 , 0.854902 , 0.725490 }},
{"PeachPuff1",  {1.000000 , 0.854902 , 0.725490 }},
{"PeachPuff2",  {0.933333 , 0.796078 , 0.678431 }},
{"PeachPuff3",  {0.803922 , 0.686275 , 0.584314 }},
{"PeachPuff4",  {0.545098 , 0.466667 , 0.396078 }},
{"PowderBlue",  {0.690196 , 0.878431 , 0.901961 }},
{"RosyBrown",  {0.737255 , 0.560784 , 0.560784 }},
{"RosyBrown1",  {1.000000 , 0.756863 , 0.756863 }},
{"RosyBrown2",  {0.933333 , 0.705882 , 0.705882 }},
{"RosyBrown3",  {0.803922 , 0.607843 , 0.607843 }},
{"RosyBrown4",  {0.545098 , 0.411765 , 0.411765 }},
{"RoyalBlue",  {0.254902 , 0.411765 , 0.882353 }},
{"RoyalBlue1",  {0.282353 , 0.462745 , 1.000000 }},
{"RoyalBlue2",  {0.262745 , 0.431373 , 0.933333 }},
{"RoyalBlue3",  {0.227451 , 0.372549 , 0.803922 }},
{"RoyalBlue4",  {0.152941 , 0.250980 , 0.545098 }},
{"SaddleBrown",  {0.545098 , 0.270588 , 0.074510 }},
{"SandyBrown",  {0.956863 , 0.643137 , 0.376471 }},
{"SeaGreen",  {0.180392 , 0.545098 , 0.341176 }},
{"SeaGreen1",  {0.329412 , 1.000000 , 0.623529 }},
{"SeaGreen2",  {0.305882 , 0.933333 , 0.580392 }},
{"SeaGreen3",  {0.262745 , 0.803922 , 0.501961 }},
{"SeaGreen4",  {0.180392 , 0.545098 , 0.341176 }},
{"SkyBlue",  {0.529412 , 0.807843 , 0.921569 }},
{"SkyBlue1",  {0.529412 , 0.807843 , 1.000000 }},
{"SkyBlue2",  {0.494118 , 0.752941 , 0.933333 }},
{"SkyBlue3",  {0.423529 , 0.650980 , 0.803922 }},
{"SkyBlue4",  {0.290196 , 0.439216 , 0.545098 }},
{"SlateBlue",  {0.415686 , 0.352941 , 0.803922 }},
{"SlateBlue1",  {0.513725 , 0.435294 , 1.000000 }},
{"SlateBlue2",  {0.478431 , 0.403922 , 0.933333 }},
{"SlateBlue3",  {0.411765 , 0.349020 , 0.803922 }},
{"SlateBlue4",  {0.278431 , 0.235294 , 0.545098 }},
{"SlateGray",  {0.439216 , 0.501961 , 0.564706 }},
{"SlateGray1",  {0.776471 , 0.886275 , 1.000000 }},
{"SlateGray2",  {0.725490 , 0.827451 , 0.933333 }},
{"SlateGray3",  {0.623529 , 0.713725 , 0.803922 }},
{"SlateGray4",  {0.423529 , 0.482353 , 0.545098 }},
{"SlateGrey",  {0.439216 , 0.501961 , 0.564706 }},
{"SpringGreen",  {0.000000 , 1.000000 , 0.498039 }},
{"SpringGreen1",  {0.000000 , 1.000000 , 0.498039 }},
{"SpringGreen2",  {0.000000 , 0.933333 , 0.462745 }},
{"SpringGreen3",  {0.000000 , 0.803922 , 0.400000 }},
{"SpringGreen4",  {0.000000 , 0.545098 , 0.270588 }},
{"SteelBlue",  {0.274510 , 0.509804 , 0.705882 }},
{"SteelBlue1",  {0.388235 , 0.721569 , 1.000000 }},
{"SteelBlue2",  {0.360784 , 0.674510 , 0.933333 }},
{"SteelBlue3",  {0.309804 , 0.580392 , 0.803922 }},
{"SteelBlue4",  {0.211765 , 0.392157 , 0.545098 }},
{"VioletRed",  {0.815686 , 0.125490 , 0.564706 }},
{"VioletRed1",  {1.000000 , 0.243137 , 0.588235 }},
{"VioletRed2",  {0.933333 , 0.227451 , 0.549020 }},
{"VioletRed3",  {0.803922 , 0.196078 , 0.470588 }},
{"VioletRed4",  {0.545098 , 0.133333 , 0.321569 }},
{"WhiteSmoke",  {0.960784 , 0.960784 , 0.960784 }},
{"YellowGreen",  {0.603922 , 0.803922 , 0.196078 }},
{"alice_blue",  {0.941176 , 0.972549 , 1.000000 }},
{"antique_white",  {0.980392 , 0.921569 , 0.843137 }},
{"aquamarine",  {0.498039 , 1.000000 , 0.831373 }},
{"aquamarine1",  {0.498039 , 1.000000 , 0.831373 }},
{"aquamarine2",  {0.462745 , 0.933333 , 0.776471 }},
{"aquamarine3",  {0.400000 , 0.803922 , 0.666667 }},
{"aquamarine4",  {0.270588 , 0.545098 , 0.454902 }},
{"azure",  {0.941176 , 1.000000 , 1.000000 }},
{"azure1",  {0.941176 , 1.000000 , 1.000000 }},
{"azure2",  {0.878431 , 0.933333 , 0.933333 }},
{"azure3",  {0.756863 , 0.803922 , 0.803922 }},
{"azure4",  {0.513725 , 0.545098 , 0.545098 }},
{"beige",  {0.960784 , 0.960784 , 0.862745 }},
{"bisque",  {1.000000 , 0.894118 , 0.768627 }},
{"bisque1",  {1.000000 , 0.894118 , 0.768627 }},
{"bisque2",  {0.933333 , 0.835294 , 0.717647 }},
{"bisque3",  {0.803922 , 0.717647 , 0.619608 }},
{"bisque4",  {0.545098 , 0.490196 , 0.419608 }},
{"black",  {0.000000 , 0.000000 , 0.000000 }},
{"blanched_almond",  {1.000000 , 0.921569 , 0.803922 }},
{"blue",  {0.000000 , 0.000000 , 1.000000 }},
{"blue1",  {0.000000 , 0.000000 , 1.000000 }},
{"blue2",  {0.000000 , 0.000000 , 0.933333 }},
{"blue3",  {0.000000 , 0.000000 , 0.803922 }},
{"blue4",  {0.000000 , 0.000000 , 0.545098 }},
{"blue_violet",  {0.541176 , 0.168627 , 0.886275 }},
{"brown",  {0.647059 , 0.164706 , 0.164706 }},
{"brown1",  {1.000000 , 0.250980 , 0.250980 }},
{"brown2",  {0.933333 , 0.231373 , 0.231373 }},
{"brown3",  {0.803922 , 0.200000 , 0.200000 }},
{"brown4",  {0.545098 , 0.137255 , 0.137255 }},
{"burlywood",  {0.870588 , 0.721569 , 0.529412 }},
{"burlywood1",  {1.000000 , 0.827451 , 0.607843 }},
{"burlywood2",  {0.933333 , 0.772549 , 0.568627 }},
{"burlywood3",  {0.803922 , 0.666667 , 0.490196 }},
{"burlywood4",  {0.545098 , 0.450980 , 0.333333 }},
{"cadet_blue",  {0.372549 , 0.619608 , 0.627451 }},
{"chartreuse",  {0.498039 , 1.000000 , 0.000000 }},
{"chartreuse1",  {0.498039 , 1.000000 , 0.000000 }},
{"chartreuse2",  {0.462745 , 0.933333 , 0.000000 }},
{"chartreuse3",  {0.400000 , 0.803922 , 0.000000 }},
{"chartreuse4",  {0.270588 , 0.545098 , 0.000000 }},
{"chocolate",  {0.823529 , 0.411765 , 0.117647 }},
{"chocolate1",  {1.000000 , 0.498039 , 0.141176 }},
{"chocolate2",  {0.933333 , 0.462745 , 0.129412 }},
{"chocolate3",  {0.803922 , 0.400000 , 0.113725 }},
{"chocolate4",  {0.545098 , 0.270588 , 0.074510 }},
{"coral",  {1.000000 , 0.498039 , 0.313725 }},
{"coral1",  {1.000000 , 0.447059 , 0.337255 }},
{"coral2",  {0.933333 , 0.415686 , 0.313725 }},
{"coral3",  {0.803922 , 0.356863 , 0.270588 }},
{"coral4",  {0.545098 , 0.243137 , 0.184314 }},
{"cornflower_blue",  {0.392157 , 0.584314 , 0.929412 }},
{"cornsilk",  {1.000000 , 0.972549 , 0.862745 }},
{"cornsilk1",  {1.000000 , 0.972549 , 0.862745 }},
{"cornsilk2",  {0.933333 , 0.909804 , 0.803922 }},
{"cornsilk3",  {0.803922 , 0.784314 , 0.694118 }},
{"cornsilk4",  {0.545098 , 0.533333 , 0.470588 }},
{"cyan",  {0.000000 , 1.000000 , 1.000000 }},
{"cyan1",  {0.000000 , 1.000000 , 1.000000 }},
{"cyan2",  {0.000000 , 0.933333 , 0.933333 }},
{"cyan3",  {0.000000 , 0.803922 , 0.803922 }},
{"cyan4",  {0.000000 , 0.545098 , 0.545098 }},
{"dark_goldenrod",  {0.721569 , 0.525490 , 0.043137 }},
{"dark_green",  {0.000000 , 0.392157 , 0.000000 }},
{"dark_khaki",  {0.741176 , 0.717647 , 0.419608 }},
{"dark_olive_green",  {0.333333 , 0.419608 , 0.184314 }},
{"dark_orange",  {1.000000 , 0.549020 , 0.000000 }},
{"dark_orchid",  {0.600000 , 0.196078 , 0.800000 }},
{"dark_salmon",  {0.913725 , 0.588235 , 0.478431 }},
{"dark_sea_green",  {0.560784 , 0.737255 , 0.560784 }},
{"dark_slate_blue",  {0.282353 , 0.239216 , 0.545098 }},
{"dark_slate_gray",  {0.184314 , 0.309804 , 0.309804 }},
{"dark_slate_grey",  {0.184314 , 0.309804 , 0.309804 }},
{"dark_turquoise",  {0.000000 , 0.807843 , 0.819608 }},
{"dark_violet",  {0.580392 , 0.000000 , 0.827451 }},
{"deep_pink",  {1.000000 , 0.078431 , 0.576471 }},
{"deep_sky_blue",  {0.000000 , 0.749020 , 1.000000 }},
{"dim_gray",  {0.411765 , 0.411765 , 0.411765 }},
{"dim_grey",  {0.411765 , 0.411765 , 0.411765 }},
{"dodger_blue",  {0.117647 , 0.564706 , 1.000000 }},
{"firebrick",  {0.698039 , 0.133333 , 0.133333 }},
{"firebrick1",  {1.000000 , 0.188235 , 0.188235 }},
{"firebrick2",  {0.933333 , 0.172549 , 0.172549 }},
{"firebrick3",  {0.803922 , 0.149020 , 0.149020 }},
{"firebrick4",  {0.545098 , 0.101961 , 0.101961 }},
{"floral_white",  {1.000000 , 0.980392 , 0.941176 }},
{"forest_green",  {0.133333 , 0.545098 , 0.133333 }},
{"gainsboro",  {0.862745 , 0.862745 , 0.862745 }},
{"ghost_white",  {0.972549 , 0.972549 , 1.000000 }},
{"gold",  {1.000000 , 0.843137 , 0.000000 }},
{"gold1",  {1.000000 , 0.843137 , 0.000000 }},
{"gold2",  {0.933333 , 0.788235 , 0.000000 }},
{"gold3",  {0.803922 , 0.678431 , 0.000000 }},
{"gold4",  {0.545098 , 0.458824 , 0.000000 }},
{"goldenrod",  {0.854902 , 0.647059 , 0.125490 }},
{"goldenrod1",  {1.000000 , 0.756863 , 0.145098 }},
{"goldenrod2",  {0.933333 , 0.705882 , 0.133333 }},
{"goldenrod3",  {0.803922 , 0.607843 , 0.113725 }},
{"goldenrod4",  {0.545098 , 0.411765 , 0.078431 }},
{"gray",  {0.752941 , 0.752941 , 0.752941 }},
{"gray0",  {0.000000 , 0.000000 , 0.000000 }},
{"gray1",  {0.011765 , 0.011765 , 0.011765 }},
{"gray10",  {0.101961 , 0.101961 , 0.101961 }},
{"gray100",  {1.000000 , 1.000000 , 1.000000 }},
{"gray11",  {0.109804 , 0.109804 , 0.109804 }},
{"gray12",  {0.121569 , 0.121569 , 0.121569 }},
{"gray13",  {0.129412 , 0.129412 , 0.129412 }},
{"gray14",  {0.141176 , 0.141176 , 0.141176 }},
{"gray15",  {0.149020 , 0.149020 , 0.149020 }},
{"gray16",  {0.160784 , 0.160784 , 0.160784 }},
{"gray17",  {0.168627 , 0.168627 , 0.168627 }},
{"gray18",  {0.180392 , 0.180392 , 0.180392 }},
{"gray19",  {0.188235 , 0.188235 , 0.188235 }},
{"gray2",  {0.019608 , 0.019608 , 0.019608 }},
{"gray20",  {0.200000 , 0.200000 , 0.200000 }},
{"gray21",  {0.211765 , 0.211765 , 0.211765 }},
{"gray22",  {0.219608 , 0.219608 , 0.219608 }},
{"gray23",  {0.231373 , 0.231373 , 0.231373 }},
{"gray24",  {0.239216 , 0.239216 , 0.239216 }},
{"gray25",  {0.250980 , 0.250980 , 0.250980 }},
{"gray26",  {0.258824 , 0.258824 , 0.258824 }},
{"gray27",  {0.270588 , 0.270588 , 0.270588 }},
{"gray28",  {0.278431 , 0.278431 , 0.278431 }},
{"gray29",  {0.290196 , 0.290196 , 0.290196 }},
{"gray3",  {0.031373 , 0.031373 , 0.031373 }},
{"gray30",  {0.301961 , 0.301961 , 0.301961 }},
{"gray31",  {0.309804 , 0.309804 , 0.309804 }},
{"gray32",  {0.321569 , 0.321569 , 0.321569 }},
{"gray33",  {0.329412 , 0.329412 , 0.329412 }},
{"gray34",  {0.341176 , 0.341176 , 0.341176 }},
{"gray35",  {0.349020 , 0.349020 , 0.349020 }},
{"gray36",  {0.360784 , 0.360784 , 0.360784 }},
{"gray37",  {0.368627 , 0.368627 , 0.368627 }},
{"gray38",  {0.380392 , 0.380392 , 0.380392 }},
{"gray39",  {0.388235 , 0.388235 , 0.388235 }},
{"gray4",  {0.039216 , 0.039216 , 0.039216 }},
{"gray40",  {0.400000 , 0.400000 , 0.400000 }},
{"gray41",  {0.411765 , 0.411765 , 0.411765 }},
{"gray42",  {0.419608 , 0.419608 , 0.419608 }},
{"gray43",  {0.431373 , 0.431373 , 0.431373 }},
{"gray44",  {0.439216 , 0.439216 , 0.439216 }},
{"gray45",  {0.450980 , 0.450980 , 0.450980 }},
{"gray46",  {0.458824 , 0.458824 , 0.458824 }},
{"gray47",  {0.470588 , 0.470588 , 0.470588 }},
{"gray48",  {0.478431 , 0.478431 , 0.478431 }},
{"gray49",  {0.490196 , 0.490196 , 0.490196 }},
{"gray5",  {0.050980 , 0.050980 , 0.050980 }},
{"gray50",  {0.498039 , 0.498039 , 0.498039 }},
{"gray51",  {0.509804 , 0.509804 , 0.509804 }},
{"gray52",  {0.521569 , 0.521569 , 0.521569 }},
{"gray53",  {0.529412 , 0.529412 , 0.529412 }},
{"gray54",  {0.541176 , 0.541176 , 0.541176 }},
{"gray55",  {0.549020 , 0.549020 , 0.549020 }},
{"gray56",  {0.560784 , 0.560784 , 0.560784 }},
{"gray57",  {0.568627 , 0.568627 , 0.568627 }},
{"gray58",  {0.580392 , 0.580392 , 0.580392 }},
{"gray59",  {0.588235 , 0.588235 , 0.588235 }},
{"gray6",  {0.058824 , 0.058824 , 0.058824 }},
{"gray60",  {0.600000 , 0.600000 , 0.600000 }},
{"gray61",  {0.611765 , 0.611765 , 0.611765 }},
{"gray62",  {0.619608 , 0.619608 , 0.619608 }},
{"gray63",  {0.631373 , 0.631373 , 0.631373 }},
{"gray64",  {0.639216 , 0.639216 , 0.639216 }},
{"gray65",  {0.650980 , 0.650980 , 0.650980 }},
{"gray66",  {0.658824 , 0.658824 , 0.658824 }},
{"gray67",  {0.670588 , 0.670588 , 0.670588 }},
{"gray68",  {0.678431 , 0.678431 , 0.678431 }},
{"gray69",  {0.690196 , 0.690196 , 0.690196 }},
{"gray7",  {0.070588 , 0.070588 , 0.070588 }},
{"gray70",  {0.701961 , 0.701961 , 0.701961 }},
{"gray71",  {0.709804 , 0.709804 , 0.709804 }},
{"gray72",  {0.721569 , 0.721569 , 0.721569 }},
{"gray73",  {0.729412 , 0.729412 , 0.729412 }},
{"gray74",  {0.741176 , 0.741176 , 0.741176 }},
{"gray75",  {0.749020 , 0.749020 , 0.749020 }},
{"gray76",  {0.760784 , 0.760784 , 0.760784 }},
{"gray77",  {0.768627 , 0.768627 , 0.768627 }},
{"gray78",  {0.780392 , 0.780392 , 0.780392 }},
{"gray79",  {0.788235 , 0.788235 , 0.788235 }},
{"gray8",  {0.078431 , 0.078431 , 0.078431 }},
{"gray80",  {0.800000 , 0.800000 , 0.800000 }},
{"gray81",  {0.811765 , 0.811765 , 0.811765 }},
{"gray82",  {0.819608 , 0.819608 , 0.819608 }},
{"gray83",  {0.831373 , 0.831373 , 0.831373 }},
{"gray84",  {0.839216 , 0.839216 , 0.839216 }},
{"gray85",  {0.850980 , 0.850980 , 0.850980 }},
{"gray86",  {0.858824 , 0.858824 , 0.858824 }},
{"gray87",  {0.870588 , 0.870588 , 0.870588 }},
{"gray88",  {0.878431 , 0.878431 , 0.878431 }},
{"gray89",  {0.890196 , 0.890196 , 0.890196 }},
{"gray9",  {0.090196 , 0.090196 , 0.090196 }},
{"gray90",  {0.898039 , 0.898039 , 0.898039 }},
{"gray91",  {0.909804 , 0.909804 , 0.909804 }},
{"gray92",  {0.921569 , 0.921569 , 0.921569 }},
{"gray93",  {0.929412 , 0.929412 , 0.929412 }},
{"gray94",  {0.941176 , 0.941176 , 0.941176 }},
{"gray95",  {0.949020 , 0.949020 , 0.949020 }},
{"gray96",  {0.960784 , 0.960784 , 0.960784 }},
{"gray97",  {0.968627 , 0.968627 , 0.968627 }},
{"gray98",  {0.980392 , 0.980392 , 0.980392 }},
{"gray99",  {0.988235 , 0.988235 , 0.988235 }},
{"green",  {0.000000 , 1.000000 , 0.000000 }},
{"green1",  {0.000000 , 1.000000 , 0.000000 }},
{"green2",  {0.000000 , 0.933333 , 0.000000 }},
{"green3",  {0.000000 , 0.803922 , 0.000000 }},
{"green4",  {0.000000 , 0.545098 , 0.000000 }},
{"green_yellow",  {0.678431 , 1.000000 , 0.184314 }},
{"grey",  {0.752941 , 0.752941 , 0.752941 }},
{"grey0",  {0.000000 , 0.000000 , 0.000000 }},
{"grey1",  {0.011765 , 0.011765 , 0.011765 }},
{"grey10",  {0.101961 , 0.101961 , 0.101961 }},
{"grey100",  {1.000000 , 1.000000 , 1.000000 }},
{"grey11",  {0.109804 , 0.109804 , 0.109804 }},
{"grey12",  {0.121569 , 0.121569 , 0.121569 }},
{"grey13",  {0.129412 , 0.129412 , 0.129412 }},
{"grey14",  {0.141176 , 0.141176 , 0.141176 }},
{"grey15",  {0.149020 , 0.149020 , 0.149020 }},
{"grey16",  {0.160784 , 0.160784 , 0.160784 }},
{"grey17",  {0.168627 , 0.168627 , 0.168627 }},
{"grey18",  {0.180392 , 0.180392 , 0.180392 }},
{"grey19",  {0.188235 , 0.188235 , 0.188235 }},
{"grey2",  {0.019608 , 0.019608 , 0.019608 }},
{"grey20",  {0.200000 , 0.200000 , 0.200000 }},
{"grey21",  {0.211765 , 0.211765 , 0.211765 }},
{"grey22",  {0.219608 , 0.219608 , 0.219608 }},
{"grey23",  {0.231373 , 0.231373 , 0.231373 }},
{"grey24",  {0.239216 , 0.239216 , 0.239216 }},
{"grey25",  {0.250980 , 0.250980 , 0.250980 }},
{"grey26",  {0.258824 , 0.258824 , 0.258824 }},
{"grey27",  {0.270588 , 0.270588 , 0.270588 }},
{"grey28",  {0.278431 , 0.278431 , 0.278431 }},
{"grey29",  {0.290196 , 0.290196 , 0.290196 }},
{"grey3",  {0.031373 , 0.031373 , 0.031373 }},
{"grey30",  {0.301961 , 0.301961 , 0.301961 }},
{"grey31",  {0.309804 , 0.309804 , 0.309804 }},
{"grey32",  {0.321569 , 0.321569 , 0.321569 }},
{"grey33",  {0.329412 , 0.329412 , 0.329412 }},
{"grey34",  {0.341176 , 0.341176 , 0.341176 }},
{"grey35",  {0.349020 , 0.349020 , 0.349020 }},
{"grey36",  {0.360784 , 0.360784 , 0.360784 }},
{"grey37",  {0.368627 , 0.368627 , 0.368627 }},
{"grey38",  {0.380392 , 0.380392 , 0.380392 }},
{"grey39",  {0.388235 , 0.388235 , 0.388235 }},
{"grey4",  {0.039216 , 0.039216 , 0.039216 }},
{"grey40",  {0.400000 , 0.400000 , 0.400000 }},
{"grey41",  {0.411765 , 0.411765 , 0.411765 }},
{"grey42",  {0.419608 , 0.419608 , 0.419608 }},
{"grey43",  {0.431373 , 0.431373 , 0.431373 }},
{"grey44",  {0.439216 , 0.439216 , 0.439216 }},
{"grey45",  {0.450980 , 0.450980 , 0.450980 }},
{"grey46",  {0.458824 , 0.458824 , 0.458824 }},
{"grey47",  {0.470588 , 0.470588 , 0.470588 }},
{"grey48",  {0.478431 , 0.478431 , 0.478431 }},
{"grey49",  {0.490196 , 0.490196 , 0.490196 }},
{"grey5",  {0.050980 , 0.050980 , 0.050980 }},
{"grey50",  {0.498039 , 0.498039 , 0.498039 }},
{"grey51",  {0.509804 , 0.509804 , 0.509804 }},
{"grey52",  {0.521569 , 0.521569 , 0.521569 }},
{"grey53",  {0.529412 , 0.529412 , 0.529412 }},
{"grey54",  {0.541176 , 0.541176 , 0.541176 }},
{"grey55",  {0.549020 , 0.549020 , 0.549020 }},
{"grey56",  {0.560784 , 0.560784 , 0.560784 }},
{"grey57",  {0.568627 , 0.568627 , 0.568627 }},
{"grey58",  {0.580392 , 0.580392 , 0.580392 }},
{"grey59",  {0.588235 , 0.588235 , 0.588235 }},
{"grey6",  {0.058824 , 0.058824 , 0.058824 }},
{"grey60",  {0.600000 , 0.600000 , 0.600000 }},
{"grey61",  {0.611765 , 0.611765 , 0.611765 }},
{"grey62",  {0.619608 , 0.619608 , 0.619608 }},
{"grey63",  {0.631373 , 0.631373 , 0.631373 }},
{"grey64",  {0.639216 , 0.639216 , 0.639216 }},
{"grey65",  {0.650980 , 0.650980 , 0.650980 }},
{"grey66",  {0.658824 , 0.658824 , 0.658824 }},
{"grey67",  {0.670588 , 0.670588 , 0.670588 }},
{"grey68",  {0.678431 , 0.678431 , 0.678431 }},
{"grey69",  {0.690196 , 0.690196 , 0.690196 }},
{"grey7",  {0.070588 , 0.070588 , 0.070588 }},
{"grey70",  {0.701961 , 0.701961 , 0.701961 }},
{"grey71",  {0.709804 , 0.709804 , 0.709804 }},
{"grey72",  {0.721569 , 0.721569 , 0.721569 }},
{"grey73",  {0.729412 , 0.729412 , 0.729412 }},
{"grey74",  {0.741176 , 0.741176 , 0.741176 }},
{"grey75",  {0.749020 , 0.749020 , 0.749020 }},
{"grey76",  {0.760784 , 0.760784 , 0.760784 }},
{"grey77",  {0.768627 , 0.768627 , 0.768627 }},
{"grey78",  {0.780392 , 0.780392 , 0.780392 }},
{"grey79",  {0.788235 , 0.788235 , 0.788235 }},
{"grey8",  {0.078431 , 0.078431 , 0.078431 }},
{"grey80",  {0.800000 , 0.800000 , 0.800000 }},
{"grey81",  {0.811765 , 0.811765 , 0.811765 }},
{"grey82",  {0.819608 , 0.819608 , 0.819608 }},
{"grey83",  {0.831373 , 0.831373 , 0.831373 }},
{"grey84",  {0.839216 , 0.839216 , 0.839216 }},
{"grey85",  {0.850980 , 0.850980 , 0.850980 }},
{"grey86",  {0.858824 , 0.858824 , 0.858824 }},
{"grey87",  {0.870588 , 0.870588 , 0.870588 }},
{"grey88",  {0.878431 , 0.878431 , 0.878431 }},
{"grey89",  {0.890196 , 0.890196 , 0.890196 }},
{"grey9",  {0.090196 , 0.090196 , 0.090196 }},
{"grey90",  {0.898039 , 0.898039 , 0.898039 }},
{"grey91",  {0.909804 , 0.909804 , 0.909804 }},
{"grey92",  {0.921569 , 0.921569 , 0.921569 }},
{"grey93",  {0.929412 , 0.929412 , 0.929412 }},
{"grey94",  {0.941176 , 0.941176 , 0.941176 }},
{"grey95",  {0.949020 , 0.949020 , 0.949020 }},
{"grey96",  {0.960784 , 0.960784 , 0.960784 }},
{"grey97",  {0.968627 , 0.968627 , 0.968627 }},
{"grey98",  {0.980392 , 0.980392 , 0.980392 }},
{"grey99",  {0.988235 , 0.988235 , 0.988235 }},
{"honeydew",  {0.941176 , 1.000000 , 0.941176 }},
{"honeydew1",  {0.941176 , 1.000000 , 0.941176 }},
{"honeydew2",  {0.878431 , 0.933333 , 0.878431 }},
{"honeydew3",  {0.756863 , 0.803922 , 0.756863 }},
{"honeydew4",  {0.513725 , 0.545098 , 0.513725 }},
{"hot_pink",  {1.000000 , 0.411765 , 0.705882 }},
{"indian_red",  {0.803922 , 0.360784 , 0.360784 }},
{"ivory",  {1.000000 , 1.000000 , 0.941176 }},
{"ivory1",  {1.000000 , 1.000000 , 0.941176 }},
{"ivory2",  {0.933333 , 0.933333 , 0.878431 }},
{"ivory3",  {0.803922 , 0.803922 , 0.756863 }},
{"ivory4",  {0.545098 , 0.545098 , 0.513725 }},
{"khaki",  {0.941176 , 0.901961 , 0.549020 }},
{"khaki1",  {1.000000 , 0.964706 , 0.560784 }},
{"khaki2",  {0.933333 , 0.901961 , 0.521569 }},
{"khaki3",  {0.803922 , 0.776471 , 0.450980 }},
{"khaki4",  {0.545098 , 0.525490 , 0.305882 }},
{"lavender",  {0.901961 , 0.901961 , 0.980392 }},
{"lavender_blush",  {1.000000 , 0.941176 , 0.960784 }},
{"lawn_green",  {0.486275 , 0.988235 , 0.000000 }},
{"lemon_chiffon",  {1.000000 , 0.980392 , 0.803922 }},
{"light_blue",  {0.678431 , 0.847059 , 0.901961 }},
{"light_coral",  {0.941176 , 0.501961 , 0.501961 }},
{"light_cyan",  {0.878431 , 1.000000 , 1.000000 }},
{"light_goldenrod",  {0.933333 , 0.866667 , 0.509804 }},
{"light_goldenrod_yellow",  {0.980392 , 0.980392 , 0.823529 }},
{"light_gray",  {0.827451 , 0.827451 , 0.827451 }},
{"light_grey",  {0.827451 , 0.827451 , 0.827451 }},
{"light_pink",  {1.000000 , 0.713725 , 0.756863 }},
{"light_salmon",  {1.000000 , 0.627451 , 0.478431 }},
{"light_sea_green",  {0.125490 , 0.698039 , 0.666667 }},
{"light_sky_blue",  {0.529412 , 0.807843 , 0.980392 }},
{"light_slate_blue",  {0.517647 , 0.439216 , 1.000000 }},
{"light_slate_gray",  {0.466667 , 0.533333 , 0.600000 }},
{"light_slate_grey",  {0.466667 , 0.533333 , 0.600000 }},
{"light_steel_blue",  {0.690196 , 0.768627 , 0.870588 }},
{"light_yellow",  {1.000000 , 1.000000 , 0.878431 }},
{"lime_green",  {0.196078 , 0.803922 , 0.196078 }},
{"linen",  {0.980392 , 0.941176 , 0.901961 }},
{"magenta",  {1.000000 , 0.000000 , 1.000000 }},
{"magenta1",  {1.000000 , 0.000000 , 1.000000 }},
{"magenta2",  {0.933333 , 0.000000 , 0.933333 }},
{"magenta3",  {0.803922 , 0.000000 , 0.803922 }},
{"magenta4",  {0.545098 , 0.000000 , 0.545098 }},
{"maroon",  {0.690196 , 0.188235 , 0.376471 }},
{"maroon1",  {1.000000 , 0.203922 , 0.701961 }},
{"maroon2",  {0.933333 , 0.188235 , 0.654902 }},
{"maroon3",  {0.803922 , 0.160784 , 0.564706 }},
{"maroon4",  {0.545098 , 0.109804 , 0.384314 }},
{"medium_aquamarine",  {0.400000 , 0.803922 , 0.666667 }},
{"medium_blue",  {0.000000 , 0.000000 , 0.803922 }},
{"medium_orchid",  {0.729412 , 0.333333 , 0.827451 }},
{"medium_purple",  {0.576471 , 0.439216 , 0.858824 }},
{"medium_sea_green",  {0.235294 , 0.701961 , 0.443137 }},
{"medium_slate_blue",  {0.482353 , 0.407843 , 0.933333 }},
{"medium_spring_green",  {0.000000 , 0.980392 , 0.603922 }},
{"medium_turquoise",  {0.282353 , 0.819608 , 0.800000 }},
{"medium_violet_red",  {0.780392 , 0.082353 , 0.521569 }},
{"midnight_blue",  {0.098039 , 0.098039 , 0.439216 }},
{"mint_cream",  {0.960784 , 1.000000 , 0.980392 }},
{"misty_rose",  {1.000000 , 0.894118 , 0.882353 }},
{"moccasin",  {1.000000 , 0.894118 , 0.709804 }},
{"navajo_white",  {1.000000 , 0.870588 , 0.678431 }},
{"navy",  {0.000000 , 0.000000 , 0.501961 }},
{"navy_blue",  {0.000000 , 0.000000 , 0.501961 }},
{"old_lace",  {0.992157 , 0.960784 , 0.901961 }},
{"olive_drab",  {0.419608 , 0.556863 , 0.137255 }},
{"orange",  {1.000000 , 0.647059 , 0.000000 }},
{"orange1",  {1.000000 , 0.647059 , 0.000000 }},
{"orange2",  {0.933333 , 0.603922 , 0.000000 }},
{"orange3",  {0.803922 , 0.521569 , 0.000000 }},
{"orange4",  {0.545098 , 0.352941 , 0.000000 }},
{"orange_red",  {1.000000 , 0.270588 , 0.000000 }},
{"orchid",  {0.854902 , 0.439216 , 0.839216 }},
{"orchid1",  {1.000000 , 0.513725 , 0.980392 }},
{"orchid2",  {0.933333 , 0.478431 , 0.913725 }},
{"orchid3",  {0.803922 , 0.411765 , 0.788235 }},
{"orchid4",  {0.545098 , 0.278431 , 0.537255 }},
{"pale_goldenrod",  {0.933333 , 0.909804 , 0.666667 }},
{"pale_green",  {0.596078 , 0.984314 , 0.596078 }},
{"pale_turquoise",  {0.686275 , 0.933333 , 0.933333 }},
{"pale_violet_red",  {0.858824 , 0.439216 , 0.576471 }},
{"papaya_whip",  {1.000000 , 0.937255 , 0.835294 }},
{"peach_puff",  {1.000000 , 0.854902 , 0.725490 }},
{"peru",  {0.803922 , 0.521569 , 0.247059 }},
{"pink",  {1.000000 , 0.752941 , 0.796078 }},
{"pink1",  {1.000000 , 0.709804 , 0.772549 }},
{"pink2",  {0.933333 , 0.662745 , 0.721569 }},
{"pink3",  {0.803922 , 0.568627 , 0.619608 }},
{"pink4",  {0.545098 , 0.388235 , 0.423529 }},
{"plum",  {0.866667 , 0.627451 , 0.866667 }},
{"plum1",  {1.000000 , 0.733333 , 1.000000 }},
{"plum2",  {0.933333 , 0.682353 , 0.933333 }},
{"plum3",  {0.803922 , 0.588235 , 0.803922 }},
{"plum4",  {0.545098 , 0.400000 , 0.545098 }},
{"powder_blue",  {0.690196 , 0.878431 , 0.901961 }},
{"purple",  {0.627451 , 0.125490 , 0.941176 }},
{"purple1",  {0.607843 , 0.188235 , 1.000000 }},
{"purple2",  {0.568627 , 0.172549 , 0.933333 }},
{"purple3",  {0.490196 , 0.149020 , 0.803922 }},
{"purple4",  {0.333333 , 0.101961 , 0.545098 }},
{"red",  {1.000000 , 0.000000 , 0.000000 }},
{"red1",  {1.000000 , 0.000000 , 0.000000 }},
{"red2",  {0.933333 , 0.000000 , 0.000000 }},
{"red3",  {0.803922 , 0.000000 , 0.000000 }},
{"red4",  {0.545098 , 0.000000 , 0.000000 }},
{"rosy_brown",  {0.737255 , 0.560784 , 0.560784 }},
{"royal_blue",  {0.254902 , 0.411765 , 0.882353 }},
{"saddle_brown",  {0.545098 , 0.270588 , 0.074510 }},
{"salmon",  {0.980392 , 0.501961 , 0.447059 }},
{"salmon1",  {1.000000 , 0.549020 , 0.411765 }},
{"salmon2",  {0.933333 , 0.509804 , 0.384314 }},
{"salmon3",  {0.803922 , 0.439216 , 0.329412 }},
{"salmon4",  {0.545098 , 0.298039 , 0.223529 }},
{"sandy_brown",  {0.956863 , 0.643137 , 0.376471 }},
{"sea_green",  {0.180392 , 0.545098 , 0.341176 }},
{"seashell",  {1.000000 , 0.960784 , 0.933333 }},
{"seashell1",  {1.000000 , 0.960784 , 0.933333 }},
{"seashell2",  {0.933333 , 0.898039 , 0.870588 }},
{"seashell3",  {0.803922 , 0.772549 , 0.749020 }},
{"seashell4",  {0.545098 , 0.525490 , 0.509804 }},
{"sienna",  {0.627451 , 0.321569 , 0.176471 }},
{"sienna1",  {1.000000 , 0.509804 , 0.278431 }},
{"sienna2",  {0.933333 , 0.474510 , 0.258824 }},
{"sienna3",  {0.803922 , 0.407843 , 0.223529 }},
{"sienna4",  {0.545098 , 0.278431 , 0.149020 }},
{"sky_blue",  {0.529412 , 0.807843 , 0.921569 }},
{"slate_blue",  {0.415686 , 0.352941 , 0.803922 }},
{"slate_gray",  {0.439216 , 0.501961 , 0.564706 }},
{"slate_grey",  {0.439216 , 0.501961 , 0.564706 }},
{"snow",  {1.000000 , 0.980392 , 0.980392 }},
{"snow1",  {1.000000 , 0.980392 , 0.980392 }},
{"snow2",  {0.933333 , 0.913725 , 0.913725 }},
{"snow3",  {0.803922 , 0.788235 , 0.788235 }},
{"snow4",  {0.545098 , 0.537255 , 0.537255 }},
{"spring_green",  {0.000000 , 1.000000 , 0.498039 }},
{"steel_blue",  {0.274510 , 0.509804 , 0.705882 }},
{"tan",  {0.823529 , 0.705882 , 0.549020 }},
{"tan1",  {1.000000 , 0.647059 , 0.309804 }},
{"tan2",  {0.933333 , 0.603922 , 0.286275 }},
{"tan3",  {0.803922 , 0.521569 , 0.247059 }},
{"tan4",  {0.545098 , 0.352941 , 0.168627 }},
{"thistle",  {0.847059 , 0.749020 , 0.847059 }},
{"thistle1",  {1.000000 , 0.882353 , 1.000000 }},
{"thistle2",  {0.933333 , 0.823529 , 0.933333 }},
{"thistle3",  {0.803922 , 0.709804 , 0.803922 }},
{"thistle4",  {0.545098 , 0.482353 , 0.545098 }},
{"tomato",  {1.000000 , 0.388235 , 0.278431 }},
{"tomato1",  {1.000000 , 0.388235 , 0.278431 }},
{"tomato2",  {0.933333 , 0.360784 , 0.258824 }},
{"tomato3",  {0.803922 , 0.309804 , 0.223529 }},
{"tomato4",  {0.545098 , 0.211765 , 0.149020 }},
{"turquoise",  {0.250980 , 0.878431 , 0.815686 }},
{"turquoise1",  {0.000000 , 0.960784 , 1.000000 }},
{"turquoise2",  {0.000000 , 0.898039 , 0.933333 }},
{"turquoise3",  {0.000000 , 0.772549 , 0.803922 }},
{"turquoise4",  {0.000000 , 0.525490 , 0.545098 }},
{"violet",  {0.933333 , 0.509804 , 0.933333 }},
{"violet_red",  {0.815686 , 0.125490 , 0.564706 }},
{"wheat",  {0.960784 , 0.870588 , 0.701961 }},
{"wheat1",  {1.000000 , 0.905882 , 0.729412 }},
{"wheat2",  {0.933333 , 0.847059 , 0.682353 }},
{"wheat3",  {0.803922 , 0.729412 , 0.588235 }},
{"wheat4",  {0.545098 , 0.494118 , 0.400000 }},
{"white",  {1.000000 , 1.000000 , 1.000000 }},
{"white_smoke",  {0.960784 , 0.960784 , 0.960784 }},
{"yellow",  {1.000000 , 1.000000 , 0.000000 }},
{"yellow1",  {1.000000 , 1.000000 , 0.000000 }},
{"yellow2",  {0.933333 , 0.933333 , 0.000000 }},
{"yellow3",  {0.803922 , 0.803922 , 0.000000 }},
{"yellow4",  {0.545098 , 0.545098 , 0.000000 }},
{"yellow_green",  {0.603922 , 0.803922 , 0.196078} }
} ;

int
LookupColorByName(name, color)
 char * name ;
 Vec color ;
{
	int rc ;
	rc = BinarySearch(name, 0, NCOLORS - 1 , Colors) ;
	if (rc < 0) {
		return(0) ;
	}

	VecCopy(Colors[rc].ce_color, color) ;
	return 1 ;
}


int 
BinarySearch(name, l, h, array)
 char * name ;
 int l, h ;
 ColorEntry array[] ;
{
	int m, rc ;
	if (l > h)
		return(-1) ;
	
	m = (l + h) / 2 ;

	rc = strcmp(name, array[m].ce_name) ;
	if (rc == 0)
		return m ;
	else if (rc < 0)
		return BinarySearch(name, l, m-1, array) ;
	else
		return BinarySearch(name, m + 1, h, array) ;
}
