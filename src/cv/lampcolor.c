#ifndef lint
static const char	RCSid[] = "$Id: lampcolor.c,v 2.6 2003/02/22 02:07:23 greg Exp $";
#endif
/*
 * Program to convert lamp color from table and compute radiance.
 */

#include <stdio.h>

#include <math.h>

#include "color.h"

#define PI	3.14159265358979323846

extern char	*gets(), *strcpy();
extern float	*matchlamp();

				/* lamp parameters */
#define LTYPE		0
#define LUNIT		1
#define LGEOM		2
#define LOUTP		3
#define NPARAMS		4

int	typecheck(), unitcheck(), geomcheck(), outpcheck();

float	*lampcolor;		/* the lamp color (RGB) */
double	unit2meter;		/* conversion from units to meters */
double	area;			/* radiating area for this geometry */
double	lumens;			/* total lamp lumens */

struct {
	char	*name;
	char	value[64];
	int	(*check)();
	char	*help;
} param[NPARAMS] = {
	{ "lamp type", "WHITE", typecheck,
"The lamp type is a string which corresponds to one of the types registered\n\
in the lamp table file.  A value of \"WHITE\" means an uncolored source,\n\
which may be preferable because it results in a color balanced image." },
	{ "length unit", "meter", unitcheck,
"Unit must be one of:  \"meter\", \"centimeter\", \"foot\", or \"inch\".\n\
These may be abbreviated as a single letter." },
	{ "lamp geometry", "polygon", geomcheck,
"The lamp geometry must be one of:  \"polygon\", \"sphere\", \"cylinder\"\n\
or \"ring\".  These may be abbreviated as a single letter." },
	{ "total lamp lumens", "0", outpcheck,
"This is the overall light output of the lamp and its fixture.  If you do\n\
not know this value explicitly, you can compute the approximate lumens\n\
by multiplying the input wattage by 14 for incandescent fixtures or 70\n\
for fluorescent fixtures." },
};


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*lamptab = "lamp.tab";
	char	buf[64];
	int	i;

	if (argc > 1) lamptab = argv[1];
	if (loadlamps(lamptab) == 0) {
		fprintf(stderr, "%s: no such lamp table\n", lamptab);
		exit(1);
	}
	printf("Program to compute lamp radiance.  Enter '?' for help.\n");
	for ( ; ; ) {
		i = 0;
		while (i < NPARAMS) {
			printf("Enter %s [%s]: ", param[i].name,
					param[i].value);
			if (gets(buf) == NULL)
				exit(0);
			if (buf[0] == '?') {
				puts(param[i].help);
				continue;
			}
			if (buf[0])
				strcpy(param[i].value, buf);
			if (!(*param[i].check)(param[i].value)) {
				fprintf(stderr, "%s: bad value for %s\n",
						param[i].value, param[i].name);
				continue;
			}
			i++;
		}
		compute();
	}
}


typecheck(s)			/* check lamp type */
char	*s;
{
	lampcolor = matchlamp(s);
	return(lampcolor != NULL);
}


unitcheck(s)			/* compute conversion to meters */
char	*s;
{
	int	len = strlen(s);

	switch (*s) {
	case 'm':
		if (strncmp(s, "meters", len))
			return(0);
		unit2meter = 1.0;
		return(1);
	case 'c':
		if (strncmp(s, "centimeters", len) && strncmp(s, "cms", len))
			return(0);
		unit2meter = 0.01;
		return(1);
	case 'f':
		if (strncmp(s, "foot", len) && strncmp(s, "feet", len))
			return(0);
		unit2meter = 0.3048;
		return(1);
	case 'i':
		if (strncmp(s, "inches", len))
			return(0);
		unit2meter = 0.0254;
		return(1);
	}
	return(0);
}


geomcheck(s)			/* check/set lamp geometry */
char	*s;
{
	int	len = strlen(s);

	switch (*s) {
	case 'p':
		if (strncmp(s, "polygonal", len))
			return(0);
		return(getpolygon());
	case 's':
		if (strncmp(s, "sphere", len) && strncmp(s, "spherical", len))
			return(0);
		return(getsphere());
	case 'c':
		if (strncmp(s,"cylinder",len) && strncmp(s,"cylindrical",len))
			return(0);
		return(getcylinder());
	case 'r':
		if (strncmp(s, "ring", len) && strncmp(s, "disk", len))
			return(0);
		return(getring());
	}
	return(0);
}


outpcheck(s)			/* check lumen output value */
register char	*s;
{
	if ((*s < '0' || *s > '9') && *s != '.')
		return(0);
	lumens = atof(s);
	return(1);
}


compute()			/* compute lamp radiance */
{
	double	whiteval;

	whiteval = lumens/area/(WHTEFFICACY*PI);

	printf("Lamp color (RGB) = %f %f %f\n",
			lampcolor[0]*whiteval,
			lampcolor[1]*whiteval,
			lampcolor[2]*whiteval);
}


getd(name, dp, help)		/* get a positive double from stdin */
char	*name;
double	*dp;
char	*help;
{
	char	buf[32];
again:
	printf("%s [%g]: ", name, *dp);
	if (gets(buf) == NULL)
		return(0);
	if (buf[0] == '?') {
		puts(help);
		goto again;
	}
	if ((buf[0] < '0' || buf[0] > '9') && buf[0] != '.')
		return(0);
	*dp = atof(buf);
	return(1);
}


getpolygon()			/* get projected area for a polygon */
{
	static double	parea = 1.0;

	getd("Polygon area", &parea,
		"Enter the total radiating area of the polygon.");
	area = unit2meter*unit2meter * parea;
	return(1);
}


getsphere()			/* get projected area for a sphere */
{
	static double	radius = 1.0;

	getd("Sphere radius", &radius,
		"Enter the distance from the sphere's center to its surface.");
	area = 4.*PI*unit2meter*unit2meter * radius*radius;
	return(1);
}


getcylinder()			/* get projected area for a cylinder */
{
	static double	length = 1.0, radius = 0.1;

	getd("Cylinder length", &length,
		"Enter the length of the cylinder.");
	getd("Cylinder radius", &radius,
		"Enter the distance from the cylinder's axis to its surface.");
	area = 2.*PI*unit2meter*unit2meter * radius*length;
	return(1);
}


getring()			/* get projected area for a ring */
{
	static double	radius = 1.0;

	getd("Disk radius", &radius,
"Enter the distance from the ring's center to its outer edge.\n\
The inner radius must be zero.");
	area = PI*unit2meter*unit2meter * radius*radius;
	return(1);
}
