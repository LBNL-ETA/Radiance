#ifndef lint
static const char	RCSid[] = "$Id: genmarble.c,v 2.6 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 *  genmarble.c - generate a marble with bubbles inside.
 *
 *     1/8/86
 */

#include  <stdio.h>

#include <stdlib.h>

#include  <math.h>

#include  "random.h"

#define  PI	3.14159265359

typedef double  FVECT[3];

static double  bubble();	/* pretty cute, huh? */
static void sphere_cart();


main(argc, argv)
int  argc;
char  **argv;
{
	char  *cmtype, *cname;
	FVECT  cent;
	double  rad;
	int  nbubbles, i;
	double  bubrad;
	FVECT  v;
	double  brad;
	
	if (argc != 9) {
		fprintf(stderr,
			"Usage: %s material name cent rad #bubbles bubrad\n",
					 argv[0]);
		exit(1);
	}
	cmtype = argv[1];
	cname = argv[2];
	cent[0] = atof(argv[3]);
	cent[1] = atof(argv[4]);
	cent[2] = atof(argv[5]);
	rad = atof(argv[6]);
	nbubbles = atoi(argv[7]);
	bubrad = atof(argv[8]);

	if (bubrad >= rad) {
		fprintf(stderr, "%s: bubbles too big for marble\n", argv[0]);
		exit(1);
	}

	printf("\n%s sphere %s\n", cmtype, cname);
	printf("0\n0\n4 %f %f %f %f\n", cent[0], cent[1], cent[2], rad);

	for (i = 0; i < nbubbles; i++) {
		brad = bubble(v, cent, rad, bubrad);
		printf("\n%s bubble %s.%d\n", cmtype, cname, i);
		printf("0\n0\n4 %f %f %f %f\n", v[0], v[1], v[2], brad);
	}
	
	return(0);
}


static double
bubble(v, cent, rad, bubrad)	/* compute location of random bubble */
FVECT  v, cent;
double  rad, bubrad;
{
	double  r, ro, theta, phi;

	r = frandom()*bubrad;
	ro = sqrt(frandom())*(rad-r);
	theta = frandom()*(2.0*PI);
	phi = frandom()*PI;
	sphere_cart(v, ro, theta, phi);
	v[0] += cent[0]; v[1] += cent[1]; v[2] += cent[2];
	return(r);
}


static void
sphere_cart(v, ro, theta, phi)	/* spherical to cartesian coord. conversion */
FVECT  v;
double  ro, theta, phi;
{
	double  d;
	
	d = sin(phi);
	v[0] = ro*d*cos(theta);
	v[1] = ro*d*sin(theta);
	v[2] = ro*cos(phi);
}
