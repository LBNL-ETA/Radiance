#ifndef lint
static const char	RCSid[] = "$Id: genbeads.c,v 2.7 2003/11/16 10:29:38 schorsch Exp $";
#endif
/*
 *  genbeads.c - generate a string of spheres using Hermite
 *		curve specification.
 *
 *     10/29/85
 */

#include  <stdio.h>
#include <stdlib.h>
#include  <math.h>

void hermite3(double hp[3], double p0[3], double p1[3],
		double r0[3], double r1[3], double t);
void htan3(double ht[3], double p0[3], double p1[3],
		double r0[3], double r1[3], double t);

static void genstring(
		char  *mtype,
		char  *name,
		double  p0[3],
		double  p1[3],
		double  r0[3],
		double  r1[3],
		double  rad,
		double  inc
		)
{
	register int  i;
	double  v[3];
	double  t;
	
	t = 0.0;
	for (i = 0; t <= 1.0; i++) {

		printf("\n%s sphere %s.%d\n", mtype, name, i);
		
		hermite3(v, p0, p1, r0, r1, t);
		
		printf("0\n0\n4 %18.12g %18.12g %18.12g %18.12g\n",
					v[0], v[1], v[2], rad);

		htan3(v, p0, p1, r0, r1, t);
		t += inc / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	}
}

int main(int  argc, char  **argv)
{
	double  p0[3], p1[3], r0[3], r1[3];
	double  rad, inc;
	char  *mtype;		/* material type */
	char  *name;		/* name */
	
	if (argc != 17) {
		fprintf(stderr, "Usage: %s material name p0 p1 r0 r1 rad inc\n",
					 argv[0]);
		exit(1);
	}
	mtype = argv[1];
	name = argv[2];
	p0[0] = atof(argv[3]);
	p0[1] = atof(argv[4]);
	p0[2] = atof(argv[5]);
	p1[0] = atof(argv[6]);
	p1[1] = atof(argv[7]);
	p1[2] = atof(argv[8]);
	r0[0] = atof(argv[9]);
	r0[1] = atof(argv[10]);
	r0[2] = atof(argv[11]);
	r1[0] = atof(argv[12]);
	r1[1] = atof(argv[13]);
	r1[2] = atof(argv[14]);
	rad = atof(argv[15]);
	inc = atof(argv[16]);
	
	genstring(mtype, name, p0, p1, r0, r1, rad, inc);
	
	return(0);
}

