#ifndef lint
static const char	RCSid[] = "$Id: genblinds.c,v 2.10 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 *  genblind2.c - make some curved or flat venetian blinds.
 *
 *	Jean-Louis Scartezzini and Greg Ward
 * 
 *  parameters: 
 *   		depth  -  depth of blinds
 *		width  -  width of slats
 *		height -  height of blinds
 *		nslats -  number of slats
 *		angle  -  blind incidence angle ( in degrees )
 *		rcurv  -  curvature radius of slats (up:>0;down:<0;flat:=0)
 */

#include  <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  <string.h>

#define  PI		3.14159265358979323846
#define  DELTA		10.  /*  MINIMAL SUSTAINED ANGLE IN DEGREES */

double  baseflat[4][3], baseblind[4][3][180];
double  A[3],X[3];
char  *material, *name;
double  height;
int  nslats,  nsurf;



makeflat(w,d,a)
double  w, d, a;
{
	double  h;

	h = d*sin(a);
	d *= cos(a);
	baseflat[0][0] = 0.0;
	baseflat[0][1] = 0.0;
	baseflat[0][2] = 0.0;
	baseflat[1][0] = 0.0;
	baseflat[1][1] = w;
	baseflat[1][2] = 0.0;
	baseflat[2][0] = d;
	baseflat[2][1] = w;
	baseflat[2][2] = h;
	baseflat[3][0] = d;
	baseflat[3][1] = 0.0;
	baseflat[3][2] = h;

}


printslat(n)			/* print slat # n */
int  n;
{
	register int  i, k;

	for (k=0; k < nsurf; k++)  {
 		printf("\n%s polygon %s.%d.%d\n", material, name, n, k);
		printf("0\n0\n12\n");
		for (i = 0; i < 4; i++)
			printf("\t%18.12g\t%18.12g\t%18.12g\n",
				baseblind[i][0][k],
				baseblind[i][1][k],
				baseblind[i][2][k] + height*(n-.5)/nslats);
	}		
}


printhead(ac, av)		/* print command header */
register int  ac;
register char  **av;
{
	putchar('#');
	while (ac--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');
}


main(argc, argv)
int  argc;
char  *argv[];
{
	double  width, delem, depth, rcurv = 0.0, angle;
	double  beta, gamma, theta, chi;
	int     i, j, k, l;


	if (argc != 8 && argc != 10)
		goto userr;
	material = argv[1];
	name = argv[2];
	depth = atof(argv[3]);
	width = atof(argv[4]);
	height = atof(argv[5]);
	nslats  = atoi(argv[6]);
	angle = atof(argv[7]);
	if (argc == 10)
		if (!strcmp(argv[8], "-r"))
			rcurv = atof(argv[9]);
		else if (!strcmp(argv[8], "+r"))
			rcurv = -atof(argv[9]);
		else
			goto userr;

/* CURVED BLIND CALCULATION */

	if (rcurv != 0) {

	/* BLINDS SUSTAINED ANGLE */

	theta = 2*asin(depth/(2*fabs(rcurv)));

 	/* HOW MANY ELEMENTARY SURFACES SHOULD BE CALCULATED ? */

	nsurf = (int)(theta / ((PI/180.)*DELTA)) + 1;

	/* WHAT IS THE DEPTH OF THE ELEMENTARY SURFACES ? */

	delem = 2*fabs(rcurv)*sin((PI/180.)*(DELTA/2.));

	beta = (PI-theta)/2.;
	gamma = beta -((PI/180.)*angle);



	if (rcurv < 0) {
		A[0]=fabs(rcurv)*cos(gamma);
		A[0] *= -1;
		A[1]=0.;
		A[2]=fabs(rcurv)*sin(gamma);
	}
	if (rcurv > 0) {
		A[0]=fabs(rcurv)*cos(gamma+theta);
		A[1]=0.;
		A[2]=fabs(rcurv)*sin(gamma+theta);
		A[2] *= -1;
	}

	for (k=0; k < nsurf; k++) {
	if (rcurv < 0) {
        chi=(PI/180.)*((180.-DELTA)/2.) - (gamma+(k*(PI/180.)*DELTA));
	}
	if (rcurv > 0) {
   	chi=(PI-(gamma+theta)+(k*(PI/180.)*DELTA))-(PI/180.)*   
        ((180.-DELTA)/2.);
	}
		makeflat(width, delem, chi);
	if (rcurv < 0.) {
		X[0]=(-fabs(rcurv))*cos(gamma+(k*(PI/180.)*DELTA))-A[0];
		X[1]=0.;
		X[2]=fabs(rcurv)*sin(gamma+(k*(PI/180.)*DELTA))-A[2];
	}
	if (rcurv > 0.) {
		X[0]=fabs(rcurv)*cos(gamma+theta-(k*(PI/180.)*DELTA))-A[0];
		X[1]=0.;
		X[2]=(-fabs(rcurv))*sin(gamma+theta-(k*(PI/180.)*DELTA))-A[2];
	}

	 	for (i=0; i < 4; i++)  {
		    for (j=0; j < 3; j++) {
			baseblind[i][j][k] = baseflat[i][j]+X[j];
		    } 
		}	
 	}
	}

 /* FLAT BLINDS CALCULATION */
	
	if (rcurv == 0.) {

		nsurf=1;
		makeflat(width,depth,angle*(PI/180.));
		for (i=0; i < 4; i++) {
		    for (j=0; j < 3; j++) {
			baseblind[i][j][0] = baseflat[i][j];
		    }
		}
	}
	
	printhead(argc, argv);


/* REPEAT THE BASIC CURVED OR FLAT SLAT TO GET THE OVERALL BLIND */

	for (l = 1; l <= nslats; l++)
		printslat(l);
	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s mat name depth width height nslats angle [-r|+r rcurv]\n",
			argv[0]);
	exit(1);
}



