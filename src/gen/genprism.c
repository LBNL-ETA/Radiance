#ifndef lint
static const char	RCSid[] = "$Id: genprism.c,v 2.12 2003/11/16 10:29:38 schorsch Exp $";
#endif
/*
 *  genprism.c - generate a prism.
 *		2D vertices in the xy plane are given on the
 *		command line or from a file.  Their order together
 *		with the extrude direction will determine surface
 *		orientation.
 */

#include  <stdio.h>
#include  <string.h>
#include <stdlib.h>
#include  <math.h>
#include  <ctype.h>

#define  MAXVERT	1024		/* maximum # vertices */

#define  FTINY		1e-6

char  *pmtype;		/* material type */
char  *pname;		/* name */

double  lvect[3] = {0.0, 0.0, 1.0};
int	lvdir = 1;
double	llen = 1.0;

double  vert[MAXVERT][2];
int  nverts = 0;

double	u[MAXVERT][2];		/* edge unit vectors */
double	a[MAXVERT];		/* corner trim sizes */

int  do_ends = 1;		/* include end caps */
int  iscomplete = 0;		/* polygon is already completed */
double  crad = 0.0;		/* radius for corners (sign from lvdir) */

static double  compute_rounding(void);


static void
readverts(fname)		/* read vertices from a file */
char  *fname;
{
	FILE  *fp;

	if (fname == NULL)
		fp = stdin;
	else if ((fp = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fname);
		exit(1);
	}
	while (fscanf(fp, "%lf %lf", &vert[nverts][0], &vert[nverts][1]) == 2)
		nverts++;
	fclose(fp);
}


static void
side(n0, n1)			/* print single side */
register int  n0, n1;
{
	printf("\n%s polygon %s.%d\n", pmtype, pname, n0+1);
	printf("0\n0\n12\n");
	printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[n0][0],
			vert[n0][1], 0.0);
	printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[n0][0]+lvect[0],
			vert[n0][1]+lvect[1], lvect[2]);
	printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[n1][0]+lvect[0],
			vert[n1][1]+lvect[1], lvect[2]);
	printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[n1][0],
			vert[n1][1], 0.0);
}


static void
rside(n0, n1)			/* print side with rounded edge */
register int  n0, n1;
{
	double  s, c, t[3];

					/* compute tanget offset vector */
	s = lvdir*(lvect[1]*u[n1][0] - lvect[0]*u[n1][1])/llen;
	if (s < -FTINY || s > FTINY) {
		c = sqrt(1. - s*s);
		t[0] = (c - 1.)*u[n1][1];
		t[1] = (1. - c)*u[n1][0];
		t[2] = s;
	} else
		t[0] = t[1] = t[2] = 0.;
					/* output side */
	printf("\n%s polygon %s.%d\n", pmtype, pname, n0+1);
	printf("0\n0\n12\n");
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
			vert[n0][0] + crad*(t[0] - a[n0]*u[n1][0]),
			vert[n0][1] + crad*(t[1] - a[n0]*u[n1][1]),
			crad*(t[2] + 1.));
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
			vert[n0][0] + lvect[0] + crad*(t[0] - a[n0]*u[n1][0]),
			vert[n0][1] + lvect[1] + crad*(t[1] - a[n0]*u[n1][1]),
			lvect[2] + crad*(t[2] - 1.));
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
			vert[n1][0] + lvect[0] + crad*(t[0] + a[n1]*u[n1][0]),
			vert[n1][1] + lvect[1] + crad*(t[1] + a[n1]*u[n1][1]),
			lvect[2] + crad*(t[2] - 1.));
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
			vert[n1][0] + crad*(t[0] + a[n1]*u[n1][0]),
			vert[n1][1] + crad*(t[1] + a[n1]*u[n1][1]),
			crad*(t[2] + 1.));
					/* output joining edge */
	if (lvdir*a[n1] < 0.)
		return;
	printf("\n%s cylinder %s.e%d\n", pmtype, pname, n0+1);
	printf("0\n0\n7\n");
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
		vert[n1][0] + crad*(a[n1]*u[n1][0] - u[n1][1]),
		vert[n1][1] + crad*(a[n1]*u[n1][1] + u[n1][0]),
		crad);
	printf("\t%18.12g\t%18.12g\t%18.12g\n",
		vert[n1][0] + lvect[0] + crad*(a[n1]*u[n1][0] - u[n1][1]),
		vert[n1][1] + lvect[1] + crad*(a[n1]*u[n1][1] + u[n1][0]),
		lvect[2] - crad);
	printf("\t%18.12g\n", lvdir*crad);
}


static double
compute_rounding()		/* compute vectors for rounding operations */
{
	register int  i;
	register double	*v0, *v1;
	double  l, asum;

	v0 = vert[nverts-1];
	for (i = 0; i < nverts; i++) {		/* compute u[*] */
		v1 = vert[i];
		u[i][0] = v0[0] - v1[0];
		u[i][1] = v0[1] - v1[1];
		l = sqrt(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
		if (l <= FTINY) {
			fprintf(stderr, "Degenerate side in prism\n");
			exit(1);
		}
		u[i][0] /= l;
		u[i][1] /= l;
		v0 = v1;
	}
	asum = 0.;
	v1 = u[0];
	for (i = nverts; i--; ) {		/* compute a[*] */
		v0 = u[i];
		l = v0[0]*v1[0] + v0[1]*v1[1];
		if (1+l <= FTINY) {
			fprintf(stderr, "Overlapping sides in prism\n");
			exit(1);
		}
		if (1-l <= 0.)
			a[i] = 0.;
		else {
			a[i] = sqrt((1-l)/(1+l));
			asum += l = v1[0]*v0[1]-v1[1]*v0[0];
			if (l < 0.)
				a[i] = -a[i];
		}
		v1 = v0;
	}
	return(asum*.5);
}


static void
printends()			/* print ends of prism */
{
	register int  i;
						/* bottom face */
	printf("\n%s polygon %s.b\n", pmtype, pname);
	printf("0\n0\n%d\n", nverts*3);
	for (i = 0; i < nverts; i++)
		printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[i][0],
				vert[i][1], 0.0);
						/* top face */
	printf("\n%s polygon %s.t\n", pmtype, pname);
	printf("0\n0\n%d\n", nverts*3);
	for (i = nverts; i--; )
		printf("\t%18.12g\t%18.12g\t%18.12g\n", vert[i][0]+lvect[0],
				vert[i][1]+lvect[1], lvect[2]);
}


static void
printrends()			/* print ends of prism with rounding */
{
	register int  i;
	double	c0[3], c1[3], cl[3];
						/* bottom face */
	printf("\n%s polygon %s.b\n", pmtype, pname);
	printf("0\n0\n%d\n", nverts*3);
	for (i = 0; i < nverts; i++) {
		printf("\t%18.12g\t%18.12g\t%18.12g\n",
			vert[i][0] + crad*(a[i]*u[i][0] - u[i][1]),
			vert[i][1] + crad*(a[i]*u[i][1] + u[i][0]),
			0.0);
	}
						/* top face */
	printf("\n%s polygon %s.t\n", pmtype, pname);
	printf("0\n0\n%d\n", nverts*3);
	for (i = nverts; i--; ) {
		printf("\t%18.12g\t%18.12g\t%18.12g\n",
		vert[i][0] + lvect[0] + crad*(a[i]*u[i][0] - u[i][1]),
		vert[i][1] + lvect[1] + crad*(a[i]*u[i][1] + u[i][0]),
		lvect[2]);
	}
						/* bottom corners and edges */
	c0[0] = cl[0] = vert[nverts-1][0] +
			crad*(a[nverts-1]*u[nverts-1][0] - u[nverts-1][1]);
	c0[1] = cl[1] = vert[nverts-1][1] +
			crad*(a[nverts-1]*u[nverts-1][1] + u[nverts-1][0]);
	c0[2] = cl[2] = crad;
	for (i = 0; i < nverts; i++) {
		if (i < nverts-1) {
			c1[0] = vert[i][0] + crad*(a[i]*u[i][0] - u[i][1]);
			c1[1] = vert[i][1] + crad*(a[i]*u[i][1] + u[i][0]);
			c1[2] = crad;
		} else {
			c1[0] = cl[0]; c1[1] = cl[1]; c1[2] = cl[2];
		}
		if (lvdir*a[i] > 0.) {
			printf("\n%s sphere %s.bc%d\n", pmtype, pname, i+1);
			printf("0\n0\n4 %18.12g %18.12g %18.12g %18.12g\n",
				c1[0], c1[1], c1[2], lvdir*crad);
		}
		printf("\n%s cylinder %s.be%d\n", pmtype, pname, i+1);
		printf("0\n0\n7\n");
		printf("\t%18.12g\t%18.12g\t%18.12g\n", c0[0], c0[1], c0[2]);
		printf("\t%18.12g\t%18.12g\t%18.12g\n", c1[0], c1[1], c1[2]);
		printf("\t%18.12g\n", lvdir*crad);
		c0[0] = c1[0]; c0[1] = c1[1]; c0[2] = c1[2];
	}
						/* top corners and edges */
	c0[0] = cl[0] = vert[nverts-1][0] + lvect[0] +
			crad*(a[nverts-1]*u[nverts-1][0] - u[nverts-1][1]);
	c0[1] = cl[1] = vert[nverts-1][1] + lvect[1] +
			crad*(a[nverts-1]*u[nverts-1][1] + u[nverts-1][0]);
	c0[2] = cl[2] = lvect[2] - crad;
	for (i = 0; i < nverts; i++) {
		if (i < nverts-1) {
			c1[0] = vert[i][0] + lvect[0] +
					crad*(a[i]*u[i][0] - u[i][1]);
			c1[1] = vert[i][1] + lvect[1] +
					crad*(a[i]*u[i][1] + u[i][0]);
			c1[2] = lvect[2] - crad;
		} else {
			c1[0] = cl[0]; c1[1] = cl[1]; c1[2] = cl[2];
		}
		if (lvdir*a[i] > 0.) {
			printf("\n%s sphere %s.tc%d\n", pmtype, pname, i+1);
			printf("0\n0\n4 %18.12g %18.12g %18.12g %18.12g\n",
				c1[0], c1[1], c1[2], lvdir*crad);
		}
		printf("\n%s cylinder %s.te%d\n", pmtype, pname, i+1);
		printf("0\n0\n7\n");
		printf("\t%18.12g\t%18.12g\t%18.12g\n", c0[0], c0[1], c0[2]);
		printf("\t%18.12g\t%18.12g\t%18.12g\n", c1[0], c1[1], c1[2]);
		printf("\t%18.12g\n", lvdir*crad);
		c0[0] = c1[0]; c0[1] = c1[1]; c0[2] = c1[2];
	}
}


static void
printsides(round)		/* print prism sides */
int  round;
{
	register int  i;

	for (i = 0; i < nverts-1; i++)
		if (round)
			rside(i, i+1);
		else
			side(i, i+1);
	if (!iscomplete) {
		if (round)
			rside(nverts-1, 0);
		else
			side(nverts-1, 0);
	}
}


static void
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


int
main(argc, argv)
int  argc;
char  **argv;
{
	int  an;
	
	if (argc < 4)
		goto userr;

	pmtype = argv[1];
	pname = argv[2];

	if (!strcmp(argv[3], "-")) {
		readverts(NULL);
		an = 4;
	} else if (isdigit(argv[3][0])) {
		nverts = atoi(argv[3]);
		if (argc-3 < 2*nverts)
			goto userr;
		for (an = 0; an < nverts; an++) {
			vert[an][0] = atof(argv[2*an+4]);
			vert[an][1] = atof(argv[2*an+5]);
		}
		an = 2*nverts+4;
	} else {
		readverts(argv[3]);
		an = 4;
	}
	if (nverts < 3) {
		fprintf(stderr, "%s: not enough vertices\n", argv[0]);
		exit(1);
	}

	for ( ; an < argc; an++) {
		if (argv[an][0] != '-')
			goto userr;
		switch (argv[an][1]) {
		case 'l':				/* length vector */
			lvect[0] = atof(argv[++an]);
			lvect[1] = atof(argv[++an]);
			lvect[2] = atof(argv[++an]);
			if (lvect[2] < -FTINY)
				lvdir = -1;
			else if (lvect[2] > FTINY)
				lvdir = 1;
			else {
				fprintf(stderr,
					"%s: illegal extrusion vector\n",
						argv[0]);
				exit(1);
			}
			llen = sqrt(lvect[0]*lvect[0] + lvect[1]*lvect[1] +
					lvect[2]*lvect[2]);
			break;
		case 'r':				/* radius */
			crad = atof(argv[++an]);
			break;
		case 'e':				/* ends */
			do_ends = !do_ends;
			break;
		case 'c':				/* complete */
			iscomplete = !iscomplete;
			break;
		default:
			goto userr;
		}
	}
	if (crad > FTINY) {
		if (crad > lvdir*lvect[2]) {
			fprintf(stderr, "%s: rounding greater than height\n",
					argv[0]);
			exit(1);
		}
		crad *= lvdir;		/* simplifies formulas */
		compute_rounding();
		printhead(argc, argv);
		if (do_ends)
			printrends();
		printsides(1);
	} else {
		printhead(argc, argv);
		if (do_ends)
			printends();
		printsides(0);
	}
	exit(0);
userr:
	fprintf(stderr, "Usage: %s material name ", argv[0]);
	fprintf(stderr, "{ - | vfile | N v1 v2 .. vN } ");
	fprintf(stderr, "[-l lvect][-r radius][-c][-e]\n");
	exit(1);
}

