#ifndef lint
static const char	RCSid[] = "$Id: geom2mgf.c,v 1.2 2003/02/28 20:11:29 greg Exp $";
#endif
/* Convert Pete Shirley's simple format to MGF */

#include <stdio.h>

#define PI		3.14159265358979323846
#define WHTEFFICACY	179.0		/* lumens/watt conversion factor */

#define FEQ(a,b)	((a)-(b) < 1e-4 && (a)-(b) > -1e-4)


main(argc, argv)		/* input files are passed parameters */
int	argc;
char	**argv;
{
	int	i;

	if (argc == 1)
		translate(NULL);
	else
		for (i = 1; i < argc; i++)
			translate(argv[i]);
	exit(0);
}


translate(fname)		/* translate file fname */
char	*fname;
{
	FILE	*fp;
	register int	c;

	if (fname == NULL) {
		fname = "<stdin>";
		fp = stdin;
	} else if ((fp = fopen(fname, "r")) == NULL) {
		perror(fname);
		exit(1);
	}
	while ((c = getc(fp)) != EOF)
		switch (c) {
		case ' ':		/* white space */
		case '\t':
		case '\n':
		case '\f':
		case '\r':
			break;
		case 'q':		/* quad */
			doquad(fname, fp);
			break;
		case 's':		/* sphere */
			dosphere(fname, fp);
			break;
		default:		/* illegal */
			fprintf(stderr, "%s: illegal character (%d)\n",
					fname, c);
			exit(1);
		}
}


doquad(fn, fp)			/* translate quadrilateral */
char	*fn;
FILE	*fp;
{
	double	vert[3];
	register int	i;

	for (i = 0; i < 4; i++) {
		if (fscanf(fp, "%lf %lf %lf", &vert[0],
				&vert[1], &vert[2]) != 3) {
			fprintf(stderr, "%s: bad quad\n", fn);
			exit(1);
		}
		printf("v qv%d =\n\tp %.9g %.9g %.9g\n",
				i, vert[0], vert[1], vert[2]);
	}
	domat(fn, fp);
	printf("f qv0 qv1 qv2 qv3\n");
}


dosphere(fn, fp)		/* translate sphere */
char	*fn;
FILE	*fp;
{
	double	cent[3], radius;

	if (fscanf(fp, "%lf %lf %lf %lf", &cent[0], &cent[1],
			&cent[2], &radius) != 4) {
		fprintf(stderr, "%s: bad sphere\n", fn);
		exit(1);
	}
	printf("v sc = \n\tp %.9g %.9g %.9g\n",
			cent[0], cent[1], cent[2]);
	domat(fn, fp);
	printf("sph sc %.9g\n", radius);
}


domat(fn, fp)			/* translate material */
char	*fn;
FILE	*fp;
{
	static short	curtype = 0;
	static double	curv1, curv2;
	double	d1, d2;
	double	f;

nextchar:
	switch(getc(fp)) {
	case ' ':
	case '\t':
	case '\n':
	case '\f':
	case '\r':
		goto nextchar;
	case 'd':			/* diffuse reflector */
		if (fscanf(fp, "%lf", &d1) != 1)
			break;
		if (curtype == 'd' && FEQ(d1, curv1))
			return;
		curtype = 'd';
		curv1 = d1;
		printf("m\n\tc\n");
		printf("\trd %.4f\n", d1);
		return;
	case 's':			/* specular reflector */
		if (fscanf(fp, "%lf", &d1) != 1)
			break;
		if (curtype == 's' && FEQ(d1, curv1))
			return;
		curtype = 's';
		curv1 = d1;
		printf("m\n\tc\n");
		printf("\trs %.4f 0\n", d1);
		return;
	case 't':			/* dielectric */
		if (fscanf(fp, "%lf", &d1) != 1)
			break;
		if (curtype == 't' && FEQ(d1, curv1))
			return;
		curtype = 't';
		curv1 = d1;
		printf("m\n\tc\n");
		printf("\tsides 1\n\tir %.3f 0\n", d1);
		f = (1.0 - d1)/(1.0 + d1);
		f *= f;
		printf("\trs %.4f 0\n", f);
		printf("\tts %.4f 0\n", 1.0-f);
		return;
	case 'l':			/* light source */
		if (fscanf(fp, "%lf %lf", &d1, &d2) != 2)
			break;
		if (curtype == 'l' && FEQ(d1, curv1) && FEQ(d2, curv2))
			return;
		curtype = 'l';
		curv1 = d1;
		curv2 = d2;
		printf("m\n\tc\n");
		printf("\tsides 1\n\trd %.4f\n\ted %.2f\n",
				d1, PI*WHTEFFICACY*d2);
		return;
	}
	fprintf(stderr, "%s: material format error\n", fn);
	exit(1);
}
