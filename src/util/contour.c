#ifndef lint
static const char	RCSid[] = "$Id: contour.c,v 1.4 2004/03/26 23:34:23 schorsch Exp $";
#endif
/*
 *  contour.c - program to make contour plots, mappings from 3 to
 *		2 dimenstions.
 *
 *     8/22/86
 */

#include  <stdio.h>

#include  <stdlib.h>

#include  <ctype.h>

#define  MAXPTS		2048		/* maximum number of input points */

typedef float  DATAPT[3];

DATAPT  xyz[MAXPTS];			/* the input data */
int  xyzsiz;

double  zmin = 1e20;			/* z minimum */
double  zmax = -1e20;			/* z maximum */
int  minset = 0;			/* user set minimum? */
int  maxset = 0;			/* user set maximum? */

int  ncurves = 6;			/* number of contours */

static void getdata(FILE  *fp);
static int xycmp(const void*  p1, const void* p2);
static void contour(int  n);
static void crossings(double  z);
static void cross(DATAPT  p0, DATAPT  p1, double  z);


int
main(
	int  argc,
	char  *argv[]
)
{
	FILE  *fp;
	int  i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'n':
			ncurves = atoi(argv[++i]);
			break;
		case 'l':
			zmin = atof(argv[++i]);
			minset = 1;
			break;
		case 'u':
			zmax = atof(argv[++i]);
			maxset = 1;
			break;
		default:
			fprintf(stderr, "%s: unknown option: %s\n",
					argv[0], argv[i]);
			exit(1);
		}

	if (i == argc)
		getdata(stdin);
	else if (i < argc-1) {
		fprintf(stderr, "%s: too many file names\n", argv[0]);
		exit(1);
	} else if ((fp = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "%s: file not found\n", argv[i]);
		exit(1);
	} else
		getdata(fp);

	qsort(xyz, xyzsiz, sizeof(DATAPT), xycmp);	/* sort data */

	for (i = 0; i < ncurves; i++)			/* do contours */
		contour(i);

	exit(0);
}


static void
getdata(			/* read input data */
	FILE  *fp
)
{
	register int  i;

	while (xyzsiz < MAXPTS) {
		for (i = 0; i < 3; i++)
			if (fscanf(fp, "%f", &xyz[xyzsiz][i]) != 1)
				break;
		if (i != 3)
			break;
		if (xyz[xyzsiz][2] < zmin) {
			if (minset) {
				continue;
			} else {
				zmin = xyz[xyzsiz][2];
			}
		}
		if (xyz[xyzsiz][2] > zmax) {
			if (maxset) {
				continue;
			} else {
				zmax = xyz[xyzsiz][2];
			}
		}
		xyzsiz++;
	}
	if (!feof(fp)) {
		fprintf(stderr, "Error reading input data\n");
		exit(1);
	}
}


static int
xycmp(			/* -1 if p1 < p2, 0 if equal, 1 if p1 > p2 */
	const void*  p1,
	const void*	p2
)
{
	const float* pp1 = p1;
	const float* pp2 = p2;
	if (pp1[0] > pp2[0])
		return(1);
	if (pp1[0] < pp2[0])
		return(-1);
	if (pp1[1] > pp2[1])
		return(1);
	if (pp1[1] < pp2[1])
		return(-1);
	return(0);
}


static void
contour(			/* make contour n */
	int  n
)
{
	double  z;

	z = (n+.5)*(zmax-zmin)/ncurves + zmin;
	printf("%clabel=\"%f\";\n", n+'A', z);
	printf("%cdata=\n", n+'A');
	crossings(z);
	printf(";\n");
}


static void
crossings(			/* find crossings for z */
	double  z
)
{
	register DATAPT  *p0, *p1;

	p0 = p1 = xyz;
	while (p0 < xyz+xyzsiz-2) {
		while (p1 < xyz+xyzsiz-2) {
			if (p0[0][0] < p1[0][0] && p0[0][1] <= p1[0][1])
				break;
			p1++;
		}
		if (p0[0][0] == p0[1][0])
			cross(p0[0], p0[1], z);
		if (p0[0][1] == p1[0][1])
			cross(p0[0], p1[0], z);
		p0++;
	}
}


static void
cross(		/* mark crossing between p0 and p1 */
	register DATAPT  p0,
	register DATAPT  p1,
	double  z
)
{
	if (p1[2] - p0[2] == 0.0)
		return;
	z = (z - p0[2])/(p1[2] - p0[2]);
	if (z < 0.0 || z >= 1.0)
		return;
	printf("%e\t", p0[0]*(1.0 - z) + p1[0]*z);
	printf("%e\n", p0[1]*(1.0 - z) + p1[1]*z);
}
