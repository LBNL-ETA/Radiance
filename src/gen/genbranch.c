#ifndef lint
static const char	RCSid[] = "$Id: genbranch.c,v 2.6 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 *  genbranch.c - program to generate 3D Christmas tree branches.
 *
 *     8/23/86
 */

#include  <stdio.h>

#include <stdlib.h>

#include  <math.h>

#include  "random.h"

#define  errf()		(var*(0.5-frandom()))

double  bstart[3] = {0.0, 0.0, 0.0};	/* start of branch */
double  bend[3] = {28.0, 8.0, 0.0};	/* end of branch */
double  bthick = .6;			/* branch radius at base */
double  bnarrow = .4;			/* ratio of tip radius to base */
double  bratio = .4;			/* ratio of limb to branch */
char  *branchmat = "m_bark";		/* bark material */

char  *leafmat = "m_leaf";		/* leaf material */

int  nshoots = 7;			/* number of offshoots */
int  rdepth = 3;			/* recursion depth */

double  var = 0.3;			/* variability */


static void
stick(mat, beg, end, rad)		/* output a branch or leaf */
char  *mat;
double  beg[3], end[3];
double  rad;
{
	static int  nsticks = 0;
	
	printf("\n%s cone s%d\n", mat, nsticks);
	printf("0\n0\n8\n");
	printf("\t%18.12g\t%18.12g\t%18.12g\n", beg[0], beg[1], beg[2]);
	printf("\t%18.12g\t%18.12g\t%18.12g\n", end[0], end[1], end[2]);
	printf("\t%18.12g\t%18.12g\n", rad, bnarrow*rad);

	printf("\n%s sphere e%d\n", mat, nsticks);
	printf("0\n0\n4");
	printf("\t%18.12g\t%18.12g\t%18.12g\t%18.12g\n",
			end[0], end[1], end[2], bnarrow*rad);

	nsticks++;
}


static void
branch(beg, end, rad, lvl)		/* generate branch recursively */
double  beg[3], end[3];
double  rad;
int  lvl;
{
	double  sqrt();
	double  newbeg[3], newend[3];
	double  t;
	int  i, j;
	
	if (lvl == 0) {
		stick(leafmat, beg, end, rad);
		return;
	}

	stick(branchmat, beg, end, rad);

	for (i = 1; i <= nshoots; i++) {
						/* right branch */
		t = (i+errf())/(nshoots+2);
		t = (t + sqrt(t))/2.0;
		for (j = 0; j < 3; j++) {
			newbeg[j] = newend[j] = (end[j]-beg[j])*t + beg[j];
			newend[j] += (end[j]-beg[j])*(1+errf())/(nshoots+2);
		}
		newend[0] += (end[2]-newbeg[2])*bratio*(1+errf());
		newend[1] += (end[1]-newbeg[1])*bratio*(1+errf());
		newend[2] -= (end[0]-newbeg[0])*bratio*(1+errf());
		branch(newbeg, newend,
				(1-(1-bnarrow)*t)*(1+2*bratio)/3*rad, lvl-1);
		
						/* left branch */
		t = (i+errf())/(nshoots+2);
		t = (t + sqrt(t))/2.0;
		for (j = 0; j < 3; j++) {
			newbeg[j] = newend[j] = (end[j]-beg[j])*t + beg[j];
			newend[j] += (end[j]-beg[j])*(1+errf())/(nshoots+2);
		}
		newend[0] -= (end[2]-newbeg[2])*bratio*(1+errf());
		newend[1] += (end[1]-newbeg[1])*bratio*(1+errf());
		newend[2] += (end[0]-newbeg[0])*bratio*(1+errf());
		branch(newbeg, newend,
				(1-(1-bnarrow)*t)*(1+2*bratio)/3*rad, lvl-1);
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
	int  i, j;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'b':			/* branch */
			switch (argv[i][2]) {
			case 's':			/* start */
				bstart[0] = atof(argv[++i]);
				bstart[1] = atof(argv[++i]);
				bstart[2] = atof(argv[++i]);
				break;
			case 'e':			/* end */
				bend[0] = atof(argv[++i]);
				bend[1] = atof(argv[++i]);
				bend[2] = atof(argv[++i]);
				break;
			case 't':			/* thickness */
				bthick = atof(argv[++i]);
				break;
			case 'n':			/* narrow */
				bnarrow = atof(argv[++i]);
				break;
			case 'r':			/* ratio */
				bratio = atof(argv[++i]);
				break;
			case 'm':			/* material */
				branchmat = argv[++i];
				break;
			default:
				goto unkopt;
			}
			break;
		case 'l':			/* leaf */
			switch (argv[i][2]) {
			case 'm':			/* material */
				leafmat = argv[++i];
				break;
			default:
				goto unkopt;
			}
			break;
		case 'n':			/* number of offshoots */
			nshoots = atoi(argv[++i]);
			break;
		case 'r':			/* recursion depth */
			rdepth = atoi(argv[++i]);
			break;
		case 's':			/* seed */
			j = atoi(argv[++i]);
			while (j-- > 0)
				frandom();
			break;
		case 'v':			/* variability */
			var = atof(argv[++i]);
			break;
		default:;
unkopt:			fprintf(stderr, "%s: unknown option: %s\n",
					argv[0], argv[i]);
			exit(1);
		}
	
	if (i != argc) {
		fprintf(stderr, "%s: bad argument\n", argv[0]);
		exit(1);
	}
	printhead(argc, argv);

	branch(bstart, bend, bthick, rdepth);

	return(0);
}

