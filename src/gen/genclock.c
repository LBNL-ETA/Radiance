#ifndef lint
static const char	RCSid[] = "$Id: genclock.c,v 2.5 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 * Generate an analog clock.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define PI		3.14159265358979323846

#define FACEBITMAP	"clockface.hex"

char	myfacemat[] = "white_plastic";
double	myfacearg[5] = {.85,.85,.85,0,0};
char	mycasemat[] = "black_plastic";
double	mycasearg[5] = {.08,.08,.08,.03,.04};
char	*facemat = myfacemat;
char	*casemat = mycasemat;
char	*name = "clock";

static void
genmats()				/* put out our materials */
{
	if (facemat == myfacemat)
		printf("\nvoid plastic %s\n0\n0\n5 %f %f %f %f %f\n",
				myfacemat, myfacearg[0], myfacearg[1],
				myfacearg[2], myfacearg[3], myfacearg[4]);
	if (casemat == mycasemat)
		printf("\nvoid plastic %s\n0\n0\n5 %f %f %f %f %f\n",
				mycasemat, mycasearg[0], mycasearg[1],
				mycasearg[2], mycasearg[3], mycasearg[4]);
	printf("\n%s brighttext clock_face_paint\n2 hexbit4x1.fnt %s\n",
			facemat, FACEBITMAP);
	printf("0\n11\n\t0\t-1\t1\n\t0\t.0185\t0\n\t0\t0\t-.00463\n");
	printf("\t.02\t1\n");
	printf("\nvoid glass clock_crystal\n0\n0\n3 .95 .95 .95\n");
	printf("\nvoid plastic hand_paint\n0\n0\n5 .03 .03 .03 0 0\n");
}


static void
genclock()				/* put out clock body */
{
	printf("\n%s ring %s.case_back\n", casemat, name);
	printf("0\n0\n8\t0\t0\t0\n\t-1\t0\t0\n\t0\t1.1\n");
	printf("\n%s cylinder %s.case_outer\n", casemat, name);
	printf("0\n0\n7\t0\t0\t0\n\t.12\t0\t0\n\t1.1\n");
	printf("\n%s ring %s.case_front\n", casemat, name);
	printf("0\n0\n8\t.12\t0\t0\n\t1\t0\t0\n\t1\t1.1\n");
	printf("\n%s cylinder %s.case_inner\n", casemat, name);
	printf("0\n0\n7\t.05\t0\t0\n\t.12\t0\t0\n\t1\n");
	printf("\nclock_crystal ring %s.crystal\n", name);
	printf("0\n0\n8\t.10\t0\t0\n\t1\t0\t0\n\t0\t1\n");
	printf("\nclock_face_paint ring %s.face\n", name);
	printf("0\n0\n8\t.05\t0\t0\n\t1\t0\t0\n\t0\t1\n");
}


static void
rvert(x, y, z, ang)			/* print rotated vertex */
double	x, y, z, ang;
{
	static double	lastang=0, sa=0, ca=1;

	if (ang != lastang) {
		sa = sin(-ang);
		ca = cos(-ang);
		lastang = ang;
	}
	printf("%15.12g %15.12g %15.12g\n", x, y*ca-z*sa, z*ca+y*sa);
}


static void
genhands(hour)				/* generate correct hand positions */
double	hour;
{
	double	hrot, mrot;

	hrot = 2.*PI/12. * hour;
	mrot = 2.*PI * (hour - floor(hour));

	printf("\nhand_paint polygon %s.hour_hand\n", name);
	printf("0\n0\n12\n");
	rvert(.06, -.03, -.06, hrot);
	rvert(.06, .03, -.06, hrot);
	rvert(.06, .025, .5, hrot);
	rvert(.06, -.025, .5, hrot);

	printf("\nhand_paint polygon %s.minute_hand\n", name);
	printf("0\n0\n12\n");
	rvert(.07, -.02, -.1, mrot);
	rvert(.07, .02, -.1, mrot);
	rvert(.07, .01, .9, mrot);
	rvert(.07, -.01, .9, mrot);
}


int
main(argc, argv)
int	argc;
char	*argv[];
{
	int	i, j;
	double	hour;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'f':
			facemat = argv[++i];
			break;
		case 'c':
			casemat = argv[++i];
			break;
		case 'n':
			name = argv[++i];
			break;
		default:
			goto userr;
		}
	if (i >= argc)
		goto userr;
	if (!isdigit(argv[i][0]))
		goto userr;
	for (j = 1; isdigit(argv[i][j]); j++)
		;
	if (argv[i][j] == ':')
		hour = atoi(argv[i]) + atoi(argv[i]+j+1)/60.0;
	else if (!argv[i][j] || argv[i][j] == '.')
		hour = atof(argv[i]);
	else
		goto userr;
	putchar('#');			/* print header */
	for (i = 0; i < argc; i++) {
		putchar(' ');
		fputs(argv[i], stdout);
	}
	putchar('\n');
	genmats();			/* print materials */
	genclock();			/* generate clock */
	genhands(hour);			/* generate hands */
	exit(0);
userr:
	fputs("Usage: ", stderr);
	fputs(argv[0], stderr);
	fputs(" [-f face_mat][-c case_mat][-n name] {HH:MM | HH.hh}\n", stderr);
	exit(1);
}

