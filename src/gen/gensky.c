#ifndef lint
static const char	RCSid[] = "$Id: gensky.c,v 2.25 2010/12/08 21:22:33 greg Exp $";
#endif
/*
 *  gensky.c - program to generate sky functions.
 *		Our zenith is along the Z-axis, the X-axis
 *		points east, and the Y-axis points north.
 *		Radiance is in watts/steradian/sq. meter.
 *
 *     3/26/86
 */

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <ctype.h>

#include  "color.h"

extern int jdate(int month, int day);
extern double stadj(int  jd);
extern double sdec(int  jd);
extern double salt(double sd, double st);
extern double sazi(double sd, double st);

#ifndef  PI
#define  PI		3.14159265358979323846
#endif

#define  DOT(v1,v2)	(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define  S_CLEAR	1
#define  S_OVER		2
#define  S_UNIF		3
#define  S_INTER	4

#define  overcast	((skytype==S_OVER)|(skytype==S_UNIF))

double  normsc();
					/* sun calculation constants */
extern double  s_latitude;
extern double  s_longitude;
extern double  s_meridian;

#undef  toupper
#define  toupper(c)	((c) & ~0x20)	/* ASCII trick to convert case */

					/* European and North American zones */
struct {
	char	zname[8];	/* time zone name (all caps) */
	float	zmer;		/* standard meridian */
} tzone[] = {
	{"YST", 135}, {"YDT", 120},
	{"PST", 120}, {"PDT", 105},
	{"MST", 105}, {"MDT", 90},
	{"CST", 90}, {"CDT", 75},
	{"EST", 75}, {"EDT", 60},
	{"AST", 60}, {"ADT", 45},
	{"NST", 52.5}, {"NDT", 37.5},
	{"GMT", 0}, {"BST", -15},
	{"CET", -15}, {"CEST", -30},
	{"EET", -30}, {"EEST", -45},
	{"AST", -45}, {"ADT", -60},
	{"GST", -60}, {"GDT", -75},
	{"IST", -82.5}, {"IDT", -97.5},
	{"JST", -135}, {"NDT", -150},
	{"NZST", -180}, {"NZDT", -195},
	{"", 0}
};
					/* required values */
int  month, day;				/* date */
double  hour;					/* time */
int  tsolar;					/* 0=standard, 1=solar */
double  altitude, azimuth;			/* or solar angles */
					/* default values */
int  skytype = S_CLEAR;				/* sky type */
int  dosun = 1;
double  zenithbr = 0.0;
int	u_zenith = 0;				/* -1=irradiance, 1=radiance */
double  turbidity = 2.45;
double  gprefl = 0.2;
					/* computed values */
double  sundir[3];
double  groundbr;
double  F2;
double  solarbr = 0.0;
int	u_solar = 0;				/* -1=irradiance, 1=radiance */

char  *progname;
char  errmsg[128];

void computesky(void);
void printsky(void);
void printdefaults(void);
void userror(char  *msg);
double normsc(void);
int cvthour(char  *hs);
void printhead(register int  ac, register char  **av);


int
main(argc, argv)
int  argc;
char  *argv[];
{
	int  got_meridian = 0;
	int  i;

	progname = argv[0];
	if (argc == 2 && !strcmp(argv[1], "-defaults")) {
		printdefaults();
		exit(0);
	}
	if (argc < 4)
		userror("arg count");
	if (!strcmp(argv[1], "-ang")) {
		altitude = atof(argv[2]) * (PI/180);
		azimuth = atof(argv[3]) * (PI/180);
		month = 0;
	} else {
		month = atoi(argv[1]);
		if (month < 1 || month > 12)
			userror("bad month");
		day = atoi(argv[2]);
		if (day < 1 || day > 31)
			userror("bad day");
		got_meridian = cvthour(argv[3]);
	}
	for (i = 4; i < argc; i++)
		if (argv[i][0] == '-' || argv[i][0] == '+')
			switch (argv[i][1]) {
			case 's':
				skytype = S_CLEAR;
				dosun = argv[i][0] == '+';
				break;
			case 'r':
			case 'R':
				u_solar = argv[i][1]=='R' ? -1 : 1;
				solarbr = atof(argv[++i]);
				break;
			case 'c':
				skytype = S_OVER;
				break;
			case 'u':
				skytype = S_UNIF;
				break;
			case 'i':
				skytype = S_INTER;
				dosun = argv[i][0] == '+';
				break;
			case 't':
				turbidity = atof(argv[++i]);
				break;
			case 'b':
			case 'B':
				u_zenith = argv[i][1]=='B' ? -1 : 1;
				zenithbr = atof(argv[++i]);
				break;
			case 'g':
				gprefl = atof(argv[++i]);
				break;
			case 'a':
				s_latitude = atof(argv[++i]) * (PI/180);
				break;
			case 'o':
				s_longitude = atof(argv[++i]) * (PI/180);
				break;
			case 'm':
				if (got_meridian) {
					++i;
					break;		/* time overrides */
				}
				s_meridian = atof(argv[++i]) * (PI/180);
				break;
			default:
				sprintf(errmsg, "unknown option: %s", argv[i]);
				userror(errmsg);
			}
		else
			userror("bad option");

	if (month && !tsolar && fabs(s_meridian-s_longitude) > 45*PI/180)
		fprintf(stderr,
	"%s: warning - %.1f hours btwn. standard meridian and longitude\n",
			progname, (s_longitude-s_meridian)*12/PI);

	printhead(argc, argv);

	computesky();
	printsky();

	exit(0);
}


void
computesky(void)			/* compute sky parameters */
{
	double	normfactor;
					/* compute solar direction */
	if (month) {			/* from date and time */
		int  jd;
		double  sd, st;

		jd = jdate(month, day);		/* Julian date */
		sd = sdec(jd);			/* solar declination */
		if (tsolar)			/* solar time */
			st = hour;
		else
			st = hour + stadj(jd);
		altitude = salt(sd, st);
		azimuth = sazi(sd, st);
		printf("# Local solar time: %.2f\n", st);
		printf("# Solar altitude and azimuth: %.1f %.1f\n",
				180./PI*altitude, 180./PI*azimuth);
	}
	if (!overcast && altitude > 87.*PI/180.) {
		fprintf(stderr,
"%s: warning - sun too close to zenith, reducing altitude to 87 degrees\n",
				progname);
		printf(
"# warning - sun too close to zenith, reducing altitude to 87 degrees\n");
		altitude = 87.*PI/180.;
	}
	sundir[0] = -sin(azimuth)*cos(altitude);
	sundir[1] = -cos(azimuth)*cos(altitude);
	sundir[2] = sin(altitude);

					/* Compute normalization factor */
	switch (skytype) {
	case S_UNIF:
		normfactor = 1.0;
		break;
	case S_OVER:
		normfactor = 0.777778;
		break;
	case S_CLEAR:
		F2 = 0.274*(0.91 + 10.0*exp(-3.0*(PI/2.0-altitude)) +
				0.45*sundir[2]*sundir[2]);
		normfactor = normsc()/F2/PI;
		break;
	case S_INTER:
		F2 = (2.739 + .9891*sin(.3119+2.6*altitude)) *
			exp(-(PI/2.0-altitude)*(.4441+1.48*altitude));
		normfactor = normsc()/F2/PI;
		break;
	}
					/* Compute zenith brightness */
	if (u_zenith == -1)
		zenithbr /= normfactor*PI;
	else if (u_zenith == 0) {
		if (overcast)
			zenithbr = 8.6*sundir[2] + .123;
		else
			zenithbr = (1.376*turbidity-1.81)*tan(altitude)+0.38;
		if (skytype == S_INTER)
			zenithbr = (zenithbr + 8.6*sundir[2] + .123)/2.0;
		if (zenithbr < 0.0)
			zenithbr = 0.0;
		else
			zenithbr *= 1000.0/SKYEFFICACY;
	}
					/* Compute horizontal radiance */
	groundbr = zenithbr*normfactor;
	printf("# Ground ambient level: %.1f\n", groundbr);
	if (!overcast && sundir[2] > 0.0 && (!u_solar || solarbr > 0.0)) {
		if (u_solar == -1)
			solarbr /= 6e-5*sundir[2];
		else if (u_solar == 0) {
			solarbr = 1.5e9/SUNEFFICACY *
			(1.147 - .147/(sundir[2]>.16?sundir[2]:.16));
			if (skytype == S_INTER)
				solarbr *= 0.15;	/* fudge factor! */
		}
		groundbr += 6e-5/PI*solarbr*sundir[2];
	} else
		dosun = 0;
	groundbr *= gprefl;
}


void
printsky(void)			/* print out sky */
{
	if (dosun) {
		printf("\nvoid light solar\n");
		printf("0\n0\n");
		printf("3 %.2e %.2e %.2e\n", solarbr, solarbr, solarbr);
		printf("\nsolar source sun\n");
		printf("0\n0\n");
		printf("4 %f %f %f 0.5\n", sundir[0], sundir[1], sundir[2]);
	}
	
	printf("\nvoid brightfunc skyfunc\n");
	printf("2 skybr skybright.cal\n");
	printf("0\n");
	if (overcast)
		printf("3 %d %.2e %.2e\n", skytype, zenithbr, groundbr);
	else
		printf("7 %d %.2e %.2e %.2e %f %f %f\n",
				skytype, zenithbr, groundbr, F2,
				sundir[0], sundir[1], sundir[2]);
}


void
printdefaults(void)			/* print default values */
{
	switch (skytype) {
	case S_OVER:
		printf("-c\t\t\t\t# Cloudy sky\n");
		break;
	case S_UNIF:
		printf("-u\t\t\t\t# Uniform cloudy sky\n");
		break;
	case S_INTER:
		if (dosun)
			printf("+i\t\t\t\t# Intermediate sky with sun\n");
		else
			printf("-i\t\t\t\t# Intermediate sky without sun\n");
		break;
	case S_CLEAR:
		if (dosun)
			printf("+s\t\t\t\t# Sunny sky with sun\n");
		else
			printf("-s\t\t\t\t# Sunny sky without sun\n");
		break;
	}
	printf("-g %f\t\t\t# Ground plane reflectance\n", gprefl);
	if (zenithbr > 0.0)
		printf("-b %f\t\t\t# Zenith radiance (watts/ster/m2\n", zenithbr);
	else
		printf("-t %f\t\t\t# Atmospheric turbidity\n", turbidity);
	printf("-a %f\t\t\t# Site latitude (degrees)\n", s_latitude*(180/PI));
	printf("-o %f\t\t\t# Site longitude (degrees)\n", s_longitude*(180/PI));
	printf("-m %f\t\t\t# Standard meridian (degrees)\n", s_meridian*(180/PI));
}


void
userror(			/* print usage error and quit */
	char  *msg
)
{
	if (msg != NULL)
		fprintf(stderr, "%s: Use error - %s\n", progname, msg);
	fprintf(stderr, "Usage: %s month day hour [options]\n", progname);
	fprintf(stderr, "   Or: %s -ang altitude azimuth [options]\n", progname);
	fprintf(stderr, "   Or: %s -defaults\n", progname);
	exit(1);
}


double
normsc(void)			/* compute normalization factor (E0*F2/L0) */
{
	static double  nfc[2][5] = {
				/* clear sky approx. */
		{2.766521, 0.547665, -0.369832, 0.009237, 0.059229},
				/* intermediate sky approx. */
		{3.5556, -2.7152, -1.3081, 1.0660, 0.60227},
	};
	register double  *nf;
	double  x, nsc;
	register int  i;
					/* polynomial approximation */
	nf = nfc[skytype==S_INTER];
	x = (altitude - PI/4.0)/(PI/4.0);
	nsc = nf[i=4];
	while (i--)
		nsc = nsc*x + nf[i];

	return(nsc);
}


int
cvthour(			/* convert hour string */
	char  *hs
)
{
	register char  *cp = hs;
	register int	i, j;

	if ( (tsolar = *cp == '+') ) cp++;		/* solar time? */
	while (isdigit(*cp)) cp++;
	if (*cp == ':')
		hour = atoi(hs) + atoi(++cp)/60.0;
	else {
		hour = atof(hs);
		if (*cp == '.') cp++;
	}
	while (isdigit(*cp)) cp++;
	if (!*cp)
		return(0);
	if (tsolar || !isalpha(*cp)) {
		fprintf(stderr, "%s: bad time format: %s\n", progname, hs);
		exit(1);
	}
	i = 0;
	do {
		for (j = 0; cp[j]; j++)
			if (toupper(cp[j]) != tzone[i].zname[j])
				break;
		if (!cp[j] && !tzone[i].zname[j]) {
			s_meridian = tzone[i].zmer * (PI/180);
			return(1);
		}
	} while (tzone[i++].zname[0]);

	fprintf(stderr, "%s: unknown time zone: %s\n", progname, cp);
	fprintf(stderr, "Known time zones:\n\t%s", tzone[0].zname);
	for (i = 1; tzone[i].zname[0]; i++)
		fprintf(stderr, " %s", tzone[i].zname);
	putc('\n', stderr);
	exit(1);
}


void
printhead(		/* print command header */
	register int  ac,
	register char  **av
)
{
	putchar('#');
	while (ac--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');
}
