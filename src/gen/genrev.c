#ifndef lint
static const char	RCSid[] = "$Id: genrev.c,v 2.6 2003/08/04 22:37:53 greg Exp $";
#endif
/*
 *  genrev.c - program to generate functions of rotation about z
 *
 *	The program takes as input the functions of t for z and r
 *	(the radius).  If z is increasing, the normal points
 *	outward, inward for decreasing z.  Negative radii are forbidden.
 *
 *	8/6/86
 */

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <math.h>
#include  "calcomp.h"

#define  ZNAME		"Z`SYS`"		/* z function name */
#define  RNAME		"R`SYS`"		/* r function name */

#define  PI		3.14159265358979323846

#define  FTINY		1e-9
					/* orientation */
#define  OUT		01
#define  IN		02
#define  UP		04
#define  DOWN		010

void	quit(), eputs(), wputs();


computen(nzp, nrp, z0, r0, z1, r1)		/* compute normal */
double  *nzp, *nrp, z0, r0, z1, r1;
{
	double  dr, dz, len;

	dz = r0 - r1;				/* right angle vector */
	dr = z1 - z0;
	len = sqrt(dr*dr + dz*dz);
	*nzp = dz/len;
	*nrp = dr/len;
}


void
eputs(msg)
char  *msg;
{
	fputs(msg, stderr);
}


void
wputs(msg)
char  *msg;
{
	eputs(msg);
}


void
quit(code)
int  code;
{
	exit(code);
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


double
l_hermite(char *nm)
{
	double  t;
	
	t = argument(5);
	return(	argument(1)*((2.0*t-3.0)*t*t+1.0) +
		argument(2)*(-2.0*t+3.0)*t*t +
		argument(3)*((t-2.0)*t+1.0)*t +
		argument(4)*(t-1.0)*t*t );
}


double
l_bezier(char *nm)
{
	double  t;

	t = argument(5);
	return(	argument(1) * (1.+t*(-3.+t*(3.-t))) +
		argument(2) * 3.*t*(1.+t*(-2.+t)) +
		argument(3) * 3.*t*t*(1.-t) +
		argument(4) * t*t*t );
}


double
l_bspline(char *nm)
{
	double  t;

	t = argument(5);
	return(	argument(1) * (1./6.+t*(-1./2.+t*(1./2.-1./6.*t))) +
		argument(2) * (2./3.+t*t*(-1.+1./2.*t)) +
		argument(3) * (1./6.+t*(1./2.+t*(1./2.-1./2.*t))) +
		argument(4) * (1./6.*t*t*t) );
}


main(argc, argv)
int  argc;
char  *argv[];
{
	char  stmp[256];
	char  *modname;
	int  smooth = 0;
	double  t, lastz, z, nextz, lastr, r, nextr;
	double  lastnz, lastnr, nz, nr, nextnz, nextnr;
	int  i, nseg;
	int  orient;

	varset("PI", ':', PI);
	funset("hermite", 5, ':', l_hermite);
	funset("bezier", 5, ':', l_bezier);
	funset("bspline", 5, ':', l_bspline);

	if (argc < 6)
		goto userror;

	for (i = 6; i < argc; i++)
		if (!strcmp(argv[i], "-e"))
			scompile(argv[++i], NULL, 0);
		else if (!strcmp(argv[i], "-f"))
			fcompile(argv[++i]);
		else if (!strcmp(argv[i], "-s"))
			smooth = 1;
		else
			goto userror;

	sprintf(stmp, "%s(t)=%s;", ZNAME, argv[3]);
	scompile(stmp, NULL, 0);
	sprintf(stmp, "%s(t)=%s;", RNAME, argv[4]);
	scompile(stmp, NULL, 0);
	nseg = atoi(argv[5]);
	if (nseg <= 0)
		goto userror;
	modname = smooth ? "Phong" : argv[1];

	printhead(argc, argv);
	eclock = 0;

	lastnz = lastnr = 0.0;
	t = 0.0;
	lastz = funvalue(ZNAME, 1, &t);
	lastr = funvalue(RNAME, 1, &t);
	t = 1.0/nseg;
	z = funvalue(ZNAME, 1, &t);
	r = funvalue(RNAME, 1, &t);
	computen(&nz, &nr, lastz, lastr, z, r);
	for (i = 1; i <= nseg; i++) {
		if (i < nseg) {
			t = (double)(i+1)/nseg;
			nextz = funvalue(ZNAME, 1, &t);
			nextr = funvalue(RNAME, 1, &t);
			computen(&nextnz, &nextnr, z, r, nextz, nextr);
		} else
			nextnz = nextnr = 0.0;
		orient = 0;
		if (z < lastz-FTINY)
			orient |= DOWN;
		else if (z > lastz+FTINY)
			orient |= UP;
		if (r < lastr-FTINY)
			orient |= IN;
		else if (r > lastr+FTINY)
			orient |= OUT;
		if (!orient)
			goto endfor;
		if (smooth) {
			printf("\n%s texfunc Phong\n", argv[1]);
			printf("4 rev_dx rev_dy rev_dz rev.cal\n");
			printf("0\n4\n");
			if (orient&(UP|DOWN)) {
				t = (nextnz - lastnz)/(z - lastz);
				printf("\t%15.9g\t%15.9g\n",
						t, lastnz - t*lastz);
			} else
				printf("\t0\t%d\n", orient&IN ? 1 : -1);
			if (orient&(OUT|IN))  {
				t = (nextnr - lastnr)/(r - lastr);
				printf("\t%15.9g\t%15.9g\n",
						t, lastnr - t*lastr);
			} else
				printf("\t0\t%d\n", orient&UP ? 1 : -1);
		}
		if (!(orient&(IN|OUT))) {
			printf("\n%s %s %s.%d\n", modname,
					orient&DOWN ? "tube" : "cylinder",
					argv[2], i);
			printf("0\n0\n7\n");
			printf("\t0\t0\t%15.9g\n", lastz);
			printf("\t0\t0\t%15.9g\n", z);
			printf("\t%15.9g\n", r);
		} else if (!(orient&(UP|DOWN))) {
			printf("\n%s ring %s.%d\n", modname, argv[2], i);
			printf("0\n0\n8\n");
			printf("\t0\t0\t%15.9g\n", z);
			printf("\t0\t0\t%15.9g\n", orient&IN ? 1.0 : -1.0);
			printf("\t%15.9g\t%15.9g\n", lastr, r);
		} else {
			printf("\n%s %s %s.%d\n", modname,
					orient&DOWN ? "cup" : "cone",
					argv[2], i);
			printf("0\n0\n8\n");
			printf("\t0\t0\t%15.9g\n", lastz);
			printf("\t0\t0\t%15.9g\n", z);
			printf("\t%15.9g\t%15.9g\n", lastr, r);
		}
	endfor:
		lastz = z; lastr = r;
		z = nextz; r = nextr;
		lastnz = nz; lastnr = nr;
		nz = nextnz; nr = nextnr;
	}
	quit(0);

userror:
	fprintf(stderr,
	"Usage: %s material name z(t) r(t) nseg [-e expr] [-f file] [-s]\n",
			argv[0]);
	quit(1);
}


