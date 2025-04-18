#ifndef lint
static const char	RCSid[] = "$Id: genworm.c,v 2.14 2025/04/18 23:02:53 greg Exp $";
#endif
/*
 *  genworm.c - program to generate worms (strings with varying thickness).
 *
 *	The program takes as input the functions of t for x, y,
 *	z and r (the radius).  Cylinders and cones will be
 *	joined by spheres.  Negative radii are forbidden.
 *
 *	9/24/87
 */

#include  <stdlib.h>
#include  <math.h>

#include  "calcomp.h"
#include  "rtio.h"
#include  "resolu.h"
#include  "rterror.h"
#include  "fvect.h"

#define  XNAME		"X`SYS`"		/* x function name */
#define  YNAME		"Y`SYS`"		/* y function name */
#define  ZNAME		"Z`SYS`"		/* z function name */
#define  RNAME		"R`SYS`"		/* r function name */

#define  PI		3.14159265358979323846

#define  max(a,b)	((a) > (b) ? (a) : (b))


/* XXX redundant, move to library */
double  l_hermite(char *), l_bezier(char *), l_bspline(char *);


int
main(argc, argv)
int  argc;
char  *argv[];
{
	char  stmp[256];
	double  t, f, lastr = 0, r;
	FVECT  lastp, p;
	int  i, nseg;

	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_OUTCHAN|E_INCHAN);
	varset("PI", ':', PI);
	funset("hermite", 5, ':', l_hermite);
	funset("bezier", 5, ':', l_bezier);
	funset("bspline", 5, ':', l_bspline);

	if (argc < 8)
		goto userror;

	for (i = 8; i < argc; i++)
		if (!strcmp(argv[i], "-e"))
			scompile(argv[++i], NULL, 0);
		else if (!strcmp(argv[i], "-f")) {
			char  *fpath = getpath(argv[++i], getrlibpath(), 0);
			if (fpath == NULL) {
				fprintf(stderr, "%s: cannot find file '%s'\n",
						argv[0], argv[i]);
				quit(1);
			}
			fcompile(fpath);
		} else
			goto userror;

	sprintf(stmp, "%s(t)=%s;", XNAME, argv[3]);
	scompile(stmp, NULL, 0);
	sprintf(stmp, "%s(t)=%s;", YNAME, argv[4]);
	scompile(stmp, NULL, 0);
	sprintf(stmp, "%s(t)=%s;", ZNAME, argv[5]);
	scompile(stmp, NULL, 0);
	sprintf(stmp, "%s(t)=%s;", RNAME, argv[6]);
	scompile(stmp, NULL, 0);
	nseg = eval(argv[7]) + .5;
	if (nseg <= 0)
		goto userror;

	fputs("# ", stdout);
	printargs(argc, argv, stdout);
	doptimize(1);
	eclock++;

	for (i = 0; i <= nseg; i++) {
		t = (double)i/nseg;
		p[0] = funvalue(XNAME, 1, &t);
		p[1] = funvalue(YNAME, 1, &t);
		p[2] = funvalue(ZNAME, 1, &t);
		r = funvalue(RNAME, 1, &t);
		if (i) {
			if (lastr <= r+FTINY && lastr >= r-FTINY) {
				printf("\n%s cylinder %s.c%d\n",
						argv[1], argv[2], i);
				printf("0\n0\n7\n");
				printf("%18.12g %18.12g %18.12g\n",
						lastp[0], lastp[1], lastp[2]);
				printf("%18.12g %18.12g %18.12g\n",
						p[0], p[1], p[2]);
				printf("%18.12g\n", r);
			} else {
				printf("\n%s cone %s.c%d\n",
						argv[1], argv[2], i);
				printf("0\n0\n8\n");
				f = (lastr - r)/dist2(lastp,p);
				printf("%18.12g %18.12g %18.12g\n",
					lastp[0] + f*lastr*(p[0] - lastp[0]),
					lastp[1] + f*lastr*(p[1] - lastp[1]),
					lastp[2] + f*lastr*(p[2] - lastp[2]));
				printf("%18.12g %18.12g %18.12g\n",
					p[0] + f*r*(p[0] - lastp[0]),
					p[1] + f*r*(p[1] - lastp[1]),
					p[2] + f*r*(p[2] - lastp[2]));
				f = 1.0 - (lastr-r)*f;
				f = f <= 0.0 ? 0.0 : sqrt(f);
				printf("%18.12g %18.12g\n", f*lastr, f*r);
			}
		}
		printf("\n%s sphere %s.s%d\n", argv[1], argv[2], i);
		printf("0\n0\n4 %18.12g %18.12g %18.12g %18.12g\n",
				p[0], p[1], p[2], r);
		VCOPY(lastp, p);
		lastr = r;
	}
	return 0;

userror:
	fprintf(stderr,
"Usage: %s material name x(t) y(t) z(t) r(t) nseg [-e expr] [-f file]\n",
			argv[0]);
	return 1;
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
