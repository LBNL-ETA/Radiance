#ifndef lint
static const char RCSid[] = "$Id: gensurf.c,v 2.12 2003/08/04 22:37:53 greg Exp $";
#endif
/*
 *  gensurf.c - program to generate functional surfaces
 *
 *	Parametric functions x(s,t), y(s,t) and z(s,t)
 *  specify the surface, which is tesselated into an m by n
 *  array of paired triangles.
 *	The surface normal is defined by the right hand
 *  rule applied to (s,t).
 *
 *	4/3/87
 *
 *	4/16/02	Added conditional vertex output
 */

#include  "standard.h"

char  XNAME[] =		"X`SYS";		/* x function name */
char  YNAME[] =		"Y`SYS";		/* y function name */
char  ZNAME[] =		"Z`SYS";		/* z function name */

char  VNAME[] = 	"valid";		/* valid vertex name */

#define  ABS(x)		((x)>=0 ? (x) : -(x))

#define  ZEROVECT(v)	(DOT(v,v) <= FTINY*FTINY)

#define  pvect(p)	printf(vformat, (p)[0], (p)[1], (p)[2])

char  vformat[] = "%15.9g %15.9g %15.9g\n";
char  tsargs[] = "4 surf_dx surf_dy surf_dz surf.cal\n";
char  texname[] = "Phong";

int  smooth = 0;		/* apply smoothing? */
int  objout = 0;		/* output .OBJ format? */

char  *modname, *surfname;

				/* recorded data flags */
#define  HASBORDER	01
#define  TRIPLETS	02
				/* a data structure */
struct {
	int	flags;			/* data type */
	short	m, n;			/* number of s and t values */
	RREAL	*data;			/* the data itself, s major sort */
} datarec;			/* our recorded data */

double  l_hermite(), l_bezier(), l_bspline(), l_dataval();
extern double  funvalue(), argument();

typedef struct {
	int  valid;	/* point is valid (vertex number) */
	FVECT  p;	/* vertex position */
	FVECT  n;	/* average normal */
	RREAL  uv[2];	/* (u,v) position */
} POINT;


main(argc, argv)
int  argc;
char  *argv[];
{
	extern long	eclock;
	POINT  *row0, *row1, *row2, *rp;
	int  i, j, m, n;
	char  stmp[256];

	varset("PI", ':', PI);
	funset("hermite", 5, ':', l_hermite);
	funset("bezier", 5, ':', l_bezier);
	funset("bspline", 5, ':', l_bspline);

	if (argc < 8)
		goto userror;

	for (i = 8; i < argc; i++)
		if (!strcmp(argv[i], "-e"))
			scompile(argv[++i], NULL, 0);
		else if (!strcmp(argv[i], "-f"))
			fcompile(argv[++i]);
		else if (!strcmp(argv[i], "-s"))
			smooth++;
		else if (!strcmp(argv[i], "-o"))
			objout++;
		else
			goto userror;

	modname = argv[1];
	surfname = argv[2];
	m = atoi(argv[6]);
	n = atoi(argv[7]);
	if (m <= 0 || n <= 0)
		goto userror;
	if (!strcmp(argv[5], "-") || access(argv[5], 4) == 0) {	/* file? */
		funset(ZNAME, 2, ':', l_dataval);
		if (!strcmp(argv[5],argv[3]) && !strcmp(argv[5],argv[4])) {
			loaddata(argv[5], m, n, 3);
			funset(XNAME, 2, ':', l_dataval);
			funset(YNAME, 2, ':', l_dataval);
		} else {
			loaddata(argv[5], m, n, 1);
			sprintf(stmp, "%s(s,t)=%s;", XNAME, argv[3]);
			scompile(stmp, NULL, 0);
			sprintf(stmp, "%s(s,t)=%s;", YNAME, argv[4]);
			scompile(stmp, NULL, 0);
		}
	} else {
		sprintf(stmp, "%s(s,t)=%s;", XNAME, argv[3]);
		scompile(stmp, NULL, 0);
		sprintf(stmp, "%s(s,t)=%s;", YNAME, argv[4]);
		scompile(stmp, NULL, 0);
		sprintf(stmp, "%s(s,t)=%s;", ZNAME, argv[5]);
		scompile(stmp, NULL, 0);
	}
	row0 = (POINT *)malloc((n+3)*sizeof(POINT));
	row1 = (POINT *)malloc((n+3)*sizeof(POINT));
	row2 = (POINT *)malloc((n+3)*sizeof(POINT));
	if (row0 == NULL || row1 == NULL || row2 == NULL) {
		fprintf(stderr, "%s: out of memory\n", argv[0]);
		quit(1);
	}
	row0++; row1++; row2++;
						/* print header */
	printhead(argc, argv);
	eclock = 0;
						/* initialize */
	comprow(-1.0/m, row0, n);
	comprow(0.0, row1, n);
	comprow(1.0/m, row2, n);
	compnorms(row0, row1, row2, n);
	if (objout) {
		printf("\nusemtl %s\n\n", modname);
		putobjrow(row1, n);
	}
						/* for each row */
	for (i = 0; i < m; i++) {
						/* compute next row */
		rp = row0;
		row0 = row1;
		row1 = row2;
		row2 = rp;
		comprow((double)(i+2)/m, row2, n);
		compnorms(row0, row1, row2, n);
		if (objout)
			putobjrow(row1, n);

		for (j = 0; j < n; j++) {
			int  orient = (j & 1);
							/* put polygons */
			if (!(row0[j].valid && row1[j+1].valid))
				orient = 1;
			else if (!(row1[j].valid && row0[j+1].valid))
				orient = 0;
			if (orient)
				putsquare(&row0[j], &row1[j],
						&row0[j+1], &row1[j+1]);
			else
				putsquare(&row1[j], &row1[j+1],
						&row0[j], &row0[j+1]);
		}
	}

	quit(0);

userror:
	fprintf(stderr, "Usage: %s material name ", argv[0]);
	fprintf(stderr, "x(s,t) y(s,t) z(s,t) m n [-s][-e expr][-f file]\n");
	quit(1);
}


loaddata(file, m, n, pointsize)		/* load point data from file */
char  *file;
int  m, n;
int  pointsize;
{
	FILE  *fp;
	char  word[64];
	register int  size;
	register RREAL  *dp;

	datarec.flags = HASBORDER;		/* assume border values */
	datarec.m = m+1;
	datarec.n = n+1;
	size = datarec.m*datarec.n*pointsize;
	if (pointsize == 3)
		datarec.flags |= TRIPLETS;
	dp = (RREAL *)malloc(size*sizeof(RREAL));
	if ((datarec.data = dp) == NULL) {
		fputs("Out of memory\n", stderr);
		exit(1);
	}
	if (!strcmp(file, "-")) {
		file = "<stdin>";
		fp = stdin;
	} else if ((fp = fopen(file, "r")) == NULL) {
		fputs(file, stderr);
		fputs(": cannot open\n", stderr);
		exit(1);
	}
	while (size > 0 && fgetword(word, sizeof(word), fp) != NULL) {
		if (!isflt(word)) {
			fprintf(stderr, "%s: garbled data value: %s\n",
					file, word);
			exit(1);
		}
		*dp++ = atof(word);
		size--;
	}
	if (size == (m+n+1)*pointsize) {	/* no border after all */
		dp = (RREAL *)realloc((void *)datarec.data,
				m*n*pointsize*sizeof(RREAL));
		if (dp != NULL)
			datarec.data = dp;
		datarec.flags &= ~HASBORDER;
		datarec.m = m;
		datarec.n = n;
		size = 0;
	}
	if (datarec.m < 2 || datarec.n < 2 || size != 0 ||
			fgetword(word, sizeof(word), fp) != NULL) {
		fputs(file, stderr);
		fputs(": bad number of data points\n", stderr);
		exit(1);
	}
	fclose(fp);
}


double
l_dataval(nam)				/* return recorded data value */
char  *nam;
{
	double  u, v;
	register int  i, j;
	register RREAL  *dp;
	double  d00, d01, d10, d11;
						/* compute coordinates */
	u = argument(1); v = argument(2);
	if (datarec.flags & HASBORDER) {
		i = u *= datarec.m-1;
		j = v *= datarec.n-1;
	} else {
		i = u = u*datarec.m - .5;
		j = v = v*datarec.n - .5;
	}
	if (i < 0) i = 0;
	else if (i > datarec.m-2) i = datarec.m-2;
	if (j < 0) j = 0;
	else if (j > datarec.n-2) j = datarec.n-2;
						/* compute value */
	if (datarec.flags & TRIPLETS) {
		dp = datarec.data + 3*(j*datarec.m + i);
		if (nam == ZNAME)
			dp += 2;
		else if (nam == YNAME)
			dp++;
		d00 = dp[0]; d01 = dp[3];
		dp += 3*datarec.m;
		d10 = dp[0]; d11 = dp[3];
	} else {
		dp = datarec.data + j*datarec.m + i;
		d00 = dp[0]; d01 = dp[1];
		dp += datarec.m;
		d10 = dp[0]; d11 = dp[1];
	}
						/* bilinear interpolation */
	return((j+1-v)*((i+1-u)*d00+(u-i)*d01)+(v-j)*((i+1-u)*d10+(u-i)*d11));
}


putobjrow(rp, n)			/* output vertex row to .OBJ */
register POINT  *rp;
int  n;
{
	static int	nverts = 0;

	for ( ; n-- >= 0; rp++) {
		if (!rp->valid)
			continue;
		fputs("v ", stdout);
		pvect(rp->p);
		if (smooth && !ZEROVECT(rp->n))
			printf("\tvn %.9g %.9g %.9g\n",
					rp->n[0], rp->n[1], rp->n[2]);
		printf("\tvt %.9g %.9g\n", rp->uv[0], rp->uv[1]);
		rp->valid = ++nverts;
	}
}


putsquare(p0, p1, p2, p3)		/* put out a square */
POINT  *p0, *p1, *p2, *p3;
{
	static int  nout = 0;
	FVECT  norm[4];
	int  axis;
	FVECT  v1, v2, vc1, vc2;
	int  ok1, ok2;
					/* compute exact normals */
	ok1 = (p0->valid && p1->valid && p2->valid);
	if (ok1) {
		VSUB(v1, p1->p, p0->p);
		VSUB(v2, p2->p, p0->p);
		fcross(vc1, v1, v2);
		ok1 = (normalize(vc1) != 0.0);
	}
	ok2 = (p1->valid && p2->valid && p3->valid);
	if (ok2) {
		VSUB(v1, p2->p, p3->p);
		VSUB(v2, p1->p, p3->p);
		fcross(vc2, v1, v2);
		ok2 = (normalize(vc2) != 0.0);
	}
	if (!(ok1 | ok2))
		return;
	if (objout) {			/* output .OBJ faces */
		int	p0n=0, p1n=0, p2n=0, p3n=0;
		if (smooth) {
			if (!ZEROVECT(p0->n))
				p0n = p0->valid;
			if (!ZEROVECT(p1->n))
				p1n = p1->valid;
			if (!ZEROVECT(p2->n))
				p2n = p2->valid;
			if (!ZEROVECT(p3->n))
				p3n = p3->valid;
		}
		if (ok1 & ok2 && fdot(vc1,vc2) >= 1.0-FTINY*FTINY) {
			printf("f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n",
					p0->valid, p0->valid, p0n,
					p1->valid, p1->valid, p1n,
					p3->valid, p3->valid, p3n,
					p2->valid, p2->valid, p2n);
			return;
		}
		if (ok1)
			printf("f %d/%d/%d %d/%d/%d %d/%d/%d\n",
					p0->valid, p0->valid, p0n,
					p1->valid, p1->valid, p1n,
					p2->valid, p2->valid, p2n);
		if (ok2)
			printf("f %d/%d/%d %d/%d/%d %d/%d/%d\n",
					p2->valid, p2->valid, p2n,
					p1->valid, p1->valid, p1n,
					p3->valid, p3->valid, p3n);
		return;
	}
					/* compute normal interpolation */
	axis = norminterp(norm, p0, p1, p2, p3);

					/* put out quadrilateral? */
	if (ok1 & ok2 && fdot(vc1,vc2) >= 1.0-FTINY*FTINY) {
		printf("\n%s ", modname);
		if (axis != -1) {
			printf("texfunc %s\n", texname);
			printf(tsargs);
			printf("0\n13\t%d\n", axis);
			pvect(norm[0]);
			pvect(norm[1]);
			pvect(norm[2]);
			fvsum(v1, norm[3], vc1, -0.5);
			fvsum(v1, v1, vc2, -0.5);
			pvect(v1);
			printf("\n%s ", texname);
		}
		printf("polygon %s.%d\n", surfname, ++nout);
		printf("0\n0\n12\n");
		pvect(p0->p);
		pvect(p1->p);
		pvect(p3->p);
		pvect(p2->p);
		return;
	}
					/* put out triangles? */
	if (ok1) {
		printf("\n%s ", modname);
		if (axis != -1) {
			printf("texfunc %s\n", texname);
			printf(tsargs);
			printf("0\n13\t%d\n", axis);
			pvect(norm[0]);
			pvect(norm[1]);
			pvect(norm[2]);
			fvsum(v1, norm[3], vc1, -1.0);
			pvect(v1);
			printf("\n%s ", texname);
		}
		printf("polygon %s.%d\n", surfname, ++nout);
		printf("0\n0\n9\n");
		pvect(p0->p);
		pvect(p1->p);
		pvect(p2->p);
	}
	if (ok2) {
		printf("\n%s ", modname);
		if (axis != -1) {
			printf("texfunc %s\n", texname);
			printf(tsargs);
			printf("0\n13\t%d\n", axis);
			pvect(norm[0]);
			pvect(norm[1]);
			pvect(norm[2]);
			fvsum(v2, norm[3], vc2, -1.0);
			pvect(v2);
			printf("\n%s ", texname);
		}
		printf("polygon %s.%d\n", surfname, ++nout);
		printf("0\n0\n9\n");
		pvect(p2->p);
		pvect(p1->p);
		pvect(p3->p);
	}
}


comprow(s, row, siz)			/* compute row of values */
double  s;
register POINT  *row;
int  siz;
{
	double  st[2];
	int  end;
	int  checkvalid;
	register int  i;
	
	if (smooth) {
		i = -1;			/* compute one past each end */
		end = siz+1;
	} else {
		if (s < -FTINY || s > 1.0+FTINY)
			return;
		i = 0;
		end = siz;
	}
	st[0] = s;
	checkvalid = (fundefined(VNAME) == 2);
	while (i <= end) {
		st[1] = (double)i/siz;
		if (checkvalid && funvalue(VNAME, 2, st) <= 0.0) {
			row[i].valid = 0;
			row[i].p[0] = row[i].p[1] = row[i].p[2] = 0.0;
			row[i].uv[0] = row[i].uv[1] = 0.0;
		} else {
			row[i].valid = 1;
			row[i].p[0] = funvalue(XNAME, 2, st);
			row[i].p[1] = funvalue(YNAME, 2, st);
			row[i].p[2] = funvalue(ZNAME, 2, st);
			row[i].uv[0] = st[0];
			row[i].uv[1] = st[1];
		}
		i++;
	}
}


compnorms(r0, r1, r2, siz)		/* compute row of averaged normals */
register POINT  *r0, *r1, *r2;
int  siz;
{
	FVECT  v1, v2;

	if (!smooth)			/* not needed if no smoothing */
		return;
					/* compute row 1 normals */
	while (siz-- >= 0) {
		if (!r1[0].valid)
			continue;
		if (!r0[0].valid) {
			if (!r2[0].valid) {
				r1[0].n[0] = r1[0].n[1] = r1[0].n[2] = 0.0;
				continue;
			}
			fvsum(v1, r2[0].p, r1[0].p, -1.0);
		} else if (!r2[0].valid)
			fvsum(v1, r1[0].p, r0[0].p, -1.0);
		else
			fvsum(v1, r2[0].p, r0[0].p, -1.0);
		if (!r1[-1].valid) {
			if (!r1[1].valid) {
				r1[0].n[0] = r1[0].n[1] = r1[0].n[2] = 0.0;
				continue;
			}
			fvsum(v2, r1[1].p, r1[0].p, -1.0);
		} else if (!r1[1].valid)
			fvsum(v2, r1[0].p, r1[-1].p, -1.0);
		else
			fvsum(v2, r1[1].p, r1[-1].p, -1.0);
		fcross(r1[0].n, v1, v2);
		normalize(r1[0].n);
		r0++; r1++; r2++;
	}
}


int
norminterp(resmat, p0, p1, p2, p3)	/* compute normal interpolation */
register FVECT  resmat[4];
POINT  *p0, *p1, *p2, *p3;
{
#define u  ((ax+1)%3)
#define v  ((ax+2)%3)

	register int  ax;
	MAT4  eqnmat;
	FVECT  v1;
	register int  i, j;

	if (!smooth)			/* no interpolation if no smoothing */
		return(-1);
					/* find dominant axis */
	VCOPY(v1, p0->n);
	fvsum(v1, v1, p1->n, 1.0);
	fvsum(v1, v1, p2->n, 1.0);
	fvsum(v1, v1, p3->n, 1.0);
	ax = ABS(v1[0]) > ABS(v1[1]) ? 0 : 1;
	ax = ABS(v1[ax]) > ABS(v1[2]) ? ax : 2;
					/* assign equation matrix */
	eqnmat[0][0] = p0->p[u]*p0->p[v];
	eqnmat[0][1] = p0->p[u];
	eqnmat[0][2] = p0->p[v];
	eqnmat[0][3] = 1.0;
	eqnmat[1][0] = p1->p[u]*p1->p[v];
	eqnmat[1][1] = p1->p[u];
	eqnmat[1][2] = p1->p[v];
	eqnmat[1][3] = 1.0;
	eqnmat[2][0] = p2->p[u]*p2->p[v];
	eqnmat[2][1] = p2->p[u];
	eqnmat[2][2] = p2->p[v];
	eqnmat[2][3] = 1.0;
	eqnmat[3][0] = p3->p[u]*p3->p[v];
	eqnmat[3][1] = p3->p[u];
	eqnmat[3][2] = p3->p[v];
	eqnmat[3][3] = 1.0;
					/* invert matrix (solve system) */
	if (!invmat4(eqnmat, eqnmat))
		return(-1);			/* no solution */
					/* compute result matrix */
	for (j = 0; j < 4; j++)
		for (i = 0; i < 3; i++)
			resmat[j][i] =	eqnmat[j][0]*p0->n[i] +
					eqnmat[j][1]*p1->n[i] +
					eqnmat[j][2]*p2->n[i] +
					eqnmat[j][3]*p3->n[i];
	return(ax);

#undef u
#undef v
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
