#ifndef lint
static const char	RCSid[] = "$Id: thf2rad.c,v 2.8 2003/11/15 17:54:06 schorsch Exp $";
#endif
/* Copyright (c) 1988 Regents of the University of California */

/*
 *  thf2rad - convert GDS things file to Radiance scene description.
 *
 *	12/21/88
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXSTR		128		/* maximum string or id length */
#define MAXPTS		2048		/* maximum points per object */

typedef struct {
	double	x, y, z;
} POINT;

double	rad = 0.0;			/* line radius */

static void thf2rad(char *name, char *file);
static int cvobject(char *name, char *file, FILE *fp);
static void header(int ac, char **av);
static void point(char *mat, char *nam, POINT *pos);
static void line(char *mat, char *nam, POINT *p1, POINT *p2);
static void loc(POINT *p);
static void readerr(char *file, char *err);


int
main(
	int	argc,
	char	*argv[]
)
{
	char	*name = NULL;
	int	an;

	for (an = 1; an < argc; an++)
		if (argv[an][0] == '-')
			switch (argv[an][1]) {
			case 'n':
				name = argv[++an];
				break;
			case 'r':
				rad = atof(argv[++an]);
				break;
			default:
				fprintf(stderr,
				"Usage: %s [-n name][-r rad] [file ..]\n",
						argv[0]);
				exit(1);
			}
		else
			break;
	header(argc, argv);
	if (an < argc) 
		for ( ; an < argc; an++)
			thf2rad(name, argv[an]);
	else
		thf2rad(name, NULL);
	exit(0);
}


void
thf2rad(		/* convert things file */
	char	*name,
	char	*file
)
{
	char	nambuf[MAXSTR];
	register char	*cp;
	FILE	*fp;
	int	on;
				/* open file and hack name */
	if (file == NULL) {
		fp = stdin;
		file = "<stdin>";
		if (name == NULL)
			name = "thing";
	} else {
		if ((fp = fopen(file, "r")) == NULL) {
			fprintf(stderr, "%s: cannot open\n", file);
			exit(1);
		}
		if (name == NULL)
			name = file;
	}
	for (cp = nambuf; (*cp = *name++); cp++)
		;
				/* get objects from file */
	on = 0;
	while (cvobject(nambuf, file, fp))
		sprintf(cp, ".%d", ++on);
}


int
cvobject(	/* convert next object in things file */
	char	*name,
	char	*file,
	FILE	*fp
)
{
	static POINT	parr[MAXPTS];
	static unsigned char	haspt[MAXPTS/8];
	char	sbuf[MAXSTR], nambuf[MAXSTR];
	int	npts, n, nv, p1, p2;
	register int	i;
						/* get points */
	if (fscanf(fp, "%s", sbuf) != 1)
		return(0);
	if (strcmp(sbuf, "YIN"))
		readerr(file, "not a GDS things file");
	if (fscanf(fp, "%d", &npts) != 1 || npts < 0 || npts > MAXPTS)
		readerr(file, "bad number of points");
	for (i = 0; i < npts; i++)
		if (fscanf(fp, "%lf %lf %lf",
				&parr[i].x, &parr[i].y, &parr[i].z) != 3)
			readerr(file, "bad point");
	for (i = (npts+7)/8; i >= 0; i--)
		haspt[i] = 0;
						/* get lines */
	if (fscanf(fp, "%d", &n) != 1 || n < 0)
		readerr(file, "bad number of lines");
	for (i = 0; i < n; i++) {
		if (fscanf(fp, "%d %d %s", &p1, &p2, sbuf) != 3)
			readerr(file, "error reading vector");
		if (--p1 < 0 || p1 >= npts || --p2 < 0 || p2 >= npts)
			readerr(file, "bad point in vector");
		if (!(haspt[p1>>3] & 1<<(p1&7))) {
			sprintf(nambuf, "%s.l%d.s0", name, i+1);
			point(sbuf, nambuf, &parr[p1]);
			haspt[p1>>3] |= 1<<(p1&7);
		}
		sprintf(nambuf, "%s.l%d.c", name, i+1);
		line(sbuf, nambuf, &parr[p1], &parr[p2]);
		if (!(haspt[p2>>3] & 1<<(p2&7))) {
			sprintf(nambuf, "%s.l%d.s1", name, i+1);
			point(sbuf, nambuf, &parr[p2]);
			haspt[p2>>3] |= 1<<(p2&7);
		}
	}
						/* get surfaces */
	if (fscanf(fp, "%d", &n) != 1 || n < 0)
		readerr(file, "bad number of surfaces");
	for (i = 0; i < n; i++) {
		if (fscanf(fp, "%d %s %d", &p1, sbuf, &nv) != 3)
			readerr(file, "error reading surface");
		if (p1 != 1 || nv < 0)
			readerr(file, "bad surface type (has hole?)");
		if (nv < 3) {
			while (nv--)
				fscanf(fp, "%*d");
			continue;
		}
		printf("\n%s polygon %s.p%d\n", sbuf, name, i+1);
		printf("0\n0\n%d\n", 3*nv);
		while (nv--) {
			if (fscanf(fp, "%d", &p1) != 1)
				readerr(file, "error reading surface vertex");
			loc(&parr[p1-1]);
		}
	}
	return(1);
}


void
header(			/* print header */
	int	ac,
	char	**av
)
{
	putchar('#');
	while (ac--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');
}


void
point(		/* print point */
	char	*mat,
	char	*nam,
	POINT	*pos
)
{
	if (rad <= 0.0)
		return;
	printf("\n%s sphere %s\n", mat, nam);
	printf("0\n0\n4\n");
	loc(pos);
	printf("\t%g\n", rad);
}


void
line(		/* print line */
	char	*mat,
	char	*nam,
	POINT	*p1,
	POINT	*p2
)
{
	if (rad <= 0.0)
		return;
	printf("\n%s cylinder %s\n", mat, nam);
	printf("0\n0\n7\n");
	loc(p1);
	loc(p2);
	printf("\t%g\n", rad);
}


void
loc(				/* print location */
	register POINT	*p
)
{
	printf("\t%14.8g\t%14.8g\t%14.8g\n", p->x, p->y, p->z);
}


void
readerr(		/* print read error and exit */
	char	*file,
	char	*err
)
{
	fprintf(stderr, "%s: %s\n", file, err);
	exit(1);
}
