#ifndef lint
static const char	RCSid[] = "$Id: genbox.c,v 2.9 2021/04/07 21:13:52 greg Exp $";
#endif
/*
 *  genbox.c - generate a parallelepiped.
 *
 *     1/8/86
 */

#include  "rtio.h"
#include  <stdlib.h>
#include  <math.h>


char  let[]="0123456789._ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

char  *cmtype;		/* ppd material type */

char  *cname;		/* ppd name */

double  size[3];	/* ppd size */

double  bevel = 0.0;	/* bevel amount */

int  rounde = 0;	/* boolean true for round edges */

int  reverse = 0;	/* boolean true for reversed normals */



static void
vertex(int v)
{
	int  i;

	for (i = 0; i < 3; i++) {
		if (v & 010)
			printf("\t%18.12g", v & 01 ? size[i]-bevel : bevel);
		else
			printf("\t%18.12g", v & 01 ? size[i] : 0.0);
		v >>= 1;
	}
	fputc('\n', stdout);
}


static void
side(int a, int b, int c, int d)		/* generate a rectangular face */
{
	printf("\n%s polygon %s.%c%c%c%c\n", cmtype, cname,
			let[a], let[b], let[c], let[d]);
	printf("0\n0\n12\n");
	if (reverse) {
		vertex(d);
		vertex(c);
		vertex(b);
		vertex(a);
	} else {
		vertex(a);
		vertex(b);
		vertex(c);
		vertex(d);
	}
}


static void
corner(int a, int b, int c)			/* generate a triangular face */
{
	printf("\n%s polygon %s.%c%c%c\n", cmtype, cname,
			let[a], let[b], let[c]);
	printf("0\n0\n9\n");
	if (reverse) {
		vertex(c);
		vertex(b);
		vertex(a);
	} else {
		vertex(a);
		vertex(b);
		vertex(c);
	}
}


static void
cylinder(int v0, int v1)		/* generate a cylinder */
{
	printf("\n%s cylinder %s.%c%c\n", cmtype, cname, v0+'0', v1+'0');
	printf("0\n0\n7\n");
	vertex(v0);
	vertex(v1);
	printf("\t%18.12g\n", bevel);
}


static void
sphere(int v0)			/* generate a sphere */
{
	printf("\n%s sphere %s.%c\n", cmtype, cname, v0+'0');
	printf("0\n0\n4\n");
	vertex(v0);
	printf("\t%18.12g\n", bevel);
}


int
main(int argc, char *argv[])
{
	int  i;
	
	if (argc < 6)
		goto userr;

	cmtype = argv[1];
	cname = argv[2];
	size[0] = atof(argv[3]);
	size[1] = atof(argv[4]);
	size[2] = atof(argv[5]);
	if ((size[0] <= 0.0) | (size[1] <= 0.0) | (size[2] <= 0.0))
		goto userr;

	for (i = 6; i < argc; i++) {
		if (argv[i][0] != '-')
			goto userr;
		switch (argv[i][1]) {
		case 'i':
			reverse = 1;
			break;
		case 'r':
			rounde = 1;
			/* fall through */
		case 'b':
			bevel = atof(argv[++i]);
			if (bevel > 0.0)
				break;
			/* fall through on error */
		default:
			goto userr;
		}
	}
	if (rounde & reverse)
		fprintf(stderr, "%s: warning - option -i ignored with -r\n",
				argv[0]);

	fputs("# ", stdout);
	printargs(argc, argv, stdout);

	if (bevel > 0.0) {
					/* minor faces */
		side(051, 055, 054, 050);
		side(064, 066, 062, 060);
		side(032, 033, 031, 030);
		side(053, 052, 056, 057);
		side(065, 061, 063, 067);
		side(036, 034, 035, 037);
	}
	if (bevel > 0.0 && !rounde) {
					/* bevel faces */
		side(031, 051, 050, 030);
		side(060, 062, 032, 030);
		side(050, 054, 064, 060);
		side(034, 036, 066, 064);
		side(037, 057, 056, 036);
		side(052, 062, 066, 056);
		side(052, 053, 033, 032);
		side(057, 067, 063, 053);
		side(061, 031, 033, 063);
		side(065, 067, 037, 035);
		side(055, 051, 061, 065);
		side(034, 054, 055, 035);
					/* bevel corners */
		corner(030, 050, 060);
		corner(051, 031, 061);
		corner(032, 062, 052);
		corner(064, 054, 034);
		corner(036, 056, 066);
		corner(065, 035, 055);
		corner(053, 063, 033);
		corner(037, 067, 057);
	}
	if (bevel > 0.0 && rounde) {
					/* round edges */
		cylinder(070, 071);
		cylinder(070, 074);
		cylinder(070, 072);
		cylinder(073, 071);
		cylinder(073, 072);
		cylinder(073, 077);
		cylinder(075, 071);
		cylinder(075, 074);
		cylinder(075, 077);
		cylinder(076, 072);
		cylinder(076, 074);
		cylinder(076, 077);
					/* round corners */
		sphere(070);
		sphere(071);
		sphere(072);
		sphere(073);
		sphere(074);
		sphere(075);
		sphere(076);
		sphere(077);
	}
	if (bevel == 0.0) {
					/* only need major faces */
		side(1, 5, 4, 0);
		side(4, 6, 2, 0);
		side(2, 3, 1, 0);
		side(3, 2, 6, 7);
		side(5, 1, 3, 7);
		side(6, 4, 5, 7);
	}
	return(0);
userr:
	fprintf(stderr, "Usage: %s ", argv[0]);
	fprintf(stderr, "material name xsize ysize zsize ");
	fprintf(stderr, "[-i] [-b bevel | -r round]\n");
	return(1);
}
