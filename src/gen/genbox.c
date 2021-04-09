#ifndef lint
static const char	RCSid[] = "$Id: genbox.c,v 2.11 2021/04/09 01:48:20 greg Exp $";
#endif
/*
 *  genbox.c - generate a parallelepiped.
 *
 *     1/8/86
 */

#include  "rtio.h"
#include  "rtmath.h"
#include  "objutil.h"
#include  <stdlib.h>


char	*progname;

int	verbose = 0;

char  let[]="0123456789._ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

char	*cmtype;		/* ppd material type */

char	*cname;		/* ppd name */

double	size[3];	/* ppd size */

int	rounde = 0;	/* round edges? (#segments = 2^rounde) */

double	bevel = 0.0;	/* bevel amount or round edge radius */

int	rev = 0;	/* boolean true for reversed normals */

Scene	*obj = NULL;	/* save as .OBJ scene if not NULL */

int	vid[0100];	/* vertex ID's for .OBJ scene */


static int
overtex(int v)					/* index a .OBJ vertex */
{
	double	vpos[3];
	int	i;

	if (vid[v] >= 0)			/* already have this vertex? */
		return(vid[v]);
						/* else create new ID */
	for (i = 0; i < 3; i++)
		if (v>>i & 010)
			vpos[i] = (v>>i & 01)^rev ? size[i]-bevel : bevel;
		else
			vpos[i] = (v>>i & 01)^rev ? size[i] : 0.0;

	return(vid[v] = addVertex(obj, vpos[0], vpos[1], vpos[2]));
}


static int
onormal(int vid1, int vid2, int vid3)		/* index a .OBJ normal */
{
	double	*p1 = obj->vert[vid1].p;
	double	*p2 = obj->vert[vid2].p;
	double	*p3 = obj->vert[vid3].p;
	FVECT	sv1, sv2, nrm;

	VSUB(sv1, p2, p1);
	VSUB(sv2, p3, p2);
	VCROSS(nrm, sv1, sv2);

	return(addNormal(obj, nrm[0], nrm[1], nrm[2]));
}


static void
vertex(int v)					/* print a Radiance vertex */
{
	int  i;

	for (i = 0; i < 3; i++) {
		if (v & 010)
			printf("\t%18.12g", (v&01)^rev ? size[i]-bevel : bevel);
		else
			printf("\t%18.12g", (v&01)^rev ? size[i] : 0.0);
		v >>= 1;
	}
	fputc('\n', stdout);
}


static void
side(int a, int b, int c, int d)		/* generate a rectangular face */
{
	if (obj != NULL) {			/* working on .OBJ? */
		VNDX	quadv[4];
		memset(quadv, 0xff, sizeof(quadv));
		quadv[0][0] = overtex(a);
		quadv[1][0] = overtex(b);
		quadv[2][0] = overtex(c);
		quadv[3][0] = overtex(d);
		if (rounde)			/* add normal if rounded */
			quadv[0][2] = quadv[1][2] = quadv[2][2] = quadv[3][2]
				= onormal(quadv[0][0], quadv[1][0], quadv[2][0]);
		addFace(obj, quadv, 4);
		return;
	}
						/* Radiance output */
	printf("\n%s polygon %s.%c%c%c%c\n", cmtype, cname,
			let[a], let[b], let[c], let[d]);
	printf("0\n0\n12\n");
	vertex(a);
	vertex(b);
	vertex(c);
	vertex(d);
}


static void
corner(int a, int b, int c)			/* generate a triangular face */
{
	if (obj != NULL) {			/* working on .OBJ? */
		VNDX	triv[3];
		memset(triv, 0xff, sizeof(triv));
		triv[0][0] = overtex(a);
		triv[1][0] = overtex(b);
		triv[2][0] = overtex(c);
		addFace(obj, triv, 3);
		return;
	}
						/* Radiance output */
	printf("\n%s polygon %s.%c%c%c\n", cmtype, cname,
			let[a], let[b], let[c]);
	printf("0\n0\n9\n");
	vertex(a);
	vertex(b);
	vertex(c);
}


static void
cylinder(int v0, int v1)			/* generate a rounded edge */
{
	if (obj != NULL) {			/* segmenting for .OBJ? */
		const int	nsgn = 1 - 2*rev;
		const int	nseg = 1<<rounde;
		const double	astep = 0.5*PI/(double)nseg;
		const int	vid0 = overtex(v0);
		const int	vid1 = overtex(v1);
		FVECT		p0, p1;
		FVECT		axis[2], voff;
		int		previd[2];
		int		prenid;
		VNDX		quadv[4];
		double		coef[2];
		int		i, j;
						/* need to copy b/c realloc */
		VCOPY(p0, obj->vert[vid0].p);
		VCOPY(p1, obj->vert[vid1].p);

		memset(axis, 0, sizeof(axis));

		switch ((v0 ^ v1) & 07) {
		case 01:			/* edge along X-axis */
			axis[1][1] = bevel*(1. - 2.*(p1[1] < .5*size[1]));
			axis[0][2] = bevel*(1. - 2.*(p1[2] < .5*size[2]));
			break;
		case 02:			/* edge along Y-axis */
			axis[0][0] = bevel*(1. - 2.*(p1[0] < .5*size[0]));
			axis[1][2] = bevel*(1. - 2.*(p1[2] < .5*size[2]));
			break;
		case 04:			/* edge along Z-axis */
			axis[0][1] = bevel*(1. - 2.*(p1[1] < .5*size[1]));
			axis[1][0] = bevel*(1. - 2.*(p1[0] < .5*size[0]));
			break;
		}
		previd[0] = addVertex(obj, p0[0]+axis[0][0],
					p0[1]+axis[0][1], p0[2]+axis[0][2]);
		previd[1] = addVertex(obj, p1[0]+axis[0][0],
					p1[1]+axis[0][1], p1[2]+axis[0][2]);
		prenid = addNormal(obj, axis[0][0]*nsgn,
					axis[0][1]*nsgn, axis[0][2]*nsgn);
		for (i = 1; i <= nseg; i++) {
			memset(quadv, 0xff, sizeof(quadv));
			quadv[0][0] = previd[0];
			quadv[1][0] = previd[1];
			quadv[0][2] = quadv[1][2] = prenid;
			coef[0] = cos(astep*i);
			coef[1] = sin(astep*i);
			for (j = 0; j < 3; j++)
				voff[j] = coef[0]*axis[0][j] + coef[1]*axis[1][j];
			previd[0] = quadv[3][0]
				= addVertex(obj, p0[0]+voff[0],
						p0[1]+voff[1], p0[2]+voff[2]);
			previd[1] = quadv[2][0]
				= addVertex(obj, p1[0]+voff[0],
						p1[1]+voff[1], p1[2]+voff[2]);
			prenid = quadv[2][2] = quadv[3][2]
				= addNormal(obj, voff[0]*nsgn,
						voff[1]*nsgn, voff[2]*nsgn);
			addFace(obj, quadv, 4);
		}
		return;
	}
						/* Radiance output */
	printf("\n%s cylinder %s.%c%c\n", cmtype, cname, let[v0], let[v1]);
	printf("0\n0\n7\n");
	vertex(v0);
	vertex(v1);
	printf("\t%18.12g\n", bevel);
}


static void					/* recursive corner subdivision */
osubcorner(const FVECT orig, const FVECT c0, const FVECT c1, const FVECT c2, int lvl)
{
	if (lvl-- <= 0) {			/* reached terminal depth? */
		const int	nsgn = 1 - 2*rev;
		FVECT		vpos;
		VNDX		triv[3];	/* output smoothed triangle */
		VSUM(vpos, orig, c0, bevel);
		triv[0][0] = addVertex(obj, vpos[0], vpos[1], vpos[2]);
		triv[0][2] = addNormal(obj, c0[0]*nsgn, c0[1]*nsgn, c0[2]*nsgn);
		VSUM(vpos, orig, c1, bevel);
		triv[1][0] = addVertex(obj, vpos[0], vpos[1], vpos[2]);
		triv[1][2] = addNormal(obj, c1[0]*nsgn, c1[1]*nsgn, c1[2]*nsgn);
		VSUM(vpos, orig, c2, bevel);
		triv[2][0] = addVertex(obj, vpos[0], vpos[1], vpos[2]);
		triv[2][2] = addNormal(obj, c2[0]*nsgn, c2[1]*nsgn, c2[2]*nsgn);
		triv[0][1] = triv[1][1] = triv[2][1] = -1;
		addFace(obj, triv, 3);
	} else {				/* else subdivide 4 subcorners */
		FVECT	m01, m12, m20;
		VADD(m01, c0, c1); normalize(m01);
		VADD(m12, c1, c2); normalize(m12);
		VADD(m20, c2, c0); normalize(m20);
		osubcorner(orig, c0, m01, m20, lvl);
		osubcorner(orig, c1, m12, m01, lvl);
		osubcorner(orig, c2, m20, m12, lvl);
		osubcorner(orig, m01, m12, m20, lvl);
	}
}


static void
sphere(int v0)					/* generate a rounded corner */
{
	if (obj != NULL) {			/* segmenting for .OBJ? */
		FVECT	orig, cdir[3];
		int	i;
		memset(cdir, 0, sizeof(cdir));
		for (i = 0; i < 3; i++)
			cdir[i][i] = 2*((v0>>i & 01)^rev) - 1;
		switch (v0 & 07) {
			case 0:
			case 3:
			case 5:
			case 6:
				VCOPY(orig, cdir[0]);
				VCOPY(cdir[0], cdir[1]);
				VCOPY(cdir[1], orig);
		}
		i = overtex(v0);
		VCOPY(orig, obj->vert[i].p);	/* realloc remedy */
		osubcorner(orig, cdir[0], cdir[1], cdir[2], rounde);
		return;
	}
						/* Radiance output */
	printf("\n%s sphere %s.%c\n", cmtype, cname, let[v0]);
	printf("0\n0\n4\n");
	vertex(v0);
	printf("\t%18.12g\n", bevel);
}


int
main(int argc, char *argv[])
{
	int  nsegs = 1;
	int  objout = 0;
	int  i;

	progname = argv[0];

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
		case 'o':			/* requesting .OBJ output */
			objout = 1;
			break;
		case 'i':			/* invert surface normals */
			rev = 1;
			break;
		case 'r':			/* rounded edges/corners */
			rounde = 1;
			/* fall through */
		case 'b':			/* beveled edges */
			bevel = atof(argv[++i]);
			if (bevel <= 0.0)
				goto userr;
			break;
		case 'n':			/* #segments for rounding */
			nsegs = atoi(argv[++i]);
			if (nsegs <= 0)
				goto userr;
			break;
		default:
			goto userr;
		}
	}
	if ((objout|rev) & (nsegs==1))		/* default to 32 segments/edge */
		nsegs = 32;
	if (rounde) {				/* rounding edges/corners? */
		--nsegs;
		while ((nsegs >>= 1))		/* segmentation requested? */
			++rounde;
	}
	if (rounde > 1 || objout) {		/* save as .OBJ scene? */
		obj = newScene();
		setMaterial(obj, cmtype);
		setGroup(obj, cname);
		memset(vid, 0xff, sizeof(vid));
	}
	fputs("# ", stdout);			/* write command as comment */
	printargs(argc, argv, stdout);

	if (bevel > 0.0) {			/* minor faces */
		side(051, 055, 054, 050);
		side(064, 066, 062, 060);
		side(032, 033, 031, 030);
		side(053, 052, 056, 057);
		side(065, 061, 063, 067);
		side(036, 034, 035, 037);
	}
	if (bevel > 0.0 && !rounde) {		/* bevel faces */
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
	if (bevel > 0.0 && rounde) {		/* round edges */
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
	if (bevel == 0.0) {			/* only need major faces */
		side(1, 5, 4, 0);
		side(4, 6, 2, 0);
		side(2, 3, 1, 0);
		side(3, 2, 6, 7);
		side(5, 1, 3, 7);
		side(6, 4, 5, 7);
	}
	if (obj != NULL) {			/* need to write output? */
		if (objout) {
			coalesceVertices(obj, 2.*FTINY);
			if (toOBJ(obj, stdout) <= 0)
				return(1);
		} else if (toRadiance(obj, stdout, 0, 0) <= 0)
			return(1);
		/* freeScene(obj); 	we're exiting, anyway... */
	}
	return(0);
userr:
	fprintf(stderr, "Usage: %s ", argv[0]);
	fprintf(stderr, "material name xsize ysize zsize ");
	fprintf(stderr, "[-i] [-b bevel | -r round [-n nsegs]] [-o]\n");
	return(1);
}
