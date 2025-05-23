#ifndef lint
static const char	RCSid[] = "$Id: mgf2rad.c,v 2.36 2025/05/23 17:02:07 greg Exp $";
#endif
/*
 * Convert MGF (Materials and Geometry Format) to Radiance
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "platform.h"
#include "mgf_parser.h"
#include "color.h"
#include "tmesh.h"
#include "lookup.h"

#define putv(v)		printf("%18.12g %18.12g %18.12g\n",(v)[0],(v)[1],(v)[2])

#define invert		(xf_context != NULL && xf_context->rev)

#define	SGEN_DEF	"spec*"
#define SGEN_RS		"rs_spec*"
#define SGEN_TD		"td_spec*"
#define SGEN_TS		"ts_spec*"

char	void_str[] = "void";		/* global VOIDID */
char	sgen_str[16] = SGEN_DEF;	/* generic specular */

double	glowdist = FHUGE;		/* glow test distance */

double  emult = 1.;			/* emitter multiplier */

FILE	*matfp;				/* material output file */

int	dospectra = 0;			/* output spectral colors? */


extern int r_comment(int ac, char **av);
extern int r_color(int ac, char **av);
extern int r_cone(int ac, char **av);
extern int r_cyl(int ac, char **av);
extern int r_sph(int ac, char **av);
extern int r_ring(int ac, char **av);
extern int r_face(int ac, char **av);
extern int r_ies(int ac, char **av);
extern void putsided(char *mname);
extern char * material(void);
extern char * object(void);
extern char * addarg(char *op, char *arg);
extern void do_tri(char *mat, C_VERTEX *cv1, C_VERTEX *cv2, C_VERTEX *cv3, int iv);
extern void cvtcolor(COLOR radrgb, C_COLOR *ciec, double intensity);
extern int color_clash(int e1, int e2);
extern int isgrey(COLOR rgb);
extern void putrgbpat(char *pnm, COLOR rgb);
extern char * specolor(COLOR radrgb, C_COLOR *ciec, double intensity);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	i;

	matfp = stdout;
				/* print out parser version */
	printf("## Translated from MGF Version %d.%d\n", MG_VMAJOR, MG_VMINOR);
				/* initialize dispatch table */
	mg_ehand[MG_E_COMMENT] = r_comment;	/* we pass comments */
	mg_ehand[MG_E_COLOR] = c_hcolor;	/* they get color */
	mg_ehand[MG_E_CONE] = r_cone;		/* we do cones */
	mg_ehand[MG_E_CMIX] = c_hcolor;		/* they mix colors */
	mg_ehand[MG_E_CXY] = c_hcolor;		/* they get chromaticities */
	mg_ehand[MG_E_CSPEC] = r_color;		/* we get spectra */
	mg_ehand[MG_E_CCT] = r_color;		/* we get color temp's */
	mg_ehand[MG_E_CYL] = r_cyl;		/* we do cylinders */
	mg_ehand[MG_E_ED] = c_hmaterial;	/* they get emission */
	mg_ehand[MG_E_FACE] = r_face;		/* we do faces */
	mg_ehand[MG_E_IES] = r_ies;		/* we do IES files */
	mg_ehand[MG_E_IR] = c_hmaterial;	/* they get refractive index */
	mg_ehand[MG_E_MATERIAL] = c_hmaterial;	/* they get materials */
	mg_ehand[MG_E_NORMAL] = c_hvertex;	/* they get normals */
	mg_ehand[MG_E_OBJECT] = obj_handler;	/* they track object names */
	mg_ehand[MG_E_POINT] = c_hvertex;	/* they get points */
	mg_ehand[MG_E_RD] = c_hmaterial;	/* they get diffuse refl. */
	mg_ehand[MG_E_RING] = r_ring;		/* we do rings */
	mg_ehand[MG_E_RS] = c_hmaterial;	/* they get specular refl. */
	mg_ehand[MG_E_SIDES] = c_hmaterial;	/* they get # sides */
	mg_ehand[MG_E_SPH] = r_sph;		/* we do spheres */
	mg_ehand[MG_E_TD] = c_hmaterial;	/* they get diffuse trans. */
	mg_ehand[MG_E_TS] = c_hmaterial;	/* they get specular trans. */
	mg_ehand[MG_E_VERTEX] = c_hvertex;	/* they get vertices */
	mg_ehand[MG_E_XF] = xf_handler;		/* they track transforms */
	mg_init();		/* initialize the parser */
					/* get our options & print header */
	printf("## %s", argv[0]);
	for (i = 1; i < argc && argv[i][0] == '-'; i++) {
		printf(" %s", argv[i]);
		switch (argv[i][1]) {
		case 'g':			/* glow distance (meters) */
			if (argv[i][2] || badarg(argc-i-1, argv+i+1, "f"))
				goto userr;
			glowdist = atof(argv[++i]);
			printf(" %s", argv[i]);
			break;
		case 'e':			/* emitter multiplier */
			if (argv[i][2] || badarg(argc-i-1, argv+i+1, "f"))
				goto userr;
			emult = atof(argv[++i]);
			printf(" %s", argv[i]);
			break;
		case 'm':			/* materials file */
			matfp = fopen(argv[++i], "a");
			if (matfp == NULL) {
				fprintf(stderr, "%s: cannot append\n", argv[i]);
				exit(1);
			}
			printf(" %s", argv[i]);
			break;
		case 's':			/* spectral color output? */
			dospectra = !dospectra;
			break;
		default:
			goto userr;
		}
	}
	putchar('\n');
	if (i == argc) {		/* convert stdin */
		if (mg_load(NULL) != MG_OK)
			exit(1);
		if (mg_nunknown)
			printf("## %s: %u unknown entities\n",
					argv[0], mg_nunknown);
	} else				/* convert each file */
		for ( ; i < argc; i++) {
			printf("## %s %s ##############################\n",
					argv[0], argv[i]);
			if (mg_load(argv[i]) != MG_OK)
				exit(1);
			if (mg_nunknown) {
				printf("## %s %s: %u unknown entities\n",
						argv[0], argv[i], mg_nunknown);
				mg_nunknown = 0;
			}
		}
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-s][-g dist][-e mult][-m matf] [file.mgf] ..\n",
			argv[0]);
	exit(1);
}


int
r_comment(		/* repeat a comment verbatim */
	int	ac,
	char	**av
)
{
	putchar('#');		/* use Radiance comment character */
	while (--ac) {			/* pass through verbatim */
		putchar(' ');
		fputs(*++av, stdout);
	}
	putchar('\n');
	return(MG_OK);
}


int
r_color(		/* call color handler & remember name */
	int	ac,
	char	**av
)
{
	int	rval = c_hcolor(ac, av);

	if (rval == MG_OK)
		c_ccolor->client_data = c_ccname;

	return(rval);
}


int
r_cone(			/* put out a cone */
	int	ac,
	char	**av
)
{
	static int	ncones;
	char	*mat;
	double	r1, r2;
	C_VERTEX	*cv1, *cv2;
	FVECT	p1, p2;
	int	inv;
					/* check argument count and type */
	if (ac != 5)
		return(MG_EARGC);
	if (!isflt(av[2]) || !isflt(av[4]))
		return(MG_ETYPE);
					/* get the endpoint vertices */
	if ((cv1 = c_getvert(av[1])) == NULL ||
			(cv2 = c_getvert(av[3])) == NULL)
		return(MG_EUNDEF);
	xf_xfmpoint(p1, cv1->p);	/* transform endpoints */
	xf_xfmpoint(p2, cv2->p);
	r1 = xf_scale(atof(av[2]));	/* scale radii */
	r2 = xf_scale(atof(av[4]));
	inv = r1 < 0.;			/* check for inverted cone */
	if (r1 == 0.) {			/* check for illegal radii */
		if (r2 == 0.)
			return(MG_EILL);
		inv = r2 < 0.;
	} else if (r2 != 0. && inv ^ (r2 < 0.))
		return(MG_EILL);
	if (inv) {
		r1 = -r1;
		r2 = -r2;
	}
	if ((mat = material()) == NULL)	/* get material */
		return(MG_EBADMAT);
					/* spit the sucker out */
	printf("\n%s %s %sc%d\n", mat, inv ? "cup" : "cone",
			object(), ++ncones);
	printf("0\n0\n8\n");
	putv(p1);
	putv(p2);
	printf("%18.12g %18.12g\n", r1, r2);
	return(MG_OK);
}


int
r_cyl(			/* put out a cylinder */
	int	ac,
	char	**av
)
{
	static int	ncyls;
	char	*mat;
	double	rad;
	C_VERTEX	*cv1, *cv2;
	FVECT	p1, p2;
	int	inv;
					/* check argument count and type */
	if (ac != 4)
		return(MG_EARGC);
	if (!isflt(av[2]))
		return(MG_ETYPE);
					/* get the endpoint vertices */
	if ((cv1 = c_getvert(av[1])) == NULL ||
			(cv2 = c_getvert(av[3])) == NULL)
		return(MG_EUNDEF);
	xf_xfmpoint(p1, cv1->p);	/* transform endpoints */
	xf_xfmpoint(p2, cv2->p);
	rad = xf_scale(atof(av[2]));	/* scale radius */
	if ((inv = rad < 0.))		/* check for inverted cylinder */
		rad = -rad;
	if ((mat = material()) == NULL)	/* get material */
		return(MG_EBADMAT);
					/* spit out the primitive */
	printf("\n%s %s %scy%d\n", mat, inv ? "tube" : "cylinder",
			object(), ++ncyls);
	printf("0\n0\n7\n");
	putv(p1);
	putv(p2);
	printf("%18.12g\n", rad);
	return(MG_OK);
}


int
r_sph(			/* put out a sphere */
	int	ac,
	char	**av
)
{
	static int	nsphs;
	char	*mat;
	double	rad;
	C_VERTEX	*cv;
	FVECT	cent;
	int	inv;
					/* check argument count and type */
	if (ac != 3)
		return(MG_EARGC);
	if (!isflt(av[2]))
		return(MG_ETYPE);
	if ((cv = c_getvert(av[1])) == NULL)	/* get center vertex */
		return(MG_EUNDEF);
	xf_xfmpoint(cent, cv->p);		/* transform center */
	rad = xf_scale(atof(av[2]));		/* scale radius */
	if ((inv = rad < 0.))			/* check for inversion */
		rad = -rad;
	if ((mat = material()) == NULL)		/* get material */
		return(MG_EBADMAT);
						/* spit out primitive */
	printf("\n%s %s %ss%d\n", mat, inv ? "bubble" : "sphere",
			object(), ++nsphs);
	printf("0\n0\n4 %18.12g %18.12g %18.12g %18.12g\n",
			cent[0], cent[1], cent[2], rad);
	return(MG_OK);
}


int
r_ring(			/* put out a ring */
	int	ac,
	char	**av
)
{
	static int	nrings;
	char	*mat;
	double	r1, r2;
	C_VERTEX	*cv;
	FVECT	cent, norm;
					/* check argument count and type */
	if (ac != 4)
		return(MG_EARGC);
	if (!isflt(av[2]) || !isflt(av[3]))
		return(MG_ETYPE);
	if ((cv = c_getvert(av[1])) == NULL)	/* get center vertex */
		return(MG_EUNDEF);
	if (is0vect(cv->n))			/* make sure we have normal */
		return(MG_EILL);
	xf_xfmpoint(cent, cv->p);		/* transform center */
	xf_rotvect(norm, cv->n);		/* rotate normal */
	r1 = xf_scale(atof(av[2]));		/* scale radii */
	r2 = xf_scale(atof(av[3]));
	if ((r1 < 0.) | (r2 <= r1))
		return(MG_EILL);
	if ((mat = material()) == NULL)		/* get material */
		return(MG_EBADMAT);
						/* spit out primitive */
	printf("\n%s ring %sr%d\n", mat, object(), ++nrings);
	printf("0\n0\n8\n");
	putv(cent);
	putv(norm);
	printf("%18.12g %18.12g\n", r1, r2);
	return(MG_OK);
}


int
r_face(			/* convert a face */
	int	ac,
	char	**av
)
{
	static int	nfaces;
	int		myi = invert;
	char	*mat;
	int	i;
	C_VERTEX	*cv;
	FVECT	v;

					/* check argument count and type */
	if (ac < 4)
		return(MG_EARGC);
	if ((mat = material()) == NULL)	/* get material */
		return(MG_EBADMAT);
	if (ac <= 5) {				/* check for smoothing */
		C_VERTEX	*cva[5];
		for (i = 1; i < ac; i++) {
			if ((cva[i-1] = c_getvert(av[i])) == NULL)
				return(MG_EUNDEF);
			if (is0vect(cva[i-1]->n))
				break;
		}
		if (i < ac)
			i = ISFLAT;
		else
			i = flat_tri(cva[0]->p, cva[1]->p, cva[2]->p,
					cva[0]->n, cva[1]->n, cva[2]->n);
		if (i == DEGEN)
			return(MG_OK);		/* degenerate (error?) */
		if (i == RVBENT) {
			myi = !myi;
			i = ISBENT;
		} else if (i == RVFLAT) {
			myi = !myi;
			i = ISFLAT;
		}
		if (i == ISBENT) {		/* smoothed triangles */
			do_tri(mat, cva[0], cva[1], cva[2], myi);
			if (ac == 5)
				do_tri(mat, cva[2], cva[3], cva[0], myi);
			return(MG_OK);
		}
	}
					/* spit out unsmoothed primitive */
	printf("\n%s polygon %sf%d\n", mat, object(), ++nfaces);
	printf("0\n0\n%d\n", 3*(ac-1));
	for (i = 1; i < ac; i++) {	/* get, transform, print each vertex */
		if ((cv = c_getvert(av[myi ? ac-i : i])) == NULL)
			return(MG_EUNDEF);
		xf_xfmpoint(v, cv->p);
		putv(v);
	}
	return(MG_OK);
}


int
r_ies(				/* convert an IES luminaire file */
	int	ac,
	char	**av
)
{
	int	xa0 = 2;
	char	combuf[128];
	char	fname[48];
	char	*oname;
	char	*op;
	int	i;
					/* check argument count */
	if (ac < 2)
		return(MG_EARGC);
					/* construct output file name */
	if ((op = strrchr(av[1], '/')) != NULL)
		op++;
	else
		op = av[1];
	(void)strcpy(fname, op);
	if ((op = strrchr(fname, '.')) == NULL)
		op = fname + strlen(fname);
	(void)strcpy(op, ".rad");
					/* see if we need to run ies2rad */
	if (access(fname, 0) == -1) {
		(void)strcpy(combuf, "ies2rad");/* build ies2rad command */
		op = combuf + 7;		/* get -m option (first) */
		if (ac-xa0 >= 2 && !strcmp(av[xa0], "-m")) {
			if (!isflt(av[xa0+1]))
				return(MG_ETYPE);
			op = addarg(addarg(op, "-m"), av[xa0+1]);
			xa0 += 2;
		}
		*op++ = ' ';			/* build IES filename */
		i = 0;
		if (mg_file != NULL &&
				(oname = strrchr(mg_file->fname,'/')) != NULL) {
			i = oname - mg_file->fname + 1;
			(void)strcpy(op, mg_file->fname);
		}
		(void)strcpy(op+i, av[1]);
		if (access(op, 0) == -1)	/* check for file existence */
			return(MG_ENOFILE);
		system(combuf);			/* run ies2rad */
		if (access(fname, 0) == -1)	/* check success */
			return(MG_EINCL);
	}
	printf("\n!xform");			/* put out xform command */
	oname = object();
	if (*oname) {
		printf(" -n ");
		for (op = oname; op[1]; op++)	/* remove trailing separator */
			putchar(*op);
	}
	for (i = xa0; i < ac; i++)
		printf(" %s", av[i]);
	if (ac > xa0 && xf_argc > 0)
		printf(" -i 1");
	for (i = 0; i < xf_argc; i++)
		printf(" %s", xf_argv[i]);
	printf(" %s\n", fname);
	return(MG_OK);
}


void
do_tri(		/* put out smoothed triangle */
	char	*mat,
	C_VERTEX	*cv1,
	C_VERTEX	*cv2,
	C_VERTEX	*cv3,
	int	iv
)
{
	static int	ntris;
	BARYCCM	bvecs;
	RREAL	bcoor[3][3];
	C_VERTEX	*cvt;
	FVECT	v1, v2, v3;
	FVECT	n1, n2, n3;
	int	i;

	if (iv) {			/* swap vertex order if inverted */
		cvt = cv1;
		cv1 = cv3;
		cv3 = cvt;
	}
	xf_xfmpoint(v1, cv1->p);
	xf_xfmpoint(v2, cv2->p);
	xf_xfmpoint(v3, cv3->p);
					/* compute barycentric coords. */
	if (comp_baryc(&bvecs, v1, v2, v3) < 0)
		return;				/* degenerate triangle! */
	printf("\n%s texfunc T-nor\n", mat);	/* put out texture */
	printf("4 dx dy dz %s\n0\n", TCALNAME);
	xf_rotvect(n1, cv1->n);
	xf_rotvect(n2, cv2->n);
	xf_rotvect(n3, cv3->n);
	for (i = 0; i < 3; i++) {
		bcoor[i][0] = n1[i];
		bcoor[i][1] = n2[i];
		bcoor[i][2] = n3[i];
	}
	fput_baryc(&bvecs, bcoor, 3, stdout);
						/* put out triangle */
	printf("\nT-nor polygon %st%d\n", object(), ++ntris);
	printf("0\n0\n9\n");
	putv(v1);
	putv(v2);
	putv(v3);
}


void
putsided(char *mname)		/* print out mixfunc for sided material */
{
	fprintf(matfp, "\nvoid mixfunc %s\n", mname);
	fprintf(matfp, "4 %s void if(Rdot,1,0) .\n0\n0\n", mname);
}


char *
material(void)			/* get (and print) current material */
{
	char	*mname = "mat";
	C_COLOR	*refclr = NULL;
	char	*pname;
	COLOR	radrgb, c2;
	double	d;

	if (c_cmname != NULL)
		mname = c_cmname;
	if (!c_cmaterial->clock)
		return(mname);		/* already current */
				/* else update output */
	c_cmaterial->clock = 0;
	if (c_cmaterial->ed > .1) {	/* emitter */
		pname = specolor(radrgb, &c_cmaterial->ed_c,
				emult*c_cmaterial->ed/(PI*WHTEFFICACY));
		if (glowdist < FHUGE) {		/* do a glow */
			fprintf(matfp, "\n%s glow %s\n0\n0\n", pname, mname);
			fprintf(matfp, "4 %f %f %f %f\n", colval(radrgb,RED),
					colval(radrgb,GRN),
					colval(radrgb,BLU), glowdist);
		} else {
			fprintf(matfp, "\n%s light %s\n0\n0\n", pname, mname);
			fprintf(matfp, "3 %f %f %f\n", colval(radrgb,RED),
					colval(radrgb,GRN),
					colval(radrgb,BLU));
		}
		return(mname);
	}
	d = c_cmaterial->rd + c_cmaterial->td +
			c_cmaterial->rs + c_cmaterial->ts;
	if ((d < 0.) | (d > 1.))
		return(NULL);
					/* check for glass/dielectric */
	if (c_cmaterial->nr > 1.1 &&
			c_cmaterial->ts > .25 && c_cmaterial->rs <= .125 &&
			c_cmaterial->td <= .01 && c_cmaterial->rd <= .01 &&
			c_cmaterial->rs_a <= .01 && c_cmaterial->ts_a <= .01) {
		cvtcolor(radrgb, &c_cmaterial->ts_c,
				c_cmaterial->ts + c_cmaterial->rs);
		if (c_cmaterial->sided) {		/* dielectric */
			colval(radrgb,RED) = pow(colval(radrgb,RED),
							1./C_1SIDEDTHICK);
			colval(radrgb,GRN) = pow(colval(radrgb,GRN),
							1./C_1SIDEDTHICK);
			colval(radrgb,BLU) = pow(colval(radrgb,BLU),
							1./C_1SIDEDTHICK);
			fprintf(matfp, "\nvoid dielectric %s\n0\n0\n", mname);
			fprintf(matfp, "5 %g %g %g %f 0\n", colval(radrgb,RED),
					colval(radrgb,GRN), colval(radrgb,BLU),
					c_cmaterial->nr);
			return(mname);
		}
							/* glass */
		fprintf(matfp, "\nvoid glass %s\n0\n0\n", mname);
		fprintf(matfp, "4 %f %f %f %f\n", colval(radrgb,RED),
				colval(radrgb,GRN), colval(radrgb,BLU),
				c_cmaterial->nr);
		return(mname);
	}
					/* check for WGMDfunc */
	if (((c_cmaterial->rs > .02) & (c_cmaterial->ts > .02) &&
			fabs(c_cmaterial->rs_a - c_cmaterial->ts_a) > .02) ||
			((c_cmaterial->rs > .05) & (c_cmaterial->rd > .05) &&
				!c_isgrey(&c_cmaterial->rs_c) &&
				!c_equiv(&c_cmaterial->rd_c, &c_cmaterial->rs_c)) ||
			color_clash(MG_E_TS, MG_E_TD) ||
			color_clash(MG_E_TD, MG_E_RD) ||
			color_clash(MG_E_RD, MG_E_TS)) {
		COLOR	rs_rgb, ts_rgb, td_rgb;	/* separate modifier paths */
		char	rs_pname[128], ts_pname[128], td_pname[128];
		strcpy(sgen_str, SGEN_RS);
		strcpy(rs_pname, specolor(rs_rgb, &c_cmaterial->rs_c, c_cmaterial->rs));
		strcpy(sgen_str, SGEN_TS);
		strcpy(ts_pname, specolor(ts_rgb, &c_cmaterial->ts_c, c_cmaterial->ts));
		strcpy(sgen_str, SGEN_TD);
		strcpy(td_pname, specolor(td_rgb, &c_cmaterial->td_c, c_cmaterial->td));
		strcpy(sgen_str, SGEN_DEF);
		pname = specolor(radrgb, &c_cmaterial->rd_c, c_cmaterial->rd);
		if (!strcmp(rs_pname, void_str) && !isgrey(rs_rgb)) {
			putrgbpat(strcpy(rs_pname,"rs_rgb*"), rs_rgb);
			colval(rs_rgb,GRN) = 1;
		}
		if (!strcmp(ts_pname, void_str) && !isgrey(ts_rgb)) {
			putrgbpat(strcpy(ts_pname,"ts_rgb*"), ts_rgb);
			colval(ts_rgb,GRN) = 1;
		}
		fprintf(matfp, "\n%s WGMDfunc %s\n", pname, mname);
		fprintf(matfp, "13\t%s %f %f %f\n", rs_pname, colval(rs_rgb,GRN),
				c_cmaterial->rs_a, c_cmaterial->rs_a);
		fprintf(matfp, "\t%s %f %f %f\n", ts_pname, colval(ts_rgb,GRN),
				c_cmaterial->ts_a, c_cmaterial->ts_a);
		fprintf(matfp, "\t%s\n\t0 0 1 .\n0\n", td_pname);
		fprintf(matfp, "9\t%f %f %f\n", colval(radrgb,RED),
				colval(radrgb,GRN), colval(radrgb,BLU));
		fprintf(matfp, "\t%f %f %f\n", colval(radrgb,RED),
				colval(radrgb,GRN), colval(radrgb,BLU));
		fprintf(matfp, "\t%f %f %f\n", colval(td_rgb,RED),
				colval(td_rgb,GRN), colval(td_rgb,BLU));
		if (c_cmaterial->sided)
			putsided(mname);
		return(mname);
	}
					/* check for trans */
	if (c_cmaterial->td + c_cmaterial->ts > .01) {
		double	a5, a6;
						/* average colors */
		d = c_cmaterial->rd + c_cmaterial->td + c_cmaterial->ts;
		cvtcolor(radrgb, &c_cmaterial->rd_c, c_cmaterial->rd/d);
		cvtcolor(c2, &c_cmaterial->td_c, c_cmaterial->td/d);
		addcolor(radrgb, c2);
		cvtcolor(c2, &c_cmaterial->ts_c, c_cmaterial->ts/d);
		addcolor(radrgb, c2);
		if (c_cmaterial->rs + c_cmaterial->ts > .0001)
			a5 = (c_cmaterial->rs*c_cmaterial->rs_a +
					c_cmaterial->ts*c_cmaterial->ts_a) /
					(c_cmaterial->rs + c_cmaterial->ts);
		a6 = (c_cmaterial->td + c_cmaterial->ts) /
				(c_cmaterial->rd + c_cmaterial->td + c_cmaterial->ts);
		if (a6 < .999)
			d = c_cmaterial->rd/(1. - c_cmaterial->rs)/(1. - a6);
		else
			d = c_cmaterial->td + c_cmaterial->ts;
		scalecolor(radrgb, d);
		fprintf(matfp, "\nvoid trans %s\n0\n0\n", mname);
		fprintf(matfp, "7 %f %f %f\n", colval(radrgb,RED),
				colval(radrgb,GRN), colval(radrgb,BLU));
		fprintf(matfp, "\t%f %f %f %f\n", c_cmaterial->rs, a5, a6,
				c_cmaterial->ts/(c_cmaterial->ts + c_cmaterial->td));
		if (c_cmaterial->sided)
			putsided(mname);
		return(mname);
	}
					/* check for plastic */
	if (c_cmaterial->rs < .1 && (c_cmaterial->rs < .1*c_cmaterial->rd ||
					c_isgrey(&c_cmaterial->rs_c))) {
		pname = specolor(radrgb, &c_cmaterial->rd_c,
					c_cmaterial->rd/(1.-c_cmaterial->rs));
		fprintf(matfp, "\n%s plastic %s\n0\n0\n", pname, mname);
		fprintf(matfp, "5 %f %f %f %f %f\n", colval(radrgb,RED),
				colval(radrgb,GRN), colval(radrgb,BLU),
				c_cmaterial->rs, c_cmaterial->rs_a);
		if (c_cmaterial->sided)
			putsided(mname);
		return(mname);
	}
					/* else it's metal */
						/* compute color */
	if (c_equiv(&c_cmaterial->rd_c, &c_cmaterial->rs_c)) {
		pname = specolor(radrgb, &c_cmaterial->rs_c, c_cmaterial->rs+c_cmaterial->rd);
	} else if (c_cmaterial->rd <= .05f) {
		pname = specolor(radrgb, &c_cmaterial->rs_c, c_cmaterial->rs);
		cvtcolor(c2, &c_cmaterial->rd_c, c_cmaterial->rd);
		addcolor(radrgb, c2);
	} else {
		pname = "void";
		cvtcolor(radrgb, &c_cmaterial->rd_c, c_cmaterial->rd);
		cvtcolor(c2, &c_cmaterial->rs_c, c_cmaterial->rs);
		addcolor(radrgb, c2);
	}
	fprintf(matfp, "\n%s metal %s\n0\n0\n", pname, mname);
	fprintf(matfp, "5 %f %f %f %f %f\n", colval(radrgb,RED),
			colval(radrgb,GRN), colval(radrgb,BLU),
			c_cmaterial->rs/(c_cmaterial->rd + c_cmaterial->rs),
			c_cmaterial->rs_a);
	if (c_cmaterial->sided)
		putsided(mname);
	return(mname);
}


void
cvtcolor(	/* convert a CIE XYZ color to RGB */
	COLOR	radrgb,
	C_COLOR	*ciec,
	double	intensity
)
{
	COLOR	ciexyz;

	if (intensity <= 0) {
		setcolor(radrgb, 0, 0, 0);
		return;
	}
	c_ccvt(ciec, C_CSXY);		/* get xy representation */
	ciexyz[1] = intensity;
	ciexyz[0] = ciec->cx/ciec->cy*ciexyz[1];
	ciexyz[2] = ciexyz[1]*(1./ciec->cy - 1.) - ciexyz[0];
	cie_rgb(radrgb, ciexyz);
}


int color_clash(	/* do non-zero material components clash? */
	int e1,
	int e2
)
{
	C_COLOR	*c1;

	if (e1 == e2)
		return(0);
	switch (e1) {
	case MG_E_RD:
		if (c_cmaterial->rd <= .01) return(0);
		c1 = &c_cmaterial->rd_c;
		break;
	case MG_E_RS:
		if (c_cmaterial->rs <= .01) return(0);
		c1 = &c_cmaterial->rs_c;
		break;
	case MG_E_TD:
		if (c_cmaterial->td <= .01) return(0);
		c1 = &c_cmaterial->td_c;
		break;
	case MG_E_TS:
		if (c_cmaterial->ts <= .01) return(0);
		c1 = &c_cmaterial->ts_c;
		break;
	default:
		return(0);
	}
	switch (e2) {
	case MG_E_RD:
		if (c_cmaterial->rd <= .01) return(0);
		return(!c_equiv(c1, &c_cmaterial->rd_c));
	case MG_E_RS:
		if (c_cmaterial->rs <= .01) return(0);
		return(!c_equiv(c1, &c_cmaterial->rs_c));
	case MG_E_TD:
		if (c_cmaterial->td <= .01) return(0);
		return(!c_equiv(c1, &c_cmaterial->td_c));
	case MG_E_TS:
		if (c_cmaterial->ts <= .01) return(0);
		return(!c_equiv(c1, &c_cmaterial->ts_c));
	}
	return(0);
}


static int	/* new spectrum definition? */
newspecdef(C_COLOR *spc)
{
	static LUTAB	spc_tab = LU_SINIT(NULL,free);
	LUENT	*lp = lu_find(&spc_tab, (const char *)spc->client_data);

	if (lp == NULL)			/* should never happen */
		return(1);
	if (lp->data == NULL) {		/* new entry */
		lp->key = (char *)spc->client_data;
		lp->data = (char *)malloc(sizeof(C_COLOR));
	} else if (c_equiv(spc, (C_COLOR *)lp->data))
		return(0);		/* unchanged */

	if (lp->data != NULL)		/* else remember if we can */
		*(C_COLOR *)lp->data = *spc;
	return(1);			/* good as new */
}


char *
specolor(	/* check if color has spectra and output accordingly */
	COLOR	radrgb,
	C_COLOR	*clr,
	double	intensity
)
{
	static char	spname[128];
	double	mult;
	int	cbeg, cend, i;

	if (!dospectra | !(clr->flags & C_CDSPEC) | (intensity <= FTINY)) {
		cvtcolor(radrgb, clr, intensity);
		return(void_str);		/* just use RGB */
	}
	setcolor(radrgb, intensity, intensity, intensity);
	for (cbeg = 0; cbeg < C_CNSS; cbeg++)	/* trim zeros off beginning */
		if (clr->ssamp[cbeg])
			break;
	if (cbeg >= C_CNSS)			/* should never happen! */
		return(void_str);
	if (clr->client_data != NULL) {		/* get name if available */
		strcpy(spname, (char *)clr->client_data);
		strcat(spname, "*");		/* make sure it's special */
		if (!newspecdef(clr))		/* output already? */
			return(spname);
	} else
		strcpy(spname, sgen_str);
	c_ccvt(clr, C_CSEFF);			/* else output spectrum prim */
	for (cend = 0; !clr->ssamp[C_CNSS-1-cend]; cend++)
		;				/* trim zeros off end */
	fprintf(matfp, "\nvoid spectrum %s\n0\n0\n", spname);
	fprintf(matfp, "%d %d %d", C_CNSS+2-cbeg-cend,
		C_CMINWL+cbeg*C_CWLI, C_CMAXWL-cend*C_CWLI);
	mult = (C_CNSS*c_dfcolor.eff)/(clr->ssum*clr->eff);
	for (i = cbeg; i < C_CNSS-cend; i++) {
		if (!((i-cbeg+1)%6)) fputc('\n', matfp);
		fprintf(matfp, "\t%.5f", clr->ssamp[i]*mult);
	}
	fputc('\n', matfp);
	return(spname);
}


int
isgrey(				/* is RGB close match to grey? */
	COLOR rgb
)
{
	const double	yv = bright(rgb);
	double		diff2 = 0;
	int		i;

	for (i = 3; i--; ) {
		double	d = yv - colval(rgb,i);
		diff2 += d*d;
	}
	return(diff2 <= yv*yv*0.0025);
}


void
putrgbpat(			/* put out RGB pattern with given name */
	char *pnm,
	COLOR rgb
)
{
	fprintf(matfp, "\nvoid colorfunc %s\n", pnm);
	fprintf(matfp, "4 %f %f %f .\n0\n0\n", colval(rgb,RED),
			colval(rgb,GRN), colval(rgb,BLU));
}


char *
object(void)			/* return current object name */
{
	static char	objbuf[64];
	int	i;
	char	*cp;
	int	len;
						/* tracked by obj_handler */
	i = obj_nnames - sizeof(objbuf)/16;
	if (i < 0)
		i = 0;
	for (cp = objbuf; i < obj_nnames &&
		cp + (len=strlen(obj_name[i])) < objbuf+sizeof(objbuf)-1;
			i++, *cp++ = '.') {
		strcpy(cp, obj_name[i]);
		cp += len;
	}
	*cp = '\0';
	return(objbuf);
}


char *
addarg(				/* add argument and advance pointer */
	char *op,
	char *arg
)
{
	*op = ' ';
	while ( (*++op = *arg++) )
		;
	return(op);
}
