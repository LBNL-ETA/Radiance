#ifndef lint
static const char RCSid[] = "$Id: pabopto2bsdf.c,v 2.41 2021/04/05 21:06:32 greg Exp $";
#endif
/*
 * Load measured BSDF data in PAB-Opto format.
 * Assumes that surface-normal (Z-axis) faces into room unless -t option given.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "rtio.h"
#include "platform.h"
#include "bsdfrep.h"
				/* global argv[0] */
char			*progname;

typedef struct {
	const char	*fname;		/* input file path */
	double		theta, phi;	/* input angles (degrees) */
	double		up_phi;		/* azimuth for "up" direction */
	int		igp[2];		/* input grid position */
	int		isDSF;		/* data is DSF (rather than BSDF)? */
	int		nspec;		/* number of spectral samples */
	long		dstart;		/* data start offset in file */
} PGINPUT;

PGINPUT		*inpfile;	/* input files sorted by incidence */

int		rev_orient = 0;	/* shall we reverse surface orientation? */

double		lim_graze = 0;	/* limit scattering near grazing (deg above) */

/* Compare incident angles */
static int
cmp_indir(const void *p1, const void *p2)
{
	const PGINPUT	*inp1 = (const PGINPUT *)p1;
	const PGINPUT	*inp2 = (const PGINPUT *)p2;
	int		ydif = inp1->igp[1] - inp2->igp[1];
	
	if (ydif)
		return(ydif);

	return(inp1->igp[0] - inp2->igp[0]);
}

/* Assign grid position from theta and phi */
static void
set_grid_pos(PGINPUT *pip)
{
	FVECT	dv;

	if (pip->theta <= FTINY) {
		pip->igp[0] = pip->igp[1] = grid_res/2 - 1;
		return;
	}
	dv[2] = sin(M_PI/180.*pip->theta);
	dv[0] = cos(M_PI/180.*pip->phi)*dv[2];
	dv[1] = sin(M_PI/180.*pip->phi)*dv[2];
	dv[2] = sqrt(1. - dv[2]*dv[2]);
	pos_from_vec(pip->igp, dv);
}

/* Prepare a PAB-Opto input file by reading its header */
static int
init_pabopto_inp(const int i, const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	char	buf[2048];
	int	c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
	inpfile[i].fname = fname;
	inpfile[i].isDSF = -1;
	inpfile[i].nspec = 0;
	inpfile[i].up_phi = 0;
	inpfile[i].theta = inpfile[i].phi = -10001.;
				/* read header information */
	while ((c = getc(fp)) == '#' || c == EOF) {
		char	typ[64];
		if (fgets(buf, sizeof(buf), fp) == NULL) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			fclose(fp);
			return(0);
		}
		if (sscanf(buf, "sample_name \"%[^\"]\"", bsdf_name) == 1)
			continue;
		if (sscanf(buf, "colorimetry: %s", typ) == 1) {
			if (!strcasecmp(typ, "CIE-XYZ"))
				inpfile[i].nspec = 3;
			else if (!strcasecmp(typ, "CIE-Y"))
				inpfile[i].nspec = 1;
			continue;
		}
		if (sscanf(buf, "format: theta phi %s", typ) == 1) {
			if (!strcasecmp(typ, "DSF"))
				inpfile[i].isDSF = 1;
			else if (!strcasecmp(typ, "BSDF") ||
					!strcasecmp(typ, "BRDF") ||
					!strcasecmp(typ, "BTDF"))
				inpfile[i].isDSF = 0;
			continue;
		}
		if (sscanf(buf, "upphi %lf", &inpfile[i].up_phi) == 1)
			continue;
		if (sscanf(buf, "intheta %lf", &inpfile[i].theta) == 1)
			continue;
		if (sscanf(buf, "inphi %lf", &inpfile[i].phi) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&inpfile[i].theta, &inpfile[i].phi) == 2)
			continue;
	}
	inpfile[i].dstart = ftell(fp) - 1;
	fclose(fp);
	if (inpfile[i].isDSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		return(0);
	}
	if ((inpfile[i].theta < -10000.) | (inpfile[i].phi < -10000.)) {
		fputs(fname, stderr);
		fputs(": unknown incident angle\n", stderr);
		return(0);
	}
	if (rev_orient) {	/* reverse Z-axis to face outside */
		inpfile[i].theta = 180. - inpfile[i].theta;
		inpfile[i].phi = 360. - inpfile[i].phi;
	}
				/* convert to Y-up orientation */
	inpfile[i].phi += 90.-inpfile[i].up_phi;
				/* convert angle to grid position */
	set_grid_pos(&inpfile[i]);
	return(1);
}

/* Load a set of measurements corresponding to a particular incident angle */
static int
add_pabopto_inp(const int i)
{
	FILE		*fp = fopen(inpfile[i].fname, "r");
	double		theta_out, phi_out, val[3];
	int		n, c;
	
	if (fp == NULL || fseek(fp, inpfile[i].dstart, 0) == EOF) {
		fputs(inpfile[i].fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
					/* prepare input grid */
	if (!i || cmp_indir(&inpfile[i-1], &inpfile[i])) {
		if (i)			/* process previous incidence */
			make_rbfrep();
#ifdef DEBUG
		fprintf(stderr, "New incident (theta,phi)=(%.1f,%.1f)\n",
					inpfile[i].theta, inpfile[i].phi);
#endif
		if (inpfile[i].nspec)
			set_spectral_samples(inpfile[i].nspec);
		new_bsdf_data(inpfile[i].theta, inpfile[i].phi);
	}
#ifdef DEBUG
	fprintf(stderr, "Loading measurements from '%s'...\n", inpfile[i].fname);
#endif
					/* read scattering data */
	while (fscanf(fp, "%lf %lf %lf", &theta_out, &phi_out, val) == 3) {
		for (n = 1; n < inpfile[i].nspec; n++)
			if (fscanf(fp, "%lf", val+n) != 1) {
				fprintf(stderr, "%s: warning: unexpected EOF\n",
						inpfile[i].fname);
				fclose(fp);
				return(1);
			}
					/* check if scatter angle is too low */
		if (fabs(theta_out - 90.) < lim_graze-FTINY)
			continue;
		if (rev_orient) {	/* reverse Z-axis to face outside? */
			theta_out = 180. - theta_out;
			phi_out = 360. - phi_out;
		}
		phi_out += 90.-inpfile[i].up_phi;
		add_bsdf_data(theta_out, phi_out, val, inpfile[i].isDSF);
	}
	n = 0;
	while ((c = getc(fp)) != EOF)
		n += !isspace(c);
	if (n) 
		fprintf(stderr,
			"%s: warning: %d unexpected characters past EOD\n",
				inpfile[i].fname, n);
	fclose(fp);
	return(1);
}

#ifndef TEST_MAIN

#define	SYM_ILL		'?'		/* illegal symmetry value */
#define	SYM_ISO		'I'		/* isotropic */
#define	SYM_QUAD	'Q'		/* quadrilateral symmetry */
#define	SYM_BILAT	'B'		/* bilateral symmetry */
#define	SYM_ANISO	'A'		/* anisotropic */
#define SYM_UP		'U'		/* "up-down" (180°) symmetry */

static const char	quadrant_rep[16][16] = {
				"in-plane","0-90","90-180","0-180",
				"180-270","0-90+180-270","90-270",
				"0-270","270-360","270-90",
				"90-180+270-360","270-180","180-360",
				"180-90","90-360","0-360"
			};
static const char	quadrant_sym[16] = {
				SYM_ISO, SYM_QUAD, SYM_QUAD, SYM_BILAT,
				SYM_QUAD, SYM_ILL, SYM_BILAT, SYM_ILL,
				SYM_QUAD, SYM_BILAT, SYM_ILL, SYM_ILL,
				SYM_BILAT, SYM_ILL, SYM_ILL, SYM_ANISO
			};

/* Read in PAB-Opto BSDF files and output RBF interpolant */
int
main(int argc, char *argv[])
{
	extern int	nprocs;
	static char	gval_buf[16];
	char *		auto_grazing = NULL;
	const char	*symmetry = "0";
	int		ninpfiles, totinc;
	int		a, i;

	progname = argv[0];			/* get options */
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 't':
			rev_orient = !rev_orient;
			break;
		case 'n':
			nprocs = atoi(argv[++a]);
			break;
		case 's':
			symmetry = argv[++a];
			break;
		case 'g':
			if (toupper(argv[a+1][0]) == 'A')
				auto_grazing = argv[a+1] = gval_buf;
			else
				lim_graze = atof(argv[a+1]);
			++a;
			break;
		default:
			goto userr;
		}
	totinc = ninpfiles = argc - a;		/* initialize & sort inputs */
	if (ninpfiles < 2)
		goto userr;
	if (toupper(symmetry[0]) == SYM_UP)	/* special case for "up" symmetry */
		totinc += ninpfiles;
	inpfile = (PGINPUT *)malloc(sizeof(PGINPUT)*totinc);
	if (inpfile == NULL)
		return(1);
	if (auto_grazing)
		lim_graze = 90.;
	for (i = 0; i < ninpfiles; i++) {
		if (!init_pabopto_inp(i, argv[a+i]))
			return(1);
		if (auto_grazing && fabs(inpfile[i].theta - 90.) < lim_graze)
			lim_graze = fabs(inpfile[i].theta - 90.);
	}
	if (auto_grazing)
		sprintf(auto_grazing, "%.2f", lim_graze);
	for (i = ninpfiles; i < totinc; i++) {	/* copy for "up" symmetry */
		inpfile[i] = inpfile[i-ninpfiles];
		inpfile[i].phi += 180.;		/* invert duplicate data */
		inpfile[i].up_phi -= 180.;
		set_grid_pos(&inpfile[i]);	/* grid location for sorting */
	}
	qsort(inpfile, totinc, sizeof(PGINPUT), cmp_indir);
						/* compile measurements */
	for (i = 0; i < totinc; i++)
		if (!add_pabopto_inp(i))
			return(1);
	make_rbfrep();				/* process last data set */
						/* check input symmetry */
	switch (toupper(symmetry[0])) {
	case '0':				/* unspecified symmetry */
		if (quadrant_sym[inp_coverage] != SYM_ILL)
			break;			/* anything legal goes */
		fprintf(stderr, "%s: unsupported phi coverage (%s)\n",
				progname, quadrant_rep[inp_coverage]);
		return(1);
	case SYM_UP:				/* faux "up" symmetry */
		if (quadrant_sym[inp_coverage] == SYM_ANISO)
			break;
		/* fall through */
	case SYM_ISO:				/* usual symmetry types */
	case SYM_QUAD:
	case SYM_BILAT:
	case SYM_ANISO:
		if (quadrant_sym[inp_coverage] == toupper(symmetry[0]))
			break;			/* matches spec */
		fprintf(stderr,
		"%s: phi coverage (%s) does not match requested '%s' symmetry\n",
				progname, quadrant_rep[inp_coverage], symmetry);
		return(1);
	default:
		fprintf(stderr,
  "%s: -s option must be Isotropic, Quadrilateral, Bilateral, Up, or Anisotropic\n",
				progname);
		return(1);
	}
#ifdef DEBUG
	fprintf(stderr, "Input phi coverage (%s) has '%c' symmetry\n",
				quadrant_rep[inp_coverage],
				quadrant_sym[inp_coverage]);
#endif
	build_mesh();				/* create interpolation */
	SET_FILE_BINARY(stdout);		/* start header */
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputnow(stdout);
	save_bsdf_rep(stdout);			/* complete header + data */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-t][-n nproc][-s symmetry][-g angle|'A'] meas1.dat meas2.dat .. > bsdf.sir\n",
					progname);
	return(1);
}

#else		/* TEST_MAIN */

/* Test main produces a Radiance model from the given input file */
int
main(int argc, char *argv[])
{
	PGINPUT	pginp;
	char	buf[128];
	FILE	*pfp;
	double	bsdf, min_log;
	FVECT	dir;
	int	i, j, n;

	progname = argv[0];
	if (argc != 2) {
		fprintf(stderr, "Usage: %s input.dat > output.rad\n", progname);
		return(1);
	}
	inpfile = &pginp;
	if (!init_pabopto_inp(0, argv[1]) || !add_pabopto_inp(0))
		return(1);
						/* reduce data set */
	if (make_rbfrep() == NULL) {
		fprintf(stderr, "%s: nothing to plot!\n", progname);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "Minimum BSDF = %.4f\n", bsdf_min);
#endif
	min_log = log(bsdf_min*.5 + 1e-5);
#if 1						/* produce spheres at meas. */
	puts("void plastic yellow\n0\n0\n5 .6 .4 .01 .04 .08\n");
	n = 0;
	for (i = 0; i < grid_res; i++)
	    for (j = 0; j < grid_res; j++)
		if (dsf_grid[i][j].sum.n > 0) {
			ovec_from_pos(dir, i, j);
			bsdf = dsf_grid[i][j].sum.v /
			   ((double)dsf_grid[i][j].sum.n*output_orient*dir[2]);
			if (bsdf <= bsdf_min*.6)
				continue;
			bsdf = log(bsdf + 1e-5) - min_log;
			ovec_from_pos(dir, i, j);
			printf("yellow sphere s%04d\n0\n0\n", ++n);
			printf("4 %.6g %.6g %.6g %.6g\n\n",
					dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf,
					.007*bsdf);
		}
#endif
#if 1						/* spheres at RBF peaks */
	puts("void plastic red\n0\n0\n5 .8 .01 .01 .04 .08\n");
	for (n = 0; n < dsf_list->nrbf; n++) {
		RBFVAL	*rbf = &dsf_list->rbfa[n];
		ovec_from_pos(dir, rbf->gx, rbf->gy);
		bsdf = eval_rbfrep(dsf_list, dir);
		bsdf = log(bsdf + 1e-5) - min_log;
		printf("red sphere p%04d\n0\n0\n", ++n);
		printf("4 %.6g %.6g %.6g %.6g\n\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf,
				.011*bsdf);
	}
#endif
#if 1						/* output continuous surface */
	puts("void trans tgreen\n0\n0\n7 .7 1 .7 .04 .04 .9 1\n");
	fflush(stdout);
	sprintf(buf, "gensurf tgreen bsdf - - - %d %d", grid_res-1, grid_res-1);
	pfp = popen(buf, "w");
	if (pfp == NULL) {
		fprintf(stderr, "%s: cannot open '| %s'\n", progname, buf);
		return(1);
	}
	for (i = 0; i < grid_res; i++)
	    for (j = 0; j < grid_res; j++) {
		ovec_from_pos(dir, i, j);
		bsdf = eval_rbfrep(dsf_list, dir);
		bsdf = log(bsdf + 1e-5) - min_log;
		fprintf(pfp, "%.8e %.8e %.8e\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
	    }
	if (pclose(pfp) != 0)
		return(1);
#endif
	return(0);
}

#endif		/* TEST_MAIN */
